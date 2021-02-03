#include <iostream>
#include <vector>
#include <array>
#include <filesystem>
#include <cmath>
#include <string>
#include <list>
#include <exception>

#include <H5Cpp.h>
#include <mpi.h>

#include <libSDEToolbox/libSDEToolbox.h>
#include <libSPhaseFile/libSPhaseFile.h>
#include <libMPIFunctions/libMPIFunctions.h>

using namespace std;
namespace fs = std::filesystem;

class NewbySchwemmer: public IsotropicPlanarSSDEAdditiveNoise{
public:
    explicit NewbySchwemmer(ParameterSet pSet):
            IsotropicPlanarSSDEAdditiveNoise(pSet["D"], pSet["D"]) {
        this->pSet = pSet;
    };

    double g(double& rho) override{
        return -pSet["gamma"]*rho*(pow(rho, 2.0) - 1.0) + pSet["D"] / rho;
    };

    double f(double& rho) override{
        return pSet["omega"]*(1 + pSet["gamma"]*pSet["c"]*pow((1.0 - rho), 2.0));
    };
};

struct DimError : public exception {
    const char* what () const throw() {
        return "vector dimensions are incompatible";
    }
};

struct TSerialCorrStats {
    vector<double> Ti; // FRTs
    vector<double> TkT; // corr. btw. FRTs
    struct Corr_t {
        double mTau = 0.0;
        double varTau = 0.0;
        vector<double> rho_k;
    };
    friend Corr_t Dist_Stats(const TSerialCorrStats& stats,
                             const Config::Simulation& simConfig);
    TSerialCorrStats(size_t noLags) {
        Ti = vector<double>(noLags + 1, 0.0);
        TkT = vector<double>(noLags + 1, 0.0);
    };
    void add(const vector<double>& kth_return_times) {
        if (kth_return_times.size() != Ti.size())
            throw DimError();
        vector<double> corr_k;
        for (auto tau_k: kth_return_times){
            corr_k.push_back(tau_k*kth_return_times[0]);
        }
        for (int i = 0; i < Ti.size(); i++){
            Ti[i] += kth_return_times[i];
            TkT[i] += corr_k[i];
        }
    };
    Corr_t get_corr(size_t ensemble_size) {
        Corr_t corr;
        corr.mTau = this->Ti[0] / ensemble_size;
        double TauSq = this->TkT[0] / ensemble_size;
        corr.varTau = TauSq - pow(corr.mTau, 2);
        for (int k = 1; k < Ti.size(); k++){
            corr.rho_k.push_back(((TkT[k] - Ti[k]*corr.mTau) / ensemble_size) / corr.varTau );
        }
        return corr;
    };
};

class CorrelationFile;

struct IsoSurfaceCorr {
private:
    friend class CorrelationFile;
    string key;
    array<vector<double>, 2> x0;
    vector<double> mT;
    vector<double> varT;
    vector<vector<double>> corr_k; // serial correlations with lag k
public:
    IsoSurfaceCorr(const string& isoSurfaceName): key(isoSurfaceName) {};
    void add(const Pos_t& x0, const TSerialCorrStats::Corr_t& corr);
    string get_key();
};

class CorrelationFile: public SPhaseFile{
private:
    string modelName;
    fs::path isoSurfaceFilePath;
    fs::path configFilePath;
    list<unique_ptr<IsoSurfaceCorr>> corrList;
protected:
    void write_body(H5::H5File& file);
    void read_body(H5::H5File& file);
public:
    CorrelationFile(string modelName): modelName(modelName), SPhaseFile("CorrelationFile") {};
    CorrelationFile(string modelName,
                    const fs::path& isoSurfaceFilePath,
                    const fs::path& configFilePath):
                        modelName(modelName),
                        isoSurfaceFilePath(isoSurfaceFilePath),
                        configFilePath(configFilePath),
                        SPhaseFile("CorrelationFile"){};
    IsoSurfaceCorr& create_isoSurfaceCorr(const string& isoSurfaceName);
    CorrelationFile& operator = (const CorrelationFile& other);
};

const size_t lags = 10;
const string modelName = "NewbySchwemmer";

int main(int argc, char* argv[]) {

    if (argc != 2) {
        std::cerr << "Argument missing! Usage: "
                  << argv[0] << " <ConfigurationFilePath>.json" << std::endl;
        return EXIT_FAILURE;
    }
    fs::path config_file_path = fs::path(argv[1]);
    /*
     * initialize MPI
     */
    int world_rank, world_size;
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    bool is_master = (world_rank == 0);

    SimConfigFile config;
    Config::Simulation simConfig;
    IsoSurfaceFile isoSurfaceFile(modelName);
    CorrelationFile corrFile(modelName);
    if (is_master){
        config.read(config_file_path);
        simConfig = config.get_simConfig();
        isoSurfaceFile.read(config.get_inPath());

        // add spoke of a wheel for reference purposes
        InterpolatedCurve& spokeCurve = isoSurfaceFile.createInterpolatedCurve("Spoke");
        IsoSurfaceFile::IsoSurfaceFileIt firstElementIt = isoSurfaceFile.begin();
        auto [rho, phi] = (*firstElementIt).get_nodes();
        phi = vector<double>(phi.size(), 0.0);
        spokeCurve.set_nodes(rho, phi);
        auto pSet = (*firstElementIt).get_parameterSet();
        for (auto e: pSet){
            spokeCurve.add_parameter(e.first.data(), e.second);
        }

        corrFile = CorrelationFile(modelName, config.get_inPath(), config_file_path);
    }
    MPI_Share(world_rank, simConfig);
    MPI_Share(world_rank, isoSurfaceFile);

    EquidistantSampler sampler;
    unsigned int subEnsembleSize = SubEnsembleSize(world_rank, world_size, simConfig.EnsembleSize);

    for (auto isoSurface: isoSurfaceFile){

        auto [rho_min, rho_max] = isoSurface.get_extensions();
        ReflectiveAnnulus domain(rho_min, rho_max);

        NewbySchwemmer model(isoSurface.get_parameterSet());
        HeunIntegrator integrator(&domain, &model);
        SDEIntegrator::config_t integratorConfig = SimConfig2IntegratorConfig(simConfig);
        integrator.configure(integratorConfig);

        vector<array<double, 2>> x0s = sampler.get_samples(simConfig.SampleSize, isoSurface);

        IsoSurfaceCorr& isoSurfaceCorr = corrFile.create_isoSurfaceCorr(isoSurface.get_name());
        for (auto x0: x0s){

            TSerialCorrStats tSerialCorrStats(lags);
            for (size_t e = 0; e < subEnsembleSize; e++){

                integrator.reset(x0, simConfig.t0);
                Pos_t x = x0;
                double tau = simConfig.t0;

                vector<double> kth_RTs (lags + 1, simConfig.t0); // kth return times
                for (auto& tau_k: kth_RTs){

                    bool has_passed = false;
                    while (!has_passed && integrator.is_in_time()){
                        if(isoSurface.is_first_return_event(x, true)){
                            has_passed = true;
                            x[1] = x[1] - 2*M_PI;
                            integrator.reset(x, simConfig.t0);
                        } else {
                            if (isoSurface.is_first_return_event(x, false)){
                                has_passed = true;
                                x[1] = x[1] + 2*M_PI;
                                integrator.reset(x, simConfig.t0);
                            } else {
                                tie(x, tau) = integrator.evolve();
                            }
                        }
                    }
                    tau_k = tau;
                    tau = simConfig.t0;
                }
                tSerialCorrStats.add(kth_RTs);
            }
            TSerialCorrStats::Corr_t corr_k = Dist_Stats(tSerialCorrStats, simConfig);
            isoSurfaceCorr.add(x0, corr_k);
        }
    }
    if (is_master)
        corrFile.write(config.get_outPath());
    MPI_Finalize();
    return 0;
}

TSerialCorrStats::Corr_t Dist_Stats(const TSerialCorrStats& stats, const Config::Simulation& simConfig){
    vector<double> Sum_Ti (stats.Ti.size());
    vector<double> Sum_TkT (stats.TkT.size());
    MPI_Allreduce(stats.Ti.data(), Sum_Ti.data(), stats.Ti.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(stats.TkT.data(), Sum_TkT.data(), stats.TkT.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    TSerialCorrStats::Corr_t corr;
    corr.mTau = Sum_Ti[0] / simConfig.EnsembleSize;
    double TauSq = Sum_TkT[0] / simConfig.EnsembleSize;
    corr.varTau = TauSq - pow(corr.mTau, 2);
    for (int k = 1; k < Sum_Ti.size(); k++){
        corr.rho_k.push_back(((Sum_TkT[k] - Sum_Ti[k]*corr.mTau) / simConfig.EnsembleSize) / corr.varTau );
    }
    return corr;
}

void IsoSurfaceCorr::add(const Pos_t& init_pos, const TSerialCorrStats::Corr_t& corr) {
    this->x0[0].push_back(init_pos[0]);
    this->x0[1].push_back(init_pos[1]);
    this->corr_k.push_back(corr.rho_k);
    this->mT.push_back(corr.mTau);
    this->varT.push_back(corr.varTau);
}

string IsoSurfaceCorr::get_key() {
    return string(key);
}

void CorrelationFile::read_body(H5::H5File &file) {

    H5std_string str_buf;

    H5::Attribute attr = file.openAttribute("model");
    H5::DataType dtype = attr.getDataType();
    attr.read(dtype, str_buf);
    if (str_buf != modelName)
        throw SPhaseFileModelConflict();
    dtype.close();
    attr.close();

    H5::DataSet dset = file.openDataSet("configuration_file");
    dtype = dset.getDataType();
    dset.read(str_buf, dtype);
    configFilePath = string(str_buf);
    dtype.close();
    dset.close();

    dset = file.openDataSet("isosurface_file");
    dtype = dset.getDataType();
    dset.read(str_buf, dtype);
    isoSurfaceFilePath = string(str_buf);
    dtype.close();
    dset.close();

}

void CorrelationFile::write_body(H5::H5File &file) {

    // define attribute types
    H5::StrType str_type(0, H5T_VARIABLE);
    str_type.setCset(H5T_CSET_UTF8);
    str_type.setStrpad(H5T_STR_NULLTERM);
    H5::DataSpace dspace(H5S_SCALAR);

    // write model name attribute
    H5::Attribute attr = file.createAttribute("model", str_type, dspace);
    attr.write(str_type, modelName);
    attr.close();

    // write configuration file path
    H5::DataSet dset = file.createDataSet("configuration_file", str_type, dspace);
    dset.write(configFilePath, str_type);
    dset.close();

    // write isosurface file path
    dset = file.createDataSet("isosurface_file", str_type, dspace);
    dset.write(isoSurfaceFilePath, str_type);
    dset.close();

    dspace.close();
    str_type.close();

    for (auto& isoSurfaceCorrPtr: this->corrList) {

        size_t size = isoSurfaceCorrPtr->x0[0].size();
        if ((isoSurfaceCorrPtr->mT.size() != size) ||
            (isoSurfaceCorrPtr->varT.size() != size) ||
            (isoSurfaceCorrPtr->x0[1].size() != size) ||
            (isoSurfaceCorrPtr->corr_k.size() != size)) throw SPhaseFileDimError();

        H5::Group surface = file.createGroup(isoSurfaceCorrPtr->get_key());

        // INITIAL POSITIONS
        H5::Group xinit_g = surface.createGroup("initial_position");
        hsize_t dim_size[] = {isoSurfaceCorrPtr->x0[0].size()};
        dspace = H5::DataSpace(1, dim_size);
        H5::DataSet rhos = xinit_g.createDataSet("rho",
                                                 H5::PredType::NATIVE_DOUBLE,
                                                 dspace);
        rhos.write(isoSurfaceCorrPtr->x0[0].data(), H5::PredType::NATIVE_DOUBLE);
        rhos.close();
        H5::DataSet phis = xinit_g.createDataSet("phi",
                                                 H5::PredType::NATIVE_DOUBLE,
                                                 dspace);
        phis.write(isoSurfaceCorrPtr->x0[1].data(), H5::PredType::NATIVE_DOUBLE);
        phis.close();
        xinit_g.close();

        // FIRST RETURN TIMES FRT
        H5::Group frt_g = surface.createGroup("first_return_time");
        H5::DataSet mFRTs = frt_g.createDataSet("mFRT",
                                                H5::PredType::NATIVE_DOUBLE,
                                                dspace);
        mFRTs.write(isoSurfaceCorrPtr->mT.data(), H5::PredType::NATIVE_DOUBLE);
        mFRTs.close();
        H5::DataSet varFRTs = frt_g.createDataSet("varFRT",
                                                  H5::PredType::NATIVE_DOUBLE,
                                                  dspace);
        varFRTs.write(isoSurfaceCorrPtr->varT.data(), H5::PredType::NATIVE_DOUBLE);
        varFRTs.close();

        // serial correlation coefficients
        H5::Group sCorr_g = surface.createGroup("serial_correlation_coefficients");
        const int rank = 2;
        hsize_t dimsf[rank];
        dimsf[0] = isoSurfaceCorrPtr->corr_k.size();
        if (dimsf[0] != 0)
            dimsf[1] = isoSurfaceCorrPtr->corr_k[0].size();
        else
            dimsf[1] = 0;
        H5::DataSpace dspace(rank, dimsf);
        H5::DataSet rho_k = sCorr_g.createDataSet("rho_k", H5::PredType::NATIVE_DOUBLE, dspace);
        double* corr_k_array = new double [dimsf[0]*dimsf[1]];
        double index = 0.0;
        for (int ii = 0; ii < dimsf[0]; ii++){
            for (int kk = 0; kk < dimsf[1]; kk++){
                corr_k_array[ii*dimsf[1] + kk] = isoSurfaceCorrPtr->corr_k[ii][kk];
            }
        }
        rho_k.write(corr_k_array, H5::PredType::NATIVE_DOUBLE);
        delete [] corr_k_array;
        rho_k.close();
        dspace.close();
        sCorr_g.close();

        surface.close();

    }

}

IsoSurfaceCorr& CorrelationFile::create_isoSurfaceCorr(const string &isoSurfaceName) {
    unique_ptr<IsoSurfaceCorr> corr_ptr (new IsoSurfaceCorr(isoSurfaceName));
    corrList.push_back(move(corr_ptr));
    return *corrList.back();
}

CorrelationFile& CorrelationFile::operator = (const CorrelationFile& other) {
    this->modelName = other.modelName;
    this->isoSurfaceFilePath = other.isoSurfaceFilePath;
    this->configFilePath = other.configFilePath;
    for (auto& corr: other.corrList) {
        unique_ptr<IsoSurfaceCorr> corr_cpy (new IsoSurfaceCorr(other.modelName));
        *corr_cpy = *corr;
        this->corrList.push_back(move(corr_cpy));
    }
    return *this;
}