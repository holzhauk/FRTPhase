#include <iostream>
#include <vector>
#include <filesystem>
#include <cmath>
#include <string>
#include <exception>

#include <libSDEToolbox/libSDEToolbox.h>
#include <libSPhaseFile/libSPhaseFile.h>

using namespace std;
namespace fs = std::filesystem;

class SerialCorrStats {
private:
    vector<double> FRTdata;
    vector<double> rho_k;
    vector<double> StDev_rho_k;
    double cv;
    double StDev_cv;
    void calc(size_t offset, size_t lag_max, size_t sub_pop_size);
public:
    SerialCorrStats() = default;
    void add(double first_return_time);
    IsoSurfaceCorr get_corr(const string& isoSurfaceName, size_t offset, size_t lag_max, size_t sub_pop_size);
};

struct DimError : public exception {
    const char* what () const throw() {
        return "vector dimensions are incompatible with lag and offset";
    }
};

const size_t OFFSET = 500;
const size_t LAGS = 20;
const size_t TOT_INTERVAL_COUNT = 10000;
const size_t SUB_POP_SIZE_STATS = 10;

int main(int argc, char* argv[]) {

    if (argc != 2) {
        std::cerr << "Argument missing! Usage: "
                  << argv[0] << " <ConfigurationFilePath>.json" << std::endl;
        return EXIT_FAILURE;
    }
    fs::path config_file_path = fs::path(argv[1]);

    SimConfigFile config;
    config.read(config_file_path);
    const string modelName = config.get_modelName();
    Config::Simulation simConfig;
    simConfig = config.get_simConfig();
    IsoSurfaceFile isoSurfaceFile(modelName);
    isoSurfaceFile.read(config.get_inPath());

    SerialCorrFile serialCorrFile = SerialCorrFile(modelName,
                                                   config.get_inPath(),
                                                   config_file_path);
    ModelFactory theModelFactory;

    for (auto isoSurface: isoSurfaceFile) {

        auto[rho_min, rho_max] = isoSurface.get_extensions();
        ReflectiveAnnulus domain(rho_min, rho_max);

        unique_ptr<IsotropicPlanarSSDE> modelPtr =
                theModelFactory.createModel(config.get_modelName(), isoSurface.get_parameterSet());
        ItoEulerIntegrator integrator(&domain, modelPtr.get());
        SDEIntegrator::config_t integratorConfig = SimConfig2IntegratorConfig(simConfig);
        integrator.configure(integratorConfig);

        IsoSurfaceCorr& isoSurfaceCorr = serialCorrFile.create_isoSurfaceCorr(isoSurface.get_name());

        double omegaBar = isoSurface.get_omegaBar();
        if (omegaBar == 0.0)
            throw FRTDetectorNoSenseOfRotation();

        bool is_positive_sense_of_rotation = (omegaBar > 0.0);
        SerialCorrStats serialCorrStats;
        Pos_t x = isoSurface.get_random_point();
        double t = 0.0;
        for (size_t count = 0; count < TOT_INTERVAL_COUNT; count++){
            integrator.reset(x, t);
            while(!isoSurface.is_first_return_event(x, is_positive_sense_of_rotation)){
                tie(x, t) = integrator.evolve();
            }
            serialCorrStats.add(t);
            if (is_positive_sense_of_rotation)
                x[1] = x[1] - 2 * M_PI;
            else
                x[1] = x[1] + 2 * M_PI;
            t = 0.0;
        }

        isoSurfaceCorr = serialCorrStats.get_corr(isoSurface.get_name(), OFFSET, LAGS, SUB_POP_SIZE_STATS);
        cout << "rho_k: ";
        for (auto rho: isoSurfaceCorr.rho_k)
            cout << rho << ", ";
        cout << endl;

    }

    serialCorrFile.write(config.get_outPath());

    return 0;
}

void SerialCorrStats::add(double first_return_time) {
    FRTdata.push_back(first_return_time);
}

void SerialCorrStats::calc(size_t offset, size_t lag_max, size_t sub_pop_size) {
    // calculate values
    if ((offset + lag_max) > FRTdata.size())
        throw DimError();
    double mean = 0.0;
    vector<double> c_k(lag_max + 1, 0.0);
    rho_k = vector<double> (lag_max + 1, 0.0);

    for (size_t i = offset; i < FRTdata.size(); i++)
        mean += FRTdata[i];
    mean /= (FRTdata.size() - offset);

    for (size_t k = 0; k < lag_max + 1; k++) {
        for (size_t i = offset; i < FRTdata.size() - k; i++)
            c_k[k] += FRTdata[i] * FRTdata[i + k];

        c_k[k] /= (FRTdata.size() - offset - k);
        c_k[k] -= mean*mean;
        rho_k[k] = c_k[k] / c_k[0];
    }
    cv = sqrt(c_k[0]) / mean;

    // calculate errors
    if (sub_pop_size > FRTdata.size())
        throw DimError();

    StDev_rho_k = vector<double>(lag_max + 1, 0.0);
    double cv_mean = 0.0; // first moment
    double cv_2 = 0.0; // second_moment
    vector<double> rho_k_mean(lag_max + 1, 0.0); // first moment
    vector<double> rho_k_2(lag_max + 1, 0.0); // second moment

    size_t NoSubInt  = floor((FRTdata.size()-offset) / sub_pop_size);
    for (size_t p = 0; p < sub_pop_size; p++){

        double mean_sub = 0.0;
        double cv_sub = 0.0;
        vector<double> c_k_sub(lag_max + 1, 0.0);
        vector<double> rho_k_sub(lag_max + 1, 0.0);

        for (size_t i = 0; i < NoSubInt; i++)
            mean_sub += FRTdata[offset + p*(NoSubInt) + i];
        mean_sub /= NoSubInt;

        for (size_t k = 0; k < lag_max + 1; k++) {
            for (size_t i = 0; i < NoSubInt - k; i++)
                c_k_sub[k] += FRTdata[offset + p*(NoSubInt) + i] * FRTdata[offset + p*(NoSubInt) + i + k];

            c_k_sub[k] /= (NoSubInt - k);
            c_k_sub[k] -= mean_sub*mean_sub;
            rho_k_sub[k] = c_k_sub[k] / c_k_sub[0];
        }
        cv_sub = sqrt(c_k_sub[0]) / mean_sub;

        // update first and second moment
        cv_mean += cv_sub;
        cv_2 += pow(cv_sub, 2.0);

        for (size_t j = 0; j < rho_k_mean.size(); j++){
            rho_k_mean[j] += rho_k_sub[j];
            rho_k_2[j] += pow(rho_k_sub[j], 2.0);
        }
    }

    cv_mean /= sub_pop_size;
    cv_2 /= sub_pop_size;
    StDev_cv = sqrt(cv_2 - pow(cv_2, 2.0)) / sqrt(sub_pop_size);

    for (size_t j = 0; j < StDev_rho_k.size(); j ++){
        rho_k_mean[j] /= sub_pop_size;
        rho_k_2[j] /= sub_pop_size;
        StDev_rho_k[j] = sqrt(rho_k_2[j] - pow(rho_k_mean[j], 2.0)) / sqrt(sub_pop_size);
    }

}

IsoSurfaceCorr SerialCorrStats::get_corr(const string& isoSurfaceName, size_t offset, size_t lag_max, size_t sub_pop_size) {
    IsoSurfaceCorr isoSurfaceCorr(isoSurfaceName);
    isoSurfaceCorr.N = FRTdata.size();
    isoSurfaceCorr.offset = offset;
    isoSurfaceCorr.sub_pop_size = sub_pop_size;
    this->calc(offset, lag_max, sub_pop_size);
    isoSurfaceCorr.cv = this->cv;
    isoSurfaceCorr.StDev_cv = this->StDev_cv;
    isoSurfaceCorr.rho_k = this->rho_k;
    isoSurfaceCorr.StDev_rho_k = this->StDev_rho_k;
    return isoSurfaceCorr;
}

