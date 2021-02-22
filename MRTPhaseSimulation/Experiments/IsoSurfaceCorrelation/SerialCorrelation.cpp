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
    void calc(size_t offset, size_t lag_max);
public:
    SerialCorrStats() = default;
    void add(double first_return_time);
    IsoSurfaceCorr get_corr(const string& isoSurfaceName, size_t offset, size_t lag_max);
};

struct DimError : public exception {
    const char* what () const throw() {
        return "vector dimensions are incompatible with lag and offset";
    }
};

const size_t OFFSET = 500;
const size_t LAGS = 15;
const size_t TOT_INTERVAL_COUNT = 10000;

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

    // add spoke of a wheel for reference purposes
    InterpolatedCurve &spokeCurve = isoSurfaceFile.createInterpolatedCurve("Spoke");
    auto firstElementIt = isoSurfaceFile.begin();
    auto[rho, phi] = (*firstElementIt).get_nodes();
    phi = vector<double>(phi.size(), 0.0);
    spokeCurve.set_nodes(rho, phi);
    auto pSet = (*firstElementIt).get_parameterSet();
    for (auto e: pSet) {
        spokeCurve.add_parameter(e.first.data(), e.second);
    }
    spokeCurve.set_omegaBar((*firstElementIt).get_omegaBar());

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
        double t = simConfig.t0;
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
            t = simConfig.t0;
        }

        isoSurfaceCorr = serialCorrStats.get_corr(isoSurface.get_name(), OFFSET, LAGS);
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

void SerialCorrStats::calc(size_t offset, size_t lag_max) {
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
}

IsoSurfaceCorr SerialCorrStats::get_corr(const string& isoSurfaceName, size_t offset, size_t lag_max) {
    IsoSurfaceCorr isoSurfaceCorr(isoSurfaceName);
    isoSurfaceCorr.N = FRTdata.size();
    isoSurfaceCorr.offset = offset;
    this->calc(offset, lag_max);
    isoSurfaceCorr.rho_k = this->rho_k;
    return isoSurfaceCorr;
}

