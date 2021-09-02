//
// Created by konstantin on 1/3/21.
//
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE libSPhaseFile test
#include <filesystem>
#include <string>
#include <random>
#include <H5Cpp.h>
#include <boost/test/unit_test.hpp>
#include <libSPhaseFile/libSPhaseFile.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace fs = std::filesystem;

const string modelName = "test_model";
void write_headerOnly_IsoSurfaceFile(fs::path p,
                                     string class_attribute,
                                     string format_attribute,
                                     string version_attribute,
                                     string modelName);

void write_headerOnly_DataFile(fs::path p,
                               string class_attribute,
                               string format_attribute,
                               string version_attribute,
                               string modelName,
                               fs::path configFilePath,
                               fs::path isoSurfaceFilePath);

fs::path TestPath = "./SPhaseFileTest.h5";

BOOST_AUTO_TEST_SUITE(InterpolatedCurveTest)

BOOST_AUTO_TEST_CASE(assignment_test){
    string name = "Surface0";
    InterpolatedCurve c1(name);
    string alternate_name = "Surface1";
    InterpolatedCurve c2(alternate_name);
    c2.add_parameter("D1", 0.1);
    array<double, 2> x = {10.0, 20.0};
    c2.add_node(x);
    c1 = c2;
    BOOST_CHECK(c1.get_name() == alternate_name);
    auto [rho1, phi1] = c1.get_nodes();
    BOOST_REQUIRE(rho1.size() == phi1.size());
    BOOST_REQUIRE(rho1.size() == 1);
    BOOST_CHECK_EQUAL(rho1[0], 10.0);
    BOOST_CHECK_EQUAL(phi1[0], 20.0);
    auto [keys, vals] = c1.get_parameters();
    BOOST_REQUIRE(keys.size() == vals.size());
    BOOST_REQUIRE(keys.size() == 1);
    BOOST_CHECK_EQUAL(keys.front(), "D1");
    BOOST_CHECK_EQUAL(vals.front(), 0.1);
}

BOOST_AUTO_TEST_CASE(sample_random_point_from_curve){

    vector<double> rho_curve = {0.5, 0.5909, 0.6818, 0.77272, 0.86363,
                                0.9545, 1.04545, 1.136363, 1.227272, 1.31818,
                                1.40909, 1.5};
    vector<double> phi_curve = {0.0, -0.1, -0.3, -0.4, -0.42, -0.41, -0.2, 0.1,
                                0.3, 0.4, 1.1, 1.7};
    InterpolatedCurve curve("IsoCurve");
    curve.set_nodes(rho_curve, phi_curve);
    auto x_rand = curve.get_random_point();
    BOOST_REQUIRE(x_rand[0] >= 0.5);
    BOOST_REQUIRE(x_rand[0] <= 1.5);
    BOOST_REQUIRE(x_rand[1] >= -0.42);
    BOOST_REQUIRE(x_rand[1] <= 1.7);

}

BOOST_AUTO_TEST_CASE(detect_first_return_event_test){

    vector<double> rho_curve = {0.5, 0.5909, 0.6818, 0.77272, 0.86363,
                                0.9545, 1.04545, 1.136363, 1.227272, 1.31818,
                                1.40909, 1.5 };
    vector<double> phi_curve = {0.0, -0.1, -0.3, -0.4, -0.42, -0.41, -0.2, 0.1,
                                0.3, 0.4, 1.1, 1.7};
    InterpolatedCurve curve("IsoCurve");
    curve.set_nodes(rho_curve, phi_curve);
    auto [rho, phi] = curve.get_nodes();
    BOOST_CHECK(rho == rho_curve);
    BOOST_CHECK(phi == phi_curve);

    array<double, 2> x = {0.7, 0.3};
    BOOST_CHECK(!curve.is_first_return_event(x, true));
    BOOST_CHECK(!curve.is_first_return_event(x, false));

    x = {1.2, 2*M_PI + 0.1};
    BOOST_CHECK(curve.is_first_return_event(x, true));
    BOOST_CHECK(!curve.is_first_return_event(x, false));

    x = {0.6, 2*M_PI - 0.09};
    BOOST_CHECK(curve.is_first_return_event(x, true));
    BOOST_CHECK(!curve.is_first_return_event(x, false));

    x = {1.046, - 2*M_PI - 0.2};
    BOOST_CHECK(!curve.is_first_return_event(x, true));
    BOOST_CHECK(curve.is_first_return_event(x, false));

    x = {1.22, - 2*M_PI + 0.09};
    BOOST_CHECK(!curve.is_first_return_event(x, true));
    BOOST_CHECK(curve.is_first_return_event(x, false));

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(IsoSurfaceFileTest)

BOOST_AUTO_TEST_CASE(assignment_test){

    IsoSurfaceFile file1(modelName);
    string alternate_name = "alternate_name";
    IsoSurfaceFile file2(alternate_name);
    file1 = file2;
    BOOST_REQUIRE(file1.get_modelName() == alternate_name);

}

BOOST_AUTO_TEST_CASE(iterate_test){
    IsoSurfaceFile isoSurfaceFile(modelName);
    // testing iterate through empty File
    int i = 0;
    for (auto surface: isoSurfaceFile) {
        i++;
    }
    BOOST_REQUIRE(i == 0);

    InterpolatedCurve& c1 = isoSurfaceFile.createInterpolatedCurve("Surface0");
    c1.add_parameter("D", 0.15);
    array<double, 2> x = {0.0, 1.0};
    c1.add_node(x);

    InterpolatedCurve& c2 = isoSurfaceFile.createInterpolatedCurve("Surface1");
    c2.add_parameter("sigma", 0.1);
    vector<double> rhos = {0.1, 0.5, 1.0};
    vector<double> phis = {-1.5, 20.5, 1000.0};
    c2.set_nodes(rhos, phis);

    i = 0;
    for (InterpolatedCurve surface: isoSurfaceFile) {
        if (i == 0) {
            // check name
            BOOST_REQUIRE(surface.get_name() == "Surface0");
            // check nodes
            auto [rho, phi] = surface.get_nodes();
            BOOST_REQUIRE(rho.size() == phi.size());
            BOOST_REQUIRE(rho.size() == 1);
            BOOST_CHECK_EQUAL(rho[0], 0.0);
            BOOST_CHECK_EQUAL(phi[0], 1.0);
            // check extensions
            auto [rho_min, rho_max] = surface.get_extensions();
            BOOST_CHECK(rho_min == rho_max);
            BOOST_CHECK_EQUAL(rho_min, 0.0);
            // check parameters
            auto [keys, vals] = surface.get_parameters();
            BOOST_REQUIRE(keys.size() == vals.size());
            BOOST_REQUIRE(keys.size() == 1);
            BOOST_CHECK_EQUAL(keys.front(), "D");
            BOOST_CHECK_EQUAL(vals.front(), 0.15);
        }
        if (i == 1) {
            // check name
            BOOST_REQUIRE(surface.get_name() == "Surface1");
            // check extensions
            auto [rho_min, rho_max] = surface.get_extensions();
            BOOST_CHECK_EQUAL(rho_min, 0.1);
            BOOST_CHECK_EQUAL(rho_max, 1.0);
        }
        i++;
    }
    BOOST_REQUIRE(i == 2);
}

BOOST_AUTO_TEST_CASE(exception_tests){
    fs::path testPath;
    // test if opening a non-existing file causes an HDF5API
    // exception
    bool exceptions_catched = false;
    try{
        fs::path wrongPath = "./wrongfilename.h5";
        IsoSurfaceFile isoSurfaceFile(modelName);
        isoSurfaceFile.read(wrongPath);
    }
    catch(SPhaseFileHDF5APIError error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);

    // test if opening a file with wrong attribute values
    // throws the correct exception
    // wrong class attribute
    exceptions_catched = false;
    testPath = fs::path ("wrong_class_attribute.h5");
        write_headerOnly_IsoSurfaceFile(testPath,
                                        "WRONG_CLASS_ATTRIBUTE",
                                        SPhaseFile_FORMAT_ID,
                                        SPhaseFile_VERSION_ID,
                                        modelName);

    try{
        IsoSurfaceFile isoSurfaceFile(modelName);
        isoSurfaceFile.read(testPath);
    }
    catch(SPhaseFileClassConflict error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);
    fs::remove(testPath); // cleanup - delete the test file

    // wrong format attribute
    exceptions_catched = false;
    testPath = fs::path ("wrong_format_attribute.h5");
        write_headerOnly_IsoSurfaceFile(testPath,
                                        SPhaseFile_IsoSurfaceFile_CLASS_ID,
                                        "WRONG_FORMAT_ATTRIBUTE",
                                        SPhaseFile_VERSION_ID,
                                        modelName);

    try{
        IsoSurfaceFile isoSurfaceFile(modelName);
        isoSurfaceFile.read(testPath);
    }
    catch(SPhaseFileWrongFormat error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);
    fs::remove(testPath); // cleanup - delete the test file

    // wrong version attribute
    exceptions_catched = false;
    testPath = fs::path ("wrong_version_attribute.h5");
        write_headerOnly_IsoSurfaceFile(testPath,
                                        SPhaseFile_IsoSurfaceFile_CLASS_ID,
                                        SPhaseFile_FORMAT_ID,
                                        "WRONG_VERSION_ATTRIBUTE",
                                        modelName);

    try{
        IsoSurfaceFile isoSurfaceFile(modelName);
        isoSurfaceFile.read(testPath);
    }
    catch(SPhaseFileVersionConflict error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);
    fs::remove(testPath); // cleanup - delete the test file

    // wrong model name
    exceptions_catched = false;
    testPath = fs::path ("wrong_model_name.h5");
        write_headerOnly_IsoSurfaceFile(testPath,
                                        SPhaseFile_IsoSurfaceFile_CLASS_ID,
                                        SPhaseFile_FORMAT_ID,
                                        SPhaseFile_VERSION_ID,
                                        "WRONG_MODEL_NAME");

    try{
        IsoSurfaceFile isoSurfaceFile(modelName);
        isoSurfaceFile.read(testPath);
    }
    catch(SPhaseFileModelConflict error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);
    fs::remove(testPath); // cleanup - delete the test file

}

BOOST_AUTO_TEST_CASE(write_attribute_test){

    IsoSurfaceFile isoSurfaceFile(modelName);
    isoSurfaceFile.write(TestPath);
    BOOST_REQUIRE(fs::exists(TestPath));
    // check if the file is in the correct format
    try{
        H5::Exception::dontPrint();
        H5std_string FILENAME(TestPath);
        H5::H5File file(FILENAME, H5F_ACC_RDONLY);

        // define attribute types
        H5::StrType attr_str_type(0, H5T_VARIABLE);
        attr_str_type.setCset(H5T_CSET_UTF8);
        attr_str_type.setStrpad(H5T_STR_NULLTERM);
        H5std_string attr_content;
        H5S_class_t attr_dspace_type;

        // check if class attribute was correctly written
        H5::Attribute attr = file.openAttribute("class");
        H5::DataType attr_dtype = attr.getDataType();
        BOOST_REQUIRE(attr_dtype == attr_str_type);
        attr_dspace_type = attr.getSpace().getSimpleExtentType();
        BOOST_REQUIRE(attr_dspace_type == H5S_SCALAR);
        attr.read(attr_dtype, attr_content);
        BOOST_REQUIRE(attr_content == SPhaseFile_IsoSurfaceFile_CLASS_ID);
        attr_dtype.close();
        attr.close();

        // check if class attribute was correctly written
        attr = file.openAttribute("format");
        attr_dtype = attr.getDataType();
        BOOST_REQUIRE(attr_dtype == attr_str_type);
        attr_dspace_type = attr.getSpace().getSimpleExtentType();
        BOOST_REQUIRE(attr_dspace_type == H5S_SCALAR);
        attr.read(attr_dtype, attr_content);
        BOOST_REQUIRE(attr_content == SPhaseFile_FORMAT_ID);
        attr_dtype.close();
        attr.close();

        // check if class attribute was correctly written
        attr = file.openAttribute("version");
        attr_dtype = attr.getDataType();
        BOOST_REQUIRE(attr_dtype == attr_str_type);
        attr_dspace_type = attr.getSpace().getSimpleExtentType();
        BOOST_REQUIRE(attr_dspace_type == H5S_SCALAR);
        attr.read(attr_dtype, attr_content);
        BOOST_REQUIRE(attr_content == SPhaseFile_VERSION_ID);
        attr_dtype.close();
        attr.close();

        // check if class attribute was correctly written
        attr = file.openAttribute("model");
        attr_dtype = attr.getDataType();
        BOOST_REQUIRE(attr_dtype == attr_str_type);
        attr_dspace_type = attr.getSpace().getSimpleExtentType();
        BOOST_REQUIRE(attr_dspace_type == H5S_SCALAR);
        attr.read(attr_dtype, attr_content);
        BOOST_REQUIRE(attr_content == modelName);
        attr_dtype.close();
        attr.close();

        attr_str_type.close();
        file.close();
    }
        // catch failure caused by the H5File operations
    catch( H5::FileIException error )
    {
        error.printErrorStack();
        BOOST_CHECK(false);
    }
        // catch failure caused by the DataSet operations
    catch( H5::DataSetIException error )
    {
        error.printErrorStack();
        BOOST_CHECK(false);
    }
        // catch failure caused by the DataSpace operations
    catch( H5::DataSpaceIException error )
    {
        error.printErrorStack();
        BOOST_CHECK(false);
    }

    //cleanup and delete file
    fs::remove(TestPath);
}

BOOST_AUTO_TEST_CASE(write_read_file_test){
    //random number generator
    std::mt19937 rn_engine;
    std::normal_distribution<double> n_dist;

    //prepare file with multiple
    IsoSurfaceFile isoSurfaceFile_write(modelName);
    InterpolatedCurve& ic_1 = isoSurfaceFile_write.createInterpolatedCurve("isosurface1");
    ic_1.add_parameter("a", 1.0);
    ic_1.add_parameter("b", 2.0);
    ic_1.set_omegaBar(0.156);
    for (int i = 0; i < 10; i++) {
        array<double, 2> x;
        x[0] = i;
        x[1] = x[0];
        ic_1.add_node(x);
    }

    InterpolatedCurve& ic_2 = isoSurfaceFile_write.createInterpolatedCurve("isosurface2");
    ic_2.add_parameter("c", 3.141592);
    ic_2.add_parameter("delta", 1.4125);
    ic_2.add_parameter("omega", 2.71);
    ic_2.set_omegaBar(-5.456);
    for (int i = 0; i < 200; i++) {
        array<double, 2> x;
        x[0] = n_dist(rn_engine);
        x[1] = n_dist(rn_engine);
        ic_2.add_node(x);
    }
    isoSurfaceFile_write.write(TestPath);

    // read in the file
    IsoSurfaceFile isoSurfaceFile_read(modelName);
    isoSurfaceFile_read.read(TestPath);
    BOOST_REQUIRE(isoSurfaceFile_write == isoSurfaceFile_read);

    //cleanup and delete file
    fs::remove(TestPath);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SerialCorrFileTest)

fs::path configFilePath = "../theConfigFile.json";
fs::path isoSurfaceFilePath = "../../SimData/theisoSurfaceFile.h5";

BOOST_AUTO_TEST_CASE(write_read_file_test){

    //random number generator
    std::mt19937 rn_engine;
    std::normal_distribution<double> n_dist;

    string modelName = "theModelName";
    fs::path testPath("../../../../theTestFile.h5");
    SerialCorrFile serialCorrFile(modelName);
    IsoSurfaceCorr& iS_0 = serialCorrFile.create_isoSurfaceCorr("Isochron0");
    IsoSurfaceCorr& iS_1 = serialCorrFile.create_isoSurfaceCorr("Isovariant0");
    iS_0.N = 10000;
    iS_1.N = 10000;
    iS_0.offset = 100;
    iS_1.offset = 10;
    iS_0.sub_pop_size = 10;
    iS_1.sub_pop_size = 20;
    iS_0.cv = 1.0;
    iS_1.cv = 15.4;
    iS_0.Err_cv = 0.01;
    iS_1.Err_cv = 5.0;
    const unsigned int lags = 10;
    for (unsigned int k = 0; k < lags; k++){
        iS_0.rho_k.push_back(n_dist(rn_engine));
        iS_0.Err_rho_k.push_back(n_dist(rn_engine));
        iS_1.rho_k.push_back(n_dist(rn_engine));
        iS_1.Err_rho_k.push_back(n_dist(rn_engine));
    }

    serialCorrFile.write(testPath);

    SerialCorrFile serialCorrFile_read(modelName);
    serialCorrFile_read.read(testPath);

    BOOST_CHECK(serialCorrFile == serialCorrFile_read);

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(FRTDataFileTest)

fs::path configFilePath = "../theConfigFile.json";
fs::path isoSurfaceFilePath = "../../SimData/theisoSurfaceFile.h5";

BOOST_AUTO_TEST_CASE(exception_tests){
    //random number generator
    std::mt19937 rn_engine;
    std::normal_distribution<double> n_dist;

    fs::path testPath;
    // test if opening a non-existing file causes an HDF5API
    // exception
    bool exceptions_catched = false;
    try{
        fs::path wrongPath = "./wrongfilename.h5";
        FRTDataFile frtDataFile(modelName);
        frtDataFile.read(wrongPath);
    }
    catch(SPhaseFileHDF5APIError error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);

    // test if opening a file with wrong attribute values
    // throws the correct exception
    // wrong class attribute
    exceptions_catched = false;
    testPath = fs::path ("wrong_class_attribute.h5");
        write_headerOnly_DataFile(testPath,
                                  "WRONG_CLASS_ATTRIBUTE",
                                  SPhaseFile_FORMAT_ID,
                                  SPhaseFile_VERSION_ID,
                                  modelName,
                                  configFilePath,
                                  isoSurfaceFilePath);

    try{
        FRTDataFile frtDataFile(modelName);
        frtDataFile.read(testPath);
    }
    catch(SPhaseFileClassConflict error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);
    fs::remove(testPath); // cleanup - delete the test file

    // wrong format attribute
    exceptions_catched = false;
    testPath = fs::path ("wrong_format_attribute.h5");
        write_headerOnly_DataFile(testPath,
                                  SPhaseFile_FRTDataFile_CLASS_ID,
                                  "WRONG_FORMAT_ATTRIBUTE",
                                  SPhaseFile_VERSION_ID,
                                  modelName,
                                  configFilePath,
                                  isoSurfaceFilePath);

    try{
        FRTDataFile frtDataFile(modelName);
        frtDataFile.read(testPath);
    }
    catch(SPhaseFileWrongFormat error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);
    fs::remove(testPath); // cleanup - delete the test file

    // wrong version attribute
    exceptions_catched = false;
    testPath = fs::path ("wrong_version_attribute.h5");
        write_headerOnly_DataFile(testPath,
                                  SPhaseFile_FRTDataFile_CLASS_ID,
                                  SPhaseFile_FORMAT_ID,
                                  "WRONG_VERSION_ATTRIBUTE",
                                  modelName,
                                  configFilePath,
                                  isoSurfaceFilePath);

    try{
        FRTDataFile frtDataFile(modelName);
        frtDataFile.read(testPath);
    }
    catch(SPhaseFileVersionConflict error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);
    fs::remove(testPath); // cleanup - delete the test file

    // wrong model name
    exceptions_catched = false;
    testPath = fs::path ("wrong_model_name.h5");
        write_headerOnly_DataFile(testPath,
                                  SPhaseFile_FRTDataFile_CLASS_ID,
                                  SPhaseFile_FORMAT_ID,
                                  SPhaseFile_VERSION_ID,
                                  "WRONG_MODEL_NAME",
                                  configFilePath,
                                  isoSurfaceFilePath);

    try{
        FRTDataFile frtDataFile(modelName);
        frtDataFile.read(testPath);
    }
    catch(SPhaseFileModelConflict error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);
    fs::remove(testPath); // cleanup - delete the test file

    // corrupt data object test
    FRTDataFile frtDataFile(modelName, configFilePath, isoSurfaceFilePath);
    FRTData& frtData = frtDataFile.createDataSet("IsoSurface1");
    for (int i = 0; i < 10; i++){
        frtData.x0[0].push_back(n_dist(rn_engine));
        frtData.x0[1].push_back(n_dist(rn_engine));
        frtData.mFRT.push_back(n_dist(rn_engine));
        frtData.varFRT.push_back(n_dist(rn_engine));
    }
    // create extra elements -> frtData becomes corrupt
    frtData.x0[0].push_back(n_dist(rn_engine));
    frtData.mFRT.push_back(n_dist(rn_engine));

    exceptions_catched = false;
    testPath = fs::path ("wrong_model_name.h5");
    try{
        frtDataFile.write(testPath);
    }
    catch(SPhaseFileDimError error){
        exceptions_catched = true;
        BOOST_CHECK(true);
    }
    if (!exceptions_catched)
        BOOST_CHECK(false);
    fs::remove(testPath); // cleanup - delete the test file

}

BOOST_AUTO_TEST_CASE(write_read_file_test){
    //random number generator
    std::mt19937 rn_engine;
    std::normal_distribution<double> n_dist;

    //prepare file with multiple
    FRTDataFile frtDataFile_write(modelName, configFilePath, isoSurfaceFilePath);

    FRTData& dat_1 = frtDataFile_write.createDataSet("isosurface1");
    for (int i = 0; i < 10; i++) {
        dat_1.x0[0].push_back(i);
        dat_1.x0[1].push_back(i + 0.2);
        dat_1.mFRT.push_back(i + 2.713);
        dat_1.varFRT.push_back(i + 3.141592);
    }
    FRTData& dat_2 = frtDataFile_write.createDataSet("isosurface2");
    for (int i = 0; i < 200; i++) {
        dat_2.x0[0].push_back(n_dist(rn_engine));
        dat_2.x0[1].push_back(n_dist(rn_engine));
        dat_2.mFRT.push_back(n_dist(rn_engine));
        dat_2.varFRT.push_back(n_dist(rn_engine));
    }
    frtDataFile_write.write(TestPath);

    // read in the file
    FRTDataFile frtDataFile_read(modelName, configFilePath, isoSurfaceFilePath);
    frtDataFile_read.read(TestPath);
    BOOST_REQUIRE(frtDataFile_write == frtDataFile_read);

    //cleanup and delete file
    fs::remove(TestPath);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TbarDataFileTest)

    fs::path configFilePath = "../theConfigFile.json";
    fs::path isoSurfaceFilePath = "../../SimData/theisoSurfaceFile.h5";

    BOOST_AUTO_TEST_CASE(exception_tests){

        fs::path testPath;
        // test if opening a non-existing file causes an HDF5API
        // exception
        bool exceptions_catched = false;
        try{
            fs::path wrongPath = "./wrongfilename.h5";
            TbarDataFile tbarDataFile(modelName);
            tbarDataFile.read(wrongPath);
        }
        catch(SPhaseFileHDF5APIError error){
            exceptions_catched = true;
            BOOST_CHECK(true);
        }
        if (!exceptions_catched)
            BOOST_CHECK(false);

        // test if opening a file with wrong attribute values
        // throws the correct exception
        // wrong class attribute
        exceptions_catched = false;
        testPath = fs::path ("wrong_class_attribute.h5");
        write_headerOnly_DataFile(testPath,
                                  "WRONG_CLASS_ATTRIBUTE",
                                  SPhaseFile_FORMAT_ID,
                                  SPhaseFile_VERSION_ID,
                                  modelName,
                                  configFilePath,
                                  isoSurfaceFilePath);

        try{
            TbarDataFile tbarDataFile(modelName);
            tbarDataFile.read(testPath);
        }
        catch(SPhaseFileClassConflict error){
            exceptions_catched = true;
            BOOST_CHECK(true);
        }
        if (!exceptions_catched)
            BOOST_CHECK(false);
        fs::remove(testPath); // cleanup - delete the test file

        // wrong format attribute
        exceptions_catched = false;
        testPath = fs::path ("wrong_format_attribute.h5");
        write_headerOnly_DataFile(testPath,
                                  SPhaseFile_TbarDataFile_CLASS_ID,
                                  "WRONG_FORMAT_ATTRIBUTE",
                                  SPhaseFile_VERSION_ID,
                                  modelName,
                                  configFilePath,
                                  isoSurfaceFilePath);

        try{
            TbarDataFile tbarDataFile(modelName);
            tbarDataFile.read(testPath);
        }
        catch(SPhaseFileWrongFormat error){
            exceptions_catched = true;
            BOOST_CHECK(true);
        }
        if (!exceptions_catched)
            BOOST_CHECK(false);
        fs::remove(testPath); // cleanup - delete the test file

        // wrong version attribute
        exceptions_catched = false;
        testPath = fs::path ("wrong_version_attribute.h5");
        write_headerOnly_DataFile(testPath,
                                  SPhaseFile_TbarDataFile_CLASS_ID,
                                  SPhaseFile_FORMAT_ID,
                                  "WRONG_VERSION_ATTRIBUTE",
                                  modelName,
                                  configFilePath,
                                  isoSurfaceFilePath);

        try{
            TbarDataFile tbarDataFile(modelName);
            tbarDataFile.read(testPath);
        }
        catch(SPhaseFileVersionConflict error){
            exceptions_catched = true;
            BOOST_CHECK(true);
        }
        if (!exceptions_catched)
            BOOST_CHECK(false);
        fs::remove(testPath); // cleanup - delete the test file

        // wrong model name
        exceptions_catched = false;
        testPath = fs::path ("wrong_model_name.h5");
        write_headerOnly_DataFile(testPath,
                                  SPhaseFile_TbarDataFile_CLASS_ID,
                                  SPhaseFile_FORMAT_ID,
                                  SPhaseFile_VERSION_ID,
                                  "WRONG_MODEL_NAME",
                                  configFilePath,
                                  isoSurfaceFilePath);

        try{
            TbarDataFile tbarDataFile(modelName);
            tbarDataFile.read(testPath);
        }
        catch(SPhaseFileModelConflict error){
            exceptions_catched = true;
            BOOST_CHECK(true);
        }
        if (!exceptions_catched)
            BOOST_CHECK(false);
        fs::remove(testPath); // cleanup - delete the test file

    }

    BOOST_AUTO_TEST_CASE(write_read_file_test){
        //random number generator
        std::mt19937 rn_engine;
        std::normal_distribution<double> n_dist;

        //prepare file with multiple
        TbarDataFile tbarDataFile_write(modelName, configFilePath, isoSurfaceFilePath);

        TbarData& dat_1 = tbarDataFile_write.createDataSet("isosurface1");
        dat_1.Tbar = 2.71;
        TbarData& dat_2 = tbarDataFile_write.createDataSet("isosurface2");
        dat_2.Tbar = n_dist(rn_engine);
        tbarDataFile_write.write(TestPath);

        // read in the file
        TbarDataFile tbarDataFile_read(modelName, configFilePath, isoSurfaceFilePath);
        tbarDataFile_read.read(TestPath);
        BOOST_REQUIRE(tbarDataFile_write == tbarDataFile_read);

        //cleanup and delete file
        fs::remove(TestPath);
    }

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SimConfigFile_test)

BOOST_AUTO_TEST_CASE(basic_functionality) {
    Config config;
    config.modelName = modelName;
    config.paths.input = fs::weakly_canonical(
            fs::path("../../Test/input.h5")
            );
    config.paths.output = fs::weakly_canonical(
            fs::path("../../Test/output.h5")
            );
    config.simulation.dt = 0.001;
    config.simulation.t0 = 10.0;
    config.simulation.T = 1000.0;
    config.simulation.EnsembleSize = 1000000;
    config.simulation.SampleSize = 20;

    SimConfigFile simConfigFile(config);
    BOOST_CHECK(simConfigFile.get_modelName() == modelName);
    Config::Simulation simConfig = simConfigFile.get_simConfig();
    BOOST_CHECK_EQUAL(simConfig.dt, 0.001);
    BOOST_CHECK_EQUAL(simConfig.t0, 10.0);
    BOOST_CHECK_EQUAL(simConfig.T, 1000.0);
    BOOST_CHECK_EQUAL(simConfig.EnsembleSize, 1000000);
    BOOST_CHECK_EQUAL(simConfig.SampleSize, 20);

    BOOST_CHECK(simConfigFile.get_inPath() ==
        fs::weakly_canonical(
                fs::path("../../Test/input.h5")
                ));
    BOOST_CHECK(simConfigFile.get_outPath() ==
        fs::weakly_canonical(
                fs::path("../../Test/output.h5")
                ));
}

BOOST_AUTO_TEST_CASE(write_read_json) {
    //fs::path testPath("./Test.json");
    fs::path testPath("Test.json");
    fs::path testestpath = fs::path("");
    if (testestpath == fs::path(""))
        testestpath = fs::path(".");
    testPath = testestpath / testPath;

    // test with relative paths
    Config config_write;
    config_write.modelName = modelName;
    ParameterSet pSet;
    pSet["D"] = 0.1;
    pSet["omega"] = 1.0;
    pSet["gamma"] = 15.0;
    pSet["c"] = -15.0;
    config_write.pSetList.push_back(pSet);
    pSet["D"] = 0.5;
    config_write.pSetList.push_back(pSet);
    config_write.domainName = string("ReflectiveAnnulus");
    ParameterSet dimSet;
    dimSet["rho_min"] = 0.5;
    dimSet["rho_max"] = 1.5;
    config_write.domainDimList.push_back(dimSet);
    dimSet["rho_min"] = 0.3;
    dimSet["rho_max"] = 4.0;
    config_write.domainDimList.push_back(dimSet);
    config_write.paths.input =
            fs::absolute(testPath).parent_path() / fs::path("../../Test/input.h5");
    config_write.paths.input = fs::weakly_canonical(config_write.paths.input);
    config_write.paths.output =
            fs::absolute(testPath).parent_path() / fs::path("../../Test/output.h5");
    config_write.paths.output = fs::weakly_canonical(config_write.paths.output);
    config_write.simulation.dt = 0.002;
    config_write.simulation.t0 = 20.0;
    config_write.simulation.T = 2000.0;
    config_write.simulation.EnsembleSize = 2000000;
    config_write.simulation.SampleSize = 40;

    SimConfigFile simConfigFile_write(config_write);
    simConfigFile_write.write(testPath);

    SimConfigFile simConfigFile_read;
    simConfigFile_read.read(testPath);

    BOOST_REQUIRE(simConfigFile_write == simConfigFile_read);
    BOOST_REQUIRE(simConfigFile_write.get_modelName() == simConfigFile_read.get_modelName());
    BOOST_REQUIRE(simConfigFile_write.get_simConfig() == simConfigFile_read.get_simConfig());
    BOOST_REQUIRE(simConfigFile_write.get_inPath() == simConfigFile_read.get_inPath());
    BOOST_REQUIRE(simConfigFile_write.get_outPath() == simConfigFile_read.get_outPath());

    // test with absolute paths
    fs::path absolute_path = fs::current_path();
    config_write.paths.input = absolute_path / fs::path("input.h5");
    config_write.paths.output = absolute_path / fs::path("output.h5");

    simConfigFile_write = SimConfigFile(config_write);
    simConfigFile_write.write(testPath, false);

    simConfigFile_read.read(testPath);

    BOOST_REQUIRE(simConfigFile_write.get_inPath() == simConfigFile_read.get_inPath());
    BOOST_REQUIRE(simConfigFile_write.get_outPath() == simConfigFile_read.get_outPath());

    fs::remove(testPath);

}

BOOST_AUTO_TEST_SUITE_END()

void write_headerOnly_IsoSurfaceFile(fs::path p,
                                     string class_attribute,
                                     string format_attribute,
                                     string version_attribute,
                                     string modelName){
    try{
        H5::Exception::dontPrint();

        H5std_string FILENAME(p);
        H5::H5File file(FILENAME, H5F_ACC_TRUNC);

        // define attribute types
        H5::StrType attr_str_type(0, H5T_VARIABLE);
        attr_str_type.setCset(H5T_CSET_UTF8);
        attr_str_type.setStrpad(H5T_STR_NULLTERM);
        H5::DataSpace attr_ds(H5S_SCALAR);

        // write file attributes
        H5::Attribute attr = file.createAttribute("class", attr_str_type, attr_ds);
        attr.write(attr_str_type, class_attribute);
        attr.close();
        attr = file.createAttribute("format", attr_str_type, attr_ds);
        attr.write(attr_str_type, format_attribute);
        attr.close();
        attr = file.createAttribute("model", attr_str_type, attr_ds);
        attr.write(attr_str_type, modelName);
        attr.close();
        attr = file.createAttribute("version", attr_str_type, attr_ds);
        attr.write(attr_str_type, version_attribute);
        attr.close();

        attr_ds.close();
        attr_str_type.close();

        file.close();
    }
        // catch failure caused by the H5File operations
    catch( H5::FileIException error )
    {
        error.printErrorStack();
        BOOST_CHECK(false);
    }
        // catch failure caused by the DataSet operations
    catch( H5::DataSetIException error )
    {
        error.printErrorStack();
        BOOST_CHECK(false);
    }
        // catch failure caused by the DataSpace operations
    catch( H5::DataSpaceIException error )
    {
        error.printErrorStack();
        BOOST_CHECK(false);
    }
}

void write_headerOnly_DataFile(fs::path p,
                               string class_attribute,
                               string format_attribute,
                               string version_attribute,
                               string modelName,
                               fs::path configFilePath,
                               fs::path isoSurfaceFilePath){
    try{
        H5::Exception::dontPrint();

        H5std_string FILENAME(p);
        H5::H5File file(FILENAME, H5F_ACC_TRUNC);

        // define attribute types
        H5::StrType str_type(0, H5T_VARIABLE);
        str_type.setCset(H5T_CSET_UTF8);
        str_type.setStrpad(H5T_STR_NULLTERM);
        H5::DataSpace dspace(H5S_SCALAR);

        // write file attributes
        H5::Attribute attr = file.createAttribute("class", str_type, dspace);
        attr.write(str_type, class_attribute);
        attr.close();
        attr = file.createAttribute("format", str_type, dspace);
        attr.write(str_type, format_attribute);
        attr.close();
        attr = file.createAttribute("model", str_type, dspace);
        attr.write(str_type, modelName);
        attr.close();
        attr = file.createAttribute("version", str_type, dspace);
        attr.write(str_type, version_attribute);
        attr.close();
        // write file paths
        H5::DataSet dset = file.createDataSet("configuration_file", str_type, dspace);
        dset.write(configFilePath, str_type);
        dset.close();
        dset = file.createDataSet("isosurface_file", str_type, dspace);
        dset.write(isoSurfaceFilePath, str_type);
        dset.close();

        dspace.close();
        str_type.close();

        file.close();
    }
        // catch failure caused by the H5File operations
    catch( H5::FileIException error )
    {
        error.printErrorStack();
        BOOST_CHECK(false);
    }
        // catch failure caused by the DataSet operations
    catch( H5::DataSetIException error )
    {
        error.printErrorStack();
        BOOST_CHECK(false);
    }
        // catch failure caused by the DataSpace operations
    catch( H5::DataSpaceIException error )
    {
        error.printErrorStack();
        BOOST_CHECK(false);
    }
}