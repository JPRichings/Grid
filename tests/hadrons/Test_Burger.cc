// Quick and dirty test of new burger implimentations

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;


    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    std::vector<std::string> flavour = {"h"}; //{"l", "s", "c1", "c2", "c3"};
    std::vector<double>      mass    = {.2}; //{.01, .04, .2  , .25 , .3  };
    std::vector<std::string> lepton_flavour    = {"mu"};
    std::vector<double>      lepton_mass    = {.2};

    unsigned int  nt    = GridDefaultLatt()[Tp];

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
    globalPar.trajCounter.step  = 20;
    globalPar.runId             = "test";
    application.setPar(globalPar);


    // gauge field
    application.createModule<MGauge::Unit>("gauge");

    // sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);

    // sink EDIT FOR TEST
    MSink::ExactPhoton::Par sinkPar;
    sinkPar.mom = "0 0 0";
    sinkPar.position = "0 0 0 0";
    sinkPar.gauge = PhotonR::Gauge::feynman;
    sinkPar.zmScheme = PhotonR::ZmScheme::qedL;
    application.createModule<MSink::ScalarExactPhoton>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    int i = 0;
    int j = 0;

    // actions
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls    = 12;
    actionPar.M5    = 1.8;
    actionPar.mass  = mass[i];
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    application.createModule<MAction::DWF>("DWF_" + flavour[i], actionPar);
    
    // solvers
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = "DWF_" + flavour[i];
    solverPar.residual     = 1.0e-8;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                solverPar);
    
    // propagators
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = "CG_" + flavour[i];
    quarkPar.source = "pt";
    application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i], quarkPar);

    MContraction::Meson::Par mesPar;
    mesPar.output  = "mesons/pt_" + flavour[i] + flavour[j];
    mesPar.q1      = "Qpt_" + flavour[i];
    mesPar.q2      = "Qpt_" + flavour[j];
    mesPar.gammas  = "GammaT GammaT"; 
    mesPar.sink    = "sink";
    application.createModule<MContraction::Meson>("burger_pt_"
                                                    + flavour[i] + flavour[j],
                                                    mesPar);
    // execution
    application.saveParameterFile("burger.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;

}