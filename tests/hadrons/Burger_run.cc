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
    std::vector<std::string> flavour = {"l","s"}; //{"l", "s", "c1", "c2", "c3"};
    std::vector<double>      mass    = {0.005,0.04}; //{.01, .04, .2  , .25 , .3  };

    unsigned int  nt    = GridDefaultLatt()[Tp];
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";
    int hits = 600 + 12*64;

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 3000;
    globalPar.trajCounter.end   = 3040;
    globalPar.trajCounter.step  = 40;
    globalPar.runId             = "mesonfield-C0";
    application.setPar(globalPar);


    // gauge field
    MIO::LoadNersc::Par confPar;
    confPar.filename = "/home/dp008/dp008/dc-rich6/Scripts/ConfigsDeflQED/ckpoint_lat";
    application.createModule<MIO::LoadNersc>("gauge",confPar);

    //single precision gauge feild
    MUtilities::GaugeSinglePrecisionCast::Par singlePar;
    singlePar.field = "gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("gaugef",singlePar);


    for(unsigned int i = 0; i< flavour.size(); i++)
    {
    // actions 
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls    = 12;
    actionPar.M5    = 1.8;
    actionPar.mass  = mass[i];
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    application.createModule<MAction::DWF>("DWF_" + flavour[i], actionPar);

    // single precision action
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gaugef";
    actionPar.Ls    = 12;
    actionPar.M5    = 1.8;
    actionPar.mass  = mass[i];
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    application.createModule<MAction::DWF>("DWFF_" + flavour[i], actionPar);

    // solvers 
    //MSolver::RBPrecCG::Par solverPar;
    //solverPar.action       = "DWF_" + flavour[i];
    //solverPar.residual     = 1.0e-8;
    //solverPar.maxIteration = 10000;
    // application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
    //                                        solverPar);

    // multi precision solver
    MSolver::MiedPrecisionRBPrecCG::Par MPsolverPar;
    MPsolverPar.inneraction = "DWFF_" + flavour[i];
    MPsolverPar.outeraction = "DWF_" + flavour[i];
    MPsolverPar.residual     = 1.0e-8;
    MPsolverPar.maxinnerIteration = 30000;
    MPsolverPar.maxOutererIteration = 100;
    MPsolverPar.eigenpack = "";
    application.createModule<MSolver::MiedPrecisionRBPrecCG>("CG_" + flavour[i],
                                            solverPar);
    }

    // from hit number work out source position
    
    unsigned int spatial_hits = 2;

    std::vector<std::string> position;
    std::vector<std::string> position_name;
    for(unsigned int tt = 0; tt < Tp; tt++)
    {
        for(unsigned int zz = 0; zz < spatial_hits; zz++)
        {
            for(unsigned int yy = 0; yy < spatial_hits; yy++)
            {
                for(unsigned int xx = 0; xx < spatial_hits; xx++)
                {
                position.push_back(std::to_string(8*xx + 7)
                 + " " + std::to_string(8*yy + 7)
                 + " " + std::to_string(8*zz + 7)
                 + " " + std::to_string(tt));
                position_name.push_back(std::to_string(8*xx + 7)
                 + "_" + std::to_string(8*yy + 7)
                 + "_" + std::to_string(8*zz + 7)
                 + "_" + std::to_string(tt));
                }
            }
        }
    }


    for(unsigned int i = 0; i < flavour.size(); i++)
    {
        for(unsigned int h = 0; h<position.size(); h++)
        {
 

        // sources
        MSource::Point::Par ptPar;
        ptPar.position = position[h];
        application.createModule<MSource::Point>("pt_" + position_name[h], ptPar);

        // sink
        MSink::ExactPhoton::Par sinkPar;
        sinkPar.mom = "0 0 0";
        sinkPar.position = position[h];
        sinkPar.gauge = PhotonR::Gauge::feynman;
        sinkPar.zmScheme = PhotonR::ZmScheme::qedL;
        application.createModule<MSink::ScalarExactPhoton>("sink_" + position_name[h], sinkPar);
    
        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt_" + position_name[h];
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i] + position_name[h], quarkPar);

        MContraction::Meson::Par mesPar;
        mesPar.output  = "mesons/pt_" + flavour[i] + flavour[i] + "_GammaT_GammaT_" + position_name[h];
        mesPar.q1      = "Qpt_" + flavour[i] + position_name[h];
        mesPar.q2      = "Qpt_" + flavour[i] + position_name[h];
        mesPar.gammas  = "GammaT GammaT"; 
        mesPar.sink    = "sink_" + position_name[h];
        application.createModule<MContraction::Meson>("burger_pt_" + flavour[i]
                                                                   + flavour[i]
                                                                   + "_GammaT_GammaT_"
                                                                   + position_name[h],
                                                        mesPar);
        MContraction::Meson::Par mesPar;
        mesPar.output  = "mesons/pt_" + flavour[i] + flavour[i] + "_GammaX_GammaX_" + position_name[h];
        mesPar.q1      = "Qpt_" + flavour[i] + position_name[h];
        mesPar.q2      = "Qpt_" + flavour[i] + position_name[h];
        mesPar.gammas  = "GammaX GammaX"; 
        mesPar.sink    = "sink_" + position_name[h];
        application.createModule<MContraction::Meson>("burger_pt_" + flavour[i]
                                                                   + flavour[i]
                                                                   + "_GammaX_GammaX_"
                                                                   + position_name[h],
                                                        mesPar);
        MContraction::Meson::Par mesPar;
        mesPar.output  = "mesons/pt_" + flavour[i] + flavour[i] + "_GammaY_GammaY_" + position_name[h];
        mesPar.q1      = "Qpt_" + flavour[i] + position_name[h];
        mesPar.q2      = "Qpt_" + flavour[i] + position_name[h];
        mesPar.gammas  = "GammaY GammaY"; 
        mesPar.sink    = "sink_" + position_name[h];
        application.createModule<MContraction::Meson>("burger_pt_" + flavour[i]
                                                                   + flavour[i]
                                                                   + "_GammaY_GammaY_"
                                                                   + position_name[h],
                                                        mesPar);
        MContraction::Meson::Par mesPar;
        mesPar.output  = "mesons/pt_" + flavour[i] + flavour[i] + "_GammaZ_GammaZ_" + position_name[h];
        mesPar.q1      = "Qpt_" + flavour[i] + position_name[h];
        mesPar.q2      = "Qpt_" + flavour[i] + position_name[h];
        mesPar.gammas  = "GammaZ GammaZ"; 
        mesPar.sink    = "sink_" + position_name[h];
        application.createModule<MContraction::Meson>("burger_pt_" + flavour[i]
                                                                   + flavour[i]
                                                                   + "_GammaZ_GammaZ_"
                                                                   + position_name[h],
                                                        mesPar);                                                        
        }
    }

    // execution
    application.saveParameterFile("burger.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;

}
