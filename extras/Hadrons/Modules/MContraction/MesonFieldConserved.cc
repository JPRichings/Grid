#include <Grid/Hadrons/Modules/MContraction/MesonFieldConserved.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MContraction;

/******************************************************************************
*                  TMesonFieldConserved implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TMesonFieldConserved::TMesonFieldConserved(const std::string name)
: Module<MesonFieldConservedPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TMesonFieldConserved::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TMesonFieldConserved::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TMesonFieldConserved::setup(void)
{

}

// execution ///////////////////////////////////////////////////////////////////
void TMesonFieldConserved::execute(void)
{

}
