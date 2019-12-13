#ifndef Hadrons_MSink_ExactPhoton_hpp_
#define Hadrons_MSink_ExactPhoton_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         ExactPhoton                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class ExactPhotonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExactPhotonPar,
                                    std::string, mom,
                                    std::string, position,
                                    PhotonR::Gauge, gauge,
                                    PhotonR::ZmScheme, zmScheme);
};

template <typename FImpl>
class TExactPhoton: public Module<ExactPhotonPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
public:
    typedef PhotonR::GaugeField         EmField;
    typedef PhotonR::GaugeLinkField     EmLinkField;
public:
    // constructor
    TExactPhoton(const std::string name);
    // destructor
    virtual ~TExactPhoton(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasPhase_{false}; 
    std::string momphName_;
};

MODULE_REGISTER_TMP(ExactPhoton, TExactPhoton<FIMPL>, MSink);
MODULE_REGISTER_TMP(ScalarExactPhoton, TExactPhoton<ScalarImplCR>, MSink);

/******************************************************************************
 *                 TExactPhoton implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TExactPhoton<FImpl>::TExactPhoton(const std::string name)
: Module<ExactPhotonPar>(name)
, momphName_ (name + "_momph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TExactPhoton<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TExactPhoton<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExactPhoton<FImpl>::setup(void)
{
    envTmpLat(LatticeComplex, "coor");
    envCacheLat(LatticeComplex, momphName_);
    envCacheLat(LatticeComplex, "out");
    envCreate(SinkFn, getName(), 1, nullptr);
    envCreateLat(EmField, "a");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExactPhoton<FImpl>::execute(void)
{
    std::vector<Real> position = strToVec<Real>(par().position);

    LOG(Message) << "Setting up Exact photon sink function for momentum ["
                 << par().mom << "]" << std::endl;

    auto &ph = envGet(LatticeComplex, momphName_);
    
    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = zero;
        for(unsigned int mu = 0; mu < p.size(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        hasPhase_ = true;
    }

    // set up enviroment
    PhotonR photon(envGetGrid(EmField), par().gauge, par().zmScheme);
    auto    &a = envGet(EmField, "a");
    photon.UnitField(a);
    //auto    &out = envGet(EmField, "out");

    // get exp(ik.x0) phase

    auto                        g      = env().getGrid();
    LatticeComplex              k(g);  
    const unsigned int          nd     = g->Nd();
    std::vector<int>            l      = g->FullDimensions();
    Complex                     i(0.0,1.0);
    LatticeComplex              res(g);
    auto &out = envGet(LatticeComplex, "out");   //LatticeComplex              out(g);

    // calculate phase exp(ik.x0)
    //k.resize(nd, g);
    res = Zero();

    for (unsigned int mu = 0; mu < nd; mu++)
    {
      LatticeCoordinate(k, mu);
      res = res + (position[mu]/env().getDim(mu))*k;
    }
    res  = exp((Real)(M_PI)*i*res);
    
    // poke phase into
    for(unsigned int mu = 0; mu < nd; mu++)
    {
      pokeLorentz(a, res, mu);
    }
    // multiply phase and khat sqr inverse
    photon.MomentumSpacePropagator(a,a);

    out = peekLorentz(a,0);

    std::cout << "out: " << std::endl;//<< out 

    std::cout << "a: " << std::endl;//a << 

    // get fft
    FFT  fft(g);
    std::cout << "check A" << std::endl;
    fft.FFT_all_dim(out, out, FFT::backward);
    std::cout << "check B" << std::endl;

    auto sink = [&ph, &out](const PropagatorField &field)
    {
        SlicedPropagator res;
        std::cout << "check 1" << std::endl;
        PropagatorField  tmp = ph*out*field;
        std::cout << "check 2" << std::endl;
        // do spacial y sum
        sliceSum(tmp, res, Tp);
        
        return res;
    };

    envGet(SinkFn, getName()) = sink;   
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_ExactPhoton_hpp_
