/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/WilsonFermion.cc

Copyright (C) 2015

Author: James Richings <j.p.richings@soton.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/


#ifndef Hadrons_MContraction_MesonFieldConserved_hpp_
#define Hadrons_MContraction_MesonFieldConserved_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential A2A vector with insertion of conserved current. 
 Additionally optional insertion of a photon field A_\mu(x).
 -----------------------------
 * src_x = sum_{mu=mu_min}^{mu_max} 
     q_x * theta(x_3 - tA) * theta(tB - x_3) * J_mu * exp(i x.mom) (* A_\mu(x))
 
 * options:
 - q: input propagator (string)
 - action: fermion action used for propagator q (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - curr_type: type of conserved current to insert (Current)
 - mu_min: begin Lorentz Index (integer)
 - mu_max: end Lorentz Index (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 - photon: optional photon field (string)
 
 */

/******************************************************************************
 *                         MesonFieldConserved                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class MesonFieldConservedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFieldConservedPar,
                                    int, Nl,
                                    int, N,
                                    std::string, A2A1,
                                    std::string, A2A2,
                                    std::string, gammas,
                                    std::string, output,
                                    std::string,  action,
                                    unsigned int tA,
                                    unsigned int tB,
                                    Current,      curr_type,
                                    unsigned int, mu_min,
                                    unsigned int, mu_max,
                                    std::string,  mom,
                                    std::string,  photon);
};

template <typename FImpl>
class TMesonFieldConserved: public Module<MesonFieldConservedPar>
{
    public:
        FERM_TYPE_ALIASES(FImpl, );

        typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat> A2ABase;
        typedef A2AVectorsReturn<typename FImpl::FermionField, FMat> A2AReturn;

        class Result : Serializable
        {
            public:
                GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                                Gamma::Algebra, gamma_snk,
                                                Gamma::Algebra, gamma_src,
                                                std::vector<Complex>, corr);
        };


public:
    typedef PhotonR::GaugeField     EmField;
public:
    // constructor
    TMesonFieldConserved(const std::string name);
    // destructor
    virtual ~TMesonFieldConserved(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void parseGammaString(atd::vector<GammaPair> & gammaList);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        SeqhasPhase_{false};
    std::string SeqmomphName_;
};

MODULE_REGISTER(MesonFieldConserved, ARG(TMesonFieldConserved<FIMPL>), MContraction);
MODULE_REGISTER(ZMesonFieldConserved, ARG(TMesonFieldConserved<ZFIMPL>), MContraction);
/******************************************************************************
 *                  MesonFieldConserved  implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMesonFieldConserved<FIpml>::TMesonFieldConserved(const std::string name)
: Module<MesonFieldConservedPar>(name)
, SeqmomphName_ (name + "_Seqmom")
{}
// dependencies/products ///////////////////////////////////////////////////////
template<type FImpl>
std::vector<std::string> TMesonFieldConserved<FImpl>::getInput(void)
{
    std::vecotr<std<std::string> in = {par().action, par().A2A1, par().A2A2};
    if (!par().photon.empty()) in.push_back(par().photon);
    in.push_back(par().A2A1 + "_ret");
    in.push_back(par().A2A2 + "_ret");

    return in;
}

template <typename FImpl>
std::vector<std::string> TMesonFieldConserved<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

template <typename FImpl>
void TMesonFieldConserved<FImpl>::parse GammaString(std::vector<GammaPair> & gammaList)
{
    gammaList.clear()
    //Parse individual contractions from input string.
    gammaList = strToVec<GammaPair>(Par().gammas);
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldConserved<FImpl>::setup(void)
{   
    // for seq conserved
    auto Ls_ = env().getObjectLs(par().action);
    envTmpLat(FermionField, getName(), Ls_);
    envTmpLat(FermionField, "src_tmp_ferm");
    envTmpLat(PropagatorField, "src_tmp");
    envCacheLat(LatticeComplex, SeqmomphName_);
    envTmpLat(LatticeComplex, "coor");
    envTmpLat(LatticeComplex, "latt_compl");
    envTempLat(PropagatorField, "q"); // delete once template

    // for A2A

    int nt = env().getDim(Tp);
    int N = par().N;
    envTmp(std::vector<FermionField>, "w1", 1, N, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "v1", 1, N, FermionField(env().getGrid(1)));
    //envTmp(std::vector<FermionField>, "v1_con", Ls_, N, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "w1_5d", Ls_, N, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "v1_5d", Ls_, N, FermionField(env().getGrid(1)));
    //envTmpLat(FermionField, "tmpv_5d", Ls_);
    //envTmpLat(FermionField, "tmpw_5d", Ls_);

    envTmp(std::vector<ComplexD>, "MF_x", 1, nt);
    envTmp(std::vector<ComplexD>, "MF_y", 1, nt);
    envTmp(std::vector<ComplexD>, "MF_z1", 1, nt);
    envTmp(std::vector<ComplexD>, "MF_z2", 1, nt);
    envTmp(std::vector<ComplexD>, "tmp", 1, nt);
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldConserved<FImpl>::execute(void)
{
    // Log message to explain process that is about to start.
    LOG(Message) << "Computing A2A Sequential conserved vector" << std::endl;

    // Define results structure
    Result result;

    // Gamma matrices
    Gamma g5(Gamma::Algebra::Gamma5);
    std::vector<GammaPair> gammaList;
    
    parseGammaString(gammaList);

    // Get number of timeslices
    nt = env().getDim(Tp);

    // Fill Results with data avalible

    result.gamma_snk = gammaList[0].first
    result.gamma_src = gammaList[0].second
    result.corr.resize(nt);
    //result.MFx.resize(nt);
    //result.MFz1.resize(nt);
    //result.MFy.resize(nt);

    int Nl = par().Nl;
    int N = par().N;

    LOG(Message) << "N for A2A cont: " << N << std::endl;

    // get temporary volumes
    envGetTmp(std::vector<ComplexD>, MF_x);
    envGetTmp(std::vector<ComplexD>, MF_y);
    envGetTmp(std::vector<ComplexD>, MF_z1);
    envGetTmp(std::vector<ComplexD>, MF_z2);
    envGetTmp(std::vector<ComplexD>, tmp);

    // Check the the initialization of tmp.
    for (unsigned int t = 0; t < nt; ++t)
    {
        tmp[t] = TensorRemove(MF_x[t] * MF_z1[t] * MF_x[t] * MF_z1[t] * 0.0);
    }

    Gamma gSnk(gammaList[0].first);
    Gamma gSrc(gammaList[0].second);


    // Get a pointer to the derived tyoe that returns the all to all vectors
    auto &a2a1_fn = envGet(A2ABase, par().A2A1 + "_ret");
    auto &a2a2_fn = envGet(A2ABase, par().A2A1 + "_ret");
    
    // Get temp fermion fields 
    envGetTmp(std::vector<FermionField>, w1);
    envGetTmp(std::vector<FermionField>, v1);
    //envGetTmp(std::vector<FermionField>, v1_con);
    envGetTmp(std::vector<FermionField>, w1_5d);
    envGetTmp(std::vector<FermionField>, v1_5d);

    // Get v and w vectors
    LOG(Message) << "Finding v and w vectors for N =  " << N << std::endl;
    for (int i = 0; i < N; i++)
    {
        a2a1_fn.return_v(i, v1_5d[i], v1[i]);
        a2a1_fn.return_w(i, w1_5d[i], w1[i]); // remove this to lower down
    }
    LOG(Message) << "Found v and w vectors for N =  " << N << std::endl;


    // sequential converved //

   /* if (par().tA == par().tB)
    {
        LOG(Message) << "Generating A2A sequential V_i with conserved "
                     << par().curr_type << " current at " 
		             << "t = " << par().tA << " summed over the indices " 
		             << par().mu_min << " <= mu <= " << par().mu_max 
		             << std::endl;
    }
    else
    {
        LOG(Message) << "Generating A2A sequential V_i with conserved "
                     << par().curr_type << " current for " 
                     << par().tA << " <= t <= " 
                     << par().tB << " summed over the indices " 
		             << par().mu_min << " <= mu <= " << par().mu_max
	                 << std::endl;
    }*/

    envGetTmp(FermionField, src);
    envGetTmp(PropagatorField, src_tmp);
    envGetTmp(FermionField, src_tmp_ferm);
    src_tmp_ferm = src; // this is here for vera found that grid wouldn't compile without it.
    envGetTmp(PropagatorField, q);
    auto &mat = envGet(FMat, par().action);
    envGetTmp(LatticeComplex, latt_compl);

    //exp(ipx)
    auto &mom_phase = envGet(LatticeComplex, SeqmomphName_);
    if (!SeqhasPhase_)
    {    
        std::vector<Real> mom = strToVec<Real>(par().mom);
        mom_phase = zero;
        Complex           i(0.0,1.0);
        envGetTmp(LatticeComplex, coor);
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
            mom_phase = mom_phase + (mom[mu]/env().getGrid()->_fdimensions[mu])*coor;
        }
        mom_phase = exp((Real)(2*M_PI)*i*mom_phase);
        SeqhasPhase_ = true;
    }
    LOG(Message) << "Inserting momentum " << strToVec<Real>(par().mom) << std::endl;

    if (!par().photon.empty())    	
    {
	 LOG(Message) << "Inserting the stochastic photon field " << par().photon << std::endl;
    }


    for(int i = 0; i < N; i++)
    {
        src = zero;

        for(unsigned int mu=par().mu_min;mu<=par().mu_max;mu++)
        {
            if (!par().photon.empty())    	
            {
            //Get the stochastic photon field, if required
                auto &stoch_photon = envGet(EmField,  par().photon);
                latt_compl =  PeekIndex<LorentzIndex>(stoch_photon, mu) * mom_phase;
            }
            else
            {
                latt_compl = mom_phase;
            }
            // old code
            // mat.SeqConservedCurrent(q, src_tmp, par().curr_type, mu, //seqConservedCurrent is called from wilsonFermion5D        
            //                     par().tA, par().tB, latt_compl);

            // new code
            // Convert the fermion field of the all to all vector into a propagator
            FermToProp(v1_5d[i], q, 0, 0); // need to check that v1 is called correctly. probably need to loop over number of modes
            // Run the Multiplication
            mat.SeqConservedCurrent(q, src_tmp, par().current, mu, par().tA, par().tB, latt_compl);
            // Convert back into a fermion
            PropToFerm(src_tmp_ferm, src_tmp, 0, 0);
        src += src_tmp_ferm;
        }
        v1_5d[i] = src;
    }

    // multiply be gamma matrices as required

    for (unsigned int i = 0; i < N; i++)
    {
        v1[i] = gSnk * v1[i];
    }
    // Define ty
    int ty;

    // start contraction loop
    for(int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            for (unsigned int k = 0; k < N; k++)
            {
                for (unsigned int l = 0; l < N; l++)
                {
                    sliceInnerProductVector(MF_x, adj(w1[l]), v1[i], Tp);
                    sliceInnerProductVector(MF_y, adj(w1[j]), v1[k], Tp);
                    sliceInnerProductVector(MF_z1_5d, adj(w1_5d[i]), v1_5d[j], Tp); // MF_z1_5d define
                    sliceInnerProductVector(MF_z2_5d, adj(W1_5d[k]), v1_5d[l], Tp); // MF_z2_5d define
                    //sum over 5th dim
                    MF_z1 = zero;
                    MF_z2 = zero;
                    for (unsigned int s = 0, s < Ls_, s++)
                    {
                        ExtractSlice(tmp_4d, MF_z1_5d,s,0);
                        MF_z1+=tmp_4d;
                        ExtractSlice(tmp_4d, MF_z2_5d,s,0);
                        MF_z2+=tmp_4d;
                    }

                    //perform the contraction.
                    for (unsigned int t = 0; t < nt; ++t)
                    {
                        for (unsigned int tx = 0; tx < nt; tx++)
                        {   
                            ty = (t + tx) % nt;

                            for (unsigned int tz1 = 0; tz1 < nt; tz1++)
                            {
                                for (unsigned int tz2 = 0; tz2 < nt; tz2++)
                                {
                                tmp_exch[t] += TensorRemove(MF_x[tx] * MF_z1[tz1] * MF_y[ty] * MF_z2[tz2]);
                                tmp_self[t] += TensorRemove(MF_x[tx] * MF_z1[tz1] * MF_z2[tz2] * MF_y[ty]);
                                }
                            }

                        }
                    }
                }

            }
            
        }
        if (i % 10 == 0)
        {
            LOG(Message) << "Mf for i: " << i << " of " << N << std::endl;
        }
    }

    double NTinv = 1.0 / static_cast<double>(nt);

    for (unsigned int t = 0; t < nt; ++t)
    {
        result.corr_exch[t] = NTinv * tmp_exch[t];
        result.corr_self[t] = NTinv * tmp_self[t];
    }
    saveResult(par().output, "mesonQEDexch", result);

}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonFieldConserved_hpp_
