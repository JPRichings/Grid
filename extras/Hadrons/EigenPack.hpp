/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/EigenPack.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_EigenPack_hpp_
#define Hadrons_EigenPack_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/algorithms/iterative/Deflation.h>
#include <Grid/algorithms/iterative/LocalCoherenceLanczos.h>

BEGIN_HADRONS_NAMESPACE

// Lanczos type
#ifndef HADRONS_DEFAULT_LANCZOS_NBASIS
#define HADRONS_DEFAULT_LANCZOS_NBASIS 60
#endif

template <typename F>
class EigenPack
{
public:
    typedef F Field;
    struct PackRecord
    {
        std::string operatorXml, solverXml;
    };
    struct VecRecord: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(VecRecord,
                                        unsigned int, index,
                                        double,       eval);
        VecRecord(void): index(0), eval(0.) {}
    };
public:
//F=Grid::Lattice<Grid::QCD::vSpinColourVector>
    std::vector<RealD> eval;
    std::vector<F>     evec;
    //std::vector<LatticeSpinColourVectorF> evectmp;
    //std::vector<Grid::Lattice<Grid::QCD::vSpinColourVectorF>> evectmp;    // Grid::Lattice<Grid::QCD::vSpinColourVector> // Grid::Lattice<Grid::QCD::vSpinColourVectorF>
    std::vector<F>     evec_result; // Grid::Lattice<Grid::QCD::vSpinColourVector>
    PackRecord         record;

public:
    EigenPack(void)          = default;
    virtual ~EigenPack(void) = default;

    EigenPack(const size_t size, GridBase *grid)
    {
        resize(size, grid);
    }

    void resize(const size_t size, GridBase *grid)
    {
        eval.resize(size);
        evec.resize(size, grid);
        evectmp.resize(size, grid);
        evec_result.resize(size, grid);

        const int Ls = 16;

        // Single precision lattice set up
        GridCartesian         * UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
        GridRedBlackCartesian * UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
        GridCartesian         * FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
        GridRedBlackCartesian * FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_f);
        LatticeFermionF tmp(FGrid_f); tmp=zero;

    }

    virtual void read(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < evec.size(); ++k)
            {
                basicReadSingle(evec[k], eval[k], evecFilename(fileStem, k, traj), k, evec_result[k]);
            }
        }
        else
        {
            basicRead(evec, eval, evecFilename(fileStem, -1, traj), evec.size(), evec_result);
        }
    }

    virtual void write(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < evec.size(); ++k)
            {
                basicWriteSingle(evecFilename(fileStem, k, traj), evec[k], eval[k], k);
            }
        }
        else
        {
            basicWrite(evecFilename(fileStem, -1, traj), evec, eval, evec.size());
        }
    }
protected:
    std::string evecFilename(const std::string stem, const int vec, const int traj)
    {
        std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

        if (vec == -1)
        {
            return stem + t + ".bin";
        }
        else
        {
            return stem + t + "/v" + std::to_string(vec) + ".bin";
        }
    }

    template <typename T>
    void basicRead(std::vector<T> &evec, std::vector<double> &eval,
                   const std::string filename, const unsigned int size, std::vector<T> &evec_result)
    {
        ScidacReader    binReader;
        //LatticeFermionF evecbuf;

        binReader.open(filename);
        binReader.skipPastObjectRecord(SCIDAC_FILE_XML);
        for(int k = 0; k < size; ++k) 
        {
            VecRecord vecRecord;

            LOG(Message) << "Reading eigenvector " << k << std::endl;
            binReader.readScidacFieldRecord(evec[k], vecRecord);
            if (vecRecord.index != k)
            {
                HADRONS_ERROR(Io, "Eigenvector " + std::to_string(k) + " has a"
                              + " wrong index (expected " + std::to_string(vecRecord.index) 
                              + ") in file '" + filename + "'");
            }
            eval[k] = vecRecord.eval;
            if(true)
            {
                // convert the eigen values to single precision
                //RealF tmp = (RealF) eval[k];
                //eval[k] = (RealD) tmp;
                // convert the eigen vectors to single precision
                LOG(Message) << "beforeCast" << std::endl;
                precisionChange(tmp, evec[k]);
                LOG(Message) << "duringCast" << std::endl;
                precisionChange(evec_result[k], tmp); // this is the issue!!!

                vSpinColourVector::scalar_object  site_evec;
                
                
                
                LOG(Message) << "double: " << norm2(evec[k]) << std::endl;
                LOG(Message) << "singletmp: " << norm2(tmp) << std::endl;
                LOG(Message) << "single: " << norm2(evec_result[k]) << std::endl;
                

                std::vector<int> lcoor = {0, 0, 0, 0};
                //peekSite(site_evec, dummy, lcoor);
                LOG(Message) << "evec site: " << evec[k]._odata[0] << std::endl;
                LOG(Message) << "evectmp site: " << tmp._odata[0] << std::endl;
                LOG(Message) << "evec_result site: " << evec_result[k]._odata[0] << std::endl;
                evec[k] = evec_result[k] - evec[k];
                LOG(Message) << "diff: " << norm2(evec[k]) << std::endl;
                LOG(Message) << "evec diff site: " << evec[k]._odata[0] << std::endl;
                evec[k] = evec_result[k];

            }
        }
        binReader.close();
    }

    template <typename T>
    void basicReadSingle(T &evec, double &eval, const std::string filename, 
                         const unsigned int index, T &evec_result)
    {
        ScidacReader binReader;
        VecRecord    vecRecord;

        binReader.open(filename);
        binReader.skipPastObjectRecord(SCIDAC_FILE_XML);
        LOG(Message) << "Reading eigenvector " << index << std::endl;
        binReader.readScidacFieldRecord(evec, vecRecord);
        if (vecRecord.index != index)
        {
            HADRONS_ERROR(Io, "Eigenvector " + std::to_string(index) + " has a"
                          + " wrong index (expected " + std::to_string(vecRecord.index) 
                          + ") in file '" + filename + "'");
        }
        LOG(Message) << "Before eval extracted" << std::endl;
        eval = vecRecord.eval;
        LOG(Message) << "After eval extracted" << std::endl;
        if(true)
        {
            //site print to log set up
            std::vector<int> lcoor = {0, 0, 0, 0};
            //vSpinColourVectorD::scalar_object  site_evec; //LatticeSpinColourVectorF
            //vSpinColourVectorD::scalar_object  site_result;

            //convert the eigen values to single precision
            //RealF tmp = (RealF) eval[k];
            //eval[k] = (RealD) tmp;

            // Convert eigenvectors to Single precision

            LOG(Message) << "Before Precision change" << std::endl; 
            
            LOG(Message) << typeid(evec).name() << ":" << typeid(evec).name() << std::endl;

            precisionChange(evectmp[0], evec);
            
            precisionChange(evec_result, evectmp[0]);

            //peekSite(site_evec, evec, lcoor);
            //peekSite(site_result, evec_result, lcoor);

            //localConvertJPR(evec, evectmp[0]);
            LOG(Message) << "After Precision change" << std::endl;
            LOG(Message) << "norm2 double: " << norm2(evec) << std::endl;
            LOG(Message) << "norm2 single: " << norm2(evec_result) << std::endl;
            evec = evec - evec_result;
            LOG(Message) << "norm2 diff: " << evec << std::endl;
            //LOG(Message) << "tmp: " << evectmp[0] << std::endl;
            //LOG(Message) << "evec site: " << site_evec << std::endl;
            //LOG(Message) << "result site: " << site_result << std::endl;
            evec = evec_result;

        }
        binReader.close();
    }

    template <typename T>
    void basicWrite(const std::string filename, std::vector<T> &evec, 
                    const std::vector<double> &eval, const unsigned int size)
    {
        ScidacWriter binWriter(evec[0]._grid->IsBoss());
        XmlWriter    xmlWriter("", "eigenPackPar");

        makeFileDir(filename, evec[0]._grid);
        xmlWriter.pushXmlString(record.operatorXml);
        xmlWriter.pushXmlString(record.solverXml);
        binWriter.open(filename);
        binWriter.writeLimeObject(1, 1, xmlWriter, "parameters", SCIDAC_FILE_XML);
        for(int k = 0; k < size; ++k) 
        {
            VecRecord vecRecord;

            vecRecord.index = k;
            vecRecord.eval  = eval[k];
            LOG(Message) << "Writing eigenvector " << k << std::endl;
            binWriter.writeScidacFieldRecord(evec[k], vecRecord, DEFAULT_ASCII_PREC);
        }
        binWriter.close();
    }

    template <typename T>
    void basicWriteSingle(const std::string filename, T &evec, 
                          const double eval, const unsigned int index)
    {
        ScidacWriter binWriter(evec._grid->IsBoss());
        XmlWriter    xmlWriter("", "eigenPackPar");
        VecRecord    vecRecord;

        makeFileDir(filename, evec._grid);
        xmlWriter.pushXmlString(record.operatorXml);
        xmlWriter.pushXmlString(record.solverXml);
        binWriter.open(filename);
        binWriter.writeLimeObject(1, 1, xmlWriter, "parameters", SCIDAC_FILE_XML);
        vecRecord.index = index;
        vecRecord.eval  = eval;
        LOG(Message) << "Writing eigenvector " << index << std::endl;
        binWriter.writeScidacFieldRecord(evec, vecRecord, DEFAULT_ASCII_PREC);
        binWriter.close();
    }
};

template <typename FineF, typename CoarseF>
class CoarseEigenPack: public EigenPack<FineF>
{
public:
    typedef CoarseF CoarseField;
public:
    std::vector<RealD>   evalCoarse;
    std::vector<CoarseF> evecCoarse;
public:
    CoarseEigenPack(void)          = default;
    virtual ~CoarseEigenPack(void) = default;

    CoarseEigenPack(const size_t sizeFine, const size_t sizeCoarse, 
                    GridBase *gridFine, GridBase *gridCoarse)
    {
        resize(sizeFine, sizeCoarse, gridFine, gridCoarse);
    }

    void resize(const size_t sizeFine, const size_t sizeCoarse, 
                GridBase *gridFine, GridBase *gridCoarse)
    {
        EigenPack<FineF>::resize(sizeFine, gridFine);
        evalCoarse.resize(sizeCoarse);
        evecCoarse.resize(sizeCoarse, gridCoarse);
    }

    void readFine(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < this->evec.size(); ++k)
            {
                this->basicReadSingle(this->evec[k], this->eval[k], this->evecFilename(fileStem + "_fine", k, traj), k, this->evec[k]);
            }
        }
        else
        {
            this->basicRead(this->evec, this->eval, this->evecFilename(fileStem + "_fine", -1, traj), this->evec.size(), this->evec);
        }
    }

    void readCoarse(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < evecCoarse.size(); ++k)
            {
                this->basicReadSingle(evecCoarse[k], evalCoarse[k], this->evecFilename(fileStem + "_coarse", k, traj), k, evecCoarse[k]);
            }
        }
        else
        {
            this->basicRead(evecCoarse, evalCoarse, this->evecFilename(fileStem + "_coarse", -1, traj), evecCoarse.size(), evecCoarse);
        }
    }

    virtual void read(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        readFine(fileStem, multiFile, traj);
        readCoarse(fileStem, multiFile, traj);
    }

    void writeFine(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < this->evec.size(); ++k)
            {
                this->basicWriteSingle(this->evecFilename(fileStem + "_fine", k, traj), this->evec[k], this->eval[k], k);
            }
        }
        else
        {
            this->basicWrite(this->evecFilename(fileStem + "_fine", -1, traj), this->evec, this->eval, this->evec.size());
        }
    }

    void writeCoarse(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        if (multiFile)
        {
            for(int k = 0; k < evecCoarse.size(); ++k)
            {
                this->basicWriteSingle(this->evecFilename(fileStem + "_coarse", k, traj), evecCoarse[k], evalCoarse[k], k);
            }
        }
        else
        {
            this->basicWrite(this->evecFilename(fileStem + "_coarse", -1, traj), evecCoarse, evalCoarse, evecCoarse.size());
        }
    }
    
    virtual void write(const std::string fileStem, const bool multiFile, const int traj = -1)
    {
        writeFine(fileStem, multiFile, traj);
        writeCoarse(fileStem, multiFile, traj);
    }
};

template <typename FImpl>
using FermionEigenPack = EigenPack<typename FImpl::FermionField>;

template <typename FImpl, int nBasis>
using CoarseFermionEigenPack = CoarseEigenPack<
    typename FImpl::FermionField,
    typename LocalCoherenceLanczos<typename FImpl::SiteSpinor, 
                                   typename FImpl::SiteComplex, 
                                   nBasis>::CoarseField>;

END_HADRONS_NAMESPACE

#endif // Hadrons_EigenPack_hpp_
