    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_cg_prec.cc

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class d>
struct scal {
  d internal;
};

  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
  //Set 5th dimension length
  const int Ls=8;

  // default precision lattice set up
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
  // Single precision lattice set up
  GridCartesian         * UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian         * FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
  GridRedBlackCartesian * FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_f);
  // double precsion lattice setup
  GridCartesian         * UGrid_d   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid_d = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_d);
  GridCartesian         * FGrid_d   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_d);
  GridRedBlackCartesian * FrbGrid_d = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  // Fermion fields
  LatticeFermionD    src(FGrid); random(RNG5,src);
  LatticeFermionF    tmp(FGrid_f); tmp=zero;
  LatticeFermionD result(FGrid_d); result=zero;

  // Increase precision of source

   //src = src * 0.000000000001;

  // creat site objects

  vSpinColourVectorD::scalar_object  site_src;
  vSpinColourVectorF::scalar_object  site_tmp;
  vSpinColourVectorD::scalar_object  site_result;

  // coordiates to poke

  std::vector<int> lcoor = {0, 0, 0, 0};
  
  // Gauge fields
  //LatticeGaugeFieldD Umu(UGrid);
  //LatticeGaugeFieldF Umu_f(UGrid_f);
  //LatticeGaugeFieldF Umu_f(UGrid_d);

  // Change precision

  precisionChange(tmp,src);

  precisionChange(result,tmp);

  // Peek Site

  peekSite(site_src,src,lcoor);
  peekSite(site_tmp,tmp,lcoor);
  peekSite(site_result,result,lcoor);

  // Logging

  std::cout  << GridLogMessage << "Source:" << site_src << std::endl;
    
  std::cout << GridLogMessage << "tmp:" << site_tmp << std::endl;

  std::cout << GridLogMessage << "result:" << site_result << std::endl;

  Grid_finalize();
}