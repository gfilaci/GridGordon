/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_hadrons_wilsonFund.cc

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
    std::vector<std::string> flavour = {"l"};
    std::vector<double>      mass    = {-0.1};
    double                   csw     = 0.0;
    
    // global parameters
    Application::GlobalPar globalPar;
  
    globalPar.trajCounter.start = 309;
    globalPar.trajCounter.end   = 310;
    globalPar.trajCounter.step  = 1;
    globalPar.runId             = "test";

    application.setPar(globalPar);
    // gauge field
    application.createModule<MIO::LoadNersc>("gauge");
   
    // sources
    //MSource::Z2::Par z2Par;
    //z2Par.tA = 0;
    //z2Par.tB = 0;
    //application.createModule<MSource::Z2>("z2", z2Par);
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::WilsonClover::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        actionPar.csw_r = csw;
        actionPar.csw_t = csw;

        // !!!!! Check if Anisotropy works  !!!!!
        actionPar.clover_anisotropy.isAnisotropic= false;
        actionPar.clover_anisotropy.t_direction  = 3    ;   // Explicit for D=4
        actionPar.clover_anisotropy.xi_0         = 1.0  ;
        actionPar.clover_anisotropy.nu           = 1.0  ;

        application.createModule<MAction::WilsonClover>("WilsonClover_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "WilsonClover_" + flavour[i];
        solverPar.residual     = 1.0e-8;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);
        
        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i], quarkPar);
 //       quarkPar.source = "z2";
 //       application.createModule<MFermion::GaugeProp>("QZ2_" + flavour[i], quarkPar);
    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        MContraction::Meson::Par mesPar;
        
        mesPar.output  = "Fund_mesons/pt_" + flavour[i] + flavour[j];
        mesPar.q1      = "Qpt_" + flavour[i];
        mesPar.q2      = "Qpt_" + flavour[j];
        mesPar.gammas  = "all";
        mesPar.sink    = "sink";
        application.createModule<MContraction::Meson>("meson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
   //     mesPar.output  = "mesons/Z2_" + flavour[i] + flavour[j];
   //     mesPar.q1      = "QZ2_" + flavour[i];
   //     mesPar.q2      = "QZ2_" + flavour[j];
   //     mesPar.gammas  = "all";
   //     mesPar.sink    = "sink";
   //     application.createModule<MContraction::Meson>("meson_Z2_"
   //                                                   + flavour[i] + flavour[j],
   //                                                   mesPar);
    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    for (unsigned int k = j; k < flavour.size(); ++k)
    {
        MContraction::Baryon::Par barPar;
        
        barPar.output = "Fund_baryons/pt_" + flavour[i] + flavour[j] + flavour[k];
        barPar.q1     = "Qpt_" + flavour[i];
        barPar.q2     = "Qpt_" + flavour[j];
        barPar.q3     = "Qpt_" + flavour[k];
        application.createModule<MContraction::Baryon>(
            "baryon_pt_" + flavour[i] + flavour[j] + flavour[k], barPar);
    }
    
    // execution
    application.saveParameterFile("WilsonClover_spectrum.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
