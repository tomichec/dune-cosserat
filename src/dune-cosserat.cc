#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <stdio.h>
#include <string.h>

#include <dune/common/parametertree.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/grid/yaspgrid.hh> // Checked Inclusion
#include <dune/grid/common/gridview.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/gridinfo.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/io.hh>
#include <dune/istl/matrixmarket.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>   
#include <dune/common/timer.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>
#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/pdelab/finiteelementmap/rannacherturekfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/common/instationaryfilenamehelper.hh>
#include <dune/pdelab/instationary/implicitonestep.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
//#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/instationary/onestepparameter.hh>

#include <dune/pdelab/backend/istl/seqistlsolverbackend.hh>

 #include <dune/grid/geometrygrid/grid.hh>

#include <dune/common/parametertree.hh>


// Header files for dune-cosserat

#include "model/model.hh"
#include "model/geometry_transformation.hh"

#include "FEM/cosserat_driver.hh"



int main(int argc, char** argv)
{
   Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);


    // Setup My Model

    //Read ini file

    Dune::ParameterTree config;

    Dune::ParameterTreeParser parser;
    parser.readINITree(argv[1],config);


    const int dim = 2;

    MODEL<dim> myModel(config);

  // Setup Periodic Grid    
    


    // Set up Yasp Grid

    typedef typename Dune::YaspGrid<dim> YGRID;
      YGRID yaspgrid(myModel.getL(),myModel.getN(),myModel.getPeriodic(),myModel.getOverlap());

    // Tranform Yasp Grid to Part Geometry

   GridTransformation<MODEL<dim>, dim> mytransformation(myModel);
   typedef typename Dune::GeometryGrid<YGRID,GridTransformation<MODEL<dim>,dim>> GRID;
    GRID grid(yaspgrid,mytransformation);

    // Cosserat Driver

    cosserat_driver<MODEL<dim>, GRID>(myModel,grid);





    return 0;
  
}


