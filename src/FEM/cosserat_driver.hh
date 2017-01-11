#include "boundaryConditions.hh"
#include "LocalOperators/Cosserat2D.hh"
#include "LocalOperators/Cosserat_Spatial2D.hh"
#include "LocalOperators/Cosserat_Time2D.hh"

template <typename MODEL, typename GRID>
void cosserat_driver(MODEL& model, const GRID& grid)
{

using Dune::PDELab::Backend::native;

double time = 0.0; // Initialise time!

// Get leafGridView from grid
    typedef typename GRID::LeafGridView GV;
      GV gv = grid.leafGridView();

// Construct Finite Element Space
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,typename GV::Grid::ctype,double,1> FEM_P1;
      FEM_P1 P1(gv);
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,typename GV::Grid::ctype,double,1> FEM_P2;
      FEM_P2 P2(gv);

// Initialise Various Types for Finite element Calculations

typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
 typedef Dune::PDELab::istl::VectorBackend<> VBE;

typedef Dune::PDELab::GridFunctionSpace<GV, FEM_P1, CON, VBE> P1_GFS; // Piecewise Linear Grid Function Space
typedef Dune::PDELab::GridFunctionSpace<GV, FEM_P2, CON, VBE> P2_GFS; // Piecewise Quadratic Grid Function Space


// Initialise Nodal Variables

    P2_GFS dispU1(gv,P2); dispU1.name("U1");
    P2_GFS dispU2(gv,P2); dispU2.name("U2");
    
    P1_GFS rotR3(gv,P1); rotR3.name("ROT3");
    P1_GFS p(gv,P1); p.name("PRESSURE");

// Wrap variables into a Vector Values Function

typedef Dune::PDELab::CompositeGridFunctionSpace <VBE,Dune::PDELab::LexicographicOrderingTag, P2_GFS, P2_GFS, P1_GFS,P1_GFS> GFS;
    GFS gfs(dispU1, dispU2, rotR3,p);

 // Make constraints map and initialize it from a function
typedef typename GFS::template ConstraintsContainer<double>::Type C;
    C cg; cg.clear();



    typedef Scalar_BC<GV,double,MODEL> BC;
    BC U1_cc(gv,model), U2_cc(gv,model), ROT3_cc(gv,model), p_cc(gv,model);
        U1_cc.setDof(1); U1_cc.setTime(time);
        U2_cc.setDof(2); U2_cc.setTime(time);
        ROT3_cc.setDof(3); ROT3_cc.setTime(time);
        p_cc.setDof(4); p_cc.setTime(time);


 typedef Dune::PDELab::CompositeConstraintsParameters<BC,BC,BC,BC>
    Constraints;

    Constraints constraints(U1_cc,U2_cc,ROT3_cc,p_cc);

  Dune::PDELab::constraints(constraints,gfs,cg);

// Linear Operators

  typedef CosseratSpatial2D<MODEL,12> LOP;
    LOP     lop(model,1);

  typedef Dune::PDELab::CosseratTime2D<MODEL,12> TLOP;
    TLOP    tlop(model,1);

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
    MBE mbe(36); // 68 if quadratic 36 linear



  // ==== Contruct Grid Operators

  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,C,C> GO0;
  GO0 go0(gfs,cg,gfs,cg,lop,mbe);

  typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,double,double,double,C,C> GO1;
    GO1 go1(gfs,cg,gfs,cg,tlop,mbe);

  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
    IGO igo(go0,go1);

// << Make FE function with initial values >>
    typedef typename IGO::Traits::Domain V;
      V xold(gfs,0.0);

     model.setPressure(time);

    BC u1(gv,model), u2(gv,model), rot3(gv,model), pressure(gv,model);
    u1.setDof(1); u1.setTime(time);
    u2.setDof(2); u2.setTime(time);
    rot3.setDof(3); rot3.setTime(time);
    pressure.setDof(4); pressure.setTime(time);


    typedef Dune::PDELab::CompositeGridFunction<BC,BC,BC,BC> InitialSolution;
      InitialSolution initial_solution(u1,u2,rot3,pressure);

    Dune::PDELab::interpolate(initial_solution,gfs,xold);


// Output Results to file


    if (model.getVerb_IG()){

        stringstream sstm;
        sstm << model.getOutputFolder() << "/InitialSolution";

        string file_output_name = sstm.str();

      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
      Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,xold);
      vtkwriter.write(file_output_name,Dune::VTK::appendedraw);

      std::cout << "Initial Solution & Mesh written to : " << file_output_name << ".vtk" << std::endl;

    }
    typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<IGO> LS;
    // typedef Dune::PDELab::ISTLBackend_SEQ_UMFPack LS;
      LS ls(100,false); // Note that first argument is 'Max Iterations' which does not apply for direct solver and is ignored

    // <<<6>>> Solver for non-linear problem per stage
    typedef Dune::PDELab::Newton<IGO,LS,V> PDESOLVER;
     PDESOLVER pdesolver(igo,ls);
      pdesolver.setReassembleThreshold(model.getNewtonAssemblyThreshold());
      pdesolver.setVerbosityLevel(model.getNewtonVerbLevel());
      pdesolver.setReduction(model.getNewtonReduction());
      pdesolver.setMinLinearReduction(model.getNewtonLineReduction());
      pdesolver.setMaxIterations(model.getNewtonMaxIterations());
      pdesolver.setLineSearchMaxIterations(model.getNewtonLineSearchMaxIterations());

    // << Define time-stepper >>

    Dune::PDELab::Alexander2Parameter<double> method;
    
    Dune::PDELab::OneStepMethod<double,IGO,PDESOLVER,V,V> osm(method,igo,pdesolver);

    osm.setVerbosityLevel(model.getOSMVerb());

    // =========== TIME - INTERGRATION LOOP ===========

    V xnew(gfs,xold), res(gfs,0.0);

    double dt = model.get_dt();

    int tstep = 0;

    model.setStep(tstep);

    model.set_temp(model.get_initialTemp());
    model.set_alpha(model.get_initialAlpha());

   /* while (time < model.get_tend() - 1e-8){

      double T0 = model.get_new_temp(time,0); // Get Temperature at Beginning of the timestep

      double T;

      double a0 = model.get_alpha(); // Get Alpha at the beginning of the timestep

      model.set_alpha_temp(a0); // Initialise alpha temp

      model.setStep(tstep); // Set the new time step


      if (model.getVerb_TI() > 0){
        std::cout << "| " << tstep << " | " << time << " | " << dt << " | " << a0 << " | " << T0 << " | " << model.get_pressure() << " | " << std::endl;
      }

      //  >>>>>>> Adaptive time Stepping Algorithm <<<<<<<

      double dt_star = 1e-20;
      int m = 1;

      if(model.isAdaptive()){m = model.getAdaptM();}

      V x_mdt(gfs,0.0),  x_dt(gfs,0.0), x_dt_old(gfs,0.0);


      // >>> Adaptive Step Loop

      while (0.8 * dt > dt_star){ // Current time step dt is too big for adaptive step criterion

        if (counter > 0){ dt = std::min(dt_star, model.dtmax); } // if not the first time in this loop
        counter +=1; // increment counter

        // ======== Solve big time step - m * dt

                double mdt = m * dt;

               // Compute constraints for this time step
               cg.clear();

               model.setPressure(time + mdt); // Set the current pressure boundary conditions

               if (model.getVerb_TI() > 1){
                std::cout << "== Computing Big Step ... ";
               }

              typedef Scalar_BC<GV,double,MODEL> BC;
                BC U1_cc(gv,model), U2_cc(gv,model), ROT3_cc(gv,model),p_cc(gv,model);
              U1_cc.setDof(1); U1_cc.setTime(time + mdt);
              U2_cc.setDof(2); U2_cc.setTime(time + mdt);
              ROT3_cc.setDof(3); ROT3_cc.setTime(time + mdt);
              p_cc.setDof(4); p_cc.setTime(time + mdt);

              typedef Dune::PDELab::CompositeConstraintsParameters<BC,BC,BC,BC>
                Constraints;

              Constraints constraints(U1_cc,U2_cc,ROT3_cc,p_cc);

              Dune::PDELab::constraints(constraints,gfs,cg);
                BC u1(gv,model), u2(gv,model), rot3(gv,model), pressure(gv,model);
              u1.setDof(1); u1.setTime(time + mdt);
              u2.setDof(2); u2.setTime(time + mdt);
              rot3.setDof(3); rot3.setTime(time + mdt);
              pressure.setDof(4); pressure.setTime(time + mdt);

              typedef Dune::PDELab::CompositeGridFunction<BC,BC,BC,BC> InitialSolution;
              InitialSolution initial_solution(u1,u2,rot3,pressure);

              T = model.get_new_temp(time,mdt); // Update temperature according to the cure cycle
              
              model.cureModel(a0,T,mdt); // Update cure according to cure model and new temperature

              osm.apply(time,mdt,xold,initial_solution,x_mdt); // Solve big time step - solution stored to x_mdt

              if (model.getVerb_TI() > 1){
                std::cout << "Solved. " << std::endl;
              }


        // ======== Solve m smalltime steps

        double a = a0;

        x_dt_old = xold;

        for (int i = 1; i <= m; i++){

           model.setPressure(time + i * dt);

           cg.clear();

           if (model.getVerb_TI() > 1){
            std::cout << "== Computing Small Step " << std::endl;
           }

          typedef Scalar_BC<GV,double,MODEL> BC;
            BC U1_cc(gv,model), U2_cc(gv,model), ROT3_cc(gv,model),p_cc(gv,model);
          U1_cc.setDof(1); U1_cc.setTime(time + i * dt);
          U2_cc.setDof(2); U2_cc.setTime(time + i * dt);
          ROT3_cc.setDof(3); ROT3_cc.setTime(time + i * dt);
          p_cc.setDof(4); p_cc.setTime(time + i * dt);

          typedef Dune::PDELab::CompositeConstraintsParameters<BC,BC,BC,BC> Constraints;

          Constraints constraints(U1_cc,U2_cc,ROT3_cc,p_cc);

          Dune::PDELab::constraints(constraints,gfs,cg);

          BC u1(gv,model), u2(gv,model), rot3(gv,model), pressure(gv,model);
          u1.setDof(1); u1.setTime(time + i * dt);
          u2.setDof(2); u2.setTime(time + i * dt);
          rot3.setDof(3); rot3.setTime(time + i * dt);
          pressure.setDof(4); pressure.setTime(time + i * dt);

        typedef Dune::PDELab::CompositeGridFunction<BC,BC,BC,BC> InitialSolution;
          InitialSolution initial_solution(u1,u2,rot3,pressure);

        T = model.get_new_temp(time,i * dt);
        if (i > 1) { a = model.get_alpha_temp(); }
        model.cureModel(a,T,dt);

        osm.apply(time + (i-1) * dt,dt,x_dt_old,initial_solution,x_dt);

        if (model.getVerb_TI() > 1){
                std::cout << "Solved. " << std::endl;
              }

        if (i < m){ x_dt_old = x_dt; } // If not the last time step

     } // end for i dt time steps



     // == Estimate the temporal error & time step

     double twoNorm = 0.0;

     for (int i = 0; i < native(x_dt).size(); i++){
       twoNorm += std::pow(native(x_dt)[i] - native(x_mdt)[i],2);
     }

     dt_star = dt;


     if (model.isAdaptive()){
        dt_star = std::sqrt(model.tolerance * (dt * dt) * (m * m - 1.0) / std::sqrt(twoNorm));
     }

     if ((model.getVerb_TI() > 1) && model.adaptive()){
      std::cout << "Suggested Time Step = " << dt_star << std::endl;
     }




      } // End Adaptive Step Loop


      if(model.isAdaptive()){
          // Update time step
         for (int i = 0; i < native(xnew).size(); i++){
           native(xnew)[i] = ((m * m) * native(x_dt)[i] - native(x_mdt)[i]) / (m * m - 1.0);
         }
      }

     time += m * dt;

     model.recordTime(time); // Record at end of time steps

     model.set_temp(T); // Set current temperature
     model.setCure(); // Set current cure = temp cure

     tstep += 1;

     dt = std::min(dt_star,model.get_dtmax()); // Suggest new time step

     // Graphics!

    // Output Results to file

     if (model.getVerb_S()){


        stringstream sstm;
        sstm << model.getOutputFolder() << "/Solution_"; << tstep;
        string file_output_name = sstm.str();
  
        Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
        Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,xnew);
        vtkwriter.write(file_output_name,Dune::VTK::appendedraw);
      }

      xold = xnew; // The new becomes old! Before going onto the next time step


    } // end time integration loop

    */


} // End Cosserat_driver
