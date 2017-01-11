/** 
 \brief A function that defines Dirichlet Boundary conditions AND the extension to the interior of the domain
 
*/

#ifndef BoundaryConditions_h
#define BoundaryConditions_h

// Define Scalar Dirichlet Boundary Conditions

template<typename GV, typename RF, class MODEL>
class Scalar_BC :
public Dune::PDELab::AnalyticGridFunctionBase<
Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
Scalar_BC<GV,RF,MODEL> >,
public Dune::PDELab::InstationaryFunctionDefaults
{

    double time;
    
public:
    
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, Scalar_BC<GV,RF,MODEL>> BaseT;
    
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    
    // Constructor
    Scalar_BC(const GV & gv, MODEL& model_) : time(0.0), BaseT(gv), model(model_)
    {
    }
    
    template<typename I>
    bool isDirichlet(const I & ig,
                     const typename Dune::FieldVector<typename I::ctype, I::dimension-1> & xlocal
                     ) const
    {   
        /* isDirichlet - function returns true if Dirchlet boundary conditions using boundary data from gmsh */
       

        Dune::FieldVector<double,2> x = ig.geometry().global(xlocal);

        bool Dirich = model.isDirichlet(x,dof);

        return Dirich;
    }
    
    
    inline void evaluateGlobal(const DomainType & x, RangeType & value) const
    {
        value = model.evaluateDirichlet(x,dof,time);

     
    } // end inline function evaluateGlobal
    
    void setDof(int degree_of_freedom){
        dof = degree_of_freedom;
    }
    
    void setTime(double t) { time = t; }
    
private:
    int dof;
    MODEL& model;
    
};




#endif /* BoundaryConditions_h */