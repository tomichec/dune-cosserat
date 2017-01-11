#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include "../helpful_functions.hh"

using namespace std;

namespace Dune {
    namespace PDELab {

template<typename PARAM, int dispdofel>
class Cosserat2D:
        public Dune::PDELab::NumericalJacobianApplyVolume<Cosserat2D<PARAM,dispdofel>>,
        public Dune::PDELab::NumericalJacobianVolume<Cosserat2D<PARAM,dispdofel>>,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
// pattern assembly flags
    enum { doPatternVolume = true };

// residual assembly flags
    enum { doAlphaVolume = true };
  //  enum { doLambdaBoundary = false };


Cosserat2D(PARAM& param_,
          unsigned int intorder_ = 0):
    param(param_), intorder(intorder_), time(0.0)
{}


// Volume Integral Depending on Test and Ansatz Functions
template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
{


    // Galerkin Finite Elements - assumes lfsu == lfsv
    const int dim = 2;

    // Unwrap shape local function spaces
    typedef typename LFSU::template Child<0>::Type LFSU_U1;
    typedef typename LFSU::template Child<1>::Type LFSU_U2;
    typedef typename LFSU::template Child<2>::Type LFSU_ROT3;

    typedef typename LFSU::template Child<3>::Type LFSU_P;

    const LFSU_U1& lfsu_u1 = lfsu.template child<0>();
    const LFSU_U2& lfsu_u2 = lfsu.template child<1>();
    const LFSU_ROT3& lfsu_rot3 = lfsu.template child<2>();

    const LFSU_P& lfsu_p = lfsu.template child<3>();

    const unsigned int p1_n = lfsu_p.size();
    const unsigned int p2_n = lfsu_u1.size();

   // domain and range field type
    typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;

    typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainType D;
    typedef typename R::value_type RF;
    typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianTypeU;
    typedef typename LFSU_P::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianTypeP;

    typedef typename LFSU_U1::Traits::SizeType size_type;

    typedef typename LFSU_P::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;


    // get geometry

  /*  auto geo = eg.geometry();

    // determine quadrature order
    const int u_order = lfsu_u1.finiteElement().localBasis().order();
    const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
    const int jac_order = geo.type().isSimplex() ? 0 : 1;
    const int qorder = 3*u_order - 1 + jac_order + det_jac_order + intorder; //

    // select quadrature rule
    GeometryType gt = geo.type();
    const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,qorder);

    // Unwrap solution at node, into vector d

    Dune::FieldVector<double,dispdofel> d(0.0);
    Dune::FieldVector<double,4> p(0.0);

    for (size_type i=0; i < p2_n; ++i){
        d[i] = x(lfsu_u1,i); // U1
        d[i + p2_n] = x(lfsu_u2,i); // U2
        if (i < p1_n){
          d[i + 2 * p2_n] = x(lfsu_rot3,i);
          p[i] = x(lfsu_p,i);
        }
    }


    // Loop over quadrature points

    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it = rule.begin(),endit = rule.end(); it != endit; ++it)
    {
       // Evaluate shape functions at Integration Point

        std::vector<Range> phi_P1(lfsu_p.size());
        lfsu_p.finiteElement().localBasis().evaluateFunction(it->position(),phi_P1);

        // Evaluate gradients of shape function at integration point
        std::vector<JacobianTypeP> js_P1(p1_n);
        lfsu_p.finiteElement().localBasis().evaluateJacobian(it->position(),js_P1);

        std::vector<JacobianTypeU> js_P2(p2_n);
        lfsu_u1.finiteElement().localBasis().evaluateJacobian(it->position(),js_P2);

        // Transform gradient to real element
        const typename EG::Geometry::JacobianInverseTransposed jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF,dim> > p2_gradphi(p2_n);
        std::vector<Dune::FieldVector<RF,dim> > p1_gradphi(p1_n);

        for (int i=0; i < p2_n; i++){
            if (i < p1_n){
                p1_gradphi[i] = 0.0;
                jac.umv(js_P1[i][0],p1_gradphi[i]);
            }
            p2_gradphi[i] = 0.0;
            jac.umv(js_P2[i][0],p2_gradphi[i]);
        }

        // (A) Compute Elastic part of material

        // Compute Element Strains from nodal displacements

        
        Dune::FieldMatrix<double,6,dispdofel> B(0.0);

        for (int i = 0; i < 4; i++){
            B[0][i] = p2_gradphi[i][0]; // E11
            B[1][i + 4] = p2_gradphi[i][1]; // E22

            B[2][i] = p2_gradphi[i][1];
            B[3][i + 4] = p2_gradphi[i][0];

            if (i < p1_n){
              B[2][i + 8] = phi_P1[i];
              B[3][i + 8] = -phi_P1[i];
              B[4][i + 8] = p1_gradphi[i][0];
              B[5][i + 8] = p1_gradphi[i][1];
            }
        }


        Dune::FieldVector<double,6> e(0.0);

        B.mv(d,e);

        //std::cout << e << std::endl;

        Dune::FieldVector<double,2> xg = geo.global(it->position());

        double ROT = param.getRotation(xg);

       // std::cout << "Rotation" << ROT << std::endl;

        Dune::FieldVector<double,6> el = rotateTensor3(e,ROT,0); // Rotate tensor from global to local rotation

        double e_d = el[0] + el[1]; // Volumetric Strain

        // Compute Constitutive Law

        double Vf = param.get_Vf0() / (1.0 + e_d); // Compute new fibre volume fraction

        // std::cout << Vf << std::endl;



        Dune::FieldMatrix<double,6,6> C = param.getElasticTensor(Vf,ROT);

        

        // Compute Stress
        Dune::FieldVector<double,6> sl(0.0);

        C.mv(el,sl);

        Dune::FieldVector<double,6> s = rotateTensor3(sl,ROT,1); // Rotate stress tensor from local to global

        Dune::FieldVector<double,dispdofel> res_u(0.0);

        B.mtv(s,res_u); // res_u = Bt * sig;


        //
        // (2) Compute Pressure Terms
        //

        Dune::FieldMatrix<double,2,4> G(0.0);

        for (size_type i = 0; i < p1_n; i++){
            G[0][i] = p1_gradphi[i][0];
            G[1][i] = p1_gradphi[i][1];
        }

        double temp = param.get_temp_temp();
        double alpha = param.get_alpha_temp();

        double mu = param.getVisocity(temp,alpha);

        Dune::FieldMatrix<double,dim,dim> K = param.getPermTensor(Vf,mu,ROT);

        Dune::FieldVector<double,2> gradp(0.0),gradp_l(0.0), flux_l(0.0), flux_g(0.0);

        G.mv(p,gradp); // compute pressure gradient  grad(p) = G * p

        gradp_l = rotateTensor2(gradp,ROT,0);

        K *= -1.0;

        K.mv(gradp_l,flux_l); // Compute flux = - Perm * grad(p)

        flux_g = rotateTensor2(flux_l,ROT,1);

        Dune::FieldVector<double,4> res_p(0.0);

        G.mtv(flux_g,res_p); // Compute residual vector res = - G' * Perm * G * p

        // -----------------------------------------------------------------
        //                        Coupling Term
        // -----------------------------------------------------------------

        Dune::FieldVector<double,6> delta(0.0); delta[0] = 1.0; delta[1] = 1.0;

        Dune::FieldMatrix<double,dispdofel,4> L(0.0);

        Dune::FieldVector<double,dispdofel> tmp(0.0);

        B.mtv(delta,tmp);

        for (size_type i = 0; i < dispdofel; i++){
            for (size_type j = 0; j < 4; j++){
                L[i][j] = tmp[i] * phi_P1[j];
            }
        }

        L.mmv(p,res_u);

        // Assemble Residuals

        // geometric weight
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());

        for (size_type i = 0; i < p2_n; i++){
           r.accumulate( lfsu_u1, i, res_u[i] * factor);
           r.accumulate( lfsu_u2, i, res_u[i + p2_n] * factor);
        }

        for (size_type i=0; i < p1_n; i++){
            r.accumulate( lfsu_rot3, i, res_u[i + 2 * p2_n] * factor);
            r.accumulate( lfsu_p, i, res_p[i] * factor);
        }

    } // end for each quadrature point

*/
    

} // end alpha_volume



 /*   // residual for boundary term
    template<typename IG, typename LFSV, typename R>
    void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
    {

        const int dim = 2;
        const int dimw = 2;

        // extract local function spaces
        typedef typename LFSV::template Child<0>::Type LFSV_U1;
        typedef typename LFSV::template Child<1>::Type LFSV_U2;

        const LFSV_U1& lfsv_u1 = lfsv.template child<0>();
        const LFSV_U2& lfsv_u2 = lfsv.template child<1>();

        const unsigned int nodes_per_element = lfsv_u1.size();

        // domain and range field type
        typedef typename LFSV_U1::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename R::value_type RF;
        typedef typename LFSV_U1::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSV_U1::Traits::SizeType size_type;
        typedef typename LFSV_U1::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeType RT_V;

        // select quadrature rule
        GeometryType gtface = ig.geometry().type();

        auto geo = ig.geometry();

        // determine quadrature order
        const int u_order = lfsv_u1.finiteElement().localBasis().order();
        const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-2);
        const int jac_order = geo.type().isSimplex() ? 0 : 1;
        const int qorder = 3*u_order - 1 + jac_order + det_jac_order + intorder; //;

        const QuadratureRule<DF,dim-1>& rule = QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        double totalForce = 0.0;

        if(param.isNeumann(ig.geometry().corner(0)))
        {

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it = rule.begin(), endit = rule.end();it != endit;++it)
        {
            // position of quadrature point in local coordinates of element
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate basis functions
            std::vector<RT_V> phi(nodes_per_element);
            lfsv_u1.finiteElement().localBasis().evaluateFunction(local,phi);

            // Compute Weight Factor
            const RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(it->position());

            // Need to add nonlinear rotation to the normal: OUTSTANDING

            // evaluate flux boundary condition - needs updating
            Dune::FieldVector<double,2> xg = ig.geometry().global(it->position());
            Dune::FieldVector<double,2> neumann_stress = param.getNeumann(xg, normal);

          //  std::cout << normal << std::endl;

            for (size_t i=0; i < nodes_per_element; i++){
                //std::cout << neumann_stress << std::endl;
                r.accumulate(lfsv_u1,i, -neumann_stress[0] * phi[i] * factor);
                r.accumulate(lfsv_u2,i, -neumann_stress[1] * phi[i] * factor);
            }

        } // For each quadrature point

    }

    */




        private:

            PARAM& param;
            int intorder;
            double time;
};

    }
}
