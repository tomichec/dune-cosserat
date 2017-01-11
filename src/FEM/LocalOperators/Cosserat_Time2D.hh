#include <dune/pdelab/localoperator/idefault.hh>


namespace Dune {
    namespace PDELab {

template <class PARAM,int dispdofel>
class CosseratTime2D:
        public Dune::PDELab::NumericalJacobianApplyVolume<CosseratTime2D<PARAM,dispdofel>>,
        public Dune::PDELab::NumericalJacobianVolume<CosseratTime2D<PARAM,dispdofel>>,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
// pattern assembly flags
    enum { doPatternVolume = true };

// residual assembly flags
    enum { doAlphaVolume = true };


CosseratTime2D(PARAM& param_,unsigned int intorder_=1):
    param(param_),intorder(intorder_), time(0.0)
{}


// set time
void setTime(double t){time = t;}

// Volume Integral Depending on Test and Ansatz Functions
template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
{

    // Galerkin Finite Elements - assumes lfsu == lfsv
    const int dim = 2;

    // extract local function spaces
    typedef typename LFSU::template Child<0>::Type LFSU_U1;
    typedef typename LFSU::template Child<1>::Type LFSU_U2;
    typedef typename LFSU::template Child<2>::Type LFSU_ROT3;
    typedef typename LFSU::template Child<3>::Type LFSU_P;

    const LFSU_U1& lfsu_u1 = lfsu.template child<0>();
    const LFSU_U2& lfsu_u2 = lfsu.template child<1>();
    const LFSU_ROT3& lfsu_rot3 = lfsu.template child<2>();
    const LFSU_P& lfsu_p = lfsu.template child<3>();

    const int p1_n = lfsu_p.size();
    const int p2_n = lfsu_u1.size();

    /*

   // domain and range field type
    typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainType D;
    typedef typename R::value_type RF;
    typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianType;
    typedef typename LFSU_U1::Traits::SizeType size_type;
    typedef typename LFSU_P::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;

    auto geo = eg.geometry();

    // determine quadrature order
    const int u_order = lfsu_u1.finiteElement().localBasis().order();
    const int det_jac_order = geo.type().isSimplex() ? 0 : (dim-1);
    const int jac_order = geo.type().isSimplex() ? 0 : 1;
    const int qorder = 3*u_order - 1 + jac_order + det_jac_order + intorder;


    // select quadrature rule
    GeometryType gt = geo.type();
    const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,qorder);

    Dune::FieldVector<double,dispdofel> d(0.0);
    Dune::FieldVector<double,4> p(0.0);

    for (size_type i=0; i < p2_n; ++i){
        d[i] = x(lfsu_u1,i); // U1
        d[i + lfsu_u1.size()] = x(lfsu_u2,i); // U2

        if (i < p1_n){
          d[i + 2 * lfsu_u1.size()] = x(lfsu_rot3,i); // U3
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
        std::vector<JacobianType> js_P1(p1_n);
        lfsu_p.finiteElement().localBasis().evaluateJacobian(it->position(),js_P1);

        std::vector<JacobianType> js_P2(p2_n);
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

        Dune::FieldMatrix<double,6,22> B(0.0);

        for (int i = 0; i < 4; i++){
            B[0][i] = p2_gradphi[i][0]; // E11
            B[1][i + 4] = p2_gradphi[i][1]; // E22

            B[2][i] = p2_gradphi[i][1];
            B[3][i + 4] = p2_gradphi[i][0];

            if (i < 4){
              B[2][i + 8] = phi_P1[i];
              B[3][i + 8] = -phi_P1[i];
              B[4][i + 8] = p1_gradphi[i][0];
              B[5][i + 8] = p1_gradphi[i][1];
            }
        }

         Add coupling term 

        Dune::FieldVector<double,3> res_p(0.0);

        Dune::FieldVector<double,6> delta(0.0); delta[0] = 1.0; delta[1] = 1.0;

        Dune::FieldMatrix<double,dispdofel,4> L(0.0);

        Dune::FieldVector<double,dispdofel> tmp(0.0);

        B.mtv(delta,tmp);

        for (size_type i = 0; i < dispdofel; i++){
            for (size_type j = 0; j < 4; j++){
                L[i][j] = tmp[i] * phi_P1[j];
            }
        }

        L.mtv(d,res_p);

        // geometric weight
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());


        //std::cout << d << std::endl;

        for (size_type i = 0; i < p1_n; i++){
            r.accumulate(lfsu_p,i,-res_p[i] * factor);
        }

    } // end for each quadrature point

    */

} // end alpha_volume


        private:
            PARAM& param;
            int intorder;
            double time;

};

    }
}
