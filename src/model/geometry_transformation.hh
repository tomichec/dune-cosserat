 #include <dune/grid/geometrygrid/coordfunction.hh>

template <class MODEL,int d>
 class GridTransformation
  : public Dune :: AnalyticalCoordFunction< double, d, d, GridTransformation <MODEL,d> >
  {
    typedef GridTransformation This;
    typedef Dune :: AnalyticalCoordFunction< double, d, d, This > Base;

  public:
    typedef typename Base :: DomainVector DomainVector;
    typedef typename Base :: RangeVector RangeVector;

    GridTransformation(MODEL& model_) : model(model_){}

    void evaluate ( const DomainVector &x, RangeVector &y ) const
    {

      y = model.evaluateGeometryTransformation(x);

    }

  private:

    MODEL& model;
  };
