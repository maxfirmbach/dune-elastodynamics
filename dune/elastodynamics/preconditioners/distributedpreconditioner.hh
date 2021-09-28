// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DISTRIBUTED_PRECONDITIONER_HH
#define DISTRIBUTED_PRECONDITIONER_HH

#include <dune/istl/preconditioner.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/elastodynamics/parallel/vectordatahandle.hh>

// implements a wrapper for sequential preconditioners applied to each rank
// and communicating the necessary values in between
// currently works for: Richardson, SeqJac (n=1!!!)

namespace Amg
{ template<class T> struct ConstructionTraits; }
  
template<class GridView, class P>
class DistributedPreconditioner : public Dune::Preconditioner<typename P::domain_type, typename P::range_type> {
    
  friend struct Amg::ConstructionTraits<DistributedPreconditioner<GridView, P>>;
  
  using X = typename P::domain_type;
  using Y = typename P::range_type;
    
  private:

    std::shared_ptr<P> preconditioner_;
    const GridView& gridView_;
    
  public:

    DistributedPreconditioner(const GridView& gridView, P& p)
      : gridView_(gridView),
        preconditioner_(stackobject_to_shared_ptr(p))
    {}

    DistributedPreconditioner(const GridView& gridView, const std::shared_ptr<P>& p)
      : gridView_(gridView),
        preconditioner_(p)
    {}

    virtual void pre (X& x, Y& b)
    {
      preconditioner_->pre(x,b);
    }

    virtual void apply (X& v, const Y& d)
    {
      preconditioner_->apply(v,d); 
      // communicate values here
      VectorExchangeAdd<GridView, X> datahandle(gridView_, v);
      gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    }

    template<bool forward>
    void apply (X& v, const Y& d)
    {
      preconditioner_->template apply<forward>(v,d);      
      // comunnicate values here
      VectorExchangeAdd<GridView, X> datahandle(gridView_, v);
      gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    }

    virtual void post (X& x)
    { preconditioner_->post(x); }

    virtual Dune::SolverCategory::Category category() const
    { return Dune::SolverCategory::sequential; }

};

#endif
