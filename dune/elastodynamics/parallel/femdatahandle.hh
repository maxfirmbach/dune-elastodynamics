// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*
This is a prototype version for a FEM data handle!
*/

#ifndef VECTOR_DATA_HANDLE_HH
#define VECTOR_DATA_HANDLE_HH

#include <dune/grid/common/datahandleif.hh>

namespace Type
{
  struct Equal
  {
    template<class A, class B>
    static void apply(A& a, const B& b)
    { a = b; }
  };
    
  struct Add
  {
    template<class A, class B>
    static void apply(A& a, const B& b)
    { a += b; }
  };
}

template<class Basis, class Vector, class Operation>	
class VectorExchange : public Dune::CommDataHandleIF<VectorExchange<Basis, Vector, Operation>, typename Vector::value_type> {

  private:
      
    const Basis& basis_;
    Vector& vector_;

  public:
      
    typedef typename Vector::value_type DataType;
    
    VectorExchange(const Basis& basis, Vector& vector)  
      : basis_(basis), vector_(vector)
    {}
    
    // we want to first operate on the element
    bool contains(int dim, int codim) const
    { return (codim == 0); }
    
    bool fixedSize(int dim, int codim) const
    { return true; }
    
    template<class Entity>
    size_t size(const Entity& entity) const
    { return 1; }

    // This doesn't work like I expected it to be ...... 

    template<class MessageBuffer, class Entity>
    void gather(MessageBuffer& buffer, const Entity& entity) const
    {
      auto localView = basis_.localView();
      //localView.bind(entity);
      const auto& localFE = localView.tree().child(0).finiteElement();
      for( size_t i=0; i<localFE.size(); i++) {
        buffer.write(vector_[localView.index(i)[0]]);
      }   
    }
    
    template<class MessageBuffer, class Entity>
    void scatter(MessageBuffer& buffer, const Entity& entity, size_t n) const
    {
      auto localView = basis_.localView();
      localView.bind(entity);
      DataType x;
      buffer.read(x);
      const auto& localFE = localView.tree().child(0).finiteElement();
      for( size_t i=0; i<localFE.size(); i++) {
        Operation::apply(vector_[localView.index(i)[0]], x);
      }
    }
};

template<class Basis, class Vector>
using VectorExchangeEqual = VectorExchange<Basis, Vector, Type::Equal>;

template<class Basis, class Vector>
using VectorExchangeAdd = VectorExchange<Basis, Vector, Type::Add>;

#endif
