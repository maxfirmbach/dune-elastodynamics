// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*
This file is taken from the Book: DUNE-The Distributed and Unified
Numerics Environment by Oliver Sander and is modified!
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

template<class GridView, class Vector, class Operation>	
class VectorExchange : public Dune::CommDataHandleIF<VectorExchange<GridView, Vector, Operation>, typename Vector::value_type> {

  private:
      
    const GridView& gridView_;
    Vector& vector_;

  public:
      
    typedef typename Vector::value_type DataType;
    
    VectorExchange(const GridView& gridView, Vector& vector)
      : gridView_(gridView), vector_(vector)
    {}
    
    // we want the vertex data
    bool contains(int dim, int codim) const
    { return  (codim == dim); }
    
    bool fixedSize(int dim, int codim) const
    { return true; }
    
    template<class Entity>
    size_t size(const Entity& entity) const
    { return 1; }

    template<class MessageBuffer, class Entity>
    void gather(MessageBuffer& buffer, const Entity& entity) const
    {
      buffer.write(vector_[gridView_.indexSet().index(entity)]);
    }
    
    template<class MessageBuffer, class Entity>
    void scatter(MessageBuffer& buffer, const Entity& entity, size_t n) const
    {
      DataType x;
      buffer.read(x);
      Operation::apply(vector_[gridView_.indexSet().index(entity)], x);
    }
};

template<class GridView, class Vector>
using VectorExchangeEqual = VectorExchange<GridView, Vector, Type::Equal>;

template<class GridView, class Vector>
using VectorExchangeAdd = VectorExchange<GridView, Vector, Type::Add>;

#endif
