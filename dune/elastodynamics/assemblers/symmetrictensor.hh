// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*
This file is taken from the dune-fufem module!
Title: symmetrictensor.hh
Authors: dune-fufem team dune-fufem@lists.spline.inf.fu-berlin.de
Date: 2021
Version: -
Availability: https://git.imp.fu-berlin.de/agnumpde/dune-fufem/-/blob/master/dune/fufem/symmetrictensor.hh
*/

#ifndef SYMMETRICTENSOR_HH
#define SYMMETRICTENSOR_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

template <int dim, class field_type=double>
class SymmetricTensor : public Dune::FieldVector<field_type,dim*(dim+1)/2> {

  public:

    enum {dimension = dim};
    enum {rows = dim};
    enum {cols = dim};

    SymmetricTensor(bool eye=false): Dune::FieldVector<field_type,dim*(dim+1)/2>(0.0) {
      
      if (eye) {
        for (int i = 0; i < dim; ++i) {
          (*this)[i] = 1.0;
        }      
      }
    }

    SymmetricTensor(field_type a): Dune::FieldVector<field_type,dim*(dim+1)/2>(a) {}

    SymmetricTensor(const Dune::FieldVector<field_type,dim*(dim+1)/2> fvector): Dune::FieldVector<field_type,dim*(dim+1)/2>(fvector) {}

    field_type operator*(const SymmetricTensor<dim,field_type>& B) const {
      
      field_type r = 0.0;
      for (int i = 0; i < dim*(dim+1)/2; ++i) {
        r += (i>=dim ? 2 : 1)*(*this)[i]*B[i];
      }
      return r;
    }

    using Dune::template FieldVector<field_type,dim*(dim+1)/2>::operator=;
    using Dune::template FieldVector<field_type,dim*(dim+1)/2>::operator*=;

    void umv(const Dune::FieldVector<field_type,dim>& x, Dune::FieldVector<field_type,dim>& y) const {
   
      for (int i = 0; i < dim; ++i) {
        y[i] += (*this)[i]*x[i];
        int j = i;
        while (++j < dim) {
          y[i] += (*this)[(i+1)*dim - i*(i+1)/2 + j-i-1]*x[j];
          y[j] += (*this)[(i+1)*dim - i*(i+1)/2 + j-i-1]*x[i];
        }
      }
    }

    void mv(const Dune::FieldVector<field_type,dim>& x, Dune::FieldVector<field_type,dim>& y) const {

      y = 0;
      umv(x,y);
    }

    field_type& operator() (int i, int j) {
      
      if (i==j)
        return (*this)[i];
      else if (i<j)
        return (*this)[(i+1)*dim - i*(i+1)/2 + (j-i) - 1];
      else
        return (*this)[(j+1)*dim - j*(j+1)/2 + (i-j) - 1];
    }

    const field_type& operator() (int i, int j) const {
      
      if (i==j)
        return (*this)[i];
      else if (i<j)
        return (*this)[(i+1)*dim - i*(i+1)/2 + (j-i) - 1];
      else
        return (*this)[(j+1)*dim - j*(j+1)/2 + (i-j) - 1];
    }

    field_type trace() const {

      field_type ret = 0.0;
      for (int i = 0; i < dim; ++i)
        ret += (*this)[i];

      return ret;
    }

    void addToDiag(field_type x) {
    
      for (int i = 0; i < dim; ++i)
        (*this)[i] += x;
    }

    void setDiag(field_type diagValue) {

      *this = 0.0;
      addToDiag(diagValue);
    }

    void setDiag(const Dune::FieldVector<field_type,dim> &diagVector) {

      *this = 0.0;
      for (int i = 0; i < dim; ++i)
        (*this)[i] = diagVector[i];
    }

    Dune::FieldMatrix<field_type,dim,dim> matrix() const {

      Dune::FieldMatrix<field_type,dim,dim> mat;
      for (int i=0; i<dim; i++) {
        mat[i][i] = (*this)[i];
        for (int j=0; j<i; j++) {
          mat[i][j] = (*this)[(j+1)*dim - j*(j+1)/2 + (i-j) - 1];
          mat[j][i] = mat[i][j];
        }
      }
      return mat;
    }
};

template <int dim, class field_type>
std::ostream& operator<< (std::ostream& s, const SymmetricTensor<dim,field_type>& E)
{
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
      s << (j>0 ? " " : "") << E(i,j);
    s << std::endl;
  }
  return s;
}

namespace Dune {

  template< class K, int dim >
  struct FieldTraits< SymmetricTensor<dim, K> >
  {
    using field_type = field_t<K>;
    using real_type = real_t<K>;
  };
}

#endif
