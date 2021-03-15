// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FUFEM_SYMMETRICTENSOR_HH
#define DUNE_FUFEM_SYMMETRICTENSOR_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

/** \brief  A class implementing a 2nd order symmetric tensor.
 *
 *  A \f$ dim\times dim \f$ tensor is stored internally as a <tt>Dune::FieldVector<field_type, dim*(dim+1)/2></tt>
 *  For a symmetric <i>3x3</i>-Tensor \f$ E \f$ the components are assumed to be stored in the order
 *  \f$ [ E(1,1),\  E(2,2),\ E(3,3),\ E(1,2),\ E(1,3),\ E(2,3) ]\f$ and analogous for other dimensions 
 *  \tparam  dim number of lines/columns of tensor
 */
template <int dim, class field_type=double>
class SymmetricTensor : public Dune::FieldVector<field_type,dim*(dim+1)/2>
{

public:
    /** \brief dimension of Tensor
     *
     */
    enum {dimension = dim};

    enum {rows = dim};

    enum {cols = dim};

    /** \brief Default constructor
     *
     *  Tensor is initialized containing zeros if no argument is given.
     *  \param eye if true tensor is initialized as identity
     */
    SymmetricTensor(bool eye=false): Dune::FieldVector<field_type,dim*(dim+1)/2>(0.0)
    {
        if (eye)
            for (int i = 0; i < dim; ++i)
                (*this)[i] = 1.0;
    }

    /** \brief Constructor for tensor with all entries being the same
     *  \param a the value of all entries
     */
    SymmetricTensor(field_type a): Dune::FieldVector<field_type,dim*(dim+1)/2>(a) {}

    /** \brief Construct from FieldVector of correct dimension
     *  \param fvector the FieldVector of the internal representation of the symmetric tensor,
     *                  For 3x3 the order in the vector is expected to be \f$[ E(1,1) E(2,2) E(3,3) E(1,2) E(1,3) E(2,3) ]\f$.
     */
    SymmetricTensor(const Dune::FieldVector<field_type,dim*(dim+1)/2> fvector): Dune::FieldVector<field_type,dim*(dim+1)/2>(fvector) {}

    /** \brief This is the contraction product of two tensors of 2nd stage \f$ E*F := E:F = \sum_{ij}(E_{ij}\cdot F_{ij}) \f$ 
     *  \param B the <tt>SymmetricTensor<dim></tt> to contract with
     */
    field_type operator*(const SymmetricTensor<dim,field_type>& B) const
    {
        field_type r = 0.0;
        for (int i = 0; i < dim*(dim+1)/2; ++i)
            r += (i>=dim ? 2 : 1)*(*this)[i]*B[i];
        return r;
    }

    /** \brief Forwarding of <tt>Dune::template FieldVector<double,dim*(dim+1)/2>::operator=</tt>
     *
     */
    using Dune::template FieldVector<field_type,dim*(dim+1)/2>::operator=;
    /** \brief Forwarding of <tt>Dune::template FieldVector<double,dim*(dim+1)/2>::operator*=</tt>
     *
     */
    using Dune::template FieldVector<field_type,dim*(dim+1)/2>::operator*=;

    /** \brief Matrix-Vector product
     *
     *  This computes \f$ y+=Ax \f$.
     *  \param x the vector to multiply with
     *  \param y the vector \f$ Ax \f$ is added to
     */
    void umv(const Dune::FieldVector<field_type,dim>& x, Dune::FieldVector<field_type,dim>& y) const
    {
        for (int i = 0; i < dim; ++i)
        {
            y[i] += (*this)[i]*x[i];
            int j = i;
            while (++j < dim)
            {
                y[i] += (*this)[(i+1)*dim - i*(i+1)/2 + j-i-1]*x[j];
                y[j] += (*this)[(i+1)*dim - i*(i+1)/2 + j-i-1]*x[i];
            }
        }
    }

    /** \brief Matrix-Vector product
     *
     *  This computes \f$ y=Ax \f$.
     *  \param x the vector to multiply with
     *  \param y the result vector
     */
    void mv(const Dune::FieldVector<field_type,dim>& x, Dune::FieldVector<field_type,dim>& y) const
    {
        y = 0;
        umv(x,y);
    }

    /** \brief Matrix style random read/write access to components 
     *  \param i line index
     *  \param j column index
     */
    field_type& operator() (int i, int j)
    {
        if (i==j)
            return (*this)[i];
        else if (i<j)
            return (*this)[(i+1)*dim - i*(i+1)/2 + (j-i) - 1];
        else
            return (*this)[(j+1)*dim - j*(j+1)/2 + (i-j) - 1];
    }

    /** \brief Matrix style random read access to components
     *  \param i line index
     *  \param j column index
     */
    const field_type& operator() (int i, int j) const
    {
        if (i==j)
            return (*this)[i];
        else if (i<j)
            return (*this)[(i+1)*dim - i*(i+1)/2 + (j-i) - 1];
        else
            return (*this)[(j+1)*dim - j*(j+1)/2 + (i-j) - 1];
    }

    field_type trace() const
    {
      field_type ret = 0.0;
      for (int i = 0; i < dim; ++i)
        ret += (*this)[i];

      return ret;
    }

    /** \brief add constant value to the diagonal, leaving
     *   off-diagonal entries unaffected.
     */
    void addToDiag(field_type x) {
        for (int i = 0; i < dim; ++i)
          (*this)[i] += x;
    }

    /** \brief set diagonal to "constant" value and zero non-diagonal entries
     *  \param diagValue the value the diagonal entries are to be
     */
    void setDiag(field_type diagValue)
    {
        *this = 0.0;
        addToDiag(diagValue);
    }

    /** \brief set diagonal according to vector entries and zero non-diagonal entries
     *  \param diagVector the vector storing the diagonal entries
     */
    void setDiag(const Dune::FieldVector<field_type,dim> &diagVector)
    {
        *this = 0.0;
        for (int i = 0; i < dim; ++i)
            (*this)[i] = diagVector[i];
    }

    /** \brief Return the FieldMatrix representation of the symmetric tensor.*/
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

    /* A matrix-product has not been implemented */
};

/** \brief formatting tensor for standard output stream
 */
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

