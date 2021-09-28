// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef BOUNDARY_INDEX_BC_ASSEMBLER_HH
#define BOUNDARY_INDEX_BC_ASSEMBLER_HH

namespace Dune::Elastodynamics {

  template<class Basis>
  class BoundaryIndexBCAssembler {
  
    private:

      const Basis& basis_;
      const std::vector<int> boundaryIndex_;
      
      template<class MatrixType>
      void addMatrixBC(MatrixType& matrix) {
      
        auto gridView = basis_.gridView();
      
        for( const auto& element : elements(gridView)) {
        
          auto localView = basis_.localView();
          const auto &indexSet = gridView.indexSet();
          static const int dim = gridView.dimension;
          
		  for( const auto& isect : intersections(gridView, element)) {
            if(isect.boundary()) {
		       
		       auto ref = referenceElement<double, dim>(element.type());  
		       
		       switch(boundaryIndex_[isect.boundarySegmentIndex()]) {
		        case 0: // no boundary condition applied
		          break;
		        case 1: // set dirichlet BC
		          for(int i=0; i<ref.size(isect.indexInInside(), 1, dim); i++) {         
                    int row = indexSet.subIndex(element, ref.subEntity(isect.indexInInside(), 1, i, dim), dim);
                    FieldMatrix<double, dim, dim> I = ScaledIdentityMatrix<double, dim>(1.0);
                    FieldMatrix<double, dim, dim> O(0.0);
                    auto rowBegin = matrix[row].begin();
                    auto rowEnd = matrix[row].end();
                    for(; rowBegin!=rowEnd; rowBegin++) {
                      *rowBegin = (row==rowBegin.index()) ? I : O;
                    }
                  }
		          break;
		        case 2: // force on area
		          break;
		        case 3: // something else
		          break;
		        // ....
		      }
		    }
		  }
	    }
      }
      
      template<class VectorType, class Force>
      void addVectorBC(VectorType& vector, Force& force) {
      
        auto gridView = basis_.gridView();
      
        for( const auto& element : elements(gridView)) {
        
          auto localView = basis_.localView();
          const auto &indexSet = gridView.indexSet();
          static const int dim = gridView.dimension;
          
		  for( const auto& isect : intersections(gridView, element)) {
            if(isect.boundary()) {
		       
		       auto ref = referenceElement<double, dim>(element.type());  
		       
		       switch(boundaryIndex_[isect.boundarySegmentIndex()]) {
		        case 0: // no boundary condition applied
		          break;
		        case 1: // set fixed wall
		          for(int i=0; i<ref.size(isect.indexInInside(), 1, dim); i++) {         
                    int row = indexSet.subIndex(element, ref.subEntity(isect.indexInInside(), 1, i, dim), dim);
                    FieldVector<double, dim> O(0.0);
                    vector[row] = O;
                  }
		          break;
		        case 2: // force on area
		          for(int i=0; i<ref.size(isect.indexInInside(), 1, dim); i++) {         
                    int row = indexSet.subIndex(element, ref.subEntity(isect.indexInInside(), 1, i, dim), dim);
                    FieldVector<double, dim> O(0.0);
                    vector[row] = force;
                  }
		          break;
		        case 3: // something else
		          break;
		        // ....
		      }
		    }
		  }
	    }
      }
      
    public:
    
      BoundaryIndexBCAssembler(const Basis& basis, const std::vector<int> boundaryIndex)
        : basis_(basis)
        , boundaryIndex_(boundaryIndex)
      {}
      
      template<class MatrixType>
      void assembleMatrix(MatrixType& matrix) {
        addMatrixBC(matrix);
      }
      
      
      template<class VectorType, class Force>
      void assembleVector(VectorType& vector, Force& force) {
        addVectorBC(vector, force);
      }
    
  }; 
}     

#endif
