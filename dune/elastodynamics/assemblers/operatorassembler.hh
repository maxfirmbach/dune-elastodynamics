// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef OPERATOR_ASSEMBLER_HH
#define OPERATOR_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>

#include <omp.h>

namespace Dune::Elastodynamics {

  template <class Basis>
  class OperatorAssembler {
  
    private:
	  
	  const Basis& basis_;
		
	  void addIndices(MatrixIndexSet& occupationPattern) {
	    
	    auto gridView  = basis_.gridView();
        auto localView = basis_.localView();
			  
        for( const auto& element : elements(gridView, Dune::Partitions::all)) {	
  				
  	      localView.bind(element);
   					
    	  for( size_t i=0; i<localView.size(); i++) {
      	    auto row = localView.index(i);
      		for (size_t j=0; j<localView.size(); j++) {
              auto col = localView.index(j);
        	  occupationPattern.add(row[0], col[0]);
      		}
    	  }
        }
      }
	
	  template <class LocalAssemblerType, class GlobalMatrixType>
	  void addEntries(LocalAssemblerType& localAssembler, GlobalMatrixType& A, bool lumping) {

        auto gridView  = basis_.gridView();
       
        typedef typename LocalAssemblerType::LocalMatrix LocalMatrix;
        
        // here we could technically parallelize
        for( const auto& element : elements(gridView, Dune::Partitions::all)) {
          
          auto localView = basis_.localView();
          localView.bind(element);
		
          LocalMatrix localMatrix;
          localAssembler.assemble(localMatrix, localView);

          for( size_t i=0; i<localMatrix.N(); i++) {
            auto row = localView.index(i);
            if(lumping)
              A[row[0]][row[0]][row[1]][row[1]] += localMatrix[i][i];
            else {
              for( size_t j=0; j<localMatrix.M(); j++) {
                auto col = localView.index(j);
                A[row[0]][col[0]][row[1]][col[1]] += localMatrix[i][j];
              }
            }
          }
        }
	  }

    public:
		
	  OperatorAssembler(const Basis& basis) 
	    : basis_(basis)
	  {}


      template <class GlobalMatrixType>
      void initialize(GlobalMatrixType& A)  
      {     
        Dune::MatrixIndexSet occupationPattern(basis_.size(), basis_.size());
		addIndices(occupationPattern);
		occupationPattern.exportIdx(A);
      }

	  template <class LocalAssemblerType, class GlobalMatrixType>
	  void assemble(LocalAssemblerType& localAssembler, GlobalMatrixType& A, bool lumping)
	  {			
		A = 0.0;
		addEntries(localAssembler, A, lumping);
	  }

  };
}

#endif
