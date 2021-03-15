// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef OPERATOR_ASSEMBLER_HH
#define OPERATOR_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>

namespace Dune::Elastodynamics {

  template <class Basis>
  class OperatorAssembler {
  
    private:
	  
	  const Basis& basis_;
			
	  void addIndices(MatrixIndexSet& occupationPattern) {
	    
	    auto gridView      = basis_.gridView();
        auto localView     = basis_.localView();
        auto localIndexSet = basis_.localIndexSet();
			  
        for( const auto& element : elements(gridView)) {	
  				
  	      localView.bind(element);
   	      localIndexSet.bind(localView);
   					
    	  for( size_t i=0; i<localIndexSet.size(); i++) {
      	    auto row = localIndexSet.index(i);
      		for (size_t j=0; j<localIndexSet.size(); j++) {
              auto col = localIndexSet.index(j);
        	  occupationPattern.add(row[0], col[0]);
      		}
    	  }
        }
      }
	
	  template <class LocalAssemblerType, class GlobalMatrixType>
	  void addEntries(LocalAssemblerType& localAssembler, GlobalMatrixType& A, bool lumping) {
      
        auto gridView      = basis_.gridView();
        auto localView     = basis_.localView();
        auto localIndexSet = basis_.localIndexSet();
        
        typedef typename LocalAssemblerType::LocalMatrix LocalMatrix;
        
		for( const auto& element : elements(gridView)) {
		
		  localView.bind(element);
		  localIndexSet.bind(localView);
		
		  LocalMatrix localMatrix;
		  localAssembler.assemble(localMatrix, localView);

    	  for( size_t i=0; i<localMatrix.N(); i++) {
      	    auto row = localIndexSet.index(i);
      	    if(lumping)
      	      A[row[0]][row[0]][row[1]][row[1]] += localMatrix[i][i];
      	    else {
      		  for( size_t j=0; j<localMatrix.M(); j++) {
         	    auto col = localIndexSet.index(j);
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
