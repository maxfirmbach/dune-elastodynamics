// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef COEFFICIENTS_HH
#define COEFFICIENTS_HH

#include <math.h>

#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>

using CoeffMatrix = Dune::Matrix<Dune::FieldMatrix<double, 1, 1>>;
using CoeffVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;

namespace Dune {

  // Runge-Kutta-Nyström coefficients
  // --------------------------------
  class RKNCoefficients {
  
    private:
        
        int stages_, order_;
        
        CoeffMatrix A_;
		CoeffVector b_, b_bar_, c_;

    public:
    
        RKNCoefficients(int stages, int order, CoeffMatrix& A,
                        CoeffVector& b, CoeffVector& b_bar, CoeffVector& c)
        : stages_(stages)
        , order_(order)
        , A_(A)
        , b_(b)
        , b_bar_(b_bar)
        , c_(c)
        {}
                        
        int order()  { return order_; }
        int stages() { return stages_; }
        
        CoeffMatrix A()     { return A_; }
		CoeffVector b()     { return b_; }
		CoeffVector b_bar() { return b_bar_; }
		CoeffVector c()     { return c_; }
		
  };
  
  // Runge-Kutta-Nyström method of order 4
  // -------------------------------------
  RKNCoefficients RKN4()
  {
	
    static CoeffMatrix A(3, 3);
    A[0][0] = 0.0; 	   A[0][1] = 0.0; A[0][2] = 0.0;
	A[1][0] = 1.0/8.0; A[1][1] = 0.0; A[1][2] = 0.0;
	A[2][0] = 0.0;	   A[2][1] = 0.5; A[2][2] = 0.0;
		
	static CoeffVector b_bar(3);	
	b_bar[0] = 1.0/6.0; b_bar[1] = 1.0/3.0; b_bar[2] = 0.0;
	
	static CoeffVector b(3);
	b[0] = 1.0/6.0; b[1] = 4.0/6.0; b[2] = 1.0/6.0;
					
    static CoeffVector c(3);					   	      		
	c[0] = 0.0; c[1] = 0.5; c[2] = 1.0;
	  
	return RKNCoefficients(3, 4, A, b, b_bar, c);			
  }

  // Runge-Kutta-Nyström method of order 5
  // ------------------------------------- 
  RKNCoefficients RKN5()
  {
    	
    static CoeffMatrix A(4,4);
	A[0][0] = 0.0; 	     A[0][1] = 0.0; 	  A[0][2] = 0.0; 	  A[0][3] = 0.0;
	A[1][0] = 1.0/50.0;  A[1][1] = 0.0; 	  A[1][2] = 0.0; 	  A[1][3] = 0.0;
	A[2][0] = -1.0/27.0; A[2][1] = 7.0/27.0;  A[2][2] = 0.0; 	  A[2][3] = 0.0;
	A[3][0] = 3.0/10.0;  A[3][1] = -2.0/35.0; A[3][2] = 9.0/35.0; A[3][3] = 0.0;

    static CoeffVector b_bar(4);
	b_bar[0] = 14.0/336.0; b_bar[1] = 100.0/336.0; b_bar[2] = 54.0/336.0; b_bar[3] = 0.0;
			
	static CoeffVector b(4);
	b[0] = 14.0/336.0; b[1] = 125.0/336.0; b[2] = 162.0/336.0; b[3] = 35.0/336.0;
			
	static CoeffVector c(4);
	c[0] = 0.0; c[1] = 1.0/5.0; c[2] = 2.0/3.0; c[3] = 1.0;
			
	return RKNCoefficients(4, 5, A, b, b_bar, c);
  }
  
  
  // Embedded Runge-Kutta-Nyström coefficients
  // -----------------------------------------
  class EmbeddedRKNCoefficients {
  
    private:
        
        int stages_, order_;
        
        CoeffMatrix A_;
		CoeffVector b_, b_bar_, b_tilde_, b_bar_tilde_, c_;

    public:
    
        EmbeddedRKNCoefficients(int stages, int order, CoeffMatrix& A,
                                CoeffVector& b, CoeffVector& b_bar, 
                                CoeffVector& b_tilde, CoeffVector& b_bar_tilde,
                                CoeffVector& c)
        : stages_(stages)
        , order_(order)
        , A_(A)
        , b_(b)
        , b_bar_(b_bar)
        , b_tilde_(b_tilde)
        , b_bar_tilde_(b_bar_tilde)
        , c_(c)
        {}
                        
        
        int order()  { return order_; }
        int stages() { return stages_; }
        
        CoeffMatrix A()           { return A_; }
		CoeffVector b()           { return b_; }
		CoeffVector b_bar()       { return b_bar_; }
		CoeffVector b_tilde()     { return b_tilde_; }
		CoeffVector b_bar_tilde() { return b_bar_tilde_; }
		CoeffVector c()           { return c_; }
		
  };
  
  // Emmbeded Runge-Kutta-Nyström method of order 4 5
  // Bettis coefficients
  // ------------------------------------------------
  EmbeddedRKNCoefficients BettisRKN45()
  {
    static CoeffMatrix A(6,6);
	A[0][0] = 0.0; 	      A[0][1] = 0.0; 	   A[0][2] = 0.0; 	   A[0][3] = 0.0;       A[0][4] = 0.0;      A[0][5] = 0.0;
	A[1][0] = 1.0/128.0;  A[1][1] = 0.0; 	   A[1][2] = 0.0; 	   A[1][3] = 0.0;       A[1][4] = 0.0;      A[1][5] = 0.0;
	A[2][0] = 1.0/96.0;   A[2][1] = 1.0/48.0;  A[2][2] = 0.0; 	   A[2][3] = 0.0;       A[2][4] = 0.0;      A[2][5] = 0.0;
	A[3][0] = 1.0/24.0;   A[3][1] = 0.0;       A[3][2] = 1.0/12.0; A[3][3] = 0.0;       A[3][4] = 0.0;      A[3][5] = 0.0;
	A[4][0] = 9.0/128.0;  A[4][1] = 0.0;       A[4][2] = 9.0/64.0; A[4][3] = 9.0/128.0; A[4][4] = 0.0;      A[4][5] = 0.0;
    A[5][0] = 7.0/90.0;   A[5][1] = 0.0;       A[5][2] = 4.0/15.0; A[5][3] = 1.0/15.0;  A[5][4] = 4.0/45.0; A[5][5] = 0.0;

    static CoeffVector b_bar_tilde(6);
    b_bar_tilde[0] = 7.0/90.0; b_bar_tilde[1] = 0.0; b_bar_tilde[2] = 4.0/15.0;
    b_bar_tilde[3] = 1.0/15.0; b_bar_tilde[4] = 4.0/45.0; b_bar_tilde[5] = 0.0;   
            
    static CoeffVector b_bar(6);
    b_bar[0] = 1.0/6.0; b_bar[1] = 0.0; b_bar[2] = 0.0; b_bar[3] = 1.0/3.0; b_bar[4] = 0.0; b_bar[5] = 0.0;
			
    static CoeffVector b_tilde(6);
    b_tilde[0] = 7.0/90.0; b_tilde[1] = 0.0; b_tilde[2] = 16.0/45.0;
    b_tilde[3] = 2.0/15.0; b_tilde[4] = 16.0/45.0; b_tilde[5] = 7.0/90.0;

    static CoeffVector b(6);			
    b[0] = 0.0; b[1] = 0.0; b[2] = 2.0/3.0; b[3] = -1.0/3.0; b[4] = 2.0/3.0; b[5] = 0.0;
	
	static CoeffVector c(6);
	c[0] = 0.0; c[1] = 1.0/8.0; c[2] = 1.0/4.0; c[3] = 1.0/2.0; c[4] = 3.0/4.0, c[5] = 1.0;
			
    return EmbeddedRKNCoefficients(6, 4, A, b, b_bar, b_tilde, b_bar_tilde, c);
  }
  
  // Emmbeded Runge-Kutta-Nyström method of order 6 5
  // Domain-Prince coefficients
  // ------------------------------------------------
  EmbeddedRKNCoefficients DPRKN64()
  {
    static CoeffMatrix A(6,6);
    A[0][0] = 0.0; 	              A[0][1] = 0.0; 	           A[0][2] = 0.0; 	           A[0][3] = 0.0;              A[0][4] = 0.0;             A[0][5] = 0.0;
	A[1][0] = 1.0/200.0;          A[1][1] = 0.0; 	           A[1][2] = 0.0; 	           A[1][3] = 0.0;              A[1][4] = 0.0;             A[1][5] = 0.0;
	A[2][0] = -1.0/2200.0;        A[2][1] = 1.0/22.0;          A[2][2] = 0.0; 	           A[2][3] = 0.0;              A[2][4] = 0.0;             A[2][5] = 0.0;
	A[3][0] = 637.0/6600.0;       A[3][1] = -7.0/110.0;        A[3][2] = 7.0/33.0;         A[3][3] = 0.0;              A[3][4] = 0.0;             A[3][5] = 0.0;
	A[4][0] = 225437.0/1968750.0; A[4][1] = -30073.0/281250.0; A[4][2] = 65569.0/281250.0; A[4][3] = -9367.0/984375.0; A[4][4] = 0.0;             A[4][5] = 0.0;
    A[5][0] = 151.0/2142.0;       A[5][1] = 5.0/116.0;         A[5][2] = 385.0/1368.0;     A[5][3] = 55.0/168.0;       A[5][4] = -6250.0/28101.0; A[5][5] = 0.0;
 
    static CoeffVector b_bar_tilde(6);
    b_bar_tilde[0] = 151.0/2142.0; b_bar_tilde[1] = 5.0/116.0; b_bar_tilde[2] = 385.0/1368.0;
    b_bar_tilde[3] = 55.0/168.0; b_bar_tilde[4] = -6250.0/28101.0; b_bar_tilde[5] = 0.0;
    
    static CoeffVector b_bar(6);
    b_bar[0] = 1349.0/157500.0; b_bar[1] = 7873.0/50000.0; b_bar[2] = 192199.0/900000.0; b_bar[3] = 521683.0/2100000.0; b_bar[4] = -16.0/125.0; b_bar[5] = 0.0;
    
    static CoeffVector b_tilde(6);
    b_tilde[0] = 151.0/2142.0; b_tilde[1] = 25.0/522.0; b_tilde[2] = 275.0/684.0;
    b_tilde[3] = 275.0/252.0; b_tilde[4] = -78125.0/112404.0; b_tilde[5] = 1.0/12.0;
    
    static CoeffVector b(6);
    b[0] = 1349.0/157500.0; b[1] = 7873.0/45000.0; b[2] = 27457.0/90000.0; b[3] = 521683.0/630000.0; b[4] = -2.0/5.0; b[5] = 1.0/12.0;

    static CoeffVector c(6);
    c[0] = 0.0; c[1] = 1.0/10.0; c[2] = 3.0/10.0; c[3] = 7.0/10.0; c[4] = 17.0/25.0; c[5] = 1.0;
  
    return EmbeddedRKNCoefficients(6, 4, A, b, b_bar, b_tilde, b_bar_tilde, c);
  }
  
  EmbeddedRKNCoefficients DPRKN86()
  {
    static CoeffMatrix A(9,9);
    A[0][0] = 0.0; 	                     A[0][1] = 0.0; 	               A[0][2] = 0.0; 	                   A[0][3] = 0.0;                    A[0][4] = 0.0;             A[0][5] = 0.0;            A[0][6] = 0.0;            A[0][7] = 0.0; A[0][8] = 0.0;
	A[1][0] = 1.0/800.0;                 A[1][1] = 0.0; 	               A[1][2] = 0.0; 	                   A[1][3] = 0.0;                    A[1][4] = 0.0;             A[1][5] = 0.0;            A[1][6] = 0.0;            A[1][7] = 0.0; A[1][8] = 0.0;
	A[2][0] = 1.0/600.0;                 A[2][1] = 1.0/300.0;              A[2][2] = 0.0; 	                   A[2][3] = 0.0;                    A[2][4] = 0.0;             A[2][5] = 0.0;            A[2][6] = 0.0;            A[2][7] = 0.0; A[2][8] = 0.0;
	A[3][0] = 9.0/200.0;                 A[3][1] = -9.0/100.0;             A[3][2] = 9.0/100.0;                A[3][3] = 0.0;                    A[3][4] = 0.0;             A[3][5] = 0.0;            A[3][6] = 0.0;            A[3][7] = 0.0; A[3][8] = 0.0;
	A[4][0] = -66701.0/197352.0;         A[4][1] = 28325.0/32892.0;        A[4][2] = -2665.0/5482.0;           A[4][3] = 2170.0/24669.0;         A[4][4] = 0.0;             A[4][5] = 0.0;            A[4][6] = 0.0;            A[4][7] = 0.0; A[4][8] = 0.0;
    A[5][0] = 227015747.0/304251000.0;   A[5][1] = -54897451.0/30425100.0; A[5][2] = 12942349.0/10141700.0;    A[5][3] = -9499.0/304251.0;         A[5][4] = 539.0/9250.0;    A[5][5] = 0.0;            A[5][6] = 0.0;            A[5][7] = 0.0; A[5][8] = 0.0;
    A[6][0] = -1131891597.0/901789000.0; A[6][1] = 41964921.0/12882700.0;  A[6][2] = -6663147.0/3320675.0;     A[6][3] = 270954.0/644135.0;      A[6][4] = -108.0/5875.0;   A[6][5] = 114.0/1645.0;   A[6][6] = 0.0;            A[6][7] = 0.0; A[6][8] = 0.0;
    A[7][0] = 13836959.0/3667458.0;      A[7][1] = -17731450.0/1833729.0;  A[7][2] = 1063919505.0/156478208.0; A[7][3] = -33213845.0/39119552.0; A[7][4] = 13335.0/28544.0; A[7][5] = -705.0/14272.0; A[7][6] = 1645.0/57088.0; A[7][7] = 0.0; A[7][8] = 0.0;
    A[8][0] = 223.0/7938.0;              A[8][1] = 0.0;                    A[8][2] = 1175.0/8064.0;            A[8][3] = 925.0/6048.0;           A[8][4] = 41.0/448.0;      A[8][5] = 925.0/14112.0;  A[8][6] = 1175.0/72576.0; A[8][7] = 0.0; A[8][8] = 0.0;
 
 
    static CoeffVector b_bar_tilde(9);
    b_bar_tilde[0] = 223.0/7938.0; b_bar_tilde[1] = 0.0; b_bar_tilde[2] = 1175.0/8064.0; b_bar_tilde[3] = 925.0/6048.0; b_bar_tilde[4] = 41.0/448.0;
    b_bar_tilde[5] = 925.0/14112.0; b_bar_tilde[6] = 1175.0/72576.0; b_bar_tilde[7] = 0.0; b_bar_tilde[8] = 0.0;
    
    static CoeffVector b_bar(9);
    b_bar[0] = 7987313.0/109941300.0; b_bar[1] = 0.0; b_bar[2] = 1610737.0/44674560.0; b_bar[3] = 10023263.0/33505920.0; 
    b_bar[4] = -497221.0/12409600.0; b_bar[5] = 10023263.0/78180480.0; b_bar[6] = 1610737.0/402071040.0; b_bar[7] = 0.0; b_bar[8] = 0.0;
    
    static CoeffVector b_tilde(9);
    b_tilde[0] = 223.0/7938.0; b_tilde[1] = 0.0; b_tilde[2] = 5875.0/36288.0; b_tilde[3] = 4625.0/21168.0; b_tilde[4] = 41.0/224.0; 
    b_tilde[5] = 4625.0/21168.0; b_tilde[6] = 5875.0/36288.0; b_tilde[7] = 223.0/7938.0; b_tilde[8] = 0.0;
    
    static CoeffVector b(9);
    b[0] = 7987313.0/109941300.0; b[1] = 0.0; b[2] = 1610737.0/40207104.0; b[3] = 10023263.0/23454144.0; b[4] = -497221.0/6204800.0; 
    b[5] = 10023263.0/23454144.0; b[6] = 1610737.0/402071040.0; b[7] = -4251941.0/54970650.0; b[8] = 3.0/20.0;

    static CoeffVector c(9);
    c[0] = 0.0; c[1] = 1.0/20.0; c[2] = 1.0/10.0; c[3] = 3.0/10.0; c[4] = 1.0/2.0; c[5] = 7.0/10.0; c[6] = 9.0/10.0; c[7] = 1.0; c[8] = 1.0;
  
    return EmbeddedRKNCoefficients(9, 6, A, b, b_bar, b_tilde, b_bar_tilde, c);
  
  }
  
    
  // Newmark coefficients
  // --------------------
  class NewmarkCoefficients {
  
    private:
    
      double beta_, gamma_;
  
    public: 
    
      NewmarkCoefficients(double beta, double gamma)
      : beta_(beta)
      , gamma_(gamma)
      {}
  
      double beta()  { return beta_; }
      double gamma() { return gamma_; }
      
  };
  
  // Störmer coefficients
  // --------------------
  NewmarkCoefficients Stoermer()
  {
    static const double beta  = 0.0;
	static const double gamma = 1.0/2.0;
	  
	return NewmarkCoefficients(beta, gamma);
  }	
	
  // Fox Goodwin coefficients
  // ------------------------
  NewmarkCoefficients FoxGoodwin()
  {
    static const double beta  = 1.0/12.0;
    static const double gamma = 1.0/2.0;
	  
    return NewmarkCoefficients(beta, gamma);
  }
	
  // Linear Acceleration coefficients
  // --------------------------------
  NewmarkCoefficients LinearAcceleration()
  {
    static const double beta  = 1.0/6.0;
    static const double gamma = 1.0/2.0;
	  
    return NewmarkCoefficients(beta, gamma);
  }
		
  // Constant Acceleration coefficients
  // ----------------------------------   
  NewmarkCoefficients ConstantAcceleration()
  {
    static const double beta  = 1.0/4.0;
    static const double gamma = 1.0/2.0;
	  
    return NewmarkCoefficients(beta, gamma);
  }
  
}

#endif
