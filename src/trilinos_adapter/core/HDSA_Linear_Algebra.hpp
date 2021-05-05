#ifndef HDSA_LINEAR_ALGEBRA_HPP
#define HDSA_LINEAR_ALGEBRA_HPP

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziRTRSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialQRDenseSolver.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"

namespace HDSA
{

namespace Linear_Algebra
{
    
  template <class RealT>
  void Iterative_GEVP_Solver(std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > & evecs, std::vector<RealT> & evals, const HDSA::Ptr<HDSA::Linear_Operator<RealT> > & A, 
			     const HDSA::Ptr<HDSA::Linear_Operator<RealT> > & B, RealT tol, int numBlocks, int blocksize, std::string solver, bool is_joint = false)
  {      
    int num_eigs = evals.size();
    HDSA::Ptr<Anasazi::MultiVec<RealT> > anasazi_vec = HDSA::makePtr<HDSA_Anasazi_Vector<RealT> >(evecs[0],blocksize,is_joint);
    HDSA::Ptr<Anasazi::Operator<RealT> > A_anasazi = HDSA::makePtr<HDSA_Anasazi_Operator<RealT> >(A);
    HDSA::Ptr<Anasazi::Operator<RealT> > B_anasazi = HDSA::makePtr<HDSA_Anasazi_Operator<RealT> >(B);
    anasazi_vec->MvRandom();
    HDSA::Ptr<Anasazi::BasicEigenproblem<RealT, Anasazi::MultiVec<RealT>, Anasazi::Operator<RealT> > > problem =
      HDSA::makePtr<Anasazi::BasicEigenproblem<RealT, Anasazi::MultiVec<RealT>, Anasazi::Operator<RealT> > >(A_anasazi,B_anasazi,anasazi_vec);
    problem->setHermitian(true);
    problem->setNEV(num_eigs);
      
    bool boolret = problem->setProblem();
    if(boolret != 1)
      {
	std::cout << "Eigenvalue solver was not set propertly" << std::endl;
      }
    // Create parameter list to pass into the solver manager                                                                                                                    
    Teuchos::ParameterList MyPL;
    int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary + Anasazi::TimingDetails + Anasazi::IterationDetails;
    std::string which("LR");
    bool insitu = false;
    MyPL.set( "Verbosity", verbosity );
    MyPL.set( "Which", which );
    MyPL.set( "Block Size", blocksize );
    MyPL.set( "Num Blocks", numBlocks );
    MyPL.set( "Convergence Tolerance", tol );
    MyPL.set( "In Situ Restarting", insitu );
    MyPL.set( "Print Number of Ritz Values", blocksize*numBlocks );
    MyPL.set("Use Locking",true);
    // Create the solver manager  
    Anasazi::ReturnType returnCode = Anasazi::Converged;
    if(solver == "RTR")
      {
	Anasazi::RTRSolMgr<RealT,Anasazi::MultiVec<RealT>,Anasazi::Operator<RealT> > MySolverMgr(problem, MyPL);
	returnCode = MySolverMgr.solve();
      }
    else if(solver == "LOBPCG")
      {
	Anasazi::LOBPCGSolMgr<RealT,Anasazi::MultiVec<RealT>,Anasazi::Operator<RealT> > MySolverMgr(problem, MyPL);
	returnCode = MySolverMgr.solve();
      }
    else if(solver=="Block Davidson")
      {
	Anasazi::BlockDavidsonSolMgr<RealT,Anasazi::MultiVec<RealT>,Anasazi::Operator<RealT> > MySolverMgr(problem, MyPL);
	returnCode = MySolverMgr.solve();
      }
    if (returnCode != Anasazi::Converged)
      {
	std::cout << "Eigenvalue solver failed to converge" << std::endl;
      }
    Anasazi::Eigensolution<RealT,Anasazi::MultiVec<RealT> > sol = problem->getSolution();
    // Get the eigenvalues and eigenvectors from the eigenproblem  
    HDSA_Anasazi_Vector<RealT>* Anasazi_Evecs = dynamic_cast<HDSA_Anasazi_Vector<RealT>* >(&const_cast<Anasazi::MultiVec<RealT> &>(*sol.Evecs));
    for(int k = 0; k < num_eigs; k++)
      {
	evals[k] = sol.Evals[k].realpart;
	evecs[k]->set(*Anasazi_Evecs->vec[k]);
      } 
  }


  template <class RealT>
  // Solve the linear system A*x = b
  void Iterative_Linear_Solve(HDSA::Vector<RealT> & x, const HDSA::Vector<RealT> & b, const HDSA::Ptr<HDSA::Linear_Operator<RealT> > & A, 
			      RealT tol, std::string solver = "CG", bool verbose = false)
  {
    // Build the problem matrix
    HDSA::Ptr<HDSA_Belos_Operator<RealT> > A_Belos = HDSA::makePtr<HDSA_Belos_Operator<RealT> >(A);
    
    int frequency = 1;  // how often residuals are printed by solver
    int blocksize = 1;
    int numrhs = 1;
      
    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("blocksize",&blocksize,"Block size used by CG .");
           
    int maxits = b.dimension();  
    Teuchos::ParameterList belosList;
    belosList.set( "Block Size", blocksize );                // Blocksize to be used by iterative solver
    belosList.set( "Num Blocks", maxits );                   // Number of blocks
    belosList.set( "Maximum Iterations", maxits );           // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );           // Relative convergence tolerance requested
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
		     Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
      belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    
    HDSA::Ptr<HDSA_Belos_Vector<RealT> > soln = HDSA::makePtr<HDSA_Belos_Vector<RealT> >(b,numrhs);
    HDSA::Ptr<HDSA_Belos_Vector<RealT> > rhs = HDSA::makePtr<HDSA_Belos_Vector<RealT> >(b,numrhs);
    rhs->vec[0]->set(b);
    RealT rhs_norm = b.norm();
    if(rhs_norm != 0.0)
      {
	rhs->vec[0]->scale(1.0/rhs_norm);
	Belos::OperatorTraits<RealT,Belos::MultiVec<RealT>,Belos::Operator<RealT> >::Apply( *A_Belos, *rhs, *soln );
	
	HDSA::Ptr<Belos::LinearProblem<RealT,Belos::MultiVec<RealT>,Belos::Operator<RealT> > > problem =	
	  HDSA::makePtr<Belos::LinearProblem<RealT,Belos::MultiVec<RealT>,Belos::Operator<RealT> > >( A_Belos, soln, rhs );
	bool set = problem->setProblem();
	if (set == false) 
	  {
	    if (verbose)
	      {
		std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
	      }
	    verbose = true;
	  }
	
	HDSA::Ptr< Belos::SolverManager<RealT,Belos::MultiVec<RealT>,Belos::Operator<RealT>> > belos_solver;
	if(solver == "CG")
	  {
	    belos_solver = HDSA::makePtr<Belos::BlockCGSolMgr<RealT,Belos::MultiVec<RealT>,Belos::Operator<RealT> > >( problem, HDSA::makePtrFromRef(belosList) );
	  }
	else if(solver == "GMRES")
	  {
	    belos_solver = HDSA::makePtr<Belos::BlockGmresSolMgr<RealT,Belos::MultiVec<RealT>,Belos::Operator<RealT> > >( problem, HDSA::makePtrFromRef(belosList) );
	  }
	else
	  {
	    std::cout << "Error specifying the linear solver" << std::endl;
	  }

	Belos::ReturnType ret = belos_solver->solve();

	if(ret != Belos::Converged)
	  {
	    std::cout << "Belos solver did not converge for linear solve" << std::endl;
	  }

	x.set(*soln->vec[0]);
	x.scale(rhs_norm);

	// Test achievedTol output
	RealT ach_tol = belos_solver->achievedTol();
	if (verbose)
	  {
	    std::cout << "Achieved tol : "<<ach_tol<<std::endl;
	  }	
      }
    else
      {
	x.zero();
      }
  }

  template <class RealT>
  // Solve the linear system A*x = b
  void Iterative_Linear_Solve(HDSA::Ptr<HDSA::Vector<RealT> > & x, const HDSA::Ptr<HDSA::Vector<RealT> > & b, const HDSA::Ptr<HDSA::Linear_Operator<RealT> > & A, 
			      RealT tol, std::string solver = "CG", bool verbose = false)
  {
    HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(*x,*b,A,tol,solver,verbose);
  }

  // Solve upper triangular system R*x=b for upper triangluar matrix R
  template <class RealT>
  void Upper_Tri_Solve(const HDSA::Ptr<HDSA::Vector<RealT> > & x, const HDSA::Ptr<HDSA::Vector<RealT> > & b, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & R)
  {
    int n = R->numRows();
    for(int k = n-1; k>=0; k--)
      {
	RealT val = (*b)(k);
	for(int j = k+1; j < n; j++)
	  {
	    val -= (*x)(j)*(*R)(k,j);
	  }
	val = val/(*R)(k,k);
	x->Replace_Element(k,val);
      }
  }

  // Compute QR factorization of A, A=Q*R, with respect to Euclidean inner products, use R to request R
  template <class RealT>
  void QR_Factorization(const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & Q, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & R = HDSA::nullPtr)
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A_tmp = A->Clone();
    A_tmp->Copy_from_Teuchos_Matrix(A->Get_Teuchos_Matrix());
    Teuchos::SerialQRDenseSolver<int, RealT> QR_Solve;
    QR_Solve.setMatrix(A_tmp->Get_Teuchos_Matrix());
    QR_Solve.factorWithEquilibration(false);
    QR_Solve.factor();
    QR_Solve.formQ();
    HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > Q_teuchos = QR_Solve.getQ();
    Q->Copy_from_Teuchos_Matrix(Q_teuchos);

    if(R != HDSA::nullPtr)
      {
	QR_Solve.formR();
	HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > R_teuchos = QR_Solve.getR();
	R->Copy_from_Teuchos_Matrix(R_teuchos);
      }
  }

  // Compute Cholesky factorization of A, A=R^T*R
  template <class RealT>
  void Cholesky_Factorization(const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & R)
  {
    int n = A->numRows();
    Teuchos::SerialSpdDenseSolver<int, RealT> Chol_Solve;
    HDSA::Ptr<Teuchos::SerialSymDenseMatrix<int, RealT> > C = HDSA::makePtr<Teuchos::SerialSymDenseMatrix<int, RealT> >(n,n); 
    for(int i = 0; i < n; i++)
      {
	for(int j = 0; j < n; j++)
	  {
	    (*C)(i,j) = (*A)(i,j);
	  }
      }
    Chol_Solve.setMatrix(C);
    Chol_Solve.factor();
    HDSA::Ptr<Teuchos::SerialSymDenseMatrix<int, RealT> > Rc = Chol_Solve.getFactoredMatrix(); // R should be upper triangular. The R returned is symmetric, its upper half is what we need.
    for(int i = 0; i < n; i++)
      {
	for(int j = i; j < n; j++)
	  {
	    R->Replace_Element(i,j,(*Rc)(i,j));
	  }
      }
  }

  // Compute the SVD of A=U*S*V^T where U and V are orthogonal matrix and S is a diagaonal matrix, stored as a std::vector containing the diagonal entries
  template <class RealT>
  void SVD(const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & U, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & VT, const HDSA::Ptr<HDSA::Vector<RealT> > & S) 
  {
    int m = A->numRows();
    int n = A->numCols();
    Teuchos::LAPACK<int, RealT> lapack;
    char JOBU = 'S';
    char JOBVT = 'S';
    HDSA::Ptr<Teuchos::SerialDenseVector<int, RealT> > S_vec = HDSA::makePtr<Teuchos::SerialDenseVector<int, RealT> >(n);
    int LWORK = std::max(1,std::max(3*std::min(m,n)+std::max(m,n),5*std::min(m,n))) + 1;
    HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > WORK = HDSA::makePtr<Teuchos::SerialDenseMatrix<int, RealT> >(LWORK,1);
    HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > RWORK = HDSA::makePtr<Teuchos::SerialDenseMatrix<int, RealT> >(LWORK,1);
    int info;
    lapack.GESVD(JOBU,JOBVT,m,n,A->Get_Teuchos_Matrix()->values(),m,(*S_vec).values(),U->Get_Teuchos_Matrix()->values(),m,VT->Get_Teuchos_Matrix()->values(),n,(*WORK).values(),LWORK,(*RWORK).values(),&info);
    // Yields the decomposition A = U*diag(S)*VT, note VT is the transpose of V
    for(int k = 0; k < n; k++)
      {
	S->Replace_Element(k,(*S_vec)(k));
      }
  }

  // Compute eigenvalue decomposition A=V*S*V^T for a symmetric matrix A, store eigenvalues in a length n std::vector S instead of nxn diagonal matrix
  template <class RealT>
  void Symmetric_Eig_Decomposition(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & V, HDSA::Ptr<HDSA::Vector<RealT> > & S) 
  {
    int n = A->numRows();
    HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > B = HDSA::makePtr<Teuchos::SerialDenseMatrix<int, RealT> >(n,n);
    for(int i = 0; i < n; i++)
      {
	for(int j = 0; j < n; j++)
	  {
	    (*B)(i,j) = (*A)(i,j); 
	  }
      }
    HDSA::Ptr<Teuchos::SerialDenseVector<int, RealT> > S_rev = HDSA::makePtr<Teuchos::SerialDenseVector<int, RealT> >(n); // SYEV outputs eigenvalues if reverse order
    Teuchos::LAPACK<int, RealT> lapack;
    char JOBZ = 'V';
    char UPLO = 'U';
    HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > WORK = HDSA::makePtr<Teuchos::SerialDenseMatrix<int, RealT> >(3*n,3*n);
    int lwork = (*WORK).stride();
    int info;
    lapack.SYEV(JOBZ,UPLO,n,(*B).values(),n,(*S_rev).values(),(*WORK).values(),lwork,&info);
    for(int j = 0; j < n; j++)
      {
	S->Replace_Element(j,(*S_rev)(n-1-j));
	for(int i = 0; i < n; i++)
	  {
	    V->Replace_Element(i,j,(*B)(i,n-1-j));
	  }
      }
  }

  // Compute QR factorization A=Q*R with respect to weighted inner products defined by W with HDSA::Vector vec compatible to compute matrvecs with W,
  // use R to request R and invert_W to specify if weight matrix should be W^{-1}
  template <class RealT>
  void CholQR(const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & Q, const HDSA::Ptr<HDSA::Linear_Operator<RealT> > & W,
	      const HDSA::Ptr<HDSA::Vector<RealT> > & vec, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & WQ = HDSA::nullPtr, 
	      const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & R = HDSA::nullPtr, bool invert_W = false)
  {   
    int m = A->numRows();
    int n = A->numCols();

    // Take QR Decomposition of Y=Z*R
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Z = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,n);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_A;
    if(R != HDSA::nullPtr)
      {
	R_A = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
	HDSA::Linear_Algebra::QR_Factorization<RealT>(A,Z,R_A);
      }
    else
      {
	HDSA::Linear_Algebra::QR_Factorization<RealT>(A,Z);
      }

    // Compute X=W*Z
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > X = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,n);
    HDSA::Ptr<HDSA::Vector<RealT> > vec_in = vec->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > vec_out = vec->Clone();
    for(int j = 0; j < n; j++)
      {
	Z->Write_Column_to_Vector(j,vec_in);
	if(invert_W)
	  {
	    Iterative_Linear_Solve(vec_out, vec_in, W, 1.e-12);
	  }
	else
	  {
	    W->matvec(vec_out,vec_in);
	  }
	X->Write_Vector_to_Column(j,vec_out);
      }

    // Compute C = Z^T*W*Z = Z^T*X
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
    Z->Multiply(C,X,true,false);

    // Compute R_C=chol(C)
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_C = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
    HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(C,R_C);

    // Compute Q=YR_C^{-1} in line 4 of algorithm 4
    
    // Need to solve the linear systems R_C*x_i = e_i for i=1,2,...,n. Here e_i is the ith standard basis vector in euclidean space, hence its ith entry is 1 and others are 0
    // Need to multiply Z*x_i and store this as the ith column of Q
    for(int i = 0; i < n; i++)
      {
	// Solve the linear system R_C*x = e_i
	HDSA::Ptr<HDSA::Vector<RealT> > ei = HDSA::makePtr<Std_Vector<RealT> >(n);
	ei->basis(i);
	HDSA::Ptr<HDSA::Vector<RealT> > x = HDSA::makePtr<Std_Vector<RealT> >(n);
	HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>( x, ei, R_C);
	
	// Set the ith column of Q to Zx
	for(int k = 0; k < m; k++)
	  {
	    RealT val = 0.0;
	    for(int j = 0; j < i+1; j++)
	      {
		val += (*Z)(k,j)*(*x)(j);
	      }
	    Q->Replace_Element(k,i,val);
	  }

	if(WQ != HDSA::nullPtr)
	  {
	    // Set the ith column of WQ to Xx
	    for(int k = 0; k < m; k++)
	      {
		RealT val = 0.0;
		for(int j = 0; j < i+1; j++)
		  {
		    val += (*X)(k,j)*(*x)(j);
		  }
		WQ->Replace_Element(k,i,val);
	      }
	  }
      }
    
    if(R != HDSA::nullPtr)
      {
	R->Get_Teuchos_Matrix()->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*R_C->Get_Teuchos_Matrix(),*R_A->Get_Teuchos_Matrix(),0.0);
      }
  }
 
  // Compute QR factorization A=Q*R with respect to weighted inner products defined by W with A*W precomputed
  template <class RealT>
  void CholQR_Pre_W(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & WA, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & Q,
		    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & R, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & WQ = HDSA::nullPtr)
  {   
    int m = A->numRows();
    int n = A->numCols();

    // Compute C = A^T*W*A
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
    A->Multiply(C,WA,true,false);
    
    // Compute R=chol(C)
    HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(C,R);

    // Compute Q=YR^{-1} in line 4 of algorithm 4
    
    // Need to solve the linear systems R*x_i = e_i for i=1,2,...,n. Here e_i is the ith standard basis vector in euclidean space, hence its ith entry is 1 and others are 0
    // Need to multiply A*x_i and store this as the ith column of Q
    for(int i = 0; i < n; i++)
      {
	// Solve the linear system R*x = e_i
	HDSA::Ptr<HDSA::Vector<RealT> > ei = HDSA::makePtr<Std_Vector<RealT> >(n);
	ei->basis(i);
	HDSA::Ptr<HDSA::Vector<RealT> > x = HDSA::makePtr<Std_Vector<RealT> >(n);
	HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>( x, ei, R);
	
	// Set the ith column of Q to Ax
	for(int k = 0; k < m; k++)
	  {
	    RealT val = 0.0;
	    for(int j = 0; j < i+1; j++)
	      {
		val += (*A)(k,j)*(*x)(j);
	      }
	    Q->Replace_Element(k,i,val);
	  }

	if(WQ != HDSA::nullPtr)
	  {
	    // Set the ith column of WQ to WAx
	    for(int k = 0; k < m; k++)
	      {
		RealT val = 0.0;
		for(int j = 0; j < i+1; j++)
		  {
		    val += (*WA)(k,j)*(*x)(j);
		  }
		WQ->Replace_Element(k,i,val);
	      }
	  }
      }
    
  }

 
}

}

#endif

