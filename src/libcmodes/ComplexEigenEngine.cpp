#include "slepceps.h"
#include "slepcst.h"
#include <complex>
#include <string>
#include <iostream>
#include "ComplexEigenEngine.h"



int NCPA::ComplexEigenEngine::doESSCalculation( std::complex<double> *diag, int Nz_grid, double dz, double tol, 
	int *nev, double *k_min, double *k_max, PetscInt *nconv, std::complex<double> *k2, std::complex<double> **v ) {

	//
	// Declarations related to Slepc computations
	//
	Mat            A;           // problem matrix
	EPS            eps;         // eigenproblem solver context
	ST             stx;
	// KSP            kspx;
	// PC             pcx;
	EPSType        type;        // CHH 191022: removed const qualifier
	PetscReal      re, im;
	PetscScalar    kr, ki, *xr_, sigma;
	Vec            xr, xi;
	PetscInt       Istart, Iend, col[3], its, maxit;
	PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
	PetscScalar    value[3];	
	PetscErrorCode ierr;
	PetscMPIInt    rank, size;
	int i, j;
	double h2 = dz * dz;
	sigma = std::pow( (0.5*(*k_min + *k_max)), 2);

	// Initialize Slepc
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);  

	// Create the matrix A to use in the eigensystem problem: Ak=kx
	ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Nz_grid,Nz_grid); CHKERRQ(ierr);
	ierr = MatSetFromOptions(A); CHKERRQ(ierr);

	// the following Preallocation call is needed in PETSc version 3.3
	ierr = MatSeqAIJSetPreallocation(A, 3, PETSC_NULL); CHKERRQ(ierr);
	// or use: ierr = MatSetUp(A); 

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Compute the operator matrix that defines the eigensystem, Ax=kx
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	// Make matrix A 
	ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
	if (Istart==0) 
		FirstBlock=PETSC_TRUE;
	if (Iend==Nz_grid) 		/* @todo check if should be Nz_grid-1 */
		LastBlock=PETSC_TRUE;    

	value[0]=1.0/h2; 
	value[2]=1.0/h2;
	for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++ ) {
		value[1] = -2.0/h2 + diag[i];
		col[0]=i-1; 
		col[1]=i; 
		col[2]=i+1;
		ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES); CHKERRQ(ierr);
	}
	if (LastBlock) {
		i=Nz_grid-1; 
		col[0]=Nz_grid-2; 
		col[1]=Nz_grid-1;
		ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES); CHKERRQ(ierr);
	}
	if (FirstBlock) {
		i=0; 
		col[0]=0; 
		col[1]=1; 
		value[0]=-2.0/h2 + diag[0]; 
		value[1]=1.0/h2;
		ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES); CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// CHH 191022: MatGetVecs() is deprecated, changed to MatCreateVecs()
	ierr = MatCreateVecs(A,PETSC_NULL,&xr); CHKERRQ(ierr);
	ierr = MatCreateVecs(A,PETSC_NULL,&xi); CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Create the eigensolver and set various options
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* 
	Create eigensolver context
	*/
	ierr = EPSCreate(PETSC_COMM_WORLD,&eps); CHKERRQ(ierr);

	/* 
	Set operators. In this case, it is a standard eigenvalue problem
	*/
	ierr = EPSSetOperators(eps,A,PETSC_NULL); CHKERRQ(ierr);
	ierr = EPSSetProblemType(eps,EPS_NHEP); CHKERRQ(ierr);

	/*
	Set solver parameters at runtime
	*/
	ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
	ierr = EPSSetType(eps,"krylovschur"); CHKERRQ(ierr);
	ierr = EPSSetDimensions(eps,*nev,PETSC_DECIDE,PETSC_DECIDE); CHKERRQ(ierr); // leaving this line in speeds up the code; better if this is computed in chunks of 10? - consult Slepc manual
	ierr = EPSSetTarget(eps,sigma); CHKERRQ(ierr);
	ierr = EPSSetTolerances(eps,tol,PETSC_DECIDE); CHKERRQ(ierr);
	ierr = EPSSetTrueResidual(eps,PETSC_TRUE); CHKERRQ(ierr);

	ierr = EPSGetST(eps,&stx); CHKERRQ(ierr);

	// uncomment 2 lines
	// ierr = STGetKSP(stx,&kspx); CHKERRQ(ierr);
	// ierr = KSPGetPC(kspx,&pcx); CHKERRQ(ierr);

	ierr = STSetType(stx,"sinvert"); CHKERRQ(ierr);

	// uncomment 3 lines
	// ierr = KSPSetType(kspx,"preonly");
	// ierr = PCSetType(pcx,"cholesky");
	// ierr = EPSSetInterval(eps,pow(*k_min,2),pow(*k_max,2)); CHKERRQ(ierr);

	ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE); CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Solve the eigensystem
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = EPSSolve(eps);CHKERRQ(ierr);
	/*
	Optional: Get some information from the solver and display it
	*/
	ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
	ierr = EPSGetDimensions(eps,nev,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Display solution and clean up
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* 
	Get number of converged approximate eigenpairs
	*/
	ierr = EPSGetConverged(eps,nconv);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",*nconv);CHKERRQ(ierr);

	PetscReal error;
    // printf("        k           ||Ax-kx||/||kx||\n");
    // printf("-------------------------------------\n");
	if ((*nconv)>0) {
		for (i=0;i<(*nconv);i++) {
			ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
			ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
			re = PetscRealPart(kr);
			im = PetscImaginaryPart(kr);
			// printf("%03d  %9.6e     %9.6e\n", i, PetscRealPart(kr), error);
#else
			re = kr;
			im = ki;
#endif 
			//k2[(*nconv)-i-1] = re; // proper count of modes
			k2[i] = kr;
			ierr = VecGetArray(xr,&xr_);CHKERRQ(ierr);
			for (j = 0; j < Nz_grid; j++) {
				v[j][i] = xr_[j]/sqrt(dz); //per Slepc the 2-norm of xr_ is=1; we need sum(v^2)*dz=1 hence the scaling xr_/sqrt(dz)
			}
			ierr = VecRestoreArray(xr,&xr_);CHKERRQ(ierr);
		}
	}

	ierr = EPSDestroy(&eps);CHKERRQ(ierr);
	ierr = MatDestroy(&A);  CHKERRQ(ierr);
	ierr = VecDestroy(&xr); CHKERRQ(ierr);
	ierr = VecDestroy(&xi); CHKERRQ(ierr); 

	

	return 0;
}
