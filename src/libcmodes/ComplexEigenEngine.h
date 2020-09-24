#ifndef NCPAPROP_COMPLEXEIGENENGINE_H_INCLUDED
#define NCPAPROP_COMPLEXEIGENENGINE_H_INCLUDED

#include "slepceps.h"
#include "slepcst.h"
#include <string>
#include <complex>

namespace NCPA {

	class ComplexEigenEngine {

	public:
		static int doESSCalculation( std::complex<double> *diag, int Nz_grid, double dz, double tol, 
			int *nev, double *k_min, double *k_max, PetscInt *nconv, std::complex<double> *k2, std::complex<double> **v );
		// static int doWideAngleCalculation( int Nz_grid, double dz, double k_min, double k_max,
		// 	double tol, int nev, double *kd, double *md, double *cd, 
		// 	PetscInt *nconv, double *kH, double **v, std::string disp_msg );


	};

}








#endif