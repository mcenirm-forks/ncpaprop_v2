#ifndef NCPAPROP_ESSCOMPLEXMODESOLVER_H_INCLUDED
#define NCPAPROP_ESSCOMPLEXMODESOLVER_H_INCLUDED

#include "slepceps.h"
#include "slepcst.h"
#include "ComplexModeSolver.h"
#include "parameterset.h"
#include "Atmosphere1D.h"
#include <string>
#include <complex>

namespace NCPA {

	class ESSComplexModeSolver : public ComplexModeSolver {

	public:
		ESSComplexModeSolver( ParameterSet *param, Atmosphere1D *atm_profile );
		~ESSComplexModeSolver();

		void setParams( ParameterSet *param, Atmosphere1D *atm_prof );
		void printParams();

		int doSelect( int nz, int n_modes, double k_min, double k_max, std::complex< double > *k2, 
			std::complex< double > **v, std::complex< double > *k_s, std::complex< double > **v_s, 
			int *select_modes );

		int getModalTrace(int nz, double z_min, double sourceheight, double receiverheight, 
			double dz, NCPA::Atmosphere1D *p, double admittance, double freq, double azi, 
			std::complex<double> *diag, double *k_min, double *k_max, 
			bool turnoff_WKB, double *c_eff);

		int solve();

		// still need to retrofit

		
		//int doESSSLEPcCalculation( double *diag, double dz, double *k_min, double *k_max, 
		//	PetscInt *nconv, double *k2, double **v );
		
		
		// void getModalStarter(int nz, int select_modes, double dz, double freq,  double z_src, 
		// 	double z_rcv, double *rho, std::complex<double> *k_pert, double **v_s, 
		// 	std::string modstartfile);
		// int writePhaseAndGroupSpeeds(int nz, double dz, int select_modes, double freq, 
		// 	std::complex<double> *k_pert, double **v_s, double *c_eff);  

	protected:
		std::string modstartfile; // store the modal starter in this file
		//bool   write_speeds;
	};
}




#endif
