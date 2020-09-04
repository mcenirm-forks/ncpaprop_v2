#ifndef NCPAPROP_PEBROADBANDPROPAGATOR_H_INCLUDED
#define NCPAPROP_PEBROADBANDPROPAGATOR_H_INCLUDED

#include "parameterset.h"
#include "BroadbandPropagator.h"
#include <complex>

namespace NCPA {

	class PEBroadbandPropagator : public BroadbandPropagator {

	public:
		PEBroadbandPropagator( ParameterSet *param );
		~PEBroadbandPropagator();

		int calculate_waveform();

	protected:

		void read_pe_output_file(std::string filename, double az, double range );
		void read_calculation_dimensions( std::string filename );

		std::complex<double> **TL_r;
		size_t az_index;
		double precision_factor;

	};

}



#endif
