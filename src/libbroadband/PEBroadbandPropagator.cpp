#include "PEBroadbandPropagator.h"
#include <complex>
#include "binaryreader.h"
#include <cstdint>
#include <fstream>

NCPA::ModalBroadbandPropagator::ModalBroadbandPropagator( NCPA::ParameterSet *param ) {

	// read in and store parameters
	NFFT 					= param->getInteger( "nfft" );
	waveform_out_file 		 = param->getString( "output_waveform_file" );
	dispersion_input_file 	= param->getString( "input_dispersion_file" );
	source_type 			= param->getString( "source" );
	source_file 			= param->getString( "source_file" );
	f_center 				= param->getFloat( "f_center" );
	max_cel 				= param->getFloat( "max_celerity" );
	single_receiver 		= param->getString("receiver").compare("single") == 0;

	if (single_receiver) {
		Nr = 1;
		r_vec = new double[ 1 ];
		r_vec[ 0 ] = param->getFloat( "range_km" ) * 1000.0;
	} else {
		double rmin, rmax, rstep;
		rmin = param->getFloat( "start_range_km" );
		rmax = param->getFloat( "end_range_km" );
		rstep = param->getFloat( "range_step_km" );
		if (rmax <= rmin) {
			throw std::runtime_error( "end_range_km must be greater than start_range_km" );
		}
		Nr = (int)floor( (rmax-rmin) / rstep ) + 1;
		r_vec = new double[ Nr ];
		for (int i = 0; i < Nr; i++) {
			r_vec[ i ] = (rmin + ((double)i) * rstep) * 1000.0;  // expected in meters
		}
	}

	// read output as a function of frequency and range for a given azimuth and height
	// this also sets Nfreq and f_vec
	read_calculation_dimensions( dispersion_input_file );
	TL_r  = cmatrix( Nr, Nfreq );
	f_step = f_vec[ 1 ] - f_vec[ 0 ];
}

NCPA::PEBroadbandPropagator::~PEBroadbandPropagator() {
	delete [] f_vec;
	delete [] r_vec;
}

void NCPA::PEBroadbandPropagator::read_calculation_dimensions( std::string filename ) {

	std::ifstream ins( filename, std::ifstream::in | std::ifstream::binary );
	uint32_t uholder[1];
	int64_t holder[1];
	size_t nAz;

	// read header
	ins.read( (char*)uholder, sizeof(uint32_t) );
	nAz = (size_t)(uholder[0]);
	ins.read( (char*)uholder, sizeof(uint32_t) );
	nFreq = (int)(uholder[0]);
	ins.read( (char*)uholder, sizeof(uint32_t) );
	precision_factor = (double)(uholder[0]);

	double *az_vec = new double[ nAz ];
	f_vec = new double[ nFreq ];

	int i;
	for (i = 0; i < nAz; i++) {
		ins.read( (char*)holder, sizeof(int64_t) );
		az_vec[ i ] = (double)(holder[0]) / precision_factor;
	}
	for (i = 0; i < nFreq; i++) {
		ins.read( (char*)holder, sizeof(int64_t) );
		f_vec[ i ] = (double)(holder[0]) / precision_factor;
	}
	ins.close();
}

void NCPA::PEBroadbandPropagator::read_response_spectrum( std::string filename, double this_az, double this_z ) {
	
}