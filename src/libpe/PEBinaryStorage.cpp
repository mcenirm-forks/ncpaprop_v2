#include "PEBinaryStorage.h"
#include "util.h"
#include <fstream>
#include <iostream>
#include <cstdint>
#include <complex>
#include <stdexcept>



NCPA::PEBinaryHeader::PEBinaryHeader() { }

NCPA::PEBinaryHeader::~PEBinaryHeader() {
	if (az_init) {
		delete [] az_vec;
	}
	if (freq_init) {
		delete [] freq_vec;
	}
}

void NCPA::PEBinaryHeader::setAzimuthVector( double *new_az_vec, size_t n_new_az ) {
	if (az_init) {
		delete [] az_vec;
		az_init = false;
	}
	az_vec = new double[ n_new_az ];
	std::memcpy( this->az_vec, new_az_vec, n_new_az * sizeof( double ) );
	n_az = n_new_az;
	az_init = true;
}

void NCPA::PEBinaryHeader::setFrequencyVector( double *new_f_vec, size_t n_new_f ) {
	if (freq_init) {
		delete [] freq_vec;
		freq_init = false;
	}
	freq_vec = new double[ n_new_f ];
	std::memcpy( this->freq_vec, new_f_vec, n_new_f * sizeof( double ) );
	n_freq = n_new_f;
	freq_init = true;
}

bool NCPA::PEBinaryHeader::write( std::ostream &ofs ) {

	size_t i = 0;
	size_t buf_size = n_az;
	if (n_freq > buf_size) {
		buf_size = n_freq;
	}
	int64_t *buffer = new int64_t[ buf_size ];
	std::memset( buffer, 0, buf_size * sizeof( int64_t ) );

	// write header starting with vector sizes and multiplicative factor
	uint32_t holder = n_az;
	ofs.write( (char*)(&holder), sizeof( uint32_t ) );
	holder = n_f;
	ofs.write( (char*)(&holder), sizeof( uint32_t ) );
	holder = precision_factor;
	ofs.write( (char*)(&holder), sizeof( uint32_t ) );

	for (i = 0; i < n_az; i++) {
		buffer[ i ] = (int64_t)std::lround( az_vec[ i ] * (double)precision_factor );
	}
	ofs.write( (char*)buffer, n_az * sizeof( int64_t ) );
	std::memset( buffer, 0, buf_size * sizeof( int64_t ) );
	for (i = 0; i < n_f; i++) {
		buffer[ i ] = (int64_t)std::lround( freq_vec[ i ] * (double)precision_factor );
	}
	ofs.write( (char*)buffer, n_f * sizeof( int64_t ) );

	delete [] buffer;
	return true;
}

bool NCPA::PEBinaryHeader::read( std::istream &ifs ) {
	uint32_t uholder[1];
	int64_t holder[1];

	// read header
	ifs.read( (char*)uholder, sizeof(uint32_t) );
	n_az = (unsigned int)(uholder[0]);
	ifs.read( (char*)uholder, sizeof(uint32_t) );
	n_freq = (unsigned int)(uholder[0]);
	ifs.read( (char*)uholder, sizeof(uint32_t) );
	precision_factor = (double)(uholder[0]);

	if (az_init) {
		delete [] az_vec;
		az_init = false;
	}
	az_vec = new double[ n_az ];
	az_init = true;
	if (freq_init) {
		delete [] freq_vec;
		freq_init = false;
	}
	f_vec = new double[ n_freq ];
	freq_init = true;

	int i;
	for (i = 0; i < nAz; i++) {
		ifs.read( (char*)holder, sizeof(int64_t) );
		az_vec[ i ] = (double)(holder[0]) / precision_factor;
	}
	for (i = 0; i < nFreq; i++) {
		ifs.read( (char*)holder, sizeof(int64_t) );
		f_vec[ i ] = (double)(holder[0]) / precision_factor;
	}
}

NCPA::PEBinaryFrame::PEBinaryFrame() { }

NCPA::PEBinaryFrame::~PEBinaryFrame() {
	if (TL_init) {
		NCPA::free_cmatrix( TL_mat, n_z, n_r );
	}
	if (z_init) {
		delete [] z_vec;
	}
	if (r_init) {
		delete [] r_vec;
	}
	
}

bool NCPA::PEBinaryFrame::write( std::ostream &ofs, const NCPA::PEBinaryHeader &header ) {
	// write az, freq, n_z, n_range
	int64_t holder = (int64_t)std::lround( az * header.precision_factor );
	ofs.write( (char*)(&holder), sizeof( int64_t ) );
	holder = (int64_t)std::lround( freq * header.precision_factor );
	ofs.write( (char*)(&holder), sizeof( int64_t ) );
	uint32_t uholder = (uint32_t)n_z;
	ofs.write( (char*)(&uholder), sizeof( uint32_t ) );
	uholder = (uint32_t)n_r;
	ofs.write( (char*)(&uholder), sizeof( uint32_t ) );

	// z and r sizes and vectors
	size_t buf_size = n_r;
	if (n_z > buf_size) {
		buf_size = n_z;
	}
	int64_t *buffer = new int64_t[ buf_size ];
	std::memset( buffer, 0, buf_size * sizeof( int64_t ) );
	size_t i, j;
	for (i = 0; i < n_z; i++) {
		buffer[ i ] = (int64_t)std::lround( z_vec[ i ] * header.precision_factor );
	}
	ofs.write( (char*)buffer, n_z * sizeof( int64_t ) );
	for (i = 0; i < n_r; i++) {
		buffer[ i ] = (int64_t)std::lround( r_vec[ i ] * header.precision_factor );
	}
	ofs.write( (char*)buffer, n_r * sizeof( int64_t ) );
	for (i = 0; i < n_z; i++) {
		for (j = 0; j < n_r; j++) {
			holder = (int64_t)std::lround( TL_mat[ i ][ j ].real() * header.precision_factor );
			ofs.write( (char *)(&holder), sizeof( int64_t ) );
			holder = (int64_t)std::lround( TL_mat[ i ][ j ].imag() * header.precision_factor );
			ofs.write( (char *)(&holder), sizeof( int64_t ) );
		}
	}
	
	delete [] buffer;
	return true;
}

bool NCPA::PEBinaryFrame::read( std::istream &ifs, const NCPA::PEBinaryHeader &header ) {
	int64_t holder[2];
	uint32_t uholder[1];

	// read in azimuth and frequency
	ifs.read( (char*)holder, sizeof( int64_t) );
	az = (double)(holder[0]) / header.precision_factor;
	ifs.read( (char*)holder, sizeof( int64_t) );
	freq = (double)(holder[0]) / header.precision_factor;

	// read in altitude and range sizes
	ifs.read( (char*)uholder, sizeof( uint32_t ) );
	n_z = (unsigned int)(uholder[0]);
	ifs.read( (char*)uholder, sizeof( uint32_t ) );
	n_r = (unsigned int)(uholder[0]);

	// read in altitude and range vectors
	if (z_init) {
		delete [] z_vec;
		z_init = false;
	}
	z_vec = new double[ n_z ];
	z_init = true;
	size_t i, j;
	for (i = 0; i < n_z; i++) {
		ifs.read( (char*)holder, sizeof( int64_t ) );
		z_vec[ i ] = (double)(holder[0]) / header.precision_factor;
	}
	if (r_init) {
		delete [] r_vec;
		r_init = false;
	}
	r_vec = new double[ n_r ];
	r_init = true;
	for (i = 0; i < n_r; i++) {
		ifs.read( (char*)holder, sizeof( int64_t ) );
		r_vec[ i ] = (double)(holder[0]) / header.precision_factor;
	}

	// read in transmission loss
	for (i = 0; i < n_z; i++) {
		for (j = 0; j < n_r; j++) {
			ifs.read( (char*)holder, 2*sizeof( int64_t ) );
			TL_mat[ i ][ j ].real( (double)holder[0] / header.precision_factor );
			TL_mat[ i ][ j ].imag( (double)holder[1] / header.precision_factor );
		}
	}

	return true;
}

void NCPA::PEBinaryFrame::setAltitudeVector( double *new_z_vec, size_t n_new_z ) {
	if (z_init) {
		delete [] z_vec;
		z_init = false;
	}
	z_vec = new double[ n_new_z ];
	std::memcpy( this->z_vec, new_z_vec, n_new_z * sizeof( double ) );
	n_z = n_new_z;
	z_init = true;
}

void NCPA::PEBinaryFrame::setRangeVector( double *new_r_vec, size_t n_new_r ) {
	if (r_init) {
		delete [] r_vec;
		r_init = false;
	}
	r_vec = new double[ n_new_r ];
	std::memcpy( this->r_vec, new_r_vec, n_new_r * sizeof( double ) );
	n_r = n_new_r;
	r_init = true;
}

void NCPA::PEBinaryFrame::setTLMatrix( std::complex< double > **new_tl_mat ) {
	if (TL_init) {
		NCPA::free_cmatrix( TL_mat, n_z, n_r );
		TL_init = false;
	}
	TL_mat = NCPA::cmatrix( n_z, n_r );
	std::memcpy( this->TL_mat, new_tl_mat, n_r * n_z * sizeof( std::complex< double > ) );
	TL_init = true;
}