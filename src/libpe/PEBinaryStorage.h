#ifndef NCPAPROP_PEBINARYSTORAGE_H_INCLUDED
#define NCPAPROP_PEBINARYSTORAGE_H_INCLUDED

#include <iostream>
#include <complex>


namespace NCPA {

	class PEBinaryHeader {
	public:
		PEBinaryHeader();
		~PEBinaryHeader();
		bool read( std::istream &ifs );
		bool write( std::ostream &ofs );

		// before writing
		void setAzimuthVector( double *new_az_vec, size_t n_new_az );
		void setFrequencyVector( double *new_f_vec, size_t n_new_f );

		unsigned int n_az, n_freq;
		double precision_factor;
		double *az_vec, *freq_vec;
		bool az_init = false, freq_init = false;
	};

	class PEBinaryFrame {
	public:
		PEBinaryFrame();
		~PEBinaryFrame();
		bool read( std::istream &ifs, const NCPA::PEBinaryHeader &header );
		bool write( std::ostream &ofs, const NCPA::PEBinaryHeader &header );

		// before writing
		void setAltitudeVector( double *new_z_vec, size_t n_new_z );
		void setRangeVector( double *new_r_vec, size_t n_new_r );
		void setTLMatrix( std::complex< double > **new_tl_mat, size_t n_new_z, size_t n_new_r );

		double az, freq;
		unsigned int n_z, n_r;
		double *z_vec, *r_vec;
		std::complex< double > **TL_mat;

		bool z_init = false, r_init = false, TL_init = false;
	};

}





#endif