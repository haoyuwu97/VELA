/// Definition of correlator classes
#ifndef __correlator_h
#define __correlator_h

#include <stdio.h>

class Correlator {

protected:
//Mode 0 
	// Where the coming values are stored
	double **shift;
	double **shift2;
	double **shift3;
	// Array containing the actual calculated correlation function
	double **correlation;
	double **correlation2;
	double **correlation3;
	// Number of values accumulated in cor
	unsigned long int **ncorrelation;
	unsigned long int **ncorrelation2;
	unsigned long int **ncorrelation3;

	// Accumulator in each correlator
	double *accumulator;
	double *accumulator2;
	double *accumulator3;
	// Index that controls accumulation in each correlator
	unsigned int *naccumulator;
	unsigned int *naccumulator2;
	unsigned int *naccumulator3;
	// Index pointing at the position at which the current value is inserted
	unsigned int *insertindex;
	unsigned int *insertindex2;
	unsigned int *insertindex3;
	
//Mode 1
	// Where the coming values are stored
	double **shift_aa;
	double **shift_bb;
	double **shift_cc;
	double **shift_ab;
	double **shift_ac;
	double **shift_bc;
	// Array containing the actual calculated correlation function
	double **correlation_aa;
	double **correlation_bb;
	double **correlation_cc;
	double **correlation_ab;
	double **correlation_ac;
	double **correlation_bc;
	// Number of values accumulated in cor
	unsigned long int **ncorrelation_aa;
	unsigned long int **ncorrelation_bb;
	unsigned long int **ncorrelation_cc;
	unsigned long int **ncorrelation_ab;
	unsigned long int **ncorrelation_ac;
	unsigned long int **ncorrelation_bc;

	// Accumulator in each correlator
	double *accumulator_aa;
	double *accumulator_bb;
	double *accumulator_cc;
	double *accumulator_ab;
	double *accumulator_ac;
	double *accumulator_bc;
	// Index that controls accumulation in each correlator
	unsigned int *naccumulator_aa;
	unsigned int *naccumulator_bb;
	unsigned int *naccumulator_cc;
	unsigned int *naccumulator_ab;
	unsigned int *naccumulator_ac;
	unsigned int *naccumulator_bc;
	// Index pointing at the position at which the current value is inserted
	unsigned int *insertindex_aa;
	unsigned int *insertindex_bb;
	unsigned int *insertindex_cc;
	unsigned int *insertindex_ab;
	unsigned int *insertindex_ac;
	unsigned int *insertindex_bc;

	// Number of Correlators
	unsigned int numcorrelators;

	// Minimum distance between points for correlators k>0; dmin = p/m
	unsigned int dmin;
	
	//

	// Lenght of result arrays
	unsigned int length;
	// Maximum correlator attained during simulation 
	unsigned int kmax;

public:
	// Points per correlator
	unsigned int p;
	// Number of points over which to average; RECOMMENDED: p mod m = 0
	unsigned int m; 
	double *t, *f1_mode0, *f2_mode0, *f3_mode0;
	double *f1_mode1, *f2_mode1, *f3_mode1, *f4_mode1, *f5_mode1, *f6_mode1;
	double *wt, *Gt, *Gw, *Gw_storage, *Gw_loss;
	unsigned int npcorr;
	unsigned int npcorr_Gw;

	//Mode 0:
	// Accumulated result of incoming values
	double accval;
	double accval2;
	double accval3;
	
	//Mode 1:
	// Accumulated result of incoming values
	double accval_aa;
	double accval_bb;
	double accval_cc;
	double accval_ab;
	double accval_ac;
	double accval_bc;

	// Constructor
	Correlator () {numcorrelators=0;} ;
	Correlator (const unsigned int numcorrin, const unsigned int pin, const unsigned int min);
	~Correlator();

	// Set size of correlator
	void setsize (const unsigned int numcorrin = 32, const unsigned int pin = 16, const unsigned int min = 2);

	// Add a scalar to the correlator number k
	void add_mode0(const double w, const double w2, const double w3, const unsigned int k = 0);
	void add_mode1(const double w1, const double w2, const double w3, const double w4, const double w5, const double w6, const unsigned int k = 0);

	// Evaluate the current state of the correlator
	void evaluate_mode0(const bool norm = false);
	void evaluate_mode1(const bool norm = false);
	
	// Evaluate the current state of the correlator
	void calculate_Gw(const double dtime, const int time, const int Gtcut, const double tail_param = 2.0);

	// Initialize all values (current and average) to zero
	void initialize_mode0();
	void initialize_mode1();
	void initialize_Gw();

	};

#endif
