//
// This file is part of a Horizon example simulation model. The code bases
// on the 'routing' example included with OMNeT++/OMNEST by Andras Varga.
//
// Copyright (C) 1992-2005 Andras Varga
// Copyright (C) 2010 Oscar Punal
// Copyright (C) 2011 Mirko Stoffers, Georg Kunz, Simon Tenbusch
//
// This file is distributed WITHOUT ANY WARRANTY. See the file
// `license' for details on this and other legal matters.
//
//
// Channel model based on
//
//   C.-X. Wang, M. P��tzold, and Q. Yao. Stochastic Modeling and Simulation of
//   Frequency-Correlated Wideband Fading Channels. IEEE Transactions on Vehicular
//   Technology, vol. 56, no. 3 (2007), pp. 1050���1063.
//
//   O. Punal, H. Escudero, and J. Gross. Performance Comparison of Loading
//   Algorithms for 80 MHz IEEE 802.11 WLANs. In Proceedings of the 73rd IEEE
//   Vehicular Technology Conference, Budapest, Hungary, 2011.
//

/** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Channel Model based on a sum of sinusoids approach by M. Paetzold.
 %% Full deterministic model, where amplitude coefficients, Doppler
 %% frequency shifts, waves delays and phases are defined during the
 %% model set-up phase.

 *  @date_of_creation - 18 March 2010
 *  @place - UMIC Mobile Network Performance

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#pragma once

#include <omnetpp.h>

/***********************************************
 *************** CONSTANTS ***************
 ***********************************************/
#define ALPHA 25e-09		    					/* ALPHA - RMS Delay Spread (the larger the more frequency selectivity) */
#define DELAY_MAX 127e-09 				  			/* DELAY_MAX - Maximal Delay */
#define L 4				                			/* L - Number of multipath waves (cluster) */
#define M 20				    			        /* M - Number of multipath waves (sub-cluster) */
#define N 20	            						/* N - Number of Waves conforming the sub-cluster waves */
#define SUBCARRIER_SPACING 312500   				/* SUBCARRIER_SPACING  - 312.5 KHz */
#define TIME_OFFSET 0	            				/* Offset time of simulation (in order not to start at time=0)*/

#define SUBBANDS 12		/* Number of payload OFDM subcarriers (In general, 48 payload subcarriers over a 20MHz bandwidth)*/
#define BPSK 			1
#define QPSK			2
#define QAM16			4
#define QAM64			6
#define QAM256			8

class MySoSFading
{
    protected:
        friend class SoSFadingMapping;

        friend class LoadingAlgorithms;

        //Simulation Seed
        int SimSeed;

        //Movement speed in the environment
        double Speed;

        /** @brief Carrier frequency to be used. */
        double carrierFrequency;

        double FadingFactors[SUBBANDS];

        /***********************************************
         ********** GLOBAL VARIABLES  *************
         ***********************************************/
        double fd_max; /* Maximal Doppler frequency shift (in Hz) */
        double pdp; /* Power Delay Profile (PDP) normalization factor */

        /*
         * Tap Delays where the delays (ns) = [0, 36, 84, 127]
         * (from Paper: "Tapped Delay Line Channel Models at 5.3 GHz in Indoor Environments" - Table 1)
         * Different parameterization can be obtained from other measurements campaigns.
         */
        double tau[4];
        double delta_tau[4]; // Obtained by obtaining: tau(cluster+1) - tau(cluster)

        // Preallocating Vectors involved in the calculation of the fading
        double a[L];
        double b[L];
        double c[L][2 * N][M];
        double f[L][2 * N];
        double phi[L][M];
        double theta[L][2 * N][M];
        double phases_vector[L * 2 * N * M]; /* Pseudo-deterministic and randomized Doppler phases */

        const double BaseWorldUtility$$speedOfLight;

    public:

        MySoSFading(double carrierFrequency = 1.8e9, double speed = 0.5);

        void generatePhases(int n);

        void setupChannel(double speed);

        double* getFadingFactors(simtime_t simTime);
};
