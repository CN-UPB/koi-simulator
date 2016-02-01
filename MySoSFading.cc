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
//   C.-X. Wang, M. Pätzold, and Q. Yao. Stochastic Modeling and Simulation of
//   Frequency-Correlated Wideband Fading Channels. IEEE Transactions on Vehicular
//   Technology, vol. 56, no. 3 (2007), pp. 1050–1063.
//
//   O. Punal, H. Escudero, and J. Gross. Performance Comparison of Loading
//   Algorithms for 80 MHz IEEE 802.11 WLANs. In Proceedings of the 73rd IEEE
//   Vehicular Technology Conference, Budapest, Hungary, 2011.
//


#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "omnetpp.h"
#include "MySoSFading.h"

/**
 * %%%%%%%%%%%%%%  SETUP CHANNEL  %%%%%%%%%%%%%%%%%
 * Support function to initialize global variables of the model.
 * Parameters: No parameters
 * Outcome: No direct output. Initialization of global variables.
 *
 */
void MySoSFading::setupChannel(double speed)
{
    // Maximal Doppler frequency shift (in Hz)
    fd_max = speed * (carrierFrequency / BaseWorldUtility$$speedOfLight);

    // Setting up Power Delay Profile (PDP)
    pdp = 1 / (1 - exp(-1 * (DELAY_MAX / ALPHA)));

    /* Filling the global arrays */
    int i, j, h, p;

    for (i = 0; i < L; i++)
    {
        a[i] = 0.0; // A=zeros(1, L);
        b[i] = 0.0; // B=zeros(1, L);

        for (j = 0; j < 2 * N; j++)
        {
            f[i][j] = 0.0; // f=zeros(L,2*N);

            for (h = 0; h < M; h++)
            {
                c[i][j][h] = 0.0; // C=zeros(L,2*N,M);
                theta[i][j][h] = 0.0; // Theta=zeros(L,2*N,M);
            }
        }

        for (p = 0; p < M; p++)
        {
            phi[i][p] = 0.0; //filling phi=zeros(L,M);
        }
    }

    /*  The Doppler phases are being constructed and randomly distributed  */
    generatePhases(L * 2 * N * M);
}

/**
 * %%%%%%%%%%%%%%  PHASES GENERATOR  %%%%%%%%%%%%%%%%%
 * Support function to generate random phases required for the channel model
 * Parameters: Number of potential waves arriving at the receiver (N*M*L*2)
 * Outcome: No direct output. Initialization of global variables.
 *
 */
void MySoSFading::generatePhases(int n)
{
    int i = 0; /* Variable for loops */
    double d1 = 0.0, d2 = 2 * M_PI; /* Angles for Phases between 0 and 2*Pi */
    double phases[n]; /* Array with linear spaced values for each entry from 0 to 2-pi */

    //Variables for randomization with indexes
    int p[n];

    // Spacing linearly the values between 0 and 2pi according with n
    for (i = 0; i <= n - 2; i++)
    {
        phases[i] = (d1 + (i * (d2 - d1))) / (n - 1);
    }

    phases[n - 1] = d2; // MatLab precision 15 decimals, this program gets up to 6

    //Filling array with the index
    for (i = 0; i < n; i++)
    {
        p[i] = i;
    }

    // Randomizing process twice the size of the array in order to avoid some minimal correlation
    for (i = 0; i < 2 * n; i++)
    {
        int rand_i1 = uniform(0, n - 1);
        int rand_i2 = uniform(0, n - 1);

        while (rand_i2 == rand_i1) //In case of collisions
            rand_i1 = uniform(0, n - 1);

        //swapping values pointed by random indexes
        int temp = p[rand_i2];
        p[rand_i2] = p[rand_i1];
        p[rand_i1] = temp;
    }

    // Once is computed values of phases and unsorted indexes, it matches them in one final vector
    for (i = 0; i < n; i++)
    {
        phases_vector[p[i]] = phases[i];
    }
}

/**
 * %%%%%%%%%%%%%%  MODELLING_CHANNEL: Sum-Of-Sinusoids (SoS) Fading Model   %%%%%%%%%%%%%%%%%
 * Main function (implemented within "SoSFadingMapping::getValue")for obtaining the fading per sub-carrier (FadingFactors[subcarrier])
 * Parameters required:
 * 1.- Simulation Time: SimTime => Since the model provides support for correlation in time, the time instance is relevant for the output
 * 2.- Simulation Seed: SimSeed => Seed that controls the random numbers involved in the model
 * 3.- Scenario Speed: ScnSpeed => Average Speed associated to transmitter/receiver and objects in the surroundings of the communication link
 * 4.- Carrier Frequency: CarrierFreq => Carrier frequency of the transmission (typically either the 2.4GHz band or the 5.2GHz band)
 * 5.- (Vector) Output of the model: FadingFactors[] => In this global variable the output of the model is going to be stored
 *
 * Outcome:
 * 1.- No direct output. However the output is indirectly stored in the variable "FadingFactors[]".
 *
 * For a proper understanding of the model employed, the careful reading of the following paper is strongly recommended:
 *
 * @InProceedings{ Paetzold07,
 author = "C.-X. Wang and M. Patzold and Q. Yao",
 booktitle = "{IEEE} Transactions on Vehicular Technology",
 title = "Stochastic Modeling and Simulation of Frequency-Correlated Wideband Fading Channels",
 year = "2007",
 month = "May",
 volume = "56",
 number = "3",
 pages = "1050--1063"}
 *
 */
double* MySoSFading::getFadingFactors(simtime_t simTime)
{
    double CarrierFreq = carrierFrequency;
    simtime_t TimeInstance;
    TimeInstance = simTime + TIME_OFFSET;

    int cluster; /* Principal Signal Clusters */
    int subcluster; /* Secondary Clusters (determine the Doppler frequency components (positive and negative frequencies) */
    int subcluster_index;
    int waves; /* Single Waves that conform the N secondary clusters (total number of waves per channel: 2*N*M*L) */
    int phase_index, i;
    int sc; /* Sub-carrier index */

    //Values for defining the electromagnetic wave associated to every OFDM sub-carrier
    double phase[SUBBANDS];
    double real_part[SUBBANDS];
    double imaginary_part[SUBBANDS];

    /* Variables Initialization */
    for (i = 0; i < SUBBANDS; i++)
    {
        //initializing the phase per sub-carrier accumulated due to propagation
        phase[i] = 0.0;

        /* For each sample it resets the real and imaginary part for each sub-carrier */
        real_part[i] = 0.0;
        imaginary_part[i] = 0.0;
    }

    //Principal Signal Clusters
    for (cluster = 0; cluster < L; cluster++)
    {
        //Secondary Clusters (determine the Doppler frequency components (positive and negative frequencies))
        for (subcluster = 0; subcluster < 2 * N; subcluster++)
        {
            //Single Waves that conform the N secondary clusters (total number of waves per channel: 2*N*M*L)
            for (waves = 0; waves < M; waves++)
            {
                phase_index = (cluster) * M * N * 2 + (subcluster) * M + waves;
                theta[cluster][subcluster][waves] = phases_vector[phase_index];

                a[cluster] = (double) ((tau[cluster] - (double) (delta_tau[cluster] / 2)) / ALPHA);
                a[cluster] = exp(-1 * a[cluster]);

                if (cluster < (L - 1))
                {
                    b[cluster] = (double) ((tau[cluster] + ((double) (delta_tau[cluster + 1] / 2))) / ALPHA);
                    b[cluster] = exp(-1 * b[cluster]);
                }

                else
                {
                    b[cluster] = (double) ((tau[cluster] + (double) (DELAY_MAX / 2)) / ALPHA);
                    b[cluster] = (double) exp(-1 * b[cluster]);
                }

                // Amplitude Coefficients: C[ l=cluster ][ n=multipath wave ][ m=single waves ]
                c[cluster][subcluster][waves] = sqrt(pdp * (double) (a[cluster] - b[cluster]) / (2 * N * M));

                // Doppler Frequency Components: f[l][n] - We should go from -N+1 to N
                subcluster_index = subcluster - N + 1; //subcluster_index [-19:20] without +1 it repeats first value twice and eats the last one.
                f[cluster][subcluster] = fd_max * sin((double) (PI / (2 * N))
                        * (subcluster_index - 0.5));

                // Single Wave Delays: Phi[l][m]
                phi[cluster][waves] = ALPHA * log((double) 1 / (a[cluster] - (((double) (waves + 1) / M) * (a[cluster] - b[cluster]))));

                //TODO: Until here, the code could be migrated to a "init" function, since it does not change!

                for (sc = 0; sc < SUBBANDS; sc++)
                {
                    // Time Selectivity Component (stems from Doppler spread)
                    phase[sc] = 2 * PI * f[cluster][subcluster] * TimeInstance.dbl();

                    // Frequency Selectivity Component (stemps from delay spread)
                    phase[sc] = phase[sc] - 2 * PI * phi[cluster][waves] * (CarrierFreq + sc * SUBCARRIER_SPACING);

                    // Phase Offset Component
                    phase[sc] = phase[sc] - theta[cluster][subcluster][waves];

                    // A * exp(j*b) = A * [cos(b) + j * sin(b)]^--
                    real_part[sc] = real_part[sc] + c[cluster][subcluster][waves] * cos(phase[sc]);
                    imaginary_part[sc] = imaginary_part[sc] + c[cluster][subcluster][waves] * sin(phase[sc]);
                }
            }
        }
    }

    //fading computation from all time samples and per sub-carrier
    for (sc = 0; sc < SUBBANDS; sc++)
    {
        FadingFactors[sc] = (real_part[sc] * real_part[sc]) + (imaginary_part[sc] * imaginary_part[sc]);
    }

    // The default implementation returns just a single attenuation value, which
    // is the resulting value over all sub-carriers. In fact this is enough to
    // have a global attenuation due to fading on the signal. A problem arises
    // here if we are interested in the individual (per sub-carrier) attenuations.
    // A workaround that respects the default model structure is to return a
    // single attenuation value while storing (globally) the vector containing
    // the fading coefficients
    return FadingFactors;
}

//This code is executed only once (per node) at the beginning of the simulation.
//Right place to initialize variables, but not the place to make function calls that should occur on runtime.
MySoSFading::MySoSFading(double carrierFrequency, double speed) :
    carrierFrequency(carrierFrequency), BaseWorldUtility$$speedOfLight(3e8)
{
    tau[0] = 0;
    tau[1] = 36e-09;
    tau[2] = 84e-09;
    tau[3] = 127e-09;

    delta_tau[0] = 0;
    delta_tau[1] = 36e-09;
    delta_tau[2] = 48e-09;
    delta_tau[3] = 43e-09;

    //This function will initialize the necessary variables that are required by the fading model throughout the simulation
    setupChannel(speed);
}
