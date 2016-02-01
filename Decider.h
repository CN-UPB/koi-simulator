//
// This file is part of a Horizon example simulation model. The code bases
// on the 'routing' example included with OMNeT++/OMNEST by Andras Varga.
//
// Copyright (C) 1992-2005 Andras Varga
// Copyright (C) 2008-2009 Karl Wessel, Wei Hong
// Copyright (C) 2010 Oscar Punal
// Copyright (C) 2011 Mirko Stoffers, Georg Kunz, Simon Tenbusch
//
// This file is distributed WITHOUT ANY WARRANTY. See the file
// `license' for details on this and other legal matters.
//
//
// Error model based on
//
//   C.-X. Wang, M. P��tzold, and Q. Yao. Stochastic Modeling and Simulation of
//   Frequency-Correlated Wideband Fading Channels. IEEE Transactions on Vehicular
//   Technology, vol. 56, no. 3 (2007), pp. 1050���1063.
//
//   O. Punal, H. Escudero, and J. Gross. Performance Comparison of Loading
//   Algorithms for 80 MHz IEEE 802.11 WLANs. In Proceedings of the 73rd IEEE
//   Vehicular Technology Conference, Budapest, Hungary, 2011.
//

#pragma once

#include "omnetpp.h"

/**
 * @brief Decider for the 802.11G modules
 *
 * Depending on the minimum of the snr included in the PhySDU this
 * module computes a bit error probability. The header (6 Mbit/s)
 *
 * @ingroup decider
 *
 * @author Marc L���bbers, David Raguin, Karl Wessel(port for MiXiM), Wei Hong(80211G modification)
 */
class Decider
{
    protected:
        /** @brief The center frequency on which the decider listens for signals */
        double centerFrequency;
        long subBands; // number of subbands
        double delayRms; // the delay of the path

    public:

        /**
         * @brief Initializes the Decider with a pointer to its PhyLayer and
         * specific values for threshold and sensitivity
         */
        Decider(
                double centerFrequency = 1.8e9,
                long subBands = 12,
                double delayRms = 20e-9,
                int myIndex = -1,
                bool debug = false) :
            centerFrequency(centerFrequency),
            subBands(subBands), delayRms(delayRms)
        {
        }

        /** @brief computes if packet is ok or has errors*/
        bool packetOK(int lengthMPDU, const double* const fadingFactors,
                const double randomNumber, double snirMin = 100000,
                double tx_data_rate = 0);
};
