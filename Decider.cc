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


#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "Decider.h"

bool Decider::packetOK(int lengthMPDU, const double* const fadingFactors,
        const double randomNumber, double snirMin, double tx_data_rate)
{
    long double Pburst_upper, Pburst_upper_out = 0, c[30], packet_error_rate, ber_averaged, ber[subBands];
    int x, u, d;
    double snr_fading_per_sc[subBands], att;
    long double x0 = 0, x1 = 0, x2 = 0, y0 = 0, y1 = 0, y2 = 0;

    /* Jakes / Bell Shape Channel Model */
    simtime_t t;

    x = 0;
    int count_bpsk = 0;
    int count_qpsk = 0;
    int count_qam16 = 0;
    int count_qam64 = 0;
    int count_qam64_lte = 0;

    for (x = 0; x < subBands; x++)
    {
        snr_fading_per_sc[x] = snirMin * fadingFactors[x];
    }

    for (x = 0; x < subBands; x++)
    {
        ber[x] = 0;

        //lte qam64
        if (tx_data_rate == 0) 
        {
            count_qam64_lte++;
            int k = 6;
            int M = 64;
            ber[x] = (4./k)*(1.-1./sqrt(M))*(1./2.)*erfc(sqrt(3.*k*snr_fading_per_sc[x]/(M-1.))/sqrt(2.)); 
        }
        if (tx_data_rate == 6000000 || tx_data_rate == 9000000)
        {
            count_bpsk++;
            //debug: EV << "bpsk modulation" << endl;
            ber[x] = 0.5 * erfc(sqrt(snr_fading_per_sc[x]));
            if (ber[x] > 0.5)
                ber[x] = 0.5;
        }

        else if (tx_data_rate == 12000000 || tx_data_rate == 18000000)
        {
            count_qpsk++;
            //debug: EV << "qpsk modulation" << endl;

            snr_fading_per_sc[x] = 10 * log10(snr_fading_per_sc[x]) - 3;
            snr_fading_per_sc[x] = pow(10.0, snr_fading_per_sc[x] / 10);
            ber[x] = 0.5 * erfc(sqrt(snr_fading_per_sc[x] * 0.5));
            ber[x] = ber[x] * 2;
            if (ber[x] > 0.5)
                ber[x] = 0.5;
        }

        else if (tx_data_rate == 24000000 || tx_data_rate == 36000000)
        {
            count_qam16++;
            //debug: Ev << "16-qam modulation" << endl;
            snr_fading_per_sc[x] = 10 * log10(snr_fading_per_sc[x]) - 6;
            snr_fading_per_sc[x] = pow(10.0, snr_fading_per_sc[x] / 10);

            ber[x] = 0.375 * erfc(sqrt(0.1 * snr_fading_per_sc[x]));
            ber[x] += 0.25 * erfc(3 * sqrt(0.1 * snr_fading_per_sc[x]));
            ber[x] = ber[x] * 4;
            if (ber[x] > 0.5)
                ber[x] = 0.5;
        }

        else if (tx_data_rate == 48000000 || tx_data_rate == 54000000)
        {
            count_qam64++;
            //debug: EV << "64-qam modulation" << endl;
            snr_fading_per_sc[x] = 10 * log10(snr_fading_per_sc[x]) - 7.78;
            snr_fading_per_sc[x] = pow(10.0, snr_fading_per_sc[x] / 10);

            ber[x] = 0.2916666 * erfc(sqrt(0.0238095 * snr_fading_per_sc[x]))
                    + 0.25 * erfc(3 * sqrt(0.0238095 * snr_fading_per_sc[x]));
            ber[x] = ber[x] * 6;
            if (ber[x] > 0.5)
                ber[x] = 0.5;
        }
    }

    ber_averaged = 0;
    for (x = 0; x < subBands; x++)
        ber_averaged += ber[x];

    if ((count_bpsk + count_qpsk + count_qam16 + count_qam64 + count_qam64_lte) == 0)
        ber_averaged = 0.5;
    else
        ber_averaged = ber_averaged / (count_qam64_lte + count_bpsk + count_qpsk * 2
                + count_qam16 * 4 + count_qam64 * 6);

    for (u = 0; u <= 29; u++)
        c[u] = 0;

    if (tx_data_rate == 6000000 || tx_data_rate == 12000000 || tx_data_rate == 24000000)
    {
        for (d = 10; d <= 29; d = d + 2)
            c[d] = pow(2 * (pow((ber_averaged * (1 - ber_averaged)), 0.5)), d);

        att = 1;
        Pburst_upper = 0;
        Pburst_upper_out = 0;
        Pburst_upper = 36 * c[10] + 211 * c[12] + 1404 * c[14] + 11633 * c[16]
                + 77433 * c[18] + 502690 * c[20] + 3322763 * c[22] + 21292910
                * c[24] + 134365911 * c[26] + 843425871 * c[28];
        Pburst_upper = Pburst_upper / att;

        if (Pburst_upper <= 12.4 && Pburst_upper > 2.585)
        {
            x0 = 2.585;
            x1 = 5.719;
            x2 = 12.4;
            y0 = 45.5667;
            y1 = 60.4418;
            y2 = 109.7345;
        }
        else if (Pburst_upper <= 2.585 && Pburst_upper > 0.5043)
        {
            x0 = 0.5043;
            x1 = 1.149;
            x2 = 2.585;
            y0 = 15.291;
            y1 = 31.4192;
            y2 = 45.5667;
        }
        else if (Pburst_upper <= 0.5043 && Pburst_upper > 1.898e-2)
        {
            x0 = 1.898e-2;
            x1 = 9.604e-2;
            x2 = 0.5043;
            y0 = 3.2831;
            y1 = 7.3593;
            y2 = 15.291;
        }
        else if (Pburst_upper <= 1.898e-2 && Pburst_upper > 4.903e-4)
        {
            x0 = 4.903e-4;
            x1 = 1.985e-3;
            x2 = 1.898e-2;
            y0 = 1.3956;
            y1 = 1.9215;
            y2 = 3.2831;
        }
        else if (Pburst_upper <= 4.903e-4 && Pburst_upper >= 1.524e-5)
        {
            x0 = 1.524e-5;
            x1 = 6.26e-5;
            x2 = 4.903e-4;
            y0 = 0.94541;
            y1 = 0.99507;
            y2 = 1.3956;
        }

        Pburst_upper_out = Pburst_upper / (y0 * ((Pburst_upper - x1)
                * (Pburst_upper - x2)) / ((x0 - x1) * (x0 - x2)) + y1
                * ((Pburst_upper - x0) * (Pburst_upper - x2)) / ((x1 - x0)
                * (x1 - x2)) + y2 * ((Pburst_upper - x0) * (Pburst_upper - x1))
                / ((x2 - x0) * (x2 - x1)));
        if (Pburst_upper > 12.4)
            Pburst_upper_out = 0.15;
        else if (Pburst_upper < 1.524e-5)
            Pburst_upper_out = Pburst_upper;
    }

    else if (tx_data_rate == 9000000 || tx_data_rate == 18000000
            || tx_data_rate == 36000000 || tx_data_rate == 54000000)
    {
        d = 0;
        for (d = 5; d < 15; d++)
            c[d] = pow(2 * (pow((ber_averaged * (1 - ber_averaged)), 0.5)), d);

        att = 3;
        Pburst_upper = 0;
        Pburst_upper_out = 0;
        Pburst_upper = 42 * c[5] + 201 * c[6] + 1492 * c[7] + 10469 * c[8]
                + 62935 * c[9] + 379644 * c[10] + 2253373 * c[11] + 13073811
                * c[12] + 75152755 * c[13] + 428005675 * c[14];
        Pburst_upper = Pburst_upper / att;

        if (Pburst_upper <= 1023 && Pburst_upper > 344.3)
        {
            x0 = 344.3;
            x1 = 599.7;
            x2 = 1023;
            y0 = 3.8074e3;
            y1 = 5.2605e3;
            y2 = 7.1091e3;
        }
        else if (Pburst_upper <= 344.3 && Pburst_upper > 15.4)
        {
            x0 = 15.4;
            x1 = 57.07;
            x2 = 344.3;
            y0 = 2.3422e3;
            y1 = 2.4234e3;
            y2 = 3.8074e3;
        }
        else if (Pburst_upper <= 15.4 && Pburst_upper > 0.4095)
        {
            x0 = 0.4095;
            x1 = 1.841;
            x2 = 15.4;
            y0 = 1.0383e3;
            y1 = 1.519e3;
            y2 = 2.3422e3;
        }
        else if (Pburst_upper <= 0.4095 && Pburst_upper > 2.917e-4)
        {
            x0 = 2.917e-4;
            x1 = 9.062e-3;
            x2 = 0.4095;
            y0 = 6.71655e2;
            y1 = 677.7861;
            y2 = 1.0383e3;
        }
        else if (Pburst_upper <= 2.917e-4 && Pburst_upper > 4.107e-5)
        {
            x0 = 4.107e-5;
            x1 = 1.52e-4;
            x2 = 2.917e-4;
            y0 = 1.045e3;
            y1 = 757.7268;
            y2 = 6.71655e2;
        }
        else if (Pburst_upper <= 4.107e-5 && Pburst_upper > 1.068e-5)
        {
            x0 = 1.068e-5;
            x1 = 2.109e-5;
            x2 = 4.107e-5;
            y0 = 1.5575e3;
            y1 = 1.2606e3;
            y2 = 1.045e3;
        }
        else if (Pburst_upper <= 1.068e-5 && Pburst_upper > 1.243e-6)
        {
            x0 = 1.243e-6;
            x1 = 5.324e-6;
            x2 = 1.068e-5;
            y0 = 3.3041e3;
            y1 = 1.9646e3;
            y2 = 1.5575e3;
        }
        else if (Pburst_upper <= 1.243e-6 && Pburst_upper > 2.628e-7)
        {
            x0 = 2.628e-7;
            x1 = 5.793e-7;
            x2 = 1.243e-6;
            y0 = 5.9863e3;
            y1 = 4.4053e3;
            y2 = 3.3041e3;
        }
        else if (Pburst_upper <= 2.628e-7 && Pburst_upper > 2.041e-8)
        {
            x0 = 2.041e-8;
            x1 = 1.158e-7;
            x2 = 2.628e-7;
            y0 = 1.6896e4;
            y1 = 8.2951e3;
            y2 = 5.9863e3;
        }
        else if (Pburst_upper <= 2.041e-8 && Pburst_upper > 3.116e-9)
        {
            x0 = 3.116e-9;
            x1 = 8.127e-9;
            x2 = 2.041e-8;
            y0 = 3.7376e4;
            y1 = 2.4853e4;
            y2 = 1.6896e4;
        }
        else if (Pburst_upper <= 3.116e-9 && Pburst_upper > 1.367e-10)
        {
            x0 = 1.367e-10;
            x1 = 2.041e-8;
            x2 = 3.116e-9;
            y0 = 1.4567e5;
            y1 = 1.6896e4;
            y2 = 3.7376e4;
        }
        else if (Pburst_upper <= 1.367e-10 && Pburst_upper > 1.347e-11)
        {
            x0 = 1.347e-11;
            x1 = 4.399e-11;
            x2 = 1.367e-10;
            y0 = 4.0781e5;
            y1 = 2.4065e5;
            y2 = 1.4567e5;
        }

        Pburst_upper_out = Pburst_upper / (y0 * ((Pburst_upper - x1)
                * (Pburst_upper - x2)) / ((x0 - x1) * (x0 - x2)) + y1
                * ((Pburst_upper - x0) * (Pburst_upper - x2)) / ((x1 - x0)
                * (x1 - x2)) + y2 * ((Pburst_upper - x0) * (Pburst_upper - x1))
                / ((x2 - x0) * (x2 - x1)));
        if (Pburst_upper > 1023)
            Pburst_upper_out = 0.1439;
        else if (Pburst_upper < 1.347e-11)
            Pburst_upper_out = 3.303e-17;
    }

    else if (tx_data_rate == 48000000)
    {
        d = 0;
        for (d = 6; d < 16; d++)
            c[d] = pow(2 * (pow((ber_averaged * (1 - ber_averaged)), 0.5)), d);

        att = 2;
        Pburst_upper = 0;
        Pburst_upper_out = 0;
        Pburst_upper = 3 * c[6] + 70 * c[7] + 285 * c[8] + 1276 * c[9] + 6160
                * c[10] + 27128 * c[11] + 117019 * c[12] + 498860 * c[13]
                + 2103891 * c[14] + 8784123 * c[15];
        Pburst_upper = Pburst_upper / att;

        if (Pburst_upper <= 49.57 && Pburst_upper > 5.68)
        {
            x0 = 5.68;
            x1 = 17.45;
            x2 = 49.57;
            y0 = 184.6554;
            y1 = 212.4939;
            y2 = 430.2951;
        }
        else if (Pburst_upper <= 5.68 && Pburst_upper > 0.2468)
        {
            x0 = 0.2468;
            x1 = 0.9124;
            x2 = 5.68;
            y0 = 124.4579;
            y1 = 149.8939;
            y2 = 184.6554;
        }
        else if (Pburst_upper <= 0.2468 && Pburst_upper > 3.147e-2)
        {
            x0 = 3.147e-2;
            x1 = 0.1256;
            x2 = 0.2468;
            y0 = 116.5556;
            y1 = 139.2770;
            y2 = 124.4579;
        }
        else if (Pburst_upper <= 3.147e-2 && Pburst_upper > 7.699e-3)
        {
            x0 = 7.699e-3;
            x1 = 0.01559;
            x2 = 3.147e-2;
            y0 = 62.8490;
            y1 = 96.1752;
            y2 = 116.5556;
        }
        else if (Pburst_upper <= 7.699e-3 && Pburst_upper > 9.312e-4)
        {
            x0 = 9.312e-4;
            x1 = 3.8e-3;
            x2 = 7.699e-3;
            y0 = 85.0411;
            y1 = 70.8031;
            y2 = 62.8490;
        }
        else if (Pburst_upper <= 9.312e-4 && Pburst_upper > 2.296e-4)
        {
            x0 = 2.296e-4;
            x1 = 4.623e-4;
            x2 = 9.312e-4;
            y0 = 99.1364;
            y1 = 91.8538;
            y2 = 85.0411;
        }
        else if (Pburst_upper <= 2.296e-4 && Pburst_upper > 2.75e-5)
        {
            x0 = 2.75e-5;
            x1 = 1.138e-4;
            x2 = 2.296e-4;
            y0 = 129.9008;
            y1 = 107.5614;
            y2 = 99.1364;
        }
        else if (Pburst_upper <= 2.75e-5 && Pburst_upper > 6.364e-6)
        {
            x0 = 6.364e-6;
            x1 = 1.3321e-5;
            x2 = 2.75e-5;
            y0 = 163.4309;
            y1 = 144.9402;
            y2 = 129.9008;
        }
        else if (Pburst_upper <= 6.364e-6 && Pburst_upper > 6.239e-7)
        {
            x0 = 6.239e-7;
            x1 = 2.992e-6;
            x2 = 6.364e-6;
            y0 = 249.1613;
            y1 = 186.1854;
            y2 = 163.4309;
        }
        else if (Pburst_upper <= 6.239e-7 && Pburst_upper > 1.188e-7)
        {
            x0 = 1.188e-7;
            x1 = 2.756e-7;
            x2 = 6.239e-7;
            y0 = 346.5578;
            y1 = 292.3828;
            y2 = 249.1613;
        }
        else if (Pburst_upper <= 1.188e-7 && Pburst_upper > 8.074e-9)
        {
            x0 = 8.074e-9;
            x1 = 2.756e-7;
            x2 = 1.188e-7;
            y0 = 613.9924;
            y1 = 292.3828;
            y2 = 346.5578;
        }
        else if (Pburst_upper <= 8.074e-9 && Pburst_upper > 1.153e-9)
        {
            x0 = 1.153e-9;
            x1 = 3.102e-9;
            x2 = 8.074e-9;
            y0 = 949.7529;
            y1 = 759.1777;
            y2 = 613.9924;
        }
        else if (Pburst_upper <= 1.153e-9 && Pburst_upper > 4.779e-11)
        {
            x0 = 4.779e-11;
            x1 = 4.142e-10;
            x2 = 1.153e-9;
            y0 = 2.0165e3;
            y1 = 1.2041e3;
            y2 = 949.7529;
        }
        else if (Pburst_upper <= 4.779e-11 && Pburst_upper > 4.687e-12)
        {
            x0 = 4.687e-12;
            x1 = 1.529e-11;
            x2 = 4.779e-11;
            y0 = 3.5916e3;
            y1 = 2.6689e3;
            y2 = 2.0165e3;
        }
        else if (Pburst_upper <= 4.687e-12 && Pburst_upper > 1.018e-13)
        {
            x0 = 1.018e-13;
            x1 = 1.373e-12;
            x2 = 4.687e-12;
            y0 = 9.7791e3;
            y1 = 4.9194e3;
            y2 = 3.5916e3;
        }
        else if (Pburst_upper <= 1.018e-13 && Pburst_upper > 6.093e-15)
        {
            x0 = 6.093e-15;
            x1 = 2.563e-14;
            x2 = 1.018e-13;
            y0 = 2.1112e4;
            y1 = 1.421e4;
            y2 = 9.7791e3;
        }

        Pburst_upper_out = Pburst_upper / (y0 * ((Pburst_upper - x1)
                * (Pburst_upper - x2)) / ((x0 - x1) * (x0 - x2)) + y1
                * ((Pburst_upper - x0) * (Pburst_upper - x2)) / ((x1 - x0)
                * (x1 - x2)) + y2 * ((Pburst_upper - x0) * (Pburst_upper - x1))
                / ((x2 - x0) * (x2 - x1)));
        if (Pburst_upper > 47.57)
            Pburst_upper_out = 0.1152;
        else if (Pburst_upper < 6.093e-15)
            Pburst_upper_out = 2.886e-19;
    }

    if(tx_data_rate == 0)
        Pburst_upper_out = ber_averaged;
    //std::cout << "Pburst: " << Pburst_upper_out << std::endl;

    /* The "size_data" we have to enter in the computation of the PER has to be the total MAC length of the frame
     (including MAC header and excluding PLCP header) */
    /* As a remember: Rts: 160 bits / Cts & Ack: 112 bits */

    packet_error_rate = 1 - pow((1 - Pburst_upper_out), (lengthMPDU));// Payload (in bits) + 224 bits (28bytes MAC overhead)
    //std::cout << "per: " << packet_error_rate << std::endl;
    if (randomNumber < packet_error_rate)
        return false;
    else
        return true;
}
