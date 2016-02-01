/*
 * mpc.cc
 * 
 * Datastructure for Multiple Multi Path Components
 *
 *  Created on: May 6, 2014
 *      Author: Thomas Prinz
 * 
 * Contains parts of the COST 2100 Implementation by Meifang Zhu (meifang.zhu@eit.lth.se)
 */
 
#include "mpc.h"

MPC::MPC(const mat &mpc_bs_pos, const mat &mpc_ms_pos, int N_MPC, vec delay_add){
    MPC_BS_pos = mpc_bs_pos;
    MPC_MS_pos = mpc_ms_pos;
    DMC_delay_add = delay_add;
    
    cvec P_S = to_cvec(randn(N_MPC), randn(N_MPC));      // complex Rayleigh distribution
    MPC_amp = P_S / abs(sum(P_S));                       // normalized attenuation of each MPC's complex Rayleigh amplitude
}

MPC::MPC(const mat &mpc_bs_pos, const mat &mpc_ms_pos, int N_MPC, cvec &c_amp){
    MPC_BS_pos = mpc_bs_pos;
    MPC_MS_pos = mpc_ms_pos;
    
    MPC_amp = c_amp;
}
