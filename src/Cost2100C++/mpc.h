/*
 * mpc.h
 * 
 * Datastructure for Multiple Multi Path Components
 *
 *  Created on: May 5, 2014
 *      Author: Thomas Prinz
 * 
 * Contains parts of the COST 2100 Implementation by Meifang Zhu (meifang.zhu@eit.lth.se)
 */

#ifndef MPC_H
#define MPC_H

#include <itpp/itbase.h>

using namespace itpp;

// Side type: BS Side, MS Side or Both BS/MS Sides
enum SIDE_TYPE {
	BS, MS, BS_MS
};

/*!
 * \brief MPC(Multiple Paths Component) Class
 *
 * MPC Class stands for a set of MPCs occur in one cluster.
 *
 * Each MPC has positions with reference to both BS/MS sides, and its amplitude is calculated as Rayleigh distributed.
 */
class MPC {
public:
	//! Empty constructor
	MPC() {}
	//! Default constructor
	MPC(const mat &mpc_bs_pos, const mat &mpc_ms_pos, int N_MPC, vec delay);
	//! Default constructor
	MPC(const mat &mpc_bs_pos, const mat &mpc_ms_pos, int N_MPC, cvec &c_amp);
	//! Destructor
	virtual ~MPC() {}

	//! Update single MPC position to BS/MS side
	void update_MPC_pos(const SIDE_TYPE sd_type, int mpc_idx, const vec pos);

	//! Get MPC positions to BS/MS side
	mat get_MPC_pos(const SIDE_TYPE sd_type) { switch (sd_type) { case BS: return MPC_BS_pos; case MS: return MPC_MS_pos; default: return mat(0); }	}
	//! Get MPC amplitudes
	cvec get_MPC_amplitude() { return MPC_amp; }
	//! Get DMC delay_add
	vec get_DMC_delay_add() { return DMC_delay_add; }

protected:
	mat MPC_BS_pos;							//!< MPCs positions to BS side
	mat MPC_MS_pos;							//!< MPCs Generatepositions to MS side
	cvec MPC_amp;							//!< MPCs amplitudes
	vec DMC_delay_add;						//!< Additional delay for dense multipaths
};

#endif // MPC_H
