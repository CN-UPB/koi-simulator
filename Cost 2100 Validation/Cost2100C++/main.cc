#include "cluster.h"
#include "VisibilityRegion.h"
#include "ChannelCalc.h"
#include <iostream>
#include <fstream>
#include <stddef.h>

using namespace std;
using namespace itpp;

int main(){
	cout << "Start!" << endl;
	vector<Cluster> remote_cluster;
	Cluster local_cluster_bs;
	Cluster local_cluster_ms;
	vector<VisibilityRegion> VRs;
	vector<MPC> MPCs;
	ChannelCalc calc;
	
	// Read in Text file here:
	ifstream myfile ("cluster.txt");
	ifstream local_bs ("local_cluster_BS.txt");
	ifstream local_ms ("local_cluster_MS.txt");
	ifstream mpc ("mpc.txt");
	
	int c_id = 0;
	if (myfile.is_open()){
		double cluster_t;
		double VR_POS_x, VR_POS_y;
		double BS_POS_x, BS_POS_y, BS_POS_z;
		double C_BS_POS_x, C_BS_POS_y, C_BS_POS_z;
		double C_MS_POS_x, C_MS_POS_y, C_MS_POS_z;
		double a_c, b_c, h_c;
		double shadow_fading;
		double tau_link_delay;
		
		while(!myfile.eof()){
			myfile >> cluster_t;
			myfile >> VR_POS_x >> VR_POS_y;
			myfile >> BS_POS_x >> BS_POS_y >> BS_POS_z;
			myfile >> C_BS_POS_x >> C_BS_POS_y >> C_BS_POS_z;
			myfile >> C_MS_POS_x >> C_MS_POS_y >> C_MS_POS_z;
			myfile >> a_c >> b_c >> h_c;
			myfile >> shadow_fading;
			myfile >> tau_link_delay;
			//cout << shadow_fading << std::endl;
			
			vec BS_position, cluster_pos_BS, cluster_pos_MS, VR_position, spread;
			VR_position.ins(0,VR_POS_x);
			VR_position.ins(1,VR_POS_y);
			VR_position.ins(2,0);
			
			BS_position.ins(0,BS_POS_x);
			BS_position.ins(1,BS_POS_y);
			BS_position.ins(2,BS_POS_z);
			
			cluster_pos_BS.ins(0,C_BS_POS_x);
			cluster_pos_BS.ins(1,C_BS_POS_y);
			cluster_pos_BS.ins(2,C_BS_POS_z);
			
			cluster_pos_MS.ins(0,C_MS_POS_x);
			cluster_pos_MS.ins(1,C_MS_POS_y);
			cluster_pos_MS.ins(2,C_MS_POS_z);
			
			spread.ins(0,a_c);
			spread.ins(1,b_c);
			spread.ins(2,h_c);
			spread.ins(3,a_c);
			spread.ins(4,b_c);
			spread.ins(5,h_c);
			
			vec BS_angle = Cart_to_Sph(BS_position - cluster_pos_BS).get(0, 1);
			vec MS_angle = Cart_to_Sph(VR_position - cluster_pos_MS).get(0, 1);
			vec angle_final = zeros(4);
			angle_final.set(0,BS_angle(0));
			angle_final.set(1,BS_angle(1));
			angle_final.set(2,MS_angle(0));
			angle_final.set(3,MS_angle(1));
			
			VisibilityRegion new_vr;
			new_vr.x = VR_POS_x;
			new_vr.y = VR_POS_y;
			new_vr.bsIds.push_back(1);
			new_vr.clusterId = c_id;
			VRs.push_back(new_vr);
			
			// Read MPCs
			mat MPC_POS_BS, MPC_POS_MS;
			cvec MPC_amp;
			vec temp = zeros(3);
			double x,y,z;
			double t;
			
			// 27 MPC Hardcoded, because Im lazy.
			for(int i = 0; i < 27; i++){
				mpc >> x >> y >> z;
				temp.set(0,x);
				temp.set(1,y);
				temp.set(2,z);
				MPC_POS_BS.append_row(temp);
			}
			
			for(int i = 0; i < 27; i++){
				mpc >> x >> y >> z;
				temp.set(0,x);
				temp.set(1,y);
				temp.set(2,z);
				MPC_POS_MS.append_row(temp);
			}
			
			double real[27], imag[27];
			
			for(int i = 0; i < 27; i++){
				mpc >> t;
				real[i] = t;
			}
			
			for(int i = 0; i < 27; i++){
				mpc >> t;
				imag[i] = t;
			}
			
			for(int i = 0; i < 27; i++){
				MPC_amp.ins(i,complex<double>(real[i],imag[i]));
			}
			
			MPC mpc(MPC_POS_BS, MPC_POS_MS, 27, MPC_amp);
			
			remote_cluster.push_back(Cluster(static_cast<CLUSTER_TYPE> (cluster_t), cluster_pos_BS, cluster_pos_MS, spread, shadow_fading, tau_link_delay, angle_final, c_id, mpc));
			c_id++;
		}
		
		myfile.close();
	}else 
	cout << "Unable to open file"; 
	
	//cout << "local_bs" << endl;
	if (local_bs.is_open()){
		double VR_POS_x, VR_POS_y, VR_POS_z;
		double BS_POS_x, BS_POS_y, BS_POS_z;
		double C_BS_POS_x, C_BS_POS_y, C_BS_POS_z;
		double C_MS_POS_x, C_MS_POS_y, C_MS_POS_z;
		double a_c, b_c, h_c;
		double shadow_fading;
		double tau_link_delay;
		vec BS_position, cluster_pos_BS, cluster_pos_MS, VR_position, spread;
		
		local_bs >> VR_POS_x >> VR_POS_y >> VR_POS_z;
		local_bs >> BS_POS_x >> BS_POS_y >> BS_POS_z;
		local_bs >> C_BS_POS_x >> C_BS_POS_y >> C_BS_POS_z;
		local_bs >> C_MS_POS_x >> C_MS_POS_y >> C_MS_POS_z;
		local_bs >> a_c >> b_c >> h_c;
		local_bs >> shadow_fading;
		local_bs >> tau_link_delay;
		
		spread.ins(0,a_c);
		spread.ins(1,b_c);
		spread.ins(2,h_c);
		spread.ins(3,a_c);
		spread.ins(4,b_c);
		spread.ins(5,h_c);
		
		vec angle = zeros(4);
		
		BS_position.ins(0,BS_POS_x);
		BS_position.ins(1,BS_POS_y);
		BS_position.ins(2,BS_POS_z);
		
		// Read MPCs
		mat MPC_POS_BS, MPC_POS_MS;
		cvec MPC_amp;
		vec temp = zeros(3);
		double x,y,z;
		double a1,a2,a3,a4,a5;
		double b1,b2,b3,b4,b5;
		
		// 5 MPC Hardcoded, because Im lazy.
		for(int i = 0; i < 27; i++){
			mpc >> x >> y >> z;
			temp.set(0,x);
			temp.set(1,y);
			temp.set(2,z);
			MPC_POS_BS.append_row(temp);
		}
		
		for(int i = 0; i < 27; i++){
			mpc >> x >> y >> z;
			temp.set(0,x);
			temp.set(1,y);
			temp.set(2,z);
			MPC_POS_MS.append_row(temp);
		}
		
		mpc >> a1 >> a2 >> a3 >> a4 >> a5;
		mpc >> b1 >> b2 >> b3 >> b4 >> b5;
		MPC_amp.ins(0,complex<double>(a1,b1));
		MPC_amp.ins(1,complex<double>(a2,b2));
		MPC_amp.ins(2,complex<double>(a3,b3));
		MPC_amp.ins(3,complex<double>(a4,b4));
		MPC_amp.ins(4,complex<double>(a5,b5));
		//cout << MPC_amp << endl;
		
		MPC mpc(MPC_POS_BS, MPC_POS_MS, 27, MPC_amp);
		
		local_cluster_bs = Cluster(LOCAL_AT_BS, BS_position, BS_position, spread, shadow_fading, tau_link_delay, angle, -1, mpc);		
	}
	
	//cout << "local_ms" << endl;
	if (local_ms.is_open()){
		double VR_POS_x, VR_POS_y, VR_POS_z;
		double BS_POS_x, BS_POS_y, BS_POS_z;
		double C_BS_POS_x, C_BS_POS_y, C_BS_POS_z;
		double C_MS_POS_x, C_MS_POS_y, C_MS_POS_z;
		double a_c, b_c, h_c;
		double shadow_fading;
		double tau_link_delay;
		vec BS_position, cluster_pos_BS, cluster_pos_MS, VR_position, spread;
		
		local_ms >> VR_POS_x >> VR_POS_y >> VR_POS_z;
		local_ms >> BS_POS_x >> BS_POS_y >> BS_POS_z;
		local_ms >> C_BS_POS_x >> C_BS_POS_y >> C_BS_POS_z;
		local_ms >> C_MS_POS_x >> C_MS_POS_y >> C_MS_POS_z;
		local_ms >> a_c >> b_c >> h_c;
		local_ms >> shadow_fading;
		local_ms >> tau_link_delay;
		
		spread.ins(0,a_c);
		spread.ins(1,b_c);
		spread.ins(2,h_c);
		spread.ins(3,a_c);
		spread.ins(4,b_c);
		spread.ins(5,h_c);
		
		vec angle = zeros(4);
		
		BS_position.ins(0,BS_POS_x);
		BS_position.ins(1,BS_POS_y);
		BS_position.ins(2,BS_POS_z);
		
		// Read MPCs
		mat MPC_POS_BS, MPC_POS_MS;
		cvec MPC_amp;
		vec temp = zeros(3);
		double x,y,z;
		double a1,a2,a3,a4,a5;
		double b1,b2,b3,b4,b5;
		
		// 5 MPC Hardcoded, because Im lazy.
		for(int i = 0; i < 27; i++){
			mpc >> x >> y >> z;
			temp.set(0,x);
			temp.set(1,y);
			temp.set(2,z);
			MPC_POS_BS.append_row(temp);
		}
		
		for(int i = 0; i < 27; i++){
			mpc >> x >> y >> z;
			temp.set(0,x);
			temp.set(1,y);
			temp.set(2,z);
			MPC_POS_MS.append_row(temp);
		}
		
		mpc >> a1 >> a2 >> a3 >> a4 >> a5;
		mpc >> b1 >> b2 >> b3 >> b4 >> b5;
		MPC_amp.ins(0,complex<double>(a1,b1));
		MPC_amp.ins(1,complex<double>(a2,b2));
		MPC_amp.ins(2,complex<double>(a3,b3));
		MPC_amp.ins(3,complex<double>(a4,b4));
		MPC_amp.ins(4,complex<double>(a5,b5));
		//cout << MPC_amp << endl;
		
		MPC mpc(MPC_POS_BS, MPC_POS_MS, 27, MPC_amp);

		local_cluster_ms = Cluster(LOCAL_AT_MS, BS_position, BS_position, spread, shadow_fading, tau_link_delay, angle, -1, mpc);
	}
	
	// Hardcoded
	VisibilityRegion VR_LOS;
	ifstream read_vr("VR.txt");
	double vr_in;
	read_vr >> vr_in;
	VR_LOS.x = vr_in;
	read_vr >> vr_in;
	VR_LOS.y = vr_in;
	
	vec bs_position = zeros(3);
	vec ms_position = zeros(3);
	
	ms_position.set(0,100.0);
	ms_position.set(1,-200.0);
	ms_position.set(2,0.0);
	
	Position bs_pos;
	bs_pos.x = 0;
	bs_pos.y = 0;
	
	std::map <int,Position> map;
	map[1] = bs_pos;
	
	vector<cmat> channel = calc.H(local_cluster_bs, local_cluster_ms, remote_cluster, VRs, 27, 0, ms_position, 1, map, VR_LOS);
	
	ofstream final_result("../Results_Validation/MIMO_Matrix.txt");
	for(uint i = 0; i < channel.size(); i++){
		final_result.precision(16);
		if(channel.at(i).get(0,0).real() > 0){
			final_result << std::fixed << "   " << channel.at(i).get(0,0);
		}else{
			final_result << std::fixed << "  " << channel.at(i).get(0,0);
		}
		if(channel.at(i).get(0,1).real() > 0){
			final_result << std::fixed << "   " << channel.at(i).get(0,1) << endl;
		}else{
			final_result << std::fixed << "  " << channel.at(i).get(0,1) << endl;
		}
		if(channel.at(i).get(1,0).real() > 0){
			final_result << std::fixed << "   " << channel.at(i).get(1,0);
		}else{
			final_result << std::fixed << "  " << channel.at(i).get(1,0);
		}
		if(channel.at(i).get(1,1).real() > 0){
			final_result << std::fixed << "   " << channel.at(i).get(1,1) << endl << endl;
		}else{
			final_result << std::fixed << "  " << channel.at(i).get(1,1) << endl << endl;
		}
	}
	
	cout << channel[0] << endl;
	
	return 0;
}
