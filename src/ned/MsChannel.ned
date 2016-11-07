// Communication channel handling packet receival for MobileStation.
//
// This submodule of MobileStations is tasked with receiving the simulation 
// traffic KoiData packets from the local BaseStation and neighbouring 
// MobileStations. It uses the channel module of the current cell to compute 
// SINR values for messages send to its MobileStation. With these SINR values,
// it then determines whether the packet was actually received successfully or 
// not. In the former case, the packet is forwarded to the MsPhy module, while 
// the message is simply dropped in the latter case. 
//
// To compute SINR values for received packets, the MsChannel also needs 
// to know which interference occured during the transmission. For this reason,
// it receives TransInfo messages from neighbouring cells, informing it 
// about BaseStations and MobileStations transmitting in the same TTI and 
// thus, if sending on the same resource block, interfering with local 
// transmissions.
//
// TransInfo messages from both, neighbouring MobileStations and BaseStations 
// are necessary, because both, downlink (BS->MS) transmissions in neighbouring 
// cells as well as D2D (MS->MS) transmissions on downlink resource blocks 
// could interfere with packet receival.
simple MsChannel  {
	parameters:
		// The number of neighbouring cells possibly interfering 
		// with local transmissions.
		int maxNumberOfNeighbours;
		// The number of resource blocks for downlink transmissions
		int downResourceBlocks;
		// The number of resource blocks for uplink transmissions
		int upResourceBlocks;
		// The id of the local cell and BaseStation.
		int bsId;
		// The ID of this MobileStation. MobileStation IDs are unique only 
		// within any given cell. Thus, to uniquely identify a MobileStation,
		// it's <tt>bsId</tt> and it's <tt>msId</tt> are needed.
		int msId;
		// A small offset used control message order in cases where ordering 
		// cannot be guaranteed by generation order.
		double epsilon @unit(s);
		// Offset from the simulation start reserved for initialization messages.
		// No traffic is generated/send before it.
		double initOffset @unit(s);
		// Length of single transmission time interval
		double tti @unit(s);
		// Table of EESM values for CQI computation
		string eesm_beta = default("");
		// Bit Level Error Rate table for packet loss computations
		string bler_table = default("");
		// Number of antennas on BaseStation
		int NumBsAntenna = default(8);
		// Number of antennas on MobileStation
		int NumMsAntenna = default(2);
		// Table with min SINR for Modulation and Coding Schemes
		string MSCTable;
		@display("i=block/rxtx");
	gates:
		// Gate for incoming simulation traffic packets from the 
		// local BaseStation.
		input fromBs;
		// Gate to the MsPhy module to forward received packets
		output toPhy;
		// Input from MsPhy module for position exchange messages only
		input fromPhy;

		// For receival of D2D simulation traffic packets from local MS,
		// the gate index corresponds to the MobileStation id.
		input fromMs[];
}
