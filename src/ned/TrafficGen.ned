// Traffic generator for periodic traffic
//
// This module is responsible for the generation of simulated packet traffic 
// as a submodule of MobileStation compound modules. To generate the traffic,
// the module receives a <em>communication table</em> file in XML format. This 
// file defines any number of streams, each stream representing a traffic 
// generation process between this module's parent MobileStation and any 
// other MobileStation in the current LteCell.
//
// @todo Add documentation page for the XML format for the comm table file 
//
// Packets are generated according to each stream's unique period and then 
// send to the MobileStation's <em>mac</em> for storage and transmission.
//
// At the beginning of each simulation, the module also sends out information 
// messages, informing all interested parties about it's communication streams.
//
// The traffic generator is also intended as the sink for any traffic packets 
// it's MobileStation receives, e.g. for logging/measurement purposes.
//
// The traffic generator only connects to the MobileStation's <em>mac</em>, 
// which routes all of it's packets and messages.
simple TrafficGen {
	parameters:
		// The id of the local cell and BaseStation.
		int bsId;
		// Path to the communication table file with stream definitions
		string commTable = default("comm_table.xml");
		// The deadline for all packets generated here. This is an 
		// offset added to each packet's creation time to arrive 
		// at the actual deadline.
		double deadline @unit(s);
		// This is an initialization offset. Generation of traffic 
		// packets only starts this many seconds after simulation start.
		double initOffset @unit(s);
		// The id of the MobileStation this generator belongs to
		int msId;
		// This is the length in Bits for all packets generated at this 
		// generator
		int packetLength @unit(b);
		// Indicates whether this generator produces strictly periodic 
		// traffic. This is only used to set the correct information 
		// flags on generated packets.
		bool periodicTraffic = default(true);
	gates:
		// For incoming messages and packets from the <em>mac</em>.
		input fromMac;
		// For outgoing messages to the <em>mac</em>.
		output toMac;
}
