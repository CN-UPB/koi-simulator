simple TrafficGen {
	parameters:
		int bsId;
		string commTable = default("comm_table.xml");
		double deadline @unit(s);
		double initOffset @unit(s);
		int msId;
		int packetLength @unit(b);
		bool periodicTraffic;
	gates:
		input fromMac;
		output toMac;
}
