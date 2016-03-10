simple TrafficGen {
	parameters:
		int bsId;
		string commPartners;
		double deadline @unit(s);
		double initOffset @unit(s);
		int msId;
		int packetLength @unit(b);
		double period @unit(s);
		bool periodicTraffic;
	gates:
		input fromMac;
		output toMac;
}
