simple App  {
	parameters:
        int packetLength = default(1500b) @unit(b);
        int msId;
        int bsId;
        volatile double packetDistribution;
        double initOffset @unit(s);
        double tti @unit(s);
        double runtime @unit(s);
        double epsilon @unit(s);
        int iconPosX = default(0); //just used in ned
        
		@display("i=block/app2");
		
	gates:
		input fromMac;
		output toMac;
}
