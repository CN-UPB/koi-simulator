module LteCell  {
	parameters:
		int maxNumberOfNeighbours;
		int numberOfMobileStations;
		double radius;
		int currentChannel;
		double xPos;
		double yPos;
		double playgroundSizeX = 2 * radius;
		double playgroundSizeY = 2 * radius;
		double streamSchedPeriod @unit(s);
		int bsId;

		@display("p=$xPos,$yPos;r=$radius,,black;i=block/network2_l;bgb=$playgroundSizeX,$playgroundSizeY");	
		
	gates:
		input fromCell[];
		output toCell[];
		
	submodules:
	    bs: BaseStation  {
			parameters:
				xPos = xPos;
				yPos = yPos;
				bsId = bsId;
				currentChannel = currentChannel;
				numberOfMobileStations = numberOfMobileStations;
				
			gates:
				fromMs[numberOfMobileStations];
				toMs[numberOfMobileStations];

				fromCell[maxNumberOfNeighbours];
				toCell[maxNumberOfNeighbours];

				fromMsMac[numberOfMobileStations];
				toMsMac[numberOfMobileStations];
		}
		ms[numberOfMobileStations]: MobileStation  {
			parameters:
				bsId = bsId;
				msId = index;
				initBsXPos = xPos;
				initBsYPos = yPos;
				radius = radius;
				currentChannel = currentChannel;
		}
	    sched: Scheduler {
	   	parameters:
			numberOfMobileStations = numberOfMobileStations;
			streamSchedPeriod = streamSchedPeriod;
		}

	connections allowunconnected:
		//connect BaseStation with MobileStations
		for i=0..numberOfMobileStations-1  {
			bs.toMs[i] --> ms[i].fromBs;
			ms[i].toBs --> bs.fromMs[i];
		}

		//connect BaseStation with other Cells
		for i=0..(maxNumberOfNeighbours-1)  {
			bs.toCell[i] --> toCell[i];
			bs.fromCell[i] <-- fromCell[i];
		}
        
		//for direct schedule exchange of MS/BS
		for i=0..numberOfMobileStations-1  {
			bs.toMsMac[i] --> ms[i].fromBsMac;
			ms[i].toBsMac --> bs.fromMsMac[i];
		}
		// Access to the Streaming scheduler module for MS
		for i=0..numberOfMobileStations-1  {
			ms[i].toScheduler --> sched.fromMs[i];
			sched.toMs[i] --> ms[i].fromScheduler;
		}
		// Access to the Streaming scheduler for the BS
		sched.toBs --> bs.fromScheduler;
		bs.toScheduler --> sched.fromBs;
}
