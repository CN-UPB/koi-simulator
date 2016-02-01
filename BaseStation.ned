module BaseStation
{
    parameters:
        double xPos;
        double yPos;
        int numberOfMobileStations;
        int maxNumberOfNeighbours;
        int bsId;
        int playgroundSizeX = 550;
        int halfPlaygroundSizeX = playgroundSizeX / 2;
        int playgroundSizeY = 300;
        int currentChannel;
        
        //TTI - Transmission Time Interval
        double tti @unit(s);

        @display("p=$radius,$radius;bgb=$playgroundSizeX,$playgroundSizeY;i=device/antennatower");

    gates:
    	//for packets through DataChannel
    	input fromMs[];
    	output toMs[];
    	
    	//MS position exchange
    	input fromCell[];
    	output toCell[];
    	
    	//schedule/position exchange
    	input fromMsMac[];
    	output toMsMac[];
    	
   	submodules:
   		app[numberOfMobileStations]: App  {
   			parameters:
   				bsId = bsId;
   				msId = index;
   				iconPosX = 50 * (index + 1);
        		@display("p=$iconPosX,50");
   		}
   		mac: BsMac  {
   			parameters:
   				xPos = xPos;
   				yPos = yPos;
   				tti = tti;
   				bsId = bsId;
                currentChannel = currentChannel;
                numberOfMobileStations = numberOfMobileStations;
   				@display("p=$halfPlaygroundSizeX,120");
   				
   			gates:
   				fromApp[numberOfMobileStations];
   				toApp[numberOfMobileStations];
   				fromCell[maxNumberOfNeighbours];
   				toCell[maxNumberOfNeighbours];
   				fromMsMac[numberOfMobileStations];
   				toMsMac[numberOfMobileStations];
   				toBsChannel[numberOfMobileStations];
   		}
   		phy: BsPhy  {
   			parameters:
				numberOfMobileStations = numberOfMobileStations;
   				@display("p=$halfPlaygroundSizeX,190");
   			gates:
   				fromChannel[numberOfMobileStations];
   				toChannel[numberOfMobileStations];
   		}
   		dataChn[numberOfMobileStations]: BsChannel  {
   			parameters:
   				xPos = xPos;
   				yPos = yPos;
   				bsId = bsId;
   				iconPosX = 50 * (index + 1);
   				numberOfMobileStations = numberOfMobileStations;
        		@display("p=$iconPosX,250");
   		}
   		
   	connections allowunconnected:
   		phy.toMac --> mac.fromPhy;
   		mac.toPhy --> phy.fromMac;
   		
   		//connect BsMac with app for every MS
   		for i=0..numberOfMobileStations-1  {
   			mac.toApp[i] --> app[i].fromMac;
   			app[i].toMac --> mac.fromApp[i];
   		}
   		
   		//connect BsPhy with the BsChannels
   		for i=0..numberOfMobileStations-1  {
   			fromMs[i] --> dataChn[i].fromMs;
   			dataChn[i].toPhy --> phy.fromChannel[i];
   			phy.toChannel[i] --> toMs[i];	
   		}
   		
   		//connect BsMac with the other cells
   		for i=0..maxNumberOfNeighbours-1  {
   			mac.toCell[i] --> toCell[i];
   			fromCell[i] --> mac.fromCell[i];
   		}
   		
   		//for direct schedule exchange with the BS
   		for i=0..numberOfMobileStations-1  {
   			mac.toMsMac[i] --> toMsMac[i];
   			fromMsMac[i] --> mac.fromMsMac[i];
   		}
   		
   		//for ms position exchange with the BsChannel
   		for i=0..numberOfMobileStations-1  {
   			mac.toBsChannel[i] --> dataChn[i].fromBsMac;
   		}
}
