// A cell's central transceiver
//
// This module is the core of each LteCell, serving as a transceiver for 
// non-D2D traffic between it's local MobileStations and providing for 
// information exchange with neighbouring cells. Any cell always has exactly 
// <strong>one</strong> BaseStation. 
//
// <h3>Submodules</h3>
// A BaseStation's submodules fullfill the following tasks:
//
// <p><strong>mac:</strong> The <em>BsMac</em> submodule is responsible for 
// managing the message and packet flow inside the BaseStation module as well 
// as inside the entire cell and with neighbouring cells. It is connected to 
// all submodules and routes their messages to their intended destinations. 
// It also holds the BaseStation's position and manages it's packet queues.
// To facilitate transmission/forwarding of local simulation traffic packets,
// it is also connected to the local cell's Scheduler. 
//
// <p><strong>phy:</strong> The <em>BsPhy</em> submodule is only concerned with 
// forwarding simulation traffic packets either from the BsMac to MobileStation 
// recipients via their MsChannel or those arriving from the BsChannel to the 
// BsMac.
//
// <p><strong>dataChn:</strong> The <em>BsChannel</em> submodule receives all 
// simulation traffic KoiData packets directed to the BaseStation. It computes 
// SINR values with the METISChannel class and determines whether a packet has 
// been received successfully, and if so, forwards it to the BsMac. There is 
// one channel per MobileStation in the local cell.
module BaseStation
{
	parameters:
		// The number of resource blocks for uplink transmissions
		int upResourceBlocks;
		// The number of resource blocks for downlink transmissions
		int downResourceBlocks;
		// X coordinate of the BaseStation's position
		double xPos;
		// Y coordinate of the BaseStation's position
		double yPos;
		// The number of MobileStations in the local cell
		int numberOfMobileStations;
		// The number of neighbouring LteCells
		int maxNumberOfNeighbours;
		// The id of the BaseStation. This id is unique.
		int bsId;
		// The transmission time interval
		double tti @unit(s);
		// Bandwidth available to each resource block
		double bandwidthPerRB @unit(Hz);

		@display("bgb=300,300;i=device/antennatower");

	gates:
		// For simulation traffic packet transmissions with local MobileStations.
		// The gate index corresponds to the MobileStation id.
		input fromMs[];
		// For simulation traffic packet transmissions with local MobileStations.
		// The gate index corresponds to the MobileStation id.
		output toMs[];

		// For information exchange with neighbouring cells. 
		// <strong>Note:</strong> Gate indices do <em>NOT</em> correspond to 
		// BaseStation ids!
		input fromCell[];
		// For information exchange with neighbouring cells. 
		// <strong>Note:</strong> Gate indices do <em>NOT</em> correspond to 
		// BaseStation ids!
		output toCell[];

		// For direct information exchange with MsMacs. The gate index 
		// corresponds to the MobileStation id.
		input fromMsMac[];
		// For direct information exchange with MsMacs. The gate index 
		// corresponds to the MobileStation id.
		output toMsMac[];

		// Communication with Scheduler for transmission schedules
		input fromScheduler;
		// Communication with Scheduler for transmission schedules
		output toScheduler;
    	
	submodules:
   		mac: BsMac  {
   			parameters:
   				xPos = xPos;
   				yPos = yPos;
   				tti = tti;
   				bsId = bsId;
					numberOfMobileStations = numberOfMobileStations;
   				@display("p=220,200");
   			gates:
   				fromCell[maxNumberOfNeighbours];
   				toCell[maxNumberOfNeighbours];
   				fromMsMac[numberOfMobileStations];
   				toMsMac[numberOfMobileStations];
   				toBsChannel[numberOfMobileStations];
   		}
   		phy: BsPhy  {
   			parameters:
					numberOfMobileStations = numberOfMobileStations;
					@display("p=220,100");
   			gates:
   				fromChannel[numberOfMobileStations];
   				toChannel[numberOfMobileStations];
   		}
   		dataChn[numberOfMobileStations]: BsChannel  {
   			parameters:
   				xPos = xPos;
   				yPos = yPos;
   				bsId = bsId;
   				numberOfMobileStations = numberOfMobileStations;
					upResourceBlocks = upResourceBlocks;
					downResourceBlocks = downResourceBlocks;
					bandwidthPerRB = bandwidthPerRB;
				@display("p=100,200");
   		}
   		
   	connections allowunconnected:
   		phy.toMac --> mac.fromPhy;
   		mac.toPhy --> phy.fromMac;

		mac.toScheduler --> toScheduler;
		fromScheduler --> mac.fromScheduler;
   		
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
