//Namespace(c=1, chnIni=True, lp=1, n='Debug', ned=True, parIni=True, r=30.0, x=100.0, y=100.0)
network Debug  {
	parameters:
		double playgroundSizeX = 100.0;
		double playgroundSizeY = 100.0;
		double radius = 15.0;

		//network parameters
		int numberOfCells = 1;

		int maxNumberOfNeighbours;
		int numberOfMobileStations;

		cell[0].xPos = 30.0;
		cell[0].yPos = 30.0;
		@display("bgb=$playgroundSizeX,$playgroundSizeY");

submodules:
		cell[numberOfCells]: LteCell  {
			parameters:
				numberOfMobileStations = numberOfMobileStations;
				bsId = index;
				radius = radius;
				playgroundSizeX = playgroundSizeX;
				playgroundSizeY = playgroundSizeY;
			gates:
				fromCell[maxNumberOfNeighbours];
				toCell[maxNumberOfNeighbours];
			}	
		connections allowunconnected:
}
