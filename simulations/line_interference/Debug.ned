//Namespace(c=1, chnIni=True, lp=1, n='Debug', ned=True, parIni=True, r=30.0, x=100.0, y=100.0)
network Debug  {
	parameters:
		double playgroundSizeX = 900.0;
		double playgroundSizeY = 900.0;
		double radius = 250.0;

		//network parameters
		int numberOfCells = 2;

		int maxNumberOfNeighbours;
		int numberOfMobileStations;

		cell[*].ms[*].initQuadrant = intuniform(0, 5);
		cell[*].ms[*].initPosAlpha = uniform(0, 1);
		cell[*].ms[*].initPosBeta = uniform(0, 1);
		cell[*].ms[*].initPosGamma = uniform(0, 1);
		cell[0].xPos = 50.0;
		cell[0].yPos = 90.0;
		cell[1].xPos = 550.0;
		cell[1].yPos = 90.0;
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
        cell[0].toCell[0] --> cell[1].fromCell[0];
        cell[1].toCell[0] --> cell[0].fromCell[1];
}
