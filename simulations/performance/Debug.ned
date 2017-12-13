//Namespace(c=1, chnIni=True, lp=1, n='Debug', ned=True, parIni=True, r=30.0, x=100.0, y=100.0)
network Debug  {
	parameters:
		double playgroundSizeX = 100.0;
		double playgroundSizeY = 100.0;
		double radius = 15.0;

		//network parameters
		int numberOfCells = 5;

		int maxNumberOfNeighbours;
		int numberOfMobileStations;

		cell[*].ms[*].initQuadrant = intuniform(0, 5);
		cell[*].ms[*].initPosAlpha = uniform(0, 1);
		cell[*].ms[*].initPosBeta = uniform(0, 1);
		cell[*].ms[*].initPosGamma = uniform(0, 1);
		cell[0].xPos = 16.0;
		cell[0].yPos = 16.0;
		cell[1].xPos = 48.0;
		cell[1].yPos = 16.0;
		cell[2].xPos = 80.0;
		cell[2].yPos = 16.0;
		cell[3].xPos = 24.0;
		cell[3].yPos = 48.0;
		cell[4].xPos = 64.0;
		cell[4].yPos = 48.0;
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
				fromCell[numberOfCells];
				toCell[numberOfCells];
			}	
		connections allowunconnected:
        // Connect Cells with each other
        for i=0..numberOfCells-1, for j=0..numberOfCells-1, if i!=j {
        	cell[i].toCell[j] --> cell[j].fromCell[i];
        }
}
