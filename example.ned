network example  {
    parameters:
	double playgroundSizeX = 8000.0;
	double playgroundSizeY = 8000.0;
	double radius = 1000.0;

	//network parameters
	int numberOfCells = 18;

	int maxNumberOfNeighbours;
	int numberOfMobileStations;

        cell[*].ms[*].initQuadrant = intuniform(0, 5);
        cell[*].ms[*].initPosAlpha = uniform(0, 1);
        cell[*].ms[*].initPosBeta = uniform(0, 1);
        cell[*].ms[*].initPosGamma = uniform(0, 1);
        cell[0].xPos = 1000.0;
        cell[0].yPos = 1000.0;
        cell[1].xPos = 4000.0;
        cell[1].yPos = 1000.0;
        cell[2].xPos = 7000.0;
        cell[2].yPos = 1000.0;
        cell[3].xPos = 1000.0;
        cell[3].yPos = 2732.05080757;
        cell[4].xPos = 4000.0;
        cell[4].yPos = 2732.05080757;
        cell[5].xPos = 7000.0;
        cell[5].yPos = 2732.05080757;
        cell[6].xPos = 1000.0;
        cell[6].yPos = 4464.10161514;
        cell[7].xPos = 4000.0;
        cell[7].yPos = 4464.10161514;
        cell[8].xPos = 7000.0;
        cell[8].yPos = 4464.10161514;
        cell[9].xPos = 1000.0;
        cell[9].yPos = 6196.15242271;
        cell[10].xPos = 4000.0;
        cell[10].yPos = 6196.15242271;
        cell[11].xPos = 7000.0;
        cell[11].yPos = 6196.15242271;
        cell[12].xPos = 2500.0;
        cell[12].yPos = 1866.02540378;
        cell[13].xPos = 5500.0;
        cell[13].yPos = 1866.02540378;
        cell[14].xPos = 2500.0;
        cell[14].yPos = 3598.07621135;
        cell[15].xPos = 5500.0;
        cell[15].yPos = 3598.07621135;
        cell[16].xPos = 2500.0;
        cell[16].yPos = 5330.12701892;
        cell[17].xPos = 5500.0;
        cell[17].yPos = 5330.12701892;
	@display("bgb=$playgroundSizeX,$playgroundSizeY");

    submodules:
	cell[numberOfCells]: LteCell  {
		parameters:
			numberOfMobileStations = numberOfMobileStations;
			bsId = index;
			radius = radius;
		gates:
			fromCell[maxNumberOfNeighbours];
			toCell[maxNumberOfNeighbours];

	}	

    connections allowunconnected:
        cell[0].toCell[0] --> cell[3].fromCell[0];
        cell[0].toCell[1] --> cell[12].fromCell[0];
        cell[1].toCell[0] --> cell[4].fromCell[0];
        cell[1].toCell[1] --> cell[12].fromCell[1];
        cell[1].toCell[2] --> cell[13].fromCell[0];
        cell[2].toCell[0] --> cell[5].fromCell[0];
        cell[2].toCell[1] --> cell[13].fromCell[1];
        cell[3].toCell[0] --> cell[0].fromCell[0];
        cell[3].toCell[1] --> cell[6].fromCell[0];
        cell[3].toCell[2] --> cell[12].fromCell[2];
        cell[3].toCell[3] --> cell[14].fromCell[0];
        cell[4].toCell[0] --> cell[1].fromCell[0];
        cell[4].toCell[1] --> cell[7].fromCell[0];
        cell[4].toCell[2] --> cell[12].fromCell[3];
        cell[4].toCell[3] --> cell[13].fromCell[2];
        cell[4].toCell[4] --> cell[14].fromCell[1];
        cell[4].toCell[5] --> cell[15].fromCell[0];
        cell[5].toCell[0] --> cell[2].fromCell[0];
        cell[5].toCell[1] --> cell[8].fromCell[0];
        cell[5].toCell[2] --> cell[13].fromCell[3];
        cell[5].toCell[3] --> cell[15].fromCell[1];
        cell[6].toCell[0] --> cell[3].fromCell[1];
        cell[6].toCell[1] --> cell[9].fromCell[0];
        cell[6].toCell[2] --> cell[14].fromCell[2];
        cell[6].toCell[3] --> cell[16].fromCell[0];
        cell[7].toCell[0] --> cell[4].fromCell[1];
        cell[7].toCell[1] --> cell[10].fromCell[0];
        cell[7].toCell[2] --> cell[14].fromCell[3];
        cell[7].toCell[3] --> cell[15].fromCell[2];
        cell[7].toCell[4] --> cell[16].fromCell[1];
        cell[7].toCell[5] --> cell[17].fromCell[0];
        cell[8].toCell[0] --> cell[5].fromCell[1];
        cell[8].toCell[1] --> cell[11].fromCell[0];
        cell[8].toCell[2] --> cell[15].fromCell[3];
        cell[8].toCell[3] --> cell[17].fromCell[1];
        cell[9].toCell[0] --> cell[6].fromCell[1];
        cell[9].toCell[1] --> cell[16].fromCell[2];
        cell[10].toCell[0] --> cell[7].fromCell[1];
        cell[10].toCell[1] --> cell[16].fromCell[3];
        cell[10].toCell[2] --> cell[17].fromCell[2];
        cell[11].toCell[0] --> cell[8].fromCell[1];
        cell[11].toCell[1] --> cell[17].fromCell[3];
        cell[12].toCell[0] --> cell[0].fromCell[1];
        cell[12].toCell[1] --> cell[1].fromCell[1];
        cell[12].toCell[2] --> cell[3].fromCell[2];
        cell[12].toCell[3] --> cell[4].fromCell[2];
        cell[12].toCell[4] --> cell[14].fromCell[4];
        cell[13].toCell[0] --> cell[1].fromCell[2];
        cell[13].toCell[1] --> cell[2].fromCell[1];
        cell[13].toCell[2] --> cell[4].fromCell[3];
        cell[13].toCell[3] --> cell[5].fromCell[2];
        cell[13].toCell[4] --> cell[15].fromCell[4];
        cell[14].toCell[0] --> cell[3].fromCell[3];
        cell[14].toCell[1] --> cell[4].fromCell[4];
        cell[14].toCell[2] --> cell[6].fromCell[2];
        cell[14].toCell[3] --> cell[7].fromCell[2];
        cell[14].toCell[4] --> cell[12].fromCell[4];
        cell[14].toCell[5] --> cell[16].fromCell[4];
        cell[15].toCell[0] --> cell[4].fromCell[5];
        cell[15].toCell[1] --> cell[5].fromCell[3];
        cell[15].toCell[2] --> cell[7].fromCell[3];
        cell[15].toCell[3] --> cell[8].fromCell[2];
        cell[15].toCell[4] --> cell[13].fromCell[4];
        cell[15].toCell[5] --> cell[17].fromCell[4];
        cell[16].toCell[0] --> cell[6].fromCell[3];
        cell[16].toCell[1] --> cell[7].fromCell[4];
        cell[16].toCell[2] --> cell[9].fromCell[1];
        cell[16].toCell[3] --> cell[10].fromCell[1];
        cell[16].toCell[4] --> cell[14].fromCell[5];
        cell[17].toCell[0] --> cell[7].fromCell[5];
        cell[17].toCell[1] --> cell[8].fromCell[3];
        cell[17].toCell[2] --> cell[10].fromCell[2];
        cell[17].toCell[3] --> cell[11].fromCell[1];
        cell[17].toCell[4] --> cell[15].fromCell[5];
}
