import sys
import argparse
#for only one cell you need to you put the size of x and y smaller than 4*radius
parser = argparse.ArgumentParser(description='Creates a ned file for a specific network of the abstractLTE model.')
parser.add_argument('-x', required=True, type=float,
                   help='the x value of the playground size for the hole network')
parser.add_argument('-y', required=True, type=float,
                   help='the y value of the playground size for the hole network')
parser.add_argument('-r', required=True, type=float,
                   help='the radius of one LTE cell')
parser.add_argument('-lp', required=True, type=int,
                   help='the number of logical partitions for the cell distribution')
parser.add_argument('-n', required=True, type=str,
                   help='name of the network also for the output files')
parser.add_argument('-c', required=True, type=int,
                   help='number of channels')
parser.add_argument('-ned', required=False, type=bool, default=True,
                   help='generate the .ned file for the network')
parser.add_argument('-parIni', required=False, type=bool, default=True,
                   help='generate the partition .ini file for the network')
parser.add_argument('-chnIni', required=False, type=bool, default=True,
                   help='generate the channel .ini file for the network')

args = parser.parse_args()

networkName = args.n
playgroundSizeX = args.x
playgroundSizeY = args.y
channels = args.c
radius = args.r
lp = args.lp

bsPos = []
bsChn = []
#rows without offsets
dX = 3 * radius
dY = (3**0.5) * radius
xOffset = radius
row = 0
column = 0
cellCount = 0
while(1):
    xPos = xOffset + dX * column
    if(xPos > (playgroundSizeX - radius)):
        column = 0
        xPos = xOffset + dX * column
        row += 1
    yPos = radius + dY * row
    if(yPos >= playgroundSizeY - radius):
        break
    bsPos.append((xPos, yPos))
    if(cellCount%2 == 0):
        tmpChn = row % channels
        bsChn.append(tmpChn)
    else:
        tmpChn = (row+1) % channels
        bsChn.append(tmpChn)
    column += 1
    cellCount += 1

#rows with offsets
dX = 3 * radius
dY = (3**0.5) * radius
xOffset = 2.5 * radius
yOffset = (3**0.5)/2 * radius + radius
row = 0
column = 0
while(1):
    xPos = xOffset + dX * column
    if(xPos > (playgroundSizeX - radius)):
        column = 0
        xPos = xOffset + dX * column
        row += 1
    yPos = yOffset + dY * row
    if(yPos >= playgroundSizeY - radius):
        break
    column += 1
    cellCount += 1
    bsPos.append((xPos, yPos))
    tmpChn = (row+3) % channels
    bsChn.append(tmpChn)

numberOfCells = cellCount
print "Number of cells: " + str(int(numberOfCells))

if(args.ned):
    f = open(networkName + ".ned", "w")
    f.write("//" + str(args) + "\n")
    f.write("network " + networkName + "  {\n"
	    "    parameters:\n" 
	    "	double playgroundSizeX = " + str(playgroundSizeX)  + ";\n" 
	    "	double playgroundSizeY = " + str(playgroundSizeY)  + ";\n" 
	    "	double radius = " + str(radius) + ";\n\n" 
	    "	//network parameters\n" 
	    "	int numberOfCells = " + str(numberOfCells) + ";\n\n" 
	    "	int maxNumberOfNeighbours;\n" 
	    "	int numberOfMobileStations;\n\n"
	    "        cell[*].ms[*].initQuadrant = intuniform(0, 5);\n"
	    "        cell[*].ms[*].initPosAlpha = uniform(0, 1);\n"
	    "        cell[*].ms[*].initPosBeta = uniform(0, 1);\n"
	    "        cell[*].ms[*].initPosGamma = uniform(0, 1);\n")

    bsId = 0
    for pos in bsPos:
	f.write("        cell[" + str(bsId) + "].xPos = " + str(pos[0]) + ";\n"
		"        cell[" + str(bsId) + "].yPos = " + str(pos[1]) + ";\n")
	bsId += 1

    f.write("	@display(\"bgb=$playgroundSizeX,$playgroundSizeY\");\n\n" 
	    "    submodules:\n" 
	"   watch: StopWatch;\n"
	    "	cell[numberOfCells]: LteCell  {\n" 
	    "		parameters:\n" 
	    "			numberOfMobileStations = numberOfMobileStations;\n" 
	    "			bsId = index;\n" 
	    "			radius = radius;\n" 
	    "		gates:\n" 
	    "			fromCell[maxNumberOfNeighbours];\n" 
	    "			toCell[maxNumberOfNeighbours];\n\n" 
	    "	}	\n\n" 
	    "    connections allowunconnected:\n")

    fromCounter = []
    for i in range(0, numberOfCells):
	fromCounter.append(0)

    for i in range(0, numberOfCells):
	counter = 0
	for j in range(0, numberOfCells):
	    if(i == j):
		continue
	    p1 = bsPos[i] 
	    p2 = bsPos[j] 
	    dist = (((p1[0]-p2[0])**2) + ((p1[1]-p2[1])**2))**0.5
	    if(dist <= 4 * radius):
		f.write("        cell[" + str(i) + "].toCell[" + str(counter) + "] --> cell[" + str(j) + "].fromCell[" + str(fromCounter[j]) + "];\n")
		counter += 1
		fromCounter[j] += 1
    f.write("}\n")

#create the corresponding networkName-par.ini and networkName-chn.ini
if(args.chnIni):
    f = open(networkName + "-chn.ini", "w")
    f.write("#.ini for the network " + networkName + "\n"
	    "#" + str(args) + "\n"
	    "network=" + networkName + "\n")
    f.write("#Channel distribution\n")
    bsId = 0
    for chn in bsChn:
	f.write("**.cell[" + str(bsId) + "].currentChannel = " + str(chn) + "\n")
	bsId += 1

if(args.parIni):
    f = open(networkName + "-par.ini", "w")
    f.write("#Partitions\n")
    f.write("**.watch.partition-id = 0\n")
    for i in range(0, numberOfCells):
	f.write("**.cell[" + str(i) + "].partition-id = " + str(i % lp) + "\n")
