Topology
^^^^^^^^
The network is comprised of **LTECells**, each of which contains exactly one **BaseStation**, a number of **MobileStations** and a **Scheduler**.

The **BaseStation** behaves as a hub, forwarding the data packets received from a MobileStation to the corresponding destination. The BaseStation is connected to all the other LTECells, even though at the moment the connection between the BaseStations in different cells is used only for traffic information between neighboring cells. 

The **MobileStations** are the source and destination of the traffic, which is generated in form of KoiData packets. Each MobileStation can communicate with any other MobileStation in the same cell forwarding the data through the BaseStation or directly to the MobileStation using the device-to-device functionality (D2D communication).

.. figure:: scenario.png
    :width: 800px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: aligh-center

    Network Topology

The **Scheduler** is responsible for scheduling the traffic in each cell, in a way that the interference can be eliminated.  The transmissions are distributed in different time slots and different frequencies which are represented by the resource blocks. 

The whole simulation time is divided into fixed transmission time intervals (TTI), given as an input parameter.

BaseStation
...........

The **BaseStation** is the core of each LTECell and is responsible for forwarding the traffic for non-D2D communication, received from the source MobileStation to the destination MobileStation. It consists of the following submodules:

1. *BsMac* - the submodule is responsible for managing the messages and the packet flow inside the BaseStation module. It is connected to all the MobileStations and the single Scheduler in its cell. The submodule makes possible to forward the traffic packets to the MobileStation destination. To complete this purpose, the submodule communicates with the Scheduler during each TTI.

2. *BsPhy* - the submodule is in charge of forwarding the simulated traffic packets either from the BsMac to MobileStation destination or those arriving from the BsChannel to the BsMac. 

3. *BsChannel* - receives all the simulated traffic packets directed to the BaseStation. It computes the SINR value, using the channel coefficients calculated in the METISChannel class, taking in consideration the transmission power and the interference by the other mobile station transmitting at the same TTI. In the case, the packet has been received successfully, it will be forwarded to the BsMac, otherwise, it will be dropped. At the moment all the packets are considered to be received successfully.  




MobileStation
.............

A MobileStation is part of only one LTECell. It is the source and the destination of the traffic packets in the simulation. Each MobileStation communicates with any other MobileStation in the same cell forwarding the traffic packets through the BaseStation or using a direct transmission. A mobile station is comprised of the following submodules: 

1. *TrafficGen* - the submodule behaves as a generator and as a sink for the simulation traffic packets. To generate the traffic, a communication table in XML format is used, which defines the number of streams, each representing a traffic process between the source MobileStation and a destination MobileStation. The traffic generated is periodic and the packet size is fixed. After generation, the packets are transmitted to the MobileStation's mac submodule, for storage and transmission. 

2. *MsMac* - it holds information regarding the current position of the MobileStation. At the moment the MobileStations are considered to be static, they do not change their position during the simulation run. The position of each MobileStation is given as a parameter in the simulation. Secondly, it holds all the packets waiting for transmission, according to the communication stream. These packets are transmitted through MobileStation's phy submodule either to the BaseStation or directly to the neighboring MobileStation depending on the type of communication. To determine which packets, from which streams will be transmitted during one TTI, the submodule is connected to the Scheduler. This submodule is also directly connected to the local BaseStation, using this connection for information exchange related to the MobileStation position. 

3. *MsPhy* - the submodule is in charge of forwarding the simulation traffic packets either from the MsMac to MobileStation destinationor Base Station,  or those arriving from the MsChannel to the MsMac. 

4. *MsChannel* - receives all the packets directed to this mobile station. It computes the SINR value, using the channel coefficients calculated in the METISChannel class, taking in consideration the transmission power and the interference by the other Mobile Station transmitting at the same TTI. In the case, the packet has been received successfully, it will be forwarded to the MsPhy, otherwise, it will be dropped. At the moment all the packets are considered to be received successfully. 

Scheduler
.........

The Scheduler is responsible for scheduling the transmissions inside a LTECell in a way, such that two transmissions cannot be scheduled in the same resource block at the same TTI. At the beginning of the simulation, each MobileStation sends to the Scheduler a StreamInfo message carrying information about the communication partners, type of communication and the period between the packets of that stream. Based on this information, the Scheduler schedules the streams into different resource blocks. The Scheduler is comprised of the following submodules: 

1. *StreamScheduler* - the submodule is in charge of assigning streams to resource blocks, each of which represents a different frequency. For this purpose, the simulator uses the Round Robin algorithm.

2. *RBScheduler* - a number of streams can be attached to each resource block, but during a TTI, only packets from one of the streams can be transmitted. To fill this purpose, the submodule uses the Round Robin algorithm.



