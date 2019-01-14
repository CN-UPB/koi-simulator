Simulation Run
^^^^^^^^^^^^^^
In addition to the simulated KoIData packets, other messages are generated, which are not consider as traffic, carrying information needed during simulation. The flow of the data into simulation is given in the following sections together with the messages generated for making possible the communication between the MobileStations. 

Initialization Phase
....................

 An initialization phase ( *initOffset*) is used to disseminate information among all modules available in the simulation.

- Firstly, the TrafficGen module in MobileStation sets up the corresponding streams and then forward *StreamInfo* messages to the MsMac, which then sends the message to the local cell's Scheduler for the assignment of resource blocks.

- Each of the MobileStation sends its own position to the local BaseStation using the *PositionExchange* message.

- The BaseStation, after receiving the *PositionExchange* messages from all the MobileStations in the cell, adds its own position and sends a *BsMsPosition* message to all the neighboring cells.  Furthermore, those position messages as well as all such messages received from other cells are forwarded to the BsChannel instance with index 0.

- The BsChannel with index 0 is allways the channel which initializes the local cell's METISChannel instance, as soon as *BsMsPositions* messages have been received from all cells in the simulation. Then, the initial coefficient tables for the METIS model can be computed.

Traffic and Information Flow
..........................................

- After the initialization phase, the TrafficGen module in each MobileStation, generates KoiData packets, based on the corresponding period for each stream and then forwards the packets to the MsMac module, in which the packets are organized into queues based on their stream id.

- At the beginning of each TTI, each of the MobileStations and the BaseStation, sends a *StreamTransReq* message, for each of the streams it handles to the local Scheduler.

- At the beginning of each TTI, each local Scheduler waits for an epsilon time, for receiving all the *StreamTransReq* messages sent to this module. Once the Scheduler receives the *StreamTransReq* messages, it organizes them into a list and sends a *TransReqList* to the RBScheduler, containing this list. Each of the lists contains the *StreamTransReq* messages for all the streams handled by each RBScheduler.

- Once each of the RBScheduler has received a *TransReqList*, it uses the round robin algorithm to choose one of the packets from all the streams it is responsible. Right after it chooses a packet, the RBScheduler notes it in a *StreamTransSched* message, which is sent back to the StreamScheduler.

- The StreamScheduler forwards the *StreamTransSched* message back to the MobileStation or BaseStation from where this packet has been generated.

- After the MobileStation or BaseStation, receives the *StreamTransSched* message, it takes KoiData packets from the corresponding stream based on the available bandwidth, packs them into a *DataPacketBundle* and forwards it to the MsPhy/BsPhy respectively.

- In addition to the *DataPacketBundle*, the corresponding MobileStation or BaseStation, generates also a *TransInfo* message. The purpose of this message is to inform all the other MobileStations for this transmission happening during this TTI. This message is received by the BsMac and is forwarded to the BsChannel and the local MobileStations, which then forwards the message from the MsMac to the MsChannel. The BsChannel and MsChannel use this information to calculate the interference during this TTI and the SINR for the received packets.

- In the case when the transmission takes place at the MobileStation, once the *DataPacketBundle* reaches the MsPhy, it is forwarded to the MsChannel of the destination in the case of a D2D communication or otherwise it is send to the BsChannel of the local BaseStation.

- In the case when the transmission takes place at the BaseStation, once the *DataPacketBundle* reaches the BsPhy, it is forwarded to the MsChannel of the destination MobileStation.

- Once the DataPacketBundle reaches the BsChannal/MsChannel, the corresponding SINR value is calculated. Based on the value obtained from the calculation, if the packet is received successfully (at the moment all the packets are received successfully), it is forwarded  to the MsPhy or BsPhy, otherwise, the packet is dropped.

- In the case of the MobileStation, the received packet is forwarded to the MsMac and then to the TrafficGen, which represents the sink for the packets in this simulator.

- In the case of the BaseStation, the received packet is forwarded to the BsMac , where all the KoiData packets are sorted into queues accrding to their stream id. In a later TTI, based on the Scheduler decision, the KoiData packets are forwarded to the destination.

