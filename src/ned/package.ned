// @titlepage
// <h1>Koi Simulation</h1>
//
// This simulation model is intended to simulate wireless commnucation 
// scheduling algorithms for the Koi project. 
//
// In this documentation, simulated packet traffic is made up of KoiData 
// packets. Only those messages which simulate actual communication traffic 
// are called <em>packet</em>. Any other messages, for control or 
// information, are called messages. They are needed for the simulation, but do 
// not represent traffic in the simulated model.
//
// <h2>Topology</h2>
// The network topology consists of at least one LteCell. Those cells, in turn,
// contain at least one MobileStation and exactly one BaseStation. 
//
// The BaseStation is the local transceiver hub for local MobileStations, 
// delivering their KoiData simulation traffic packets to the intended 
// recipient. 
//
// The MobileStation modules in each LteCell are the communication partners 
// of the model. They produce and consume traffic packets and send them 
// to each other as defined by streams. At the moment, their position is 
// defined at the beginning of each simulation run and does not change over 
// the course of a run.
//
// The Scheduler module in each LteCell is tasked with ensuring that no 
// local transmissions interfere with each other. This goal is achieved by 
// distributing transmissions over different time slots and resource blocks. 
// In our model, resource blocks represent different frequencies.
//
// <h2>\Traffic</h2>
// \Traffic, in the form of KoiData packets, is generated in the MobileStation's 
// TrafficGen submodule. The traffic is defined by streams, which consist of 
// a sender and a recipient, both MobileStations, and a fixed period, 
// which defines the time between two packets of that stream. Additionally, 
// every stream is either send via the local BaseStation, or directly to 
// the recipient via D2D transmission. Streams are defined in a XML file. The 
// file format description can be found <a href="streams.html">here</a>.
//
// <strong>Note:</strong> \Traffic is <em>only</em> generated in the 
// TrafficGen module.
//
// \Traffic only occurs between pairs of MobileStations as origins and 
// sinks. 
//
// <h2>Simulation Run</h2>
// This section will give a short overview over a typical simulation run, the 
// messages exchanged and the modules involved.
//
// <h3>Initialization Phase</h3>
// At the start of each simulation, an initialization phase is used to 
// disseminate information among all modules which is only available after 
// a simulation has been set up. 
//
// Right after the TrafficGen modules in each MobileStation have set up 
// their streams, they forward StreamInfo messages to the MsMac, which then 
// sends them to the local cell's Scheduler for assignment of resource blocks.
//
// Also right at simulation start, all MobileStation modules send their 
// positions to their local BaseStation via a PositionExchange message. The 
// BaseStation, as soon as positions for all local MobileStations have been 
// received, then adds it's own position and sends a BsMsPositions message 
// to all neighbouring cells. Furthermore, those position messages as 
// well as all such messages received from other cells are forwarded 
// to the BsChannel instance with index 0. 
//
// That BsChannel with index 0 is always the channel which initializes 
// the local cell's METISChannel instance, as soon as BsMsPositions messages 
// have been received from all cells in the simulation. Then, the initial 
// coefficient tables for the METIS model can be computed.
//
// <h3>\Traffic and Information Flow</h3>
// After the <tt>initOffset</tt> phase has run, the TrafficGen modules start 
// generating traffic according to their stream definitions. With the period 
// defined for each stream, KoiData packets are continually generated. 
// These are then forwarded to the MsMac module, where they are put into 
// queues sorted by their unique stream id. 
//
// At the beginning of each TTI, the MobileStations as well as the 
// BaseStations take all of their stream queues and send them to the 
// Scheduler via StreamTransReq messages, one for each stream. After a 
// small delay of length <tt>epsilon</tt>, the StreamScheduler sends lists 
// with StreamTransReq messages, in the form of TransReqList messages, to each 
// RBScheduler. Each of those lists contain all StreamTransReq messages 
// for all streams handled by each RBScheduler. 
// The small <tt>epsilon</tt> delay is necessary because the StreamScheduler 
// module has no way of knowing when it has received all requests for the 
// current TTI. But, as all requests are generated and send without 
// delay right at the start of each TTI, the StreamScheduler only needs to 
// wait a small time to be certain that all requests have been received.
//
// After all TransReqList messages have been received, the RBSchedulers 
// choose a single packet from all the the packets for all the streams they 
// are responsible for and note it in the StreamTransSched message they then 
// send back to the StreamScheduler.
//
// The StreamScheduler then sends all StreamTransSched messages to the 
// MobileStation or BaseStation where the corresponding StreamTransReq message 
// originated.
//
// When the MsMac or BsMac receive a StreamTransSched message, they search for 
// the indicated packet in their stream queues and remove it. Then, they 
// wrap in in a DataPacketBundle and forward it to the MsPhy/BsPhy respectively.
//
// In addition, both produce TransInfo messages, which are forwarded to 
// neighbouring cells. Via the BaseStation, in the case of the MsMac. These 
// messages are received by all BsMacs and then forwarded to their own BsChannels
// as well as their local MobileStations, where they are received by the MsMac 
// and forwarded to the MsChannel. 
//
// The TransInfo messages are directed to the MsChannel/BsChannel modules because 
// they are the modules where SINR values for packet transmissions are calculated.
//
// Once a packet reaches the MsPhy/BsPhy, they are directly forwarded to their 
// destination. The Phy modules to not implement any behaviour other than 
// rudimentary routing. 
//
// In the case of packets send by a MobileStation, the MsPhy forwards the 
// packet to one of two modules. If the packet belongs to a D2D stream, it is 
// forwarded to the MsChannel of the destination MobileStation. Otherwise, 
// the packet is send to the BsChannel of the local BaseStation.
//
// In the case of packets forwarded by the BaseStation, the BsPhy sends 
// the packet to the MsChannel of the intended MobileStation destination.
//
// Once a packet reaches the MsChannel/BsChannel, those modules compute 
// SINR values for the transmission, using the METISChannel class and 
// taking into account the previously received TransInfo messages for interference.
//
// If the packet is received successfully (currently, all packets are), it 
// is then forwarded to the MsPhy or BsPhy and from their to the respective 
// Mac module. 
//
// When a packet arrives at the MsMac, it is directly forwarded to the TrafficGen,
// which is the <em>sink</em> for packets in this model.
//
// When a DataPacketBundle arrives at the BsMac, all KoiData packets within 
// are sorted into queues according to their streams, later to be forwarded 
// to their destination once the schedule allows it.
// @page streams.html, Stream Definition Syntax
// Streams are defined for all mobile stations in a XML file. That file 
// can be provided to each MobileStation in it's <tt>commTable</tt> parameter.
// The syntax is best shown via the following \example:
//
// <pre>
//<root>
//	<cell id="0">
//		<ms id="0">
//			<stream cell="0" destid="1" period="0.5" d2d="0" deadline="1.0" />
//		</ms>
//		<ms id="1">
//			<stream cell="0" destid="0" period="0.5" d2d="0" deadline="1.0" />
//		</ms>
//	</cell>
//	<cell id="1">
//		<ms id="0">
//			<stream cell="1" destid="1" period="0.5" d2d="0" deadline="1.0" />
//		</ms>
//		<ms id="1">
//			<stream cell="1" destid="0" period="0.5" d2d="0" deadline="1.0" />
//		</ms>
//	</cell>
//</root>
// </pre>
// This file defines streams for two cells with two MobileStations each. For 
// each local pair, two streams are set up, one for each direction. 
// All streams have similar periods of 0.5 seconds and similar deadlines of 1 second. 
// Both parameters, <tt>deadline</tt> and <tt>period</tt>, are expected in seconds.
// There can be an arbitrary 
// number of <tt><stream></tt> elements in each <tt><ms></tt> element. 
// Care needs to be taken to make sure that all ids are correct and MobileStation 
// modules corresponding to them exist. As inter cell communication is 
// currently not implemented, each <tt>cell</tt> attribute must correspond 
// to the cell of the sending MobileStation. 
// If the <tt>d2d</tt> attribute is set to "1", packets of that stream will 
// be transmitted directly to the recipient. Otherwise, the packet is first 
// send to the local BaseStation and then transmitted by it to the recipient.
