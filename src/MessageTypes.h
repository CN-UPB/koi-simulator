/**
 * \file MessageTypes.h
 * Contains enums used to distinguish message types
 */

#pragma once

enum TrafficType {
	periodic, /**< This is a message from a process with a fixed period*/
	event /**< This is a message from a process with event based generation*/
};

enum MessageType: short {
	traffic,
	koidata,
	streamInfo,
	streamSched,
	streamTransReq,
	transInfoBs,
	transInfoMs,
	transReqList,
	scheduleStreams,
	scheduleRBs
};

enum MessageDirection: int {
	up,
	down,
	d2d
};
