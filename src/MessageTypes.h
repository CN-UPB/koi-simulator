/**
 * \file MessageTypes.h
 * Contains enums used to distinguish message types
 */

#pragma once

enum TrafficType {
	periodic, /**< This is a message from a process with a fixed period*/
	gaussian /**< This is a message from a process with normal distributed intervals*/
};

enum MessageType: short {
	traffic,
	koidata,
	streamInfo,
	streamSched,
	streamTransReq,
	transInfo,
	transReqList,
	scheduleStreams,
	scheduleRBs,
	sendSchedules,
	sinrEst,
	longTermSinrEst,
	cleanupTTI,
	scheduleInfo,
	staticSchedule,
	genStaticSchedule
};

enum MessageDirection: int {
	up,
	down,
	d2d,
	d2dDown,
	d2dUp
};
