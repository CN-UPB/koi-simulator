#pragma once

#include "includes.h"

// uintptr_t is an unsigned int and can therefor be serialized by the pack/unpack functions
struct PtrExchange  {
    uintptr_t ptr;
};

void doPacking(cCommBuffer *buffer, PtrExchange &pointer);

void doUnpacking(cCommBuffer *buffer, PtrExchange &pointer);
