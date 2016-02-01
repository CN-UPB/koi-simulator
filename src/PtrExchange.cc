#include "PtrExchange.h"

void doPacking(cCommBuffer *buffer, PtrExchange &pointer) {
    buffer->pack(pointer.ptr);
}

void doUnpacking(cCommBuffer *buffer, PtrExchange &pointer) {
    buffer->unpack(pointer.ptr);
}
