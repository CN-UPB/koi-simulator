//
// Generated file, do not edit! Created by opp_msgc 4.3 from DataPacketBundle.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "DataPacketBundle_m.h"

// Template rule which fires if a struct or class doesn't have operator<<
template<typename T>
std::ostream& operator<<(std::ostream& out,const T&) {return out;}

// Another default rule (prevents compiler from choosing base class' doPacking())
template<typename T>
void doPacking(cCommBuffer *, T& t) {
    throw cRuntimeError("Parsim error: no doPacking() function for type %s or its base class (check .msg and _m.cc/h files!)",opp_typename(typeid(t)));
}

template<typename T>
void doUnpacking(cCommBuffer *, T& t) {
    throw cRuntimeError("Parsim error: no doUnpacking() function for type %s or its base class (check .msg and _m.cc/h files!)",opp_typename(typeid(t)));
}




Register_Class(DataPacketBundle);

DataPacketBundle::DataPacketBundle(const char *name, int kind) : cMessage(name,kind)
{
    this->msId_var = 0;
    this->bsId_var = 0;
    packets_arraysize = 0;
    this->packets_var = 0;
    RBs_arraysize = 0;
    this->RBs_var = 0;
    this->cqi_var = 0;
}

DataPacketBundle::DataPacketBundle(const DataPacketBundle& other) : cMessage(other)
{
    packets_arraysize = 0;
    this->packets_var = 0;
    RBs_arraysize = 0;
    this->RBs_var = 0;
    copy(other);
}

DataPacketBundle::~DataPacketBundle()
{
    for (unsigned int i=0; i<packets_arraysize; i++)
        drop(&(this->packets_var[i]));
    delete [] packets_var;
    delete [] RBs_var;
}

DataPacketBundle& DataPacketBundle::operator=(const DataPacketBundle& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void DataPacketBundle::copy(const DataPacketBundle& other)
{
    this->msId_var = other.msId_var;
    this->bsId_var = other.bsId_var;
    delete [] this->packets_var;
    this->packets_var = (other.packets_arraysize==0) ? NULL : new DataPacket[other.packets_arraysize];
    packets_arraysize = other.packets_arraysize;
    for (unsigned int i=0; i<packets_arraysize; i++)
    {
        take(&(this->packets_var[i]));
        this->packets_var[i] = other.packets_var[i];
        this->packets_var[i].setName(other.packets_var[i].getName());
    }
    delete [] this->RBs_var;
    this->RBs_var = (other.RBs_arraysize==0) ? NULL : new int[other.RBs_arraysize];
    RBs_arraysize = other.RBs_arraysize;
    for (unsigned int i=0; i<RBs_arraysize; i++)
        this->RBs_var[i] = other.RBs_var[i];
    this->cqi_var = other.cqi_var;
}

void DataPacketBundle::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->msId_var);
    doPacking(b,this->bsId_var);
    b->pack(packets_arraysize);
    doPacking(b,this->packets_var,packets_arraysize);
    b->pack(RBs_arraysize);
    doPacking(b,this->RBs_var,RBs_arraysize);
    doPacking(b,this->cqi_var);
}

void DataPacketBundle::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->msId_var);
    doUnpacking(b,this->bsId_var);
    delete [] this->packets_var;
    b->unpack(packets_arraysize);
    if (packets_arraysize==0) {
        this->packets_var = 0;
    } else {
        this->packets_var = new DataPacket[packets_arraysize];
        doUnpacking(b,this->packets_var,packets_arraysize);
    }
    delete [] this->RBs_var;
    b->unpack(RBs_arraysize);
    if (RBs_arraysize==0) {
        this->RBs_var = 0;
    } else {
        this->RBs_var = new int[RBs_arraysize];
        doUnpacking(b,this->RBs_var,RBs_arraysize);
    }
    doUnpacking(b,this->cqi_var);
}

int DataPacketBundle::getMsId() const
{
    return msId_var;
}

void DataPacketBundle::setMsId(int msId)
{
    this->msId_var = msId;
}

int DataPacketBundle::getBsId() const
{
    return bsId_var;
}

void DataPacketBundle::setBsId(int bsId)
{
    this->bsId_var = bsId;
}

void DataPacketBundle::setPacketsArraySize(unsigned int size)
{
    DataPacket *packets_var2 = (size==0) ? NULL : new DataPacket[size];
    unsigned int sz = packets_arraysize < size ? packets_arraysize : size;
    for (unsigned int i=0; i<sz; i++)
        packets_var2[i] = this->packets_var[i];
    for (unsigned int i=sz; i<size; i++)
        take(&(packets_var2[i]));
    packets_arraysize = size;
    delete [] this->packets_var;
    this->packets_var = packets_var2;
}

unsigned int DataPacketBundle::getPacketsArraySize() const
{
    return packets_arraysize;
}

DataPacket& DataPacketBundle::getPackets(unsigned int k)
{
    if (k>=packets_arraysize) throw cRuntimeError("Array of size %d indexed by %d", packets_arraysize, k);
    return packets_var[k];
}

void DataPacketBundle::setPackets(unsigned int k, const DataPacket& packets)
{
    if (k>=packets_arraysize) throw cRuntimeError("Array of size %d indexed by %d", packets_arraysize, k);
    this->packets_var[k] = packets;
}

void DataPacketBundle::setRBsArraySize(unsigned int size)
{
    int *RBs_var2 = (size==0) ? NULL : new int[size];
    unsigned int sz = RBs_arraysize < size ? RBs_arraysize : size;
    for (unsigned int i=0; i<sz; i++)
        RBs_var2[i] = this->RBs_var[i];
    for (unsigned int i=sz; i<size; i++)
        RBs_var2[i] = 0;
    RBs_arraysize = size;
    delete [] this->RBs_var;
    this->RBs_var = RBs_var2;
}

unsigned int DataPacketBundle::getRBsArraySize() const
{
    return RBs_arraysize;
}

int DataPacketBundle::getRBs(unsigned int k) const
{
    if (k>=RBs_arraysize) throw cRuntimeError("Array of size %d indexed by %d", RBs_arraysize, k);
    return RBs_var[k];
}

void DataPacketBundle::setRBs(unsigned int k, int RBs)
{
    if (k>=RBs_arraysize) throw cRuntimeError("Array of size %d indexed by %d", RBs_arraysize, k);
    this->RBs_var[k] = RBs;
}

int DataPacketBundle::getCqi() const
{
    return cqi_var;
}

void DataPacketBundle::setCqi(int cqi)
{
    this->cqi_var = cqi;
}

class DataPacketBundleDescriptor : public cClassDescriptor
{
  public:
    DataPacketBundleDescriptor();
    virtual ~DataPacketBundleDescriptor();

    virtual bool doesSupport(cObject *obj) const;
    virtual const char *getProperty(const char *propertyname) const;
    virtual int getFieldCount(void *object) const;
    virtual const char *getFieldName(void *object, int field) const;
    virtual int findField(void *object, const char *fieldName) const;
    virtual unsigned int getFieldTypeFlags(void *object, int field) const;
    virtual const char *getFieldTypeString(void *object, int field) const;
    virtual const char *getFieldProperty(void *object, int field, const char *propertyname) const;
    virtual int getArraySize(void *object, int field) const;

    virtual std::string getFieldAsString(void *object, int field, int i) const;
    virtual bool setFieldAsString(void *object, int field, int i, const char *value) const;

    virtual const char *getFieldStructName(void *object, int field) const;
    virtual void *getFieldStructPointer(void *object, int field, int i) const;
};

Register_ClassDescriptor(DataPacketBundleDescriptor);

DataPacketBundleDescriptor::DataPacketBundleDescriptor() : cClassDescriptor("DataPacketBundle", "cMessage")
{
}

DataPacketBundleDescriptor::~DataPacketBundleDescriptor()
{
}

bool DataPacketBundleDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<DataPacketBundle *>(obj)!=NULL;
}

const char *DataPacketBundleDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int DataPacketBundleDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 5+basedesc->getFieldCount(object) : 5;
}

unsigned int DataPacketBundleDescriptor::getFieldTypeFlags(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeFlags(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISARRAY | FD_ISCOMPOUND | FD_ISCOBJECT | FD_ISCOWNEDOBJECT,
        FD_ISARRAY | FD_ISEDITABLE,
        FD_ISEDITABLE,
    };
    return (field>=0 && field<5) ? fieldTypeFlags[field] : 0;
}

const char *DataPacketBundleDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "msId",
        "bsId",
        "packets",
        "RBs",
        "cqi",
    };
    return (field>=0 && field<5) ? fieldNames[field] : NULL;
}

int DataPacketBundleDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='m' && strcmp(fieldName, "msId")==0) return base+0;
    if (fieldName[0]=='b' && strcmp(fieldName, "bsId")==0) return base+1;
    if (fieldName[0]=='p' && strcmp(fieldName, "packets")==0) return base+2;
    if (fieldName[0]=='R' && strcmp(fieldName, "RBs")==0) return base+3;
    if (fieldName[0]=='c' && strcmp(fieldName, "cqi")==0) return base+4;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *DataPacketBundleDescriptor::getFieldTypeString(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeString(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldTypeStrings[] = {
        "int",
        "int",
        "DataPacket",
        "int",
        "int",
    };
    return (field>=0 && field<5) ? fieldTypeStrings[field] : NULL;
}

const char *DataPacketBundleDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldProperty(object, field, propertyname);
        field -= basedesc->getFieldCount(object);
    }
    switch (field) {
        default: return NULL;
    }
}

int DataPacketBundleDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    DataPacketBundle *pp = (DataPacketBundle *)object; (void)pp;
    switch (field) {
        case 2: return pp->getPacketsArraySize();
        case 3: return pp->getRBsArraySize();
        default: return 0;
    }
}

std::string DataPacketBundleDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    DataPacketBundle *pp = (DataPacketBundle *)object; (void)pp;
    switch (field) {
        case 0: return long2string(pp->getMsId());
        case 1: return long2string(pp->getBsId());
        case 2: {std::stringstream out; out << pp->getPackets(i); return out.str();}
        case 3: return long2string(pp->getRBs(i));
        case 4: return long2string(pp->getCqi());
        default: return "";
    }
}

bool DataPacketBundleDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    DataPacketBundle *pp = (DataPacketBundle *)object; (void)pp;
    switch (field) {
        case 0: pp->setMsId(string2long(value)); return true;
        case 1: pp->setBsId(string2long(value)); return true;
        case 3: pp->setRBs(i,string2long(value)); return true;
        case 4: pp->setCqi(string2long(value)); return true;
        default: return false;
    }
}

const char *DataPacketBundleDescriptor::getFieldStructName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldStructNames[] = {
        NULL,
        NULL,
        "DataPacket",
        NULL,
        NULL,
    };
    return (field>=0 && field<5) ? fieldStructNames[field] : NULL;
}

void *DataPacketBundleDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    DataPacketBundle *pp = (DataPacketBundle *)object; (void)pp;
    switch (field) {
        case 2: return (void *)static_cast<cObject *>(&pp->getPackets(i)); break;
        default: return NULL;
    }
}


