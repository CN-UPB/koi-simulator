//
// Generated file, do not edit! Created by opp_msgc 4.3 from Schedule.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "Schedule_m.h"

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




Register_Class(Schedule);

Schedule::Schedule(const char *name, int kind) : cMessage(name,kind)
{
    this->ScheduleDirection_var = 0;
    this->maxPower_var = 0;
    powerAdaptation_arraysize = 0;
    this->powerAdaptation_var = 0;
    this->bsId_var = 0;
    upSchedule_arraysize = 0;
    this->upSchedule_var = 0;
    this->channel_var = 0;
}

Schedule::Schedule(const Schedule& other) : cMessage(other)
{
    powerAdaptation_arraysize = 0;
    this->powerAdaptation_var = 0;
    upSchedule_arraysize = 0;
    this->upSchedule_var = 0;
    copy(other);
}

Schedule::~Schedule()
{
    delete [] powerAdaptation_var;
    delete [] upSchedule_var;
}

Schedule& Schedule::operator=(const Schedule& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void Schedule::copy(const Schedule& other)
{
    this->ScheduleDirection_var = other.ScheduleDirection_var;
    this->maxPower_var = other.maxPower_var;
    delete [] this->powerAdaptation_var;
    this->powerAdaptation_var = (other.powerAdaptation_arraysize==0) ? NULL : new double[other.powerAdaptation_arraysize];
    powerAdaptation_arraysize = other.powerAdaptation_arraysize;
    for (unsigned int i=0; i<powerAdaptation_arraysize; i++)
        this->powerAdaptation_var[i] = other.powerAdaptation_var[i];
    this->bsId_var = other.bsId_var;
    delete [] this->upSchedule_var;
    this->upSchedule_var = (other.upSchedule_arraysize==0) ? NULL : new int[other.upSchedule_arraysize];
    upSchedule_arraysize = other.upSchedule_arraysize;
    for (unsigned int i=0; i<upSchedule_arraysize; i++)
        this->upSchedule_var[i] = other.upSchedule_var[i];
    this->channel_var = other.channel_var;
}

void Schedule::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->ScheduleDirection_var);
    doPacking(b,this->maxPower_var);
    b->pack(powerAdaptation_arraysize);
    doPacking(b,this->powerAdaptation_var,powerAdaptation_arraysize);
    doPacking(b,this->bsId_var);
    b->pack(upSchedule_arraysize);
    doPacking(b,this->upSchedule_var,upSchedule_arraysize);
    doPacking(b,this->channel_var);
}

void Schedule::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->ScheduleDirection_var);
    doUnpacking(b,this->maxPower_var);
    delete [] this->powerAdaptation_var;
    b->unpack(powerAdaptation_arraysize);
    if (powerAdaptation_arraysize==0) {
        this->powerAdaptation_var = 0;
    } else {
        this->powerAdaptation_var = new double[powerAdaptation_arraysize];
        doUnpacking(b,this->powerAdaptation_var,powerAdaptation_arraysize);
    }
    doUnpacking(b,this->bsId_var);
    delete [] this->upSchedule_var;
    b->unpack(upSchedule_arraysize);
    if (upSchedule_arraysize==0) {
        this->upSchedule_var = 0;
    } else {
        this->upSchedule_var = new int[upSchedule_arraysize];
        doUnpacking(b,this->upSchedule_var,upSchedule_arraysize);
    }
    doUnpacking(b,this->channel_var);
}

int Schedule::getScheduleDirection() const
{
    return ScheduleDirection_var;
}

void Schedule::setScheduleDirection(int ScheduleDirection)
{
    this->ScheduleDirection_var = ScheduleDirection;
}

double Schedule::getMaxPower() const
{
    return maxPower_var;
}

void Schedule::setMaxPower(double maxPower)
{
    this->maxPower_var = maxPower;
}

void Schedule::setPowerAdaptationArraySize(unsigned int size)
{
    double *powerAdaptation_var2 = (size==0) ? NULL : new double[size];
    unsigned int sz = powerAdaptation_arraysize < size ? powerAdaptation_arraysize : size;
    for (unsigned int i=0; i<sz; i++)
        powerAdaptation_var2[i] = this->powerAdaptation_var[i];
    for (unsigned int i=sz; i<size; i++)
        powerAdaptation_var2[i] = 0;
    powerAdaptation_arraysize = size;
    delete [] this->powerAdaptation_var;
    this->powerAdaptation_var = powerAdaptation_var2;
}

unsigned int Schedule::getPowerAdaptationArraySize() const
{
    return powerAdaptation_arraysize;
}

double Schedule::getPowerAdaptation(unsigned int k) const
{
    if (k>=powerAdaptation_arraysize) throw cRuntimeError("Array of size %d indexed by %d", powerAdaptation_arraysize, k);
    return powerAdaptation_var[k];
}

void Schedule::setPowerAdaptation(unsigned int k, double powerAdaptation)
{
    if (k>=powerAdaptation_arraysize) throw cRuntimeError("Array of size %d indexed by %d", powerAdaptation_arraysize, k);
    this->powerAdaptation_var[k] = powerAdaptation;
}

int Schedule::getBsId() const
{
    return bsId_var;
}

void Schedule::setBsId(int bsId)
{
    this->bsId_var = bsId;
}

void Schedule::setUpScheduleArraySize(unsigned int size)
{
    int *upSchedule_var2 = (size==0) ? NULL : new int[size];
    unsigned int sz = upSchedule_arraysize < size ? upSchedule_arraysize : size;
    for (unsigned int i=0; i<sz; i++)
        upSchedule_var2[i] = this->upSchedule_var[i];
    for (unsigned int i=sz; i<size; i++)
        upSchedule_var2[i] = 0;
    upSchedule_arraysize = size;
    delete [] this->upSchedule_var;
    this->upSchedule_var = upSchedule_var2;
}

unsigned int Schedule::getUpScheduleArraySize() const
{
    return upSchedule_arraysize;
}

int Schedule::getUpSchedule(unsigned int k) const
{
    if (k>=upSchedule_arraysize) throw cRuntimeError("Array of size %d indexed by %d", upSchedule_arraysize, k);
    return upSchedule_var[k];
}

void Schedule::setUpSchedule(unsigned int k, int upSchedule)
{
    if (k>=upSchedule_arraysize) throw cRuntimeError("Array of size %d indexed by %d", upSchedule_arraysize, k);
    this->upSchedule_var[k] = upSchedule;
}

int Schedule::getChannel() const
{
    return channel_var;
}

void Schedule::setChannel(int channel)
{
    this->channel_var = channel;
}

class ScheduleDescriptor : public cClassDescriptor
{
  public:
    ScheduleDescriptor();
    virtual ~ScheduleDescriptor();

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

Register_ClassDescriptor(ScheduleDescriptor);

ScheduleDescriptor::ScheduleDescriptor() : cClassDescriptor("Schedule", "cMessage")
{
}

ScheduleDescriptor::~ScheduleDescriptor()
{
}

bool ScheduleDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<Schedule *>(obj)!=NULL;
}

const char *ScheduleDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int ScheduleDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 6+basedesc->getFieldCount(object) : 6;
}

unsigned int ScheduleDescriptor::getFieldTypeFlags(void *object, int field) const
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
        FD_ISARRAY | FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISARRAY | FD_ISEDITABLE,
        FD_ISEDITABLE,
    };
    return (field>=0 && field<6) ? fieldTypeFlags[field] : 0;
}

const char *ScheduleDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "ScheduleDirection",
        "maxPower",
        "powerAdaptation",
        "bsId",
        "upSchedule",
        "channel",
    };
    return (field>=0 && field<6) ? fieldNames[field] : NULL;
}

int ScheduleDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='S' && strcmp(fieldName, "ScheduleDirection")==0) return base+0;
    if (fieldName[0]=='m' && strcmp(fieldName, "maxPower")==0) return base+1;
    if (fieldName[0]=='p' && strcmp(fieldName, "powerAdaptation")==0) return base+2;
    if (fieldName[0]=='b' && strcmp(fieldName, "bsId")==0) return base+3;
    if (fieldName[0]=='u' && strcmp(fieldName, "upSchedule")==0) return base+4;
    if (fieldName[0]=='c' && strcmp(fieldName, "channel")==0) return base+5;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *ScheduleDescriptor::getFieldTypeString(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeString(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldTypeStrings[] = {
        "int",
        "double",
        "double",
        "int",
        "int",
        "int",
    };
    return (field>=0 && field<6) ? fieldTypeStrings[field] : NULL;
}

const char *ScheduleDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
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

int ScheduleDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    Schedule *pp = (Schedule *)object; (void)pp;
    switch (field) {
        case 2: return pp->getPowerAdaptationArraySize();
        case 4: return pp->getUpScheduleArraySize();
        default: return 0;
    }
}

std::string ScheduleDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    Schedule *pp = (Schedule *)object; (void)pp;
    switch (field) {
        case 0: return long2string(pp->getScheduleDirection());
        case 1: return double2string(pp->getMaxPower());
        case 2: return double2string(pp->getPowerAdaptation(i));
        case 3: return long2string(pp->getBsId());
        case 4: return long2string(pp->getUpSchedule(i));
        case 5: return long2string(pp->getChannel());
        default: return "";
    }
}

bool ScheduleDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    Schedule *pp = (Schedule *)object; (void)pp;
    switch (field) {
        case 0: pp->setScheduleDirection(string2long(value)); return true;
        case 1: pp->setMaxPower(string2double(value)); return true;
        case 2: pp->setPowerAdaptation(i,string2double(value)); return true;
        case 3: pp->setBsId(string2long(value)); return true;
        case 4: pp->setUpSchedule(i,string2long(value)); return true;
        case 5: pp->setChannel(string2long(value)); return true;
        default: return false;
    }
}

const char *ScheduleDescriptor::getFieldStructName(void *object, int field) const
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
        NULL,
        NULL,
        NULL,
        NULL,
    };
    return (field>=0 && field<6) ? fieldStructNames[field] : NULL;
}

void *ScheduleDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    Schedule *pp = (Schedule *)object; (void)pp;
    switch (field) {
        default: return NULL;
    }
}


