//
// Generated file, do not edit! Created by opp_msgc 4.3 from SINR.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "SINR_m.h"

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




Register_Class(SINR);

SINR::SINR(const char *name, int kind) : cMessage(name,kind)
{
    this->msId_var = 0;
    SINR_arraysize = 0;
    this->SINR_var = 0;
}

SINR::SINR(const SINR& other) : cMessage(other)
{
    SINR_arraysize = 0;
    this->SINR_var = 0;
    copy(other);
}

SINR::~SINR()
{
    delete [] SINR_var;
}

SINR& SINR::operator=(const SINR& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void SINR::copy(const SINR& other)
{
    this->msId_var = other.msId_var;
    delete [] this->SINR_var;
    this->SINR_var = (other.SINR_arraysize==0) ? NULL : new double[other.SINR_arraysize];
    SINR_arraysize = other.SINR_arraysize;
    for (unsigned int i=0; i<SINR_arraysize; i++)
        this->SINR_var[i] = other.SINR_var[i];
}

void SINR::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->msId_var);
    b->pack(SINR_arraysize);
    doPacking(b,this->SINR_var,SINR_arraysize);
}

void SINR::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->msId_var);
    delete [] this->SINR_var;
    b->unpack(SINR_arraysize);
    if (SINR_arraysize==0) {
        this->SINR_var = 0;
    } else {
        this->SINR_var = new double[SINR_arraysize];
        doUnpacking(b,this->SINR_var,SINR_arraysize);
    }
}

int SINR::getMsId() const
{
    return msId_var;
}

void SINR::setMsId(int msId)
{
    this->msId_var = msId;
}

void SINR::setSINRArraySize(unsigned int size)
{
    double *SINR_var2 = (size==0) ? NULL : new double[size];
    unsigned int sz = SINR_arraysize < size ? SINR_arraysize : size;
    for (unsigned int i=0; i<sz; i++)
        SINR_var2[i] = this->SINR_var[i];
    for (unsigned int i=sz; i<size; i++)
        SINR_var2[i] = 0;
    SINR_arraysize = size;
    delete [] this->SINR_var;
    this->SINR_var = SINR_var2;
}

unsigned int SINR::getSINRArraySize() const
{
    return SINR_arraysize;
}

double SINR::getSINR(unsigned int k) const
{
    if (k>=SINR_arraysize) throw cRuntimeError("Array of size %d indexed by %d", SINR_arraysize, k);
    return SINR_var[k];
}

void SINR::setSINR(unsigned int k, double SINR)
{
    if (k>=SINR_arraysize) throw cRuntimeError("Array of size %d indexed by %d", SINR_arraysize, k);
    this->SINR_var[k] = SINR;
}

class SINRDescriptor : public cClassDescriptor
{
  public:
    SINRDescriptor();
    virtual ~SINRDescriptor();

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

Register_ClassDescriptor(SINRDescriptor);

SINRDescriptor::SINRDescriptor() : cClassDescriptor("SINR", "cMessage")
{
}

SINRDescriptor::~SINRDescriptor()
{
}

bool SINRDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<SINR *>(obj)!=NULL;
}

const char *SINRDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int SINRDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 2+basedesc->getFieldCount(object) : 2;
}

unsigned int SINRDescriptor::getFieldTypeFlags(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeFlags(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,
        FD_ISARRAY | FD_ISEDITABLE,
    };
    return (field>=0 && field<2) ? fieldTypeFlags[field] : 0;
}

const char *SINRDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "msId",
        "SINR",
    };
    return (field>=0 && field<2) ? fieldNames[field] : NULL;
}

int SINRDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='m' && strcmp(fieldName, "msId")==0) return base+0;
    if (fieldName[0]=='S' && strcmp(fieldName, "SINR")==0) return base+1;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *SINRDescriptor::getFieldTypeString(void *object, int field) const
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
    };
    return (field>=0 && field<2) ? fieldTypeStrings[field] : NULL;
}

const char *SINRDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
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

int SINRDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    SINR *pp = (SINR *)object; (void)pp;
    switch (field) {
        case 1: return pp->getSINRArraySize();
        default: return 0;
    }
}

std::string SINRDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    SINR *pp = (SINR *)object; (void)pp;
    switch (field) {
        case 0: return long2string(pp->getMsId());
        case 1: return double2string(pp->getSINR(i));
        default: return "";
    }
}

bool SINRDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    SINR *pp = (SINR *)object; (void)pp;
    switch (field) {
        case 0: pp->setMsId(string2long(value)); return true;
        case 1: pp->setSINR(i,string2double(value)); return true;
        default: return false;
    }
}

const char *SINRDescriptor::getFieldStructName(void *object, int field) const
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
    };
    return (field>=0 && field<2) ? fieldStructNames[field] : NULL;
}

void *SINRDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    SINR *pp = (SINR *)object; (void)pp;
    switch (field) {
        default: return NULL;
    }
}


