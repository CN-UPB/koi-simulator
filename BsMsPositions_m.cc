//
// Generated file, do not edit! Created by opp_msgc 4.3 from BsMsPositions.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "BsMsPositions_m.h"

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




Register_Class(BsMsPositions);

BsMsPositions::BsMsPositions(const char *name, int kind) : cMessage(name,kind)
{
    this->bsId_var = 0;
    positions_arraysize = 0;
    this->positions_var = 0;
}

BsMsPositions::BsMsPositions(const BsMsPositions& other) : cMessage(other)
{
    positions_arraysize = 0;
    this->positions_var = 0;
    copy(other);
}

BsMsPositions::~BsMsPositions()
{
    delete [] positions_var;
}

BsMsPositions& BsMsPositions::operator=(const BsMsPositions& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void BsMsPositions::copy(const BsMsPositions& other)
{
    this->bsId_var = other.bsId_var;
    delete [] this->positions_var;
    this->positions_var = (other.positions_arraysize==0) ? NULL : new Position[other.positions_arraysize];
    positions_arraysize = other.positions_arraysize;
    for (unsigned int i=0; i<positions_arraysize; i++)
        this->positions_var[i] = other.positions_var[i];
}

void BsMsPositions::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->bsId_var);
    b->pack(positions_arraysize);
    doPacking(b,this->positions_var,positions_arraysize);
}

void BsMsPositions::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->bsId_var);
    delete [] this->positions_var;
    b->unpack(positions_arraysize);
    if (positions_arraysize==0) {
        this->positions_var = 0;
    } else {
        this->positions_var = new Position[positions_arraysize];
        doUnpacking(b,this->positions_var,positions_arraysize);
    }
}

int BsMsPositions::getBsId() const
{
    return bsId_var;
}

void BsMsPositions::setBsId(int bsId)
{
    this->bsId_var = bsId;
}

void BsMsPositions::setPositionsArraySize(unsigned int size)
{
    Position *positions_var2 = (size==0) ? NULL : new Position[size];
    unsigned int sz = positions_arraysize < size ? positions_arraysize : size;
    for (unsigned int i=0; i<sz; i++)
        positions_var2[i] = this->positions_var[i];
    positions_arraysize = size;
    delete [] this->positions_var;
    this->positions_var = positions_var2;
}

unsigned int BsMsPositions::getPositionsArraySize() const
{
    return positions_arraysize;
}

Position& BsMsPositions::getPositions(unsigned int k)
{
    if (k>=positions_arraysize) throw cRuntimeError("Array of size %d indexed by %d", positions_arraysize, k);
    return positions_var[k];
}

void BsMsPositions::setPositions(unsigned int k, const Position& positions)
{
    if (k>=positions_arraysize) throw cRuntimeError("Array of size %d indexed by %d", positions_arraysize, k);
    this->positions_var[k] = positions;
}

class BsMsPositionsDescriptor : public cClassDescriptor
{
  public:
    BsMsPositionsDescriptor();
    virtual ~BsMsPositionsDescriptor();

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

Register_ClassDescriptor(BsMsPositionsDescriptor);

BsMsPositionsDescriptor::BsMsPositionsDescriptor() : cClassDescriptor("BsMsPositions", "cMessage")
{
}

BsMsPositionsDescriptor::~BsMsPositionsDescriptor()
{
}

bool BsMsPositionsDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<BsMsPositions *>(obj)!=NULL;
}

const char *BsMsPositionsDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int BsMsPositionsDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 2+basedesc->getFieldCount(object) : 2;
}

unsigned int BsMsPositionsDescriptor::getFieldTypeFlags(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeFlags(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,
        FD_ISARRAY | FD_ISCOMPOUND,
    };
    return (field>=0 && field<2) ? fieldTypeFlags[field] : 0;
}

const char *BsMsPositionsDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "bsId",
        "positions",
    };
    return (field>=0 && field<2) ? fieldNames[field] : NULL;
}

int BsMsPositionsDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='b' && strcmp(fieldName, "bsId")==0) return base+0;
    if (fieldName[0]=='p' && strcmp(fieldName, "positions")==0) return base+1;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *BsMsPositionsDescriptor::getFieldTypeString(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeString(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldTypeStrings[] = {
        "int",
        "Position",
    };
    return (field>=0 && field<2) ? fieldTypeStrings[field] : NULL;
}

const char *BsMsPositionsDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
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

int BsMsPositionsDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    BsMsPositions *pp = (BsMsPositions *)object; (void)pp;
    switch (field) {
        case 1: return pp->getPositionsArraySize();
        default: return 0;
    }
}

std::string BsMsPositionsDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    BsMsPositions *pp = (BsMsPositions *)object; (void)pp;
    switch (field) {
        case 0: return long2string(pp->getBsId());
        case 1: {std::stringstream out; out << pp->getPositions(i); return out.str();}
        default: return "";
    }
}

bool BsMsPositionsDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    BsMsPositions *pp = (BsMsPositions *)object; (void)pp;
    switch (field) {
        case 0: pp->setBsId(string2long(value)); return true;
        default: return false;
    }
}

const char *BsMsPositionsDescriptor::getFieldStructName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldStructNames[] = {
        NULL,
        "Position",
    };
    return (field>=0 && field<2) ? fieldStructNames[field] : NULL;
}

void *BsMsPositionsDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    BsMsPositions *pp = (BsMsPositions *)object; (void)pp;
    switch (field) {
        case 1: return (void *)(&pp->getPositions(i)); break;
        default: return NULL;
    }
}


