//
// Generated file, do not edit! Created by opp_msgc 4.3 from VisibilityRegionMessage.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "VisibilityRegionMessage_m.h"

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




Register_Class(VisibilityRegionMessage);

VisibilityRegionMessage::VisibilityRegionMessage(const char *name, int kind) : cMessage(name,kind)
{
}

VisibilityRegionMessage::VisibilityRegionMessage(const VisibilityRegionMessage& other) : cMessage(other)
{
    copy(other);
}

VisibilityRegionMessage::~VisibilityRegionMessage()
{
}

VisibilityRegionMessage& VisibilityRegionMessage::operator=(const VisibilityRegionMessage& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void VisibilityRegionMessage::copy(const VisibilityRegionMessage& other)
{
    this->vr_var = other.vr_var;
}

void VisibilityRegionMessage::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->vr_var);
}

void VisibilityRegionMessage::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->vr_var);
}

VisibilityRegion& VisibilityRegionMessage::getVr()
{
    return vr_var;
}

void VisibilityRegionMessage::setVr(const VisibilityRegion& vr)
{
    this->vr_var = vr;
}

class VisibilityRegionMessageDescriptor : public cClassDescriptor
{
  public:
    VisibilityRegionMessageDescriptor();
    virtual ~VisibilityRegionMessageDescriptor();

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

Register_ClassDescriptor(VisibilityRegionMessageDescriptor);

VisibilityRegionMessageDescriptor::VisibilityRegionMessageDescriptor() : cClassDescriptor("VisibilityRegionMessage", "cMessage")
{
}

VisibilityRegionMessageDescriptor::~VisibilityRegionMessageDescriptor()
{
}

bool VisibilityRegionMessageDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<VisibilityRegionMessage *>(obj)!=NULL;
}

const char *VisibilityRegionMessageDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int VisibilityRegionMessageDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 1+basedesc->getFieldCount(object) : 1;
}

unsigned int VisibilityRegionMessageDescriptor::getFieldTypeFlags(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeFlags(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISCOMPOUND,
    };
    return (field>=0 && field<1) ? fieldTypeFlags[field] : 0;
}

const char *VisibilityRegionMessageDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "vr",
    };
    return (field>=0 && field<1) ? fieldNames[field] : NULL;
}

int VisibilityRegionMessageDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='v' && strcmp(fieldName, "vr")==0) return base+0;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *VisibilityRegionMessageDescriptor::getFieldTypeString(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeString(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldTypeStrings[] = {
        "VisibilityRegion",
    };
    return (field>=0 && field<1) ? fieldTypeStrings[field] : NULL;
}

const char *VisibilityRegionMessageDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
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

int VisibilityRegionMessageDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    VisibilityRegionMessage *pp = (VisibilityRegionMessage *)object; (void)pp;
    switch (field) {
        default: return 0;
    }
}

std::string VisibilityRegionMessageDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    VisibilityRegionMessage *pp = (VisibilityRegionMessage *)object; (void)pp;
    switch (field) {
        case 0: {std::stringstream out; out << pp->getVr(); return out.str();}
        default: return "";
    }
}

bool VisibilityRegionMessageDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    VisibilityRegionMessage *pp = (VisibilityRegionMessage *)object; (void)pp;
    switch (field) {
        default: return false;
    }
}

const char *VisibilityRegionMessageDescriptor::getFieldStructName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldStructNames[] = {
        "VisibilityRegion",
    };
    return (field>=0 && field<1) ? fieldStructNames[field] : NULL;
}

void *VisibilityRegionMessageDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    VisibilityRegionMessage *pp = (VisibilityRegionMessage *)object; (void)pp;
    switch (field) {
        case 0: return (void *)(&pp->getVr()); break;
        default: return NULL;
    }
}


