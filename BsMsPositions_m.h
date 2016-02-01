//
// Generated file, do not edit! Created by opp_msgc 4.3 from BsMsPositions.msg.
//

#ifndef _BSMSPOSITIONS_M_H_
#define _BSMSPOSITIONS_M_H_

#include <omnetpp.h>

// opp_msgc version check
#define MSGC_VERSION 0x0403
#if (MSGC_VERSION!=OMNETPP_VERSION)
#    error Version mismatch! Probably this file was generated by an earlier version of opp_msgc: 'make clean' should help.
#endif

// cplusplus {{
#include "Position.h"
// }}



/**
 * Class generated from <tt>BsMsPositions.msg</tt> by opp_msgc.
 * <pre>
 * message BsMsPositions  {
 * 	int bsId;
 * 	Position positions[];
 * }
 * </pre>
 */
class BsMsPositions : public ::cMessage
{
  protected:
    int bsId_var;
    Position *positions_var; // array ptr
    unsigned int positions_arraysize;

  private:
    void copy(const BsMsPositions& other);

  protected:
    // protected and unimplemented operator==(), to prevent accidental usage
    bool operator==(const BsMsPositions&);

  public:
    BsMsPositions(const char *name=NULL, int kind=0);
    BsMsPositions(const BsMsPositions& other);
    virtual ~BsMsPositions();
    BsMsPositions& operator=(const BsMsPositions& other);
    virtual BsMsPositions *dup() const {return new BsMsPositions(*this);}
    virtual void parsimPack(cCommBuffer *b);
    virtual void parsimUnpack(cCommBuffer *b);

    // field getter/setter methods
    virtual int getBsId() const;
    virtual void setBsId(int bsId);
    virtual void setPositionsArraySize(unsigned int size);
    virtual unsigned int getPositionsArraySize() const;
    virtual Position& getPositions(unsigned int k);
    virtual const Position& getPositions(unsigned int k) const {return const_cast<BsMsPositions*>(this)->getPositions(k);}
    virtual void setPositions(unsigned int k, const Position& positions);
};

inline void doPacking(cCommBuffer *b, BsMsPositions& obj) {obj.parsimPack(b);}
inline void doUnpacking(cCommBuffer *b, BsMsPositions& obj) {obj.parsimUnpack(b);}


#endif // _BSMSPOSITIONS_M_H_
