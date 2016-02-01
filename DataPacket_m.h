//
// Generated file, do not edit! Created by opp_msgc 4.3 from DataPacket.msg.
//

#ifndef _DATAPACKET_M_H_
#define _DATAPACKET_M_H_

#include <omnetpp.h>

// opp_msgc version check
#define MSGC_VERSION 0x0403
#if (MSGC_VERSION!=OMNETPP_VERSION)
#    error Version mismatch! Probably this file was generated by an earlier version of opp_msgc: 'make clean' should help.
#endif



/**
 * Class generated from <tt>DataPacket.msg</tt> by opp_msgc.
 * <pre>
 * packet DataPacket  {
 *     int msId;
 *     int bsId;
 *     int resourceBlock;
 * }
 * </pre>
 */
class DataPacket : public ::cPacket
{
  protected:
    int msId_var;
    int bsId_var;
    int resourceBlock_var;

  private:
    void copy(const DataPacket& other);

  protected:
    // protected and unimplemented operator==(), to prevent accidental usage
    bool operator==(const DataPacket&);

  public:
    DataPacket(const char *name=NULL, int kind=0);
    DataPacket(const DataPacket& other);
    virtual ~DataPacket();
    DataPacket& operator=(const DataPacket& other);
    virtual DataPacket *dup() const {return new DataPacket(*this);}
    virtual void parsimPack(cCommBuffer *b);
    virtual void parsimUnpack(cCommBuffer *b);

    // field getter/setter methods
    virtual int getMsId() const;
    virtual void setMsId(int msId);
    virtual int getBsId() const;
    virtual void setBsId(int bsId);
    virtual int getResourceBlock() const;
    virtual void setResourceBlock(int resourceBlock);
};

inline void doPacking(cCommBuffer *b, DataPacket& obj) {obj.parsimPack(b);}
inline void doUnpacking(cCommBuffer *b, DataPacket& obj) {obj.parsimUnpack(b);}


#endif // _DATAPACKET_M_H_
