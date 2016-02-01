#
# OMNeT++/OMNEST Makefile for abstractLTEChannelModel
#
# This file was generated with the command:
#  opp_makemake -f
#

# Name of target to be created (-o option)
TARGET = abstractLTEChannelModel$(EXE_SUFFIX)

# User interface (uncomment one) (-u option)
USERIF_LIBS = $(ALL_ENV_LIBS) # that is, $(TKENV_LIBS) $(CMDENV_LIBS)
#USERIF_LIBS = $(CMDENV_LIBS)
#USERIF_LIBS = $(TKENV_LIBS)

# C++ include paths (with -I)
INCLUDE_PATH = -I.

# Additional object and library files to link with
EXTRA_OBJS =

# Additional libraries (-L, -l options)
LIBS =

# Output directory
PROJECT_OUTPUT_DIR = out
PROJECTRELATIVE_PATH =
O = $(PROJECT_OUTPUT_DIR)/$(CONFIGNAME)/$(PROJECTRELATIVE_PATH)

# Object files for local .cc and .msg files
OBJS = \
    $O/App.o \
    $O/bestCQIScheduler.o \
    $O/BsChannel.o \
    $O/BsMac.o \
    $O/BsPhy.o \
    $O/Channel.o \
    $O/ChannelAlternative.o \
    $O/cluster.o \
    $O/Cost2100Channel.o \
    $O/Decider.o \
    $O/METISChannel.o \
    $O/mpc.o \
    $O/MsChannel.o \
    $O/MsMac.o \
    $O/MsPhy.o \
    $O/MySoSFading.o \
    $O/NeighbourIdMatching.o \
    $O/Partition.o \
    $O/Position.o \
    $O/proportionalFairScheduler.o \
    $O/PtrExchange.o \
    $O/randomScheduler.o \
    $O/roundRobinScheduler.o \
    $O/Scheduler.o \
    $O/SimpleChannelCalc.o \
    $O/StopWatch.o \
    $O/util.o \
    $O/VisibilityRegion.o \
    $O/BsMsPositions_m.o \
    $O/ChannelExchange_m.o \
    $O/ClusterMessage_m.o \
    $O/DataPacket_m.o \
    $O/DataPacketBundle_m.o \
    $O/PointerExchange_m.o \
    $O/PositionExchange_m.o \
    $O/Schedule_m.o \
    $O/SINR_m.o \
    $O/TransmitRequest_m.o \
    $O/VisibilityRegionMessage_m.o

# Message files
MSGFILES = \
    BsMsPositions.msg \
    ChannelExchange.msg \
    ClusterMessage.msg \
    DataPacket.msg \
    DataPacketBundle.msg \
    PointerExchange.msg \
    PositionExchange.msg \
    Schedule.msg \
    SINR.msg \
    TransmitRequest.msg \
    VisibilityRegionMessage.msg

#------------------------------------------------------------------------------

# Pull in OMNeT++ configuration (Makefile.inc or configuser.vc)

ifneq ("$(OMNETPP_CONFIGFILE)","")
CONFIGFILE = $(OMNETPP_CONFIGFILE)
else
ifneq ("$(OMNETPP_ROOT)","")
CONFIGFILE = $(OMNETPP_ROOT)/Makefile.inc
else
CONFIGFILE = $(shell opp_configfilepath)
endif
endif

ifeq ("$(wildcard $(CONFIGFILE))","")
$(error Config file '$(CONFIGFILE)' does not exist -- add the OMNeT++ bin directory to the path so that opp_configfilepath can be found, or set the OMNETPP_CONFIGFILE variable to point to Makefile.inc)
endif

include $(CONFIGFILE)

# Simulation kernel and user interface libraries
OMNETPP_LIB_SUBDIR = $(OMNETPP_LIB_DIR)/$(TOOLCHAIN_NAME)
OMNETPP_LIBS = -L"$(OMNETPP_LIB_SUBDIR)" -L"$(OMNETPP_LIB_DIR)" $(USERIF_LIBS) $(KERNEL_LIBS) $(SYS_LIBS)

COPTS = $(CFLAGS)  $(INCLUDE_PATH) -I$(OMNETPP_INCL_DIR) -std=c++11
MSGCOPTS = $(INCLUDE_PATH)

#------------------------------------------------------------------------------
# User-supplied makefile fragment(s)
# >>>
# <<<
#------------------------------------------------------------------------------

# Main target
all: $(TARGET)

$(TARGET) : $O/$(TARGET)
	$(LN) $O/$(TARGET) .

$O/$(TARGET): $(OBJS)  $(wildcard $(EXTRA_OBJS)) Makefile
	@$(MKPATH) $O
	$(CXX) $(LDFLAGS) -o $O/$(TARGET)  $(OBJS) $(EXTRA_OBJS) $(WHOLE_ARCHIVE_ON) $(LIBS) $(WHOLE_ARCHIVE_OFF) $(OMNETPP_LIBS)

.PHONY:

.SUFFIXES: .cc

$O/%.o: %.cc
	@$(MKPATH) $(dir $@)
	$(CXX) -c $(COPTS) -o $@ $<

%_m.cc %_m.h: %.msg
	$(MSGC) -s _m.cc $(MSGCOPTS) $?

msgheaders: $(MSGFILES:.msg=_m.h)

clean:
	-rm -rf $O
	-rm -f abstractLTEChannelModel abstractLTEChannelModel.exe libabstractLTEChannelModel.so libabstractLTEChannelModel.a libabstractLTEChannelModel.dll libabstractLTEChannelModel.dylib
	-rm -f ./*_m.cc ./*_m.h

cleanall: clean
	-rm -rf $(PROJECT_OUTPUT_DIR)

depend:
	$(MAKEDEPEND) $(INCLUDE_PATH) -f Makefile -P\$$O/ -- $(MSG_CC_FILES)  ./*.cc

# DO NOT DELETE THIS LINE -- make depend depends on it.
$O/App.o: App.cc \
  ./DataPacket_m.h \
  ./DataPacketBundle_m.h \
  ./App.h
$O/bestCQIScheduler.o: bestCQIScheduler.cc \
  ./bestCQIScheduler.h \
  ./Scheduler.h \
  ./Schedule_m.h
$O/BsChannel.o: BsChannel.cc \
  ./ChannelAlternative.h \
  ./NeighbourIdMatching.h \
  ./PositionExchange_m.h \
  ./DataPacketBundle_m.h \
  ./ClusterMessage_m.h \
  ./Schedule_m.h \
  ./Cost2100Channel.h \
  ./DataPacket_m.h \
  ./Position.h \
  ./mpc.h \
  ./PtrExchange.h \
  ./Channel.h \
  ./SimpleChannelCalc.h \
  ./PointerExchange_m.h \
  ./VisibilityRegionMessage_m.h \
  ./SINR_m.h \
  ./cluster.h \
  ./BsMsPositions_m.h \
  ./util.h \
  ./BsChannel.h \
  ./VisibilityRegion.h \
  ./METISChannel.h
$O/BsMac.o: BsMac.cc \
  ./VisibilityRegion.h \
  ./TransmitRequest_m.h \
  ./ChannelExchange_m.h \
  ./proportionalFairScheduler.h \
  ./BsMac.h \
  ./VisibilityRegionMessage_m.h \
  ./PointerExchange_m.h \
  ./BsMsPositions_m.h \
  ./cluster.h \
  ./SINR_m.h \
  ./bestCQIScheduler.h \
  ./util.h \
  ./randomScheduler.h \
  ./roundRobinScheduler.h \
  ./mpc.h \
  ./Position.h \
  ./PtrExchange.h \
  ./NeighbourIdMatching.h \
  ./PositionExchange_m.h \
  ./Schedule_m.h \
  ./ClusterMessage_m.h \
  ./DataPacketBundle_m.h \
  ./DataPacket_m.h \
  ./Scheduler.h
$O/BsPhy.o: BsPhy.cc \
  ./PtrExchange.h \
  ./BsPhy.h \
  ./DataPacket_m.h \
  ./DataPacketBundle_m.h \
  ./PointerExchange_m.h
$O/ChannelAlternative.o: ChannelAlternative.cc \
  ./VisibilityRegion.h \
  ./VisibilityRegionMessage_m.h \
  ./PointerExchange_m.h \
  ./cluster.h \
  ./mpc.h \
  ./Position.h \
  ./PtrExchange.h \
  ./Channel.h \
  ./NeighbourIdMatching.h \
  ./ChannelAlternative.h \
  ./ClusterMessage_m.h
$O/Channel.o: Channel.cc \
  ./Position.h \
  ./Channel.h
$O/cluster.o: cluster.cc \
  ./mpc.h \
  ./cluster.h
$O/Cost2100Channel.o: Cost2100Channel.cc \
  ./cluster.h \
  ./VisibilityRegionMessage_m.h \
  ./PointerExchange_m.h \
  ./VisibilityRegion.h \
  ./Cost2100Channel.h \
  ./ClusterMessage_m.h \
  ./NeighbourIdMatching.h \
  ./PtrExchange.h \
  ./Channel.h \
  ./mpc.h \
  ./Position.h
$O/Decider.o: Decider.cc \
  ./Decider.h
$O/METISChannel.o: METISChannel.cc \
  ./NeighbourIdMatching.h \
  ./Position.h \
  ./Channel.h \
  ./METISChannel.h
$O/mpc.o: mpc.cc \
  ./mpc.h
$O/MsChannel.o: MsChannel.cc \
  ./DataPacket_m.h \
  ./DataPacketBundle_m.h \
  ./PositionExchange_m.h \
  ./NeighbourIdMatching.h \
  ./SimpleChannelCalc.h \
  ./PtrExchange.h \
  ./Channel.h \
  ./MsChannel.h \
  ./Position.h \
  ./util.h \
  ./SINR_m.h \
  ./PointerExchange_m.h \
  ./ChannelExchange_m.h
$O/MsMac.o: MsMac.cc \
  ./util.h \
  ./MsMac.h \
  ./SINR_m.h \
  ./TransmitRequest_m.h \
  ./PositionExchange_m.h \
  ./DataPacket_m.h \
  ./DataPacketBundle_m.h \
  ./Schedule_m.h \
  ./Position.h
$O/MsPhy.o: MsPhy.cc \
  ./MsPhy.h \
  ./SINR_m.h
$O/MySoSFading.o: MySoSFading.cc \
  ./MySoSFading.h
$O/NeighbourIdMatching.o: NeighbourIdMatching.cc \
  ./NeighbourIdMatching.h
$O/Partition.o: Partition.cc \
  ./Partition.h
$O/Position.o: Position.cc \
  ./Position.h
$O/proportionalFairScheduler.o: proportionalFairScheduler.cc \
  ./proportionalFairScheduler.h \
  ./util.h \
  ./Schedule_m.h \
  ./Scheduler.h
$O/PtrExchange.o: PtrExchange.cc \
  ./PtrExchange.h
$O/randomScheduler.o: randomScheduler.cc \
  ./randomScheduler.h \
  ./Scheduler.h \
  ./Schedule_m.h
$O/roundRobinScheduler.o: roundRobinScheduler.cc \
  ./Schedule_m.h \
  ./Scheduler.h \
  ./roundRobinScheduler.h
$O/Scheduler.o: Scheduler.cc \
  ./Schedule_m.h \
  ./Scheduler.h
$O/SimpleChannelCalc.o: SimpleChannelCalc.cc \
  ./SimpleChannelCalc.h
$O/StopWatch.o: StopWatch.cc \
  ./StopWatch.h
$O/util.o: util.cc \
  ./util.h
$O/VisibilityRegion.o: VisibilityRegion.cc \
  ./VisibilityRegion.h
