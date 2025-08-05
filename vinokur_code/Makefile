
SHELL=/bin/bash

RM = rm -f

BASEDIR  = $(PWD)

TARGET   = test

CXXSRCS = main.cpp stretch.cpp

SRCS  = $(CXXSRCS)

.SUFFIXES: .o .cpp
CXXOBJS = $(CXXSRCS:.cpp=.o)
OBJS  = $(CXXOBJS)

CXX = pgc++
CXXFLAGS = -O3 -D_DEBUG
CXXLIBS =

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(@) $(OBJS) $(CXXLIBS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -I. -I../include -Y *.cpp *.h

# DO NOT DELETE

main.o: stretch.h
stretch.o: stretch.h
