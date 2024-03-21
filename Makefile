CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -I. -Wno-conversion-null -Wno-deprecated-declarations 

EXEC     = main1
OBJS     = main.o gradient_descent.o

all: $(EXEC)

%.o: %.cpp gradient_descent.h
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $<

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

clean:
	$(RM) *.o $(EXEC)

distclean: clean
	$(RM) $(EXEC)

