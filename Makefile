CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)ls

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR) -lMathMore

CPPFLAGS += -g

TARGET1 = multipleDRS4
SRC1 = multipleDRS4.cc
OBJ1 = $(SRC1:.cc=.o)

TARGET2 = singleDRS4
SRC2 = singleDRS4.cc
OBJ2 = $(SRC2:.cc=.o)


all : $(TARGET1) $(TARGET2)

$(TARGET1) : $(OBJ1)
	$(LD) $(CPPFLAGS) -o $(TARGET1) $(OBJ1) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET2) : $(OBJ2)
	$(LD) $(CPPFLAGS) -o $(TARGET2) $(OBJ2) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^


%.o : %.cc	
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@	
	@echo $<

clean :
	rm -f *.o src/*.o $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) *~

