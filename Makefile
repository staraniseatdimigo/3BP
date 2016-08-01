COPT = -O2

OBJS = main.o phi_math.o phi_type.o
OUT = out

MULTI_OBJS = multi.o phi_math.o phi_type.o
MULTI_OUT = multi

all: $(MULTI_OBJS) $(OBJS)
	g++ $(COPT) $(MULTI_OBJS) -o $(MULTI_OUT)
	g++ $(COPT) $(OBJS) -o $(OUT)

%.o: %.cpp
	g++ $(COPT) -c $^ -o $@ 
