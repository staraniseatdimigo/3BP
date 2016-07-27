COPT = -O2
OBJS = main.o phi_math.o phi_type.o
OUT = out

all: $(OBJS)
	g++ $(COPT) $(OBJS) -o $(OUT)
	g++ $(COPT) multi.cpp -o multi

%.o: %.cpp
	g++ $(COPT) -c $^ -o $@ 
