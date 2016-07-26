OBJS = main.o phi_math.o phi_type.o
OUT = out

all: $(OBJS)
	g++ $(OBJS) -o $(OUT)
	g++ multi.cpp -o multi

%.o: %.cpp
	g++ -c $^ -o $@

