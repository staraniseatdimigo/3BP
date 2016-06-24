OBJS = main.o phi_math.o phi_type.o
OUT = out

all: $(OBJS)
	g++ $(OBJS) -o $(OUT)

%.o: %.cpp
	g++ -c $^ -o $@

