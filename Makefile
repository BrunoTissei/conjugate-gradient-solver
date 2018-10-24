FLAGS = -lm -O3 -mavx -march=native -Wall -Wextra
OBJS = obj/main.o obj/helper.o obj/solver.o obj/algorithm.o
PROG = cgSolver
CXX = gcc

all: $(PROG)

obj/%.o: src/%.c
	@mkdir -p obj
	$(CXX) $(CFLAGS) -c -s $< $(FLAGS)
	@mv *.o obj/

$(PROG): $(OBJS)
	$(CXX) $(OBJS) -o $(PROG) $(FLAGS) 

clean:
	rm -f $(PROG)
	rm -rf obj/
