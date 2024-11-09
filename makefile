all: cubyboxy test4dalloc
mem.o: mem.cpp
	g++ -c -g -O3 -o $@ $<
utilities.o: utilities.cpp
	g++ -c -g -O3 -o $@ $<
cubyboxy.o: cubyboxy.cpp mem.o utilities.o
	g++ -c -g -O3 -o $@ $<
cubyboxy: cubyboxy.o mem.o utilities.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
test4dalloc.o: test4dalloc.cpp mem.o utilities.o
	g++ -c -g -O3 -o $@ $<
test4dalloc: test4dalloc.o mem.o utilities.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
poissson: poisson.cpp
	g++ -O3 -o $@ $^
simple: simple.cpp
	g++ -O3 -o $@ $^
matrix_solver: matrix_solver.cpp
	g++ -O3 -o $@ $^
clean:
	rm cubyboxy cubyboxy.o mem.o utilities.o test4dalloc.o
