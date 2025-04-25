$(wildcard $(SOURCEDIR)/*.tex)
all: cubyboxy test4dalloc ns2d_solver ns3d_solver testeigen iterative
.cpp.o:
	g++ -I/usr/include/eigen3/Eigen -std=c++11 -c -g -O3 -o $@ $<
cubyboxy: cubyboxy.o mem.o utilities.o vtk.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
test4dalloc: test4dalloc.o mem.o utilities.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
testeigen: testeigen.o
	g++ -O3 -o $@ $^
ns2d_solver: ns2d_solver.o ns2d.o vtk.o
	g++ -O3 -o $@ $^
ns3d_solver: ns3d_solver.o ns3d.o vtk.o
	g++ -O3 -o $@ $^
clean:
	rm cubyboxy cubyboxy.o mem.o utilities.o test4dalloc.o iterative iterative.o 
iterative: iterative.o utilities.o mem.o
	g++ -O3 -o $@ $^