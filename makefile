$(wildcard $(SOURCEDIR)/*.tex)
all: cubyboxy test4dalloc ns2d_solver ns3d_solver
.cpp.o:
	g++ -I/usr/include/eigen3/Eigen -c -g -O3 -o $@ $<
cubyboxy: cubyboxy.o mem.o utilities.o vtk.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
test4dalloc: test4dalloc.o mem.o utilities.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
ns2d_solver: ns2d_solver.o ns2d.o vtk.o
	g++ -O3 -o $@ $^
ns3d_solver: ns3d_solver.o ns3d.o vtk.o
	g++ -O3 -o $@ $^
clean:
	rm cubyboxy cubyboxy.o mem.o utilities.o test4dalloc.o
