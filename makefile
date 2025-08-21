$(wildcard $(SOURCEDIR)/*.tex)
all: test4dalloc ns2d_solver ns3d_solver demo staggerednavierdraft
.cpp.o:
	g++ -I/usr/include/eigen3/Eigen -c -g -O3 -o $@ $<
midpointsphere: midpointsphere.o 
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
copycube: copycube.o mem.o utilities.o vtk.o midpointsphere.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
test4dalloc: test4dalloc.o mem.o utilities.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
ns2d_solver: ns2d_solver.o ns2d.o vtk.o
	g++ -O3 -o $@ $^
ns3d_solver: ns3d_solver.o ns3d.o vtk.o
	g++ -O3 -o $@ $^
demo: demo.o
	g++ -O3 -o $@ $^
clean:
	rm copycube copycube.o mem.o utilities.o test4dalloc.o
staggerednavierdraft: staggerednavierdraft.o copycube.o mem.o utilities.o vtk.o midpointsphere.o
	g++ -g -O3 -o $@ $^
	objdump -d $@ > assembly.txt