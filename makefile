$(wildcard $(SOURCEDIR)/*.tex)
#CXXFLAGS=-mtune=barcelona -fpermissive -fprofile-generate -march=native
CXXFLAGS=-mtune=barcelona -fpermissive -march=native -O0
LDFLAGS=-O0
all: cubyboxy test4dalloc
.cpp.o:
	g++ $(CXXFLAGS) -c -g -o $@ $<
cubyboxy: cubyboxy.o mem.o utilities.o vtk.o
	g++ $(LDFLAGS) -o $@ $^
	objdump -d $@ > assembly.txt
test4dalloc: test4dalloc.o mem.o utilities.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
clean:
	rm cubyboxy cubyboxy.o mem.o utilities.o test4dalloc.o
