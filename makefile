$(wildcard $(SOURCEDIR)/*.tex)
all: cubyboxy test4dalloc
.cpp.o:
	g++ -c -g -O0 -o $@ $<
cubyboxy: cubyboxy.o mem.o utilities.o vtk.o
	g++ -O0 -o $@ $^
	objdump -d $@ > assembly.txt
test4dalloc: test4dalloc.o mem.o utilities.o
	g++ -O0 -o $@ $^
	objdump -d $@ > assembly.txt
clean:
	rm cubyboxy cubyboxy.o mem.o utilities.o test4dalloc.o
