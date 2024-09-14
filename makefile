all: cubyboxy test4dalloc
mem.o: mem.cpp
	g++ -c -g -O3 -o $@ $<
cubyboxy.o: cubyboxy.cpp mem.o
	g++ -c -g -O3 -o $@ $<
cubyboxy: cubyboxy.o mem.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
test4dalloc.o: test4dalloc.cpp mem.o
	g++ -c -g -O3 -o $@ $<
test4dalloc: test4dalloc.o mem.o
	g++ -O3 -o $@ $^
	objdump -d $@ > assembly.txt
clean:
	rm cubyboxy cubyboxy.o mem.o
