GCC= g++
OMPFLAGS = -g -fopenmp

all: 4DScanPoi

io.o: io.c
	$(GCC) -o $@ -c $<
scan.o: scan.c
	$(GCC) $(OMPFLAGS) -o $@ -c $<
mc.o: mc.c
	$(GCC) $(OMPFLAGS) -o $@ -c $< -std=c++11
4DScanPoi.o: 4DScanPoi.c
	$(GCC) -o $@ -c $<
4DScanPoi: 4DScanPoi.o io.o scan.o mc.o
	$(GCC) $(OMPFLAGS) -o ../$@ $+
clean:
	rm -f *.o ../4DScanPoi
