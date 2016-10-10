GCC= g++
OMPFLAGS = -g -fopenmp

all: 4DScanBerAgg

io.o: io.c
	$(GCC) -o $@ -c $<
scan.o: scan.c
	$(GCC) $(OMPFLAGS) -o $@ -c $<
mc.o: mc.c
	$(GCC) $(OMPFLAGS) -o $@ -c $< -std=c++11
4DScanBerAgg.o: 4DScanBerAgg.c
	$(GCC) -o $@ -c $<
4DScanBerAgg: 4DScanBerAgg.o io.o scan.o mc.o
	$(GCC) $(OMPFLAGS) -o ../$@ $+
clean:
	rm -f *.o ../4DScanBerAgg
