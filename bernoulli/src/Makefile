GCC= g++
OMPFLAGS = -g -fopenmp

all: 4DScanBer

io.o: io.c
	$(GCC) -o $@ -c $<
scan.o: scan.c
	$(GCC) $(OMPFLAGS) -o $@ -c $<
mc.o: mc.c
	$(GCC) $(OMPFLAGS) -o $@ -c $< -std=c++11
4DScanBer.o: 4DScanBer.c
	$(GCC) -o $@ -c $<
4DScanBer: 4DScanBer.o io.o scan.o mc.o
	$(GCC) $(OMPFLAGS) -o ../$@ $+
clean:
	rm -f *.o ../4DScanBer
