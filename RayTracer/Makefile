CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -lfreeimage -L./lib/mac/
endif

RM = /bin/rm -f
all: RayTracer.o
RayTracer.o: RayTracer.cpp
	$(CC) $(CFLAGS) RayTracer.cpp -o RayTracer.o
clean:
	$RM *.o RayTracer