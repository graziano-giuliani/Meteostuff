all::edge2nc

ETOR=etor-v1.3
NETCDF=/usr

CC=gcc
CFLAGS=-Wall -pedantic -O3 -mtune=native -fomit-frame-pointer -I$(ETOR) -I$(NETCDF)/include

$(ETOR)/libetor.a:
	cd $(ETOR) && make

edge2nc: edge2nc.c $(ETOR)/libetor.a
	$(CC) $(CFLAGS) $< -o $@ -lm -L$(ETOR) -letor \
               -L$(NETCDF)/lib64 -lnetcdf

clean:
	cd $(ETOR) && make clean
	rm -f edge2nc
