LIBS = -lm
CC = gcc
CFLAGS = -Wall -Wextra -ggdb

NUMERICAL = ../numerical.h

all: demo_trapz demo_ode demo_interp demo_linear_system

demo_trapz: $(NUMERICAL) demo_trapz.c
	$(CC) $(CFLAGS) demo_trapz.c -o demo_trapz $(LIBS)

demo_ode: $(NUMERICAL) demo_ode.c
	$(CC) $(CFLAGS) demo_ode.c -o demo_ode $(LIBS)

demo_interp: $(NUMERICAL) demo_interp.c
	$(CC) $(CFLAGS) demo_interp.c -o demo_interp $(LIBS)

demo_linear_system: $(NUMERICAL) demo_linear_system.c
	$(CC) $(CFLAGS) demo_linear_system.c -o demo_linear_system $(LIBS)

clean:
	rm demo_trapz demo_ode demo_interp demo_linear_system