LIBS = -lm
CC = clang
CFLAGS = -Wall -Wextra

NUMERICAL = ../numerical.h

all: demo_trapz demo_ode

demo_trapz: $(NUMERICAL) demo_trapz.c
	$(CC) $(CFLAGS) demo_trapz.c -o demo_trapz $(LIBS)

demo_ode: $(NUMERICAL) demo_ode.c
	$(CC) $(CFLAGS) demo_ode.c -o demo_ode $(LIBS)

clean:
	rm demo