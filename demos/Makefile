LIBS = -lm
CC = gcc
CFLAGS = -Wall -Wextra -ggdb

NUMERICAL = ../numerical.h

DEMOS := demo_trapz demo_ode demo_ode2d demo_interp demo_linear_system demo_inverse_matrix demo_root_find demo_fft demo_opt
all: $(DEMOS)

$(DEMOS): %: %.c $(NUMERICAL)
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

.PHONY: clean
clean:
	rm $(DEMOS)