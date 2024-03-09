LIBS = -lm
CC = gcc
CFLAGS = -Wall -Wextra -ggdb

NUMERICAL = ../numerical.h

DEMOS := demo_trapz demo_ode demo_interp demo_linear_system demo_inverse_matrix
all: $(DEMOS)

$(DEMOS): %: %.c $(NUMERICAL)
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

.PHONY: clean
clean:
	rm $(DEMOS)