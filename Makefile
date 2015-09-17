SOURCES = ring_rate.c conf_file.c parser.c parameters.c ext_inputs.c colormap.c \
	  circ_buffer.c utils.c cubature-20101018/cubature.c
OBJS = $(SOURCES:.c=.o)
CFLAGS = -std=gnu99 -O3 -DHAVE_INLINE=1 --pedantic -W -Winline \
	 -Wstrict-prototypes -Wno-sign-conversion -Wshadow -Wpointer-arith -Wcast-qual \
	 -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums -fno-common

LIBS = -lm -lgsl -lgslcblas -lring_network

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

MAIN = simulate_trajectory

all: libring_network.a $(MAIN) 

libring_network.a: $(OBJS)
	$(AR) rcs $@ $^

$(MAIN): simulate_trajectory.o
	$(CC) -o $@ $^ -L. $(LIBS)

clean:
	rm -f simulate_trajectory.o $(OBJS) libring_network.a
