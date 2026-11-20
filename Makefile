# Exact hypervolume: bitmatrix and wfg-bitmatrix.
#
#   make            build the program ./bitmatrix and the shared library
#   make test       check the program on the examples of the paper
#   make test-fronts  check the program on fronts generated as in the experiments
#   make run        run the program on the examples
#   make clean      remove the build products
#
# Usage: ./bitmatrix [-m bitmatrix|wfg-bitmatrix] [-r r1,...,rd] [-k leaf] [-s switch] [-t] file

CC      ?= cc
CFLAGS  ?= -O3 -march=native
CFLAGS  += -std=c99 -Wall -Wextra
LDLIBS  += -lm

SRC_DIR  := src
SOURCES  := $(SRC_DIR)/hypervolume_c.c $(SRC_DIR)/hypervolume_c_ie_pruned.c \
            $(SRC_DIR)/hypervolume_c_imatrix.c $(SRC_DIR)/hypervolume_c_wfg_imatrix.c

ifeq ($(shell uname -s),Darwin)
    LIBRARY := libbitmatrix.dylib
    SHARED  := -dynamiclib
else
    LIBRARY := libbitmatrix.so
    SHARED  := -shared
endif

.PHONY: all test test-fronts run clean

all: bitmatrix $(LIBRARY)

bitmatrix: $(SRC_DIR)/main.c $(SOURCES)
	$(CC) $(CFLAGS) $(SRC_DIR)/main.c -o $@ $(LDLIBS)

# Exports hypervolume_c_imatrix and hypervolume_c_wfg_imatrix (for example, for ctypes).
$(LIBRARY): $(SOURCES)
	$(CC) $(CFLAGS) -fPIC $(SHARED) $(SRC_DIR)/hypervolume_c_wfg_imatrix.c -o $@ $(LDLIBS)

test: bitmatrix
	@check() { \
	    result=$$(./bitmatrix -m $$1 -r $$2 $$3); \
	    if [ "$$result" = "$$4" ]; then echo "ok    $$1 $$3: $$result"; \
	    else echo "FAIL  $$1 $$3: $$result (expected $$4)"; exit 1; fi; \
	}; \
	check bitmatrix     4,4   examples/example_2d.txt 6 && \
	check wfg-bitmatrix 4,4   examples/example_2d.txt 6 && \
	check bitmatrix     4,4,4 examples/example_3d.txt 14.125 && \
	check wfg-bitmatrix 4,4,4 examples/example_3d.txt 14.125

# Fronts generated as in the experiments; each file stores its reference point and
# the hypervolume obtained by all the methods (relative tolerance 1e-9).
test-fronts: bitmatrix
	@for file in examples/fronts/*.txt; do \
	    reference=$$(sed -n 's/^# reference: //p' $$file); \
	    expected=$$(sed -n 's/^# hypervolume: //p' $$file); \
	    for method in bitmatrix wfg-bitmatrix; do \
	        result=$$(./bitmatrix -m $$method -r $$reference $$file); \
	        if awk -v a="$$result" -v b="$$expected" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d <= 1e-9 * b) }'; \
	        then echo "ok    $$method $$file: $$result"; \
	        else echo "FAIL  $$method $$file: $$result (expected $$expected)"; exit 1; fi; \
	    done; \
	done

run: bitmatrix
	./bitmatrix -t -r 4,4 examples/example_2d.txt
	./bitmatrix -t -r 4,4,4 examples/example_3d.txt
	./bitmatrix -t -m wfg-bitmatrix -r 4,4,4 examples/example_3d.txt

clean:
	rm -f bitmatrix libbitmatrix.so libbitmatrix.dylib
