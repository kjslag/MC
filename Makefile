CXX		= g++
CXXFLAGS	= -pipe -g -std=c++11 -Isrc \
-Wall -Wextra -Wstrict-overflow=5 -Wstrict-aliasing=1 -Wunsafe-loop-optimizations \
-Wpedantic -Wshadow -Wno-conversion -Wdisabled-optimization -Wsuggest-attribute=pure -Wsuggest-attribute=const \
-O4 -march=native -ffast-math -fwhole-program
# -fopt-info-optimized-missed=optinfo
# -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include
LDFLAGS		= -lm -lboost_program_options -rdynamic
# -lglib-2.0
# echo "" | gcc -march=native -v -E - 2>&1 | grep cc1

all: MC

MC: src/MC.cpp src/util.hh
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS) -DDEBUG

callgrind: MC
	valgrind --tool=callgrind ./MC

cachegrind: MC
	valgrind --tool=cachegrind ./MC

clean:
	rm -f MC

run: MC
	echo > jobs
	all=""; \
	for L in 1 10 100 1000; \
	do	for J in 0 0.25 0.5 1 2 4; \
		do	f="results/ising_$${L}_$${J}"; \
			echo -e "$$f:\n\t./MC --L $$L --L $$L --n 6 --J $$J --sweep 100 --file $$f\n" >> jobs; \
			all="$$all $$f"; \
		done; \
	done; \
	echo "all:$$all" >> jobs
	mkdir -p results
	make -j4 -f jobs all
