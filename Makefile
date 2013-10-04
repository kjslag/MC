CXX		= g++
CXXFLAGS	= -pipe -g -std=c++11 -Isrc \
-Wall -Wextra -Wstrict-aliasing=1 \
-Wpedantic -Wshadow -Wdisabled-optimization -Wsuggest-attribute=pure -Wsuggest-attribute=const \
-Wno-conversion -Wno-missing-braces \
-O4 -march=native -ffast-math -fwhole-program -funsafe-loop-optimizations
# -Wunsafe-loop-optimizations -Wstrict-overflow=5
# -fopt-info-optimized-missed=optinfo
# -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include
LDFLAGS		= -lm -lboost_program_options -rdynamic
# -lglib-2.0
# echo "" | gcc -march=native -v -E - 2>&1 | grep cc1

all: MC

MC: src/MC.cpp src/util.hh
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS) -DDEBUG

# valgrind --tool=callgrind  ./MC
# valgrind --tool=cachegrind ./MC

clean:
	rm -f MC

run: MC
	@echo > jobs
	@all=""; \
	   for method in smart; \
	do for L in 20; \
	do for J in `seq -1 .5 +4`; \
	do for u in `seq -4 .5 +4`; \
	do	f="results/vison_hexagon_sVBS_$${L}_$${J}_$${u}_$${method}"; \
		echo -e "$$f:\n\t./MC --L $$L $$L $$L --J $$J --J $$u --potential 'vison hexagon s-VBS' --sweeps 1000 --update-method $$method --file $$f\n" >> jobs; \
		all="$$all $$f"; \
	done; done; done; done; \
	echo "all:$$all" >> jobs
	@mkdir -p results
	make -j4 -f jobs all
