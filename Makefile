CXX		= g++
CXXFLAGS	= -pipe -g -std=gnu++14 -Isrc \
-Wall -Wextra -Wstrict-aliasing=1 \
-Wpedantic -Wshadow -Wdisabled-optimization \
-Wno-missing-braces -Wno-missing-field-initializers -Wno-unused-parameter \
-O4 -march=native -ffast-math -fwhole-program -funsafe-loop-optimizations
# -Wunsafe-loop-optimizations -Wstrict-overflow=5 -Wno-conversion
# -fopt-info-optimized-missed=optinfo
# -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include
LDFLAGS		= -lm -lboost_program_options
# -rdynamic
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
	make -j4 -f jobs all
