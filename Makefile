CXX		= clang++
CXXFLAGS	= -pipe -g -std=c++11 \
-Isrc \
-Wall -Wextra -Wstrict-overflow=5 -Winit-self -Wcast-align \
-Wno-missing-field-initializers \
-O3 -march=native -mfpmath=sse -ffast-math
# -Wno-ignored-qualifiers -Wno-comment
# -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include \
# -Wunsafe-loop-optimizations
# -fstrict-aliasing -fstrict-overflow \
# -fvisibility-inlines-hidden
# -I. -I/usr/include \
# -g3 -fno-rtti -ftree-vectorizer-verbose=7 -ftree-vectorize -falign-functions -fno-implement-inlines
# -fomit-frame-pointer -flto -fwhole-program -Ofast
LDFLAGS		= -lm -lboost_program_options -rdynamic
# -lglib-2.0
# echo "" | gcc -march=native -v -E - 2>&1 | grep cc1

all: MC

MC: src/MC.cpp src/util.hh Makefile
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS) -DDEBUG

callgrind: MC
	valgrind --tool=callgrind ./MC

cachegrind: MC
	valgrind --tool=cachegrind ./MC

clean:
	rm MC
