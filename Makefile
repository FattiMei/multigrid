CXX        = g++
WARNINGS   = -Wall -Wextra -Wpedantic
INCLUDE    = -I ./include


all: hello


hello: build/main.o build/utils.o build/poisson1D.o build/smoothers.o
	$(CXX) -o $@ $^


build/%.o: src/%.cpp
	$(CXX) -c $(WARNINGS) $(INCLUDE) -o $@ $^


.PHONY clean:
clean:
	rm -f build/*.o hello
