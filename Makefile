CXX        = g++
WARNINGS   = -Wall -Wextra -Wpedantic
INCLUDE    = -I ./include


all: hello


hello: build/main.o
	$(CXX) -o $@ $^


build/%.o: src/%.cpp
	$(CXX) -c $(WARNINGS) $(INCLUDE) -o $@ $^
