CXX        = g++ -std=c++17
WARNINGS   = -Wall -Wextra -Wpedantic
INCLUDE    = -I ./include


all: hello convergence


hello: build/main.o build/utils.o build/poisson1D.o build/smoothers.o build/solvers.o
	$(CXX) -o $@ $^


convergence: build/convergence_history.o build/utils.o build/poisson1D.o build/smoothers.o build/solvers.o build/multigrid.o
	$(CXX) -o $@ $^


mgperf: build/multilevel_performance.o build/utils.o build/poisson1D.o build/smoothers.o build/solvers.o build/multigrid.o
	$(CXX) -o $@ $^


test_convergence: convergence
	./$^ > convergence.out
	python plot/convergence_history.py convergence.out


test_mgperf: mgperf
	./$^ > mgperf.out
	python plot/convergence_history.py mgperf.out


build/%.o: src/%.cpp
	$(CXX) -c $(WARNINGS) $(INCLUDE) -o $@ $^


build/%.o: test/%.cpp
	$(CXX) -c $(WARNINGS) $(INCLUDE) -o $@ $^


.PHONY clean:
clean:
	rm -f build/*.o hello
