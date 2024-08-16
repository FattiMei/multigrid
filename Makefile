CXX        = g++-13 -std=c++20
WARNINGS   = -Wall -Wextra -Wpedantic
INCLUDE    = -I ./include


all: hello convergence sparse


hello:
	$(CXX) $(INCLUDE) -o $@ src/main.cpp src/poisson.cpp src/utils.cpp src/operator.cpp


convergence: build/convergence_history.o build/utils.o build/poisson1D.o build/smoothers.o build/solvers.o build/multigrid.o
	$(CXX) -o $@ $^


mgperf: build/multilevel_performance.o build/utils.o build/poisson1D.o build/smoothers.o build/solvers.o build/multigrid.o
	$(CXX) -o $@ $^


operators: build/critical_operators.o build/utils.o build/poisson1D.o build/smoothers.o build/solvers.o build/multigrid.o
	$(CXX) -o $@ $^


sparse: build/test_sparse.o build/stencil.o
	$(CXX) -o $@ $^


test_convergence: convergence
	./$^ > convergence.out
	python plot/convergence_history.py convergence.out


test_mgperf: mgperf
	./$^ > mgperf.out
	python plot/convergence_history.py mgperf.out


test_operators: operators
	./$^ > operators.out
	python plot/convergence_history.py operators.out


test_sparse: sparse
	./$^



build/%.o: src/%.cpp
	$(CXX) -c $(WARNINGS) $(INCLUDE) -o $@ $^


build/%.o: test/%.cpp
	$(CXX) -c $(WARNINGS) $(INCLUDE) -o $@ $^


.PHONY clean:
clean:
	rm -f build/*.o hello
