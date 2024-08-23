CXX         = g++-13 -std=c++20
WARNINGS    = -Wall -Wextra -Wpedantic
INCLUDE     = -I ./include
OPT         = -O2



test_src    = test_compute_residual.cpp test_eigen_interop.cpp test_sparse_repr.cpp test_stress_direct.cpp test_cycle_spec.cpp
example_src = example_smoother_comparison.cpp example_forced_circuit.cpp example_multilevel_performance.cpp example_cycle_performance.cpp
targets    += $(patsubst %.cpp,%,$(test_src))
targets    += $(patsubst %.cpp,%,$(example_src))


all: $(targets)


example_smoother_comparison: build/example_smoother_comparison.o build/poisson.o build/stencil.o build/solvers.o build/utils.o
	$(CXX) -o $@ $^


example_forced_circuit: build/example_forced_circuit.o build/circuit.o build/solvers.o build/stencil.o build/utils.o
	$(CXX) -o $@ $^


example_multilevel_performance: build/example_multilevel_performance.o build/poisson.o build/solvers.o build/stencil.o build/multigrid.o build/utils.o
	$(CXX) -o $@ $^


example_cycle_performance: build/example_cycle_performance.o build/poisson.o build/solvers.o build/stencil.o build/multigrid.o build/utils.o
	$(CXX) -o $@ $^


test_compute_residual: build/test_compute_residual.o build/poisson.o build/stencil.o build/utils.o
	$(CXX) -o $@ $^
	./$@


test_eigen_interop: build/test_eigen_interop.o build/utils.o
	$(CXX) -o $@ $^
	./$@


test_cycle_spec: build/test_cycle_spec.o build/multigrid.o build/solvers.o
	$(CXX) -o $@ $^
	./$@


test_sparse_repr: build/test_sparse_repr.o build/poisson.o build/stencil.o build/utils.o
	$(CXX) -o $@ $^


test_stress_direct: build/test_stress_direct.o build/poisson.o build/stencil.o build/solvers.o build/utils.o
	$(CXX) -o $@ $^


stress: test_stress_direct
	./$< > stress.out
	python3.8 plot/stress_direct.py stress.out


smoother: example_smoother_comparison
	./$< > smoother.out
	python3.8 plot/convergence_history.py smoother.out


circuit: example_forced_circuit
	./$^ > circuit.out
	python3.8 plot/solution_plot_1D.py circuit.out


multilevel: example_multilevel_performance
	./$^ > multilevel.out
	python3.8 plot/convergence_history.py multilevel.out


cycle: example_cycle_performance
	./$^ > cycle.out
	python3.8 plot/convergence_history.py cycle.out


build/%.o: src/%.cpp
	$(CXX) -c $(WARNINGS) $(OPT) $(INCLUDE) -o $@ $^


build/%.o: test/%.cpp
	$(CXX) -c $(WARNINGS) $(OPT) $(INCLUDE) -o $@ $^


build/%.o: examples/%.cpp
	$(CXX) -c $(WARNINGS) $(OPT) $(INCLUDE) -o $@ $^


.PHONY clean:
clean:
	rm -f build/*.o hello
