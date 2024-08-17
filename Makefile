CXX        = g++-13 -std=c++20
WARNINGS   = -Wall -Wextra -Wpedantic
INCLUDE    = -I ./include


test_src   = test_compute_residual.cpp test_eigen_interop.cpp test_sparse_repr.cpp test_stress_direct.cpp
targets   += $(patsubst %.cpp,%,$(test_src))


all: $(targets)


test_compute_residual: build/test_compute_residual.o build/poisson.o build/operator.o build/utils.o
	$(CXX) -o $@ $^


test_eigen_interop: build/test_eigen_interop.o build/utils.o
	$(CXX) -o $@ $^
	./$@


test_sparse_repr: build/test_sparse_repr.o build/poisson.o build/operator.o build/utils.o
	$(CXX) -o $@ $^


test_stress_direct: build/test_stress_direct.o build/poisson.o build/operator.o build/solvers.o build/utils.o
	$(CXX) -o $@ $^


build/%.o: src/%.cpp
	$(CXX) -c $(WARNINGS) $(INCLUDE) -o $@ $^


build/%.o: test/%.cpp
	$(CXX) -c $(WARNINGS) $(INCLUDE) -o $@ $^


.PHONY clean:
clean:
	rm -f build/*.o hello
