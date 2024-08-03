# See README.md for more information on these automated stella tests
.PHONY: benchmarks benchmarks-slow
	
# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
benchmarks: check-automated-stella-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella benchmarks"
	@pytest -v -rF -W ignore::DeprecationWarning $(PROJECT_DIR)/BENCHMARKS/quick_benchmarks/
	
# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
benchmarks-slow: check-automated-stella-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella benchmarks"
	@pytest -v -rF -W ignore::DeprecationWarning $(PROJECT_DIR)/BENCHMARKS/slow_benchmarks/
	
 
