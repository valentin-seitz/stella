# See README.md for more information on these automated stella tests
# The fortran tests have their own make file in the 'test_fortran_routines' subfolder

# Test directory
TEST_DIR = $(PROJECT_DIR)/AUTOMATIC_TESTS
.PHONY: create-test-virtualenv check-tests-dependencies numerical-tests numerical-tests-verbose benchmarks benchmarks-slow

# Create a virtualenv and install the test dependencies there
create-test-virtualenv:
	bash -c "cd $(TEST_DIR); \
		python3 -m venv venv; \
		source venv/bin/activate; \
		pip install -r requirements.txt; \
		deactivate; \
	"
	
# Check that all the Python dependencies are importable
check-tests-dependencies:
	@echo "Running: check-tests-dependencies"
	@python3 $(TEST_DIR)/check_automated_stella_tests_requirements.py
	
# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
numerical-tests: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests"
	@pytest -v -rF -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests/test_1_whether_stella_runs
	
# Run the automated stella tests (the -r flag shows the full stella output)
# -rA shows all outputs, -rf shows only failed outputs 
numerical-tests-verbose: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests printing all the output"
	@pytest -v -rA -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests
	
	# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
benchmarks: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella benchmarks"
	@pytest -v -rF -W ignore::DeprecationWarning $(TEST_DIR)/quick_benchmarks
	
# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
benchmarks-slow: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella benchmarks"
	@pytest -v -rF -W ignore::DeprecationWarning $(TEST_DIR)/slow_benchmarks
