#
# Copyright (c) 2024 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of LIBSGM
#
#     https://github.com/CNES/Pandora_libsgm
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# Autodocumented Makefile for python and C++ dev, see two sections.
# see: https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
# Dependencies : python3 venv g++ gcc ar
# Some Makefile global variables can be set in make command line: VENV, PYTHON, ...
# Recall: .PHONY  defines special targets not associated with files

############### GLOBAL VARIABLES ######################

.DEFAULT_GOAL := help
# Set shell to BASH
SHELL := /bin/bash

########################## Python project dev targets ##################

############### Python GLOBAL VARIABLES ######################


# Set Virtualenv directory name
# Example: VENV="other-venv/" make install
ifndef VENV
	VENV = "venv"
endif

# Browser definition for sphinx and coverage
define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"

# Python global variables definition
PYTHON_VERSION_MIN = 3.8
# Set PYTHON if not defined in command line
# Example: PYTHON="python3.10" make venv to use python 3.10 for the venv
# By default the default python3 of the system.
ifndef PYTHON
	PYTHON = "python3"
endif
PYTHON_CMD=$(shell command -v $(PYTHON))

PYTHON_VERSION_CUR=$(shell $(PYTHON_CMD) -c 'import sys; print("%d.%d"% sys.version_info[0:2])')
PYTHON_VERSION_OK=$(shell $(PYTHON_CMD) -c 'import sys; cur_ver = sys.version_info[0:2]; min_ver = tuple(map(int, "$(PYTHON_VERSION_MIN)".split("."))); print(int(cur_ver >= min_ver))')

############### Check python version supported ############

ifeq (, $(PYTHON_CMD))
    $(error "PYTHON_CMD=$(PYTHON_CMD) not found in $(PATH)")
endif

ifeq ($(PYTHON_VERSION_OK), 0)
    $(error "Requires python version >= $(PYTHON_VERSION_MIN). Current version is $(PYTHON_VERSION_CUR)")
endif

################ MAKE targets by sections ######################

help: ## this help
	@echo "                   LIBSGM MAKE HELP"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

.PHONY: venv
venv: ## create virtualenv in "venv" dir if not exists
	@test -d ${VENV} || $(PYTHON_CMD) -m venv ${VENV}
	@touch ${VENV}/bin/activate
	@${VENV}/bin/python -m pip install --upgrade wheel setuptools pip # no check to upgrade each time


.PHONY: install
install: venv  ## install environment for development target (depends venv)
	@test -f ${VENV}/bin/pylibsgm || echo "Install LibSGM package from local directory"
	@test -f ${VENV}/bin/pylibsgm || ${VENV}/bin/pip install -e .[dev,docs]
	@test -f .git/hooks/pre-commit || echo "Install pre-commit"
	@test -f .git/hooks/pre-commit || ${VENV}/bin/pre-commit install -t pre-commit
	@test -f .git/hooks/pre-push || ${VENV}/bin/pre-commit install -t pre-push
	@echo "Libsgm installed in dev mode in virtualenv ${VENV}"
	@echo "Libsgm venv usage : source ${VENV}/bin/activate; python3 -c 'import libSGM'"

## python Test section

.PHONY: test-python
test-python: ## run only python tests and coverage quickly with the default Python
	@${VENV}/bin/pytest -o log_cli=true --cov-config=.coveragerc --cov --cov-report=term-missing

.PHONY: coverage-python
coverage-python: ## check code coverage quickly with the default Python
	@${VENV}/bin/coverage run -m pytest
	@${VENV}/bin/coverage report -m
	@${VENV}/bin/coverage html
	$(BROWSER) htmlcov/index.html

## Code quality, linting section

### Format with isort and black

.PHONY: format
format: format/isort format/black  ## run black and isort formatting

.PHONY: format/isort
format/isort: ## run isort formatting 
	@echo "+ $@"
	@${VENV}/bin/isort src/libsgm_python src/libSGM tests

.PHONY: format/black
format/black: ## run black formatting
	@echo "+ $@"
	@${VENV}/bin/black src/libsgm_python src/libSGM tests

### Check code quality and linting : isort, black, pylint

.PHONY: lint
lint: lint/isort lint/black lint/pylint lint/mypy ## check code quality and linting

.PHONY: lint/isort
lint/isort: ## check imports style with isort
	@echo "+ $@"
	@${VENV}/bin/isort --check src/libsgm_python src/libSGM tests

.PHONY: lint/black
lint/black: ## check global style with black
	@echo "+ $@"
	@${VENV}/bin/black --check src/libsgm_python src/libSGM tests

.PHONY: lint/pylint
lint/pylint: ## check linting with pylint
	@echo "+ $@"
	@set -o pipefail; ${VENV}/bin/pylint src/libsgm_python src/libSGM tests --rcfile=.pylintrc --output-format=parseable | tee pylint-report.txt # pipefail to propagate pylint exit code in bash

.PHONY: lint/mypy
lint/mypy: ## check linting type hints with mypy
	@echo "+ $@"
	@${VENV}/bin/mypy src/libsgm_python tests

## Documentation section

.PHONY: docs-python
docs-python:  ## generate Sphinx HTML documentation, including API docs
	@${VENV}/bin/sphinx-build -M clean docs/source/ docs/build
	@${VENV}/bin/sphinx-build -M html docs/source/ docs/build -W --keep-going
	$(BROWSER) docs/build/html/index.html


## Release section

.PHONY: dist
dist: clean install ## clean, install, builds source and wheel package
	@${VENV}/bin/python -m pip install --upgrade build
	@${VENV}/bin/python -m build
	ls -l dist

.PHONY: release
release: dist ## package and upload a release
	@${VENV}/bin/twine check dist/*
	@${VENV}/bin/twine upload dist/* --verbose ##  update your .pypirc accordingly

## Clean section

.PHONY: clean-python
clean-python: clean-venv clean-build clean-precommit clean-pyc clean-cython clean-test clean-lint clean-docs ## clean all python dev env

.PHONY: clean-venv
clean-venv: ## clean venv
	@echo "+ $@"
	@rm -rf ${VENV}

.PHONY: clean-build
clean-build: ## clean build artifacts
	@echo "+ $@"
	@rm -fr build/
	@rm -fr dist/
	@rm -fr .eggs/
	@find . -name '*.egg-info' -exec rm -fr {} +
	@find . -name '*.egg' -exec rm -f {} +

.PHONY: clean-precommit
clean-precommit: ## clean precommit hooks in .git/hooks
	@rm -f .git/hooks/pre-commit
	@rm -f .git/hooks/pre-push

.PHONY: clean-pyc
clean-pyc: ## clean Python file artifacts
	@echo "+ $@"
	@find . -type f -name "*.py[co]" -exec rm -fr {} +
	@find . -type d -name "__pycache__" -exec rm -fr {} +
	@find . -name '*~' -exec rm -fr {} +

.PHONY: clean-cython
clean-cython: ## clean Python file artifacts
	@echo "+ $@"
	@rm -f src/libSGM/sgm_wrapper.cpython-38-x86_64-linux-gnu.so

.PHONY: clean-test
clean-test: ## clean test and coverage artifacts
	@echo "+ $@"
	@rm -fr .tox/
	@rm -f .coverage
	@rm -rf .coverage.*
	@rm -rf coverage.xml
	@rm -fr htmlcov/
	@rm -fr .pytest_cache
	@rm -f pytest-report.xml
	@find . -type f -name "debug.log" -exec rm -fr {} +

.PHONY: clean-lint
clean-lint: ## clean linting artifacts
	@echo "+ $@"
	@rm -f pylint-report.txt
	@rm -f pylint-report.xml
	@rm -rf .mypy_cache/

.PHONY: clean-docs
clean-docs: ## clean builded documentations
	@echo "+ $@"
	@rm -rf docs/build/
	@rm -rf docs/source/api_reference/
	@rm -rf docs/doxyoutput/

########################## C++ libSGM compilation targets ##################

# Variables for compilation
CC=gcc
CFLAGS=-Wall -ansi -pedantic -Wextra -g -fdiagnostics-show-option -o2
CXX=g++
# Project's variables
SOURCES=src/libSGM/lib/sgm.cpp
EXEC=sgm

# Points to the root of Google Test, relative to where this file is.
GTEST_DIR=./tests

# Where to find user code.
SRC_DIR=./src/libSGM/lib
TEST_DIR=./tests

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra -pthread -fprofile-arcs -ftest-coverage -g -O0 --coverage -std=c++11

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS=functions_unittest

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS=$(GTEST_DIR)/include/gtest/*.h \
              $(GTEST_DIR)/include/gtest/internal/*.h

# Tests
TESTS_EXEC=tests
TESTS_REPORT=tests-report.xml

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_=$(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# Build the project
# TODO: build a testing main for dev, no compilation for now
#$(EXEC): $(SOURCES)
#	$(CXX) $^ -o $(EXEC)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o: $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o: $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest_main.a: gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

.PHONY: test-cpp 
test-cpp: run_functions_unittest ## Run libSGM C++ unit tests

# Build tests
sgm.o: $(SRC_DIR)/sgm.cpp $(SRC_DIR)/sgm.hpp $(GTEST_HEADERS) ## Generate libsgm C++ library
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(SRC_DIR)/sgm.cpp

functions_unittest.o: $(TEST_DIR)/functions_unittest.cpp \
                     $(SRC_DIR)/sgm.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(TEST_DIR)/functions_unittest.cpp

functions_unittest: sgm.o functions_unittest.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@
	
run_functions_unittest: functions_unittest
	./$^ --gtest_output=xml:$^.gtest.xml

.PHONY: coverage-cpp 
coverage-cpp: install  ## Gcovr (depends on gcovr in venv)
	@gcovr --sonarqube -r . -f src > gcovr-report.xml

.PHONY: cppcheck
cppcheck: ## C++ cppcheck for CI (depends cppcheck)
	@cppcheck -v --enable=all --xml -Isrc/libSGM/lib src/libSGM/lib/*.cpp 2> cppcheck-report.xml

.PHONY: docs-cpp
docs-cpp: ## C++ doxygen doc generation (depends doxygen)
	@doxygen docs/Doxyfile


.PHONY: clean-cpp
clean-cpp: ## Clean C++ libSGM project
	@echo "+ $@"
	@rm -f $(TESTS) $(OBJECTS) $(EXEC)  gtest_main.a *.o *.gtest.xml *.gcno *.gcda


########################## GLOBAL Targets (python + cpp) ######################

.PHONY: test
test: test-cpp test-python ## run all tests: cpp tests (libSGM) then python tests (libsgm_python)

.PHONY: coverage
coverage: coverage-python coverage-cpp ## check code coverage C++ + python


.PHONY: coverage
docs: docs-cpp docs-python ## cpp doxygen and python sphinx doc generation

.PHONY: clean
clean: clean-python clean-cpp ## clean all python + C++
