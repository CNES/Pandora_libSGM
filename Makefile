#
# Copyright (c) 2020 Centre National d'Etudes Spatiales (CNES).
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
# Variables for compilation
CC=gcc
CFLAGS=-Wall -ansi -pedantic -Wextra -g -fdiagnostics-show-option -o2
CXX=g++
# Project's variables
SOURCES=sources/lib/sgm.cpp
EXEC=sgm

# Points to the root of Google Test, relative to where this file is.
GTEST_DIR=./test

# Where to find user code.
SRC_DIR=./sources/lib
TEST_DIR=./test

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
$(EXEC): $(SOURCES)
	$(CXX) $^ -o $(EXEC)

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

gtest.a: gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a: gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Run unit tests
tests: run_functions_unittest

# Build tests
sgm.o: $(SRC_DIR)/sgm.cpp $(SRC_DIR)/sgm.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(SRC_DIR)/sgm.cpp

functions_unittest.o: $(TEST_DIR)/functions_unittest.cpp \
                     $(SRC_DIR)/sgm.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(TEST_DIR)/functions_unittest.cpp

functions_unittest: sgm.o functions_unittest.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@
	
run_functions_unittest: functions_unittest
	./$^ --gtest_output=xml:$^.gtest.xml

# Clean the project
clean:
	rm -f $(TESTS) $(OBJECTS) $(EXEC) gtest.a gtest_main.a *.o *.gtest.xml *.gcno *.gcda

