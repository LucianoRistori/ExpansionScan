#==============================================================================
# Makefile for ExpansionScan
#==============================================================================

# Compiler and flags
CXX      = clang++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -Wno-cpp \
           -stdlib=libc++ -pthread -m64 -mmacosx-version-min=13.0

# ROOT paths (uses root-config)
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

# Include dir for shared common module
COMMON_DIR = ../common

INCLUDES = -I$(COMMON_DIR)

# Sources
SRCS = ExpansionScan.cpp \
       $(COMMON_DIR)/Points.cpp

# Objects
OBJS = $(SRCS:.cpp=.o)

# Output binary
TARGET = ExpansionScan

# Default rule
all: $(TARGET)

# Link rule
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(ROOTLIBS) -o $(TARGET)

# Compile rule
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	@echo "Cleaning up..."
	rm -f $(OBJS) $(TARGET)
	find . -name "*.dSYM" -type d -exec rm -rf {} +

.PHONY: all clean
