# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall
LDFLAGS = -lm -pthread

# Source files and target
SRCS = experiment.cpp iofunc.cpp logfunc.cpp genfunc.cpp filter.cpp fourier.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = experiment

# Default rule to build the program
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# Rule to compile source files into object files
%.o: %.cpp
	$(CXX) $(DEFINES) $(CXXFLAGS) -c $< -o $@

# Clean rule to remove generated files
clean:
	rm -f $(OBJS) $(TARGET)
