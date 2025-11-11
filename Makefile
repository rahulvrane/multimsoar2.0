# MultiMSOAR 2.0 - Multi-threaded Edition
# Makefile for building the parallel ortholog identification tool

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -O3 -march=native -pthread
LDFLAGS = -pthread

# Directories
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

# Source files
SOURCES = $(SRC_DIR)/MultiMSOARSoftware.cpp
HEADERS = $(SRC_DIR)/Hungarian.h \
          $(SRC_DIR)/TreeCentric.h \
          $(SRC_DIR)/NodeCentric.h \
          $(SRC_DIR)/TreeAnalysis.h \
          $(SRC_DIR)/ThreadPool.h

# Output binary
TARGET = $(BIN_DIR)/MultiMSOAR2.0

# Default target
all: directories $(TARGET)

# Create necessary directories
directories:
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(BUILD_DIR)

# Build the main executable
$(TARGET): $(SOURCES) $(HEADERS)
	@echo "Compiling MultiMSOAR 2.0 (Multi-threaded Edition)..."
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET) $(LDFLAGS)
	@echo "Build complete: $(TARGET)"
	@echo ""
	@echo "Hardware threads available: $$(nproc)"
	@echo ""

# Debug build with sanitizers (for development/testing)
debug: CXXFLAGS = -std=c++11 -Wall -g -O0 -pthread -fsanitize=thread -fsanitize=undefined
debug: LDFLAGS = -pthread -fsanitize=thread -fsanitize=undefined
debug: directories
	@echo "Compiling MultiMSOAR 2.0 (Debug with ThreadSanitizer)..."
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(BIN_DIR)/MultiMSOAR2.0_debug $(LDFLAGS)
	@echo "Debug build complete: $(BIN_DIR)/MultiMSOAR2.0_debug"

# Profile build (with profiling symbols)
profile: CXXFLAGS = -std=c++11 -Wall -O3 -march=native -pthread -pg
profile: LDFLAGS = -pthread -pg
profile: directories
	@echo "Compiling MultiMSOAR 2.0 (Profiling build)..."
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(BIN_DIR)/MultiMSOAR2.0_profile $(LDFLAGS)
	@echo "Profile build complete: $(BIN_DIR)/MultiMSOAR2.0_profile"

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	@echo "Clean complete"

# Install (copy to /usr/local/bin)
install: $(TARGET)
	@echo "Installing MultiMSOAR 2.0 to /usr/local/bin..."
	install -m 755 $(TARGET) /usr/local/bin/
	@echo "Installation complete"

# Uninstall
uninstall:
	@echo "Uninstalling MultiMSOAR 2.0..."
	rm -f /usr/local/bin/MultiMSOAR2.0
	@echo "Uninstallation complete"

# Help target
help:
	@echo "MultiMSOAR 2.0 - Multi-threaded Edition Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  make          - Build optimized release version (default)"
	@echo "  make debug    - Build debug version with ThreadSanitizer"
	@echo "  make profile  - Build profiling version with gprof support"
	@echo "  make clean    - Remove all build artifacts"
	@echo "  make install  - Install to /usr/local/bin (requires sudo)"
	@echo "  make uninstall- Remove from /usr/local/bin (requires sudo)"
	@echo "  make help     - Show this help message"
	@echo ""
	@echo "Compiler flags:"
	@echo "  -std=c++11    : Enable C++11 features (required for threading)"
	@echo "  -O3           : Maximum optimization"
	@echo "  -march=native : Optimize for current CPU architecture"
	@echo "  -pthread      : Enable POSIX threads"
	@echo ""
	@echo "Usage after building:"
	@echo "  $(TARGET) <#species> <speciesTree> <GeneFamily> <GeneInfo> <OrthoGroups>"
	@echo ""

# Phony targets
.PHONY: all directories debug profile clean install uninstall help
