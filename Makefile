# === Compiler and Flags ===
CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall

# === Include Paths ===
INCLUDES := -Iperfgaslib -Igridlib -Ilinalglib -Iwritefilelib

# === Source Files ===
LIB_SOURCES := \
    perfgaslib/perfgas.cpp \
    gridlib/grid.cpp \
    linalglib/linalg.cpp \
    writefilelib/writefile.cpp

MAIN_SOURCE := perfgassolver/main.cpp

# === Object Files ===
OBJ_DIR := build
LIB_OBJECTS := $(LIB_SOURCES:%.cpp=$(OBJ_DIR)/%.o)
MAIN_OBJECT := $(OBJ_DIR)/perfgassolver_main.o

# === Target Executable ===
TARGET := solve_perf_gas

# === Default Rule ===
all: $(TARGET)

# === Compile Main ===
$(MAIN_OBJECT): $(MAIN_SOURCE)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# === Compile Libraries ===
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# === Link Executable ===
$(TARGET): $(LIB_OBJECTS) $(MAIN_OBJECT)
	$(CXX) $^ -o $@

# === Clean Rule ===
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean
