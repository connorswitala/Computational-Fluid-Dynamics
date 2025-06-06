# === Compiler and Flags ===
CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall

BIN_DIR := programs

# === Include Paths ===
COMMON_INCLUDES := -Igridlib -Ilinalglib -Iwritefilelib
PERF_INCLUDES := -Iperfgaslib
REAL_INCLUDES := -Irealgaslib

# === Source Files ===
COMMON_LIB_SOURCES := \
    gridlib/grid.cpp \
    linalglib/linalg.cpp \
    writefilelib/writefile.cpp

PERF_SOURCE := perfgassolver/main.cpp
REAL_SOURCE := realgassolver/main.cpp
GIBBS_SOURCE := gibbs/main.cpp
TESTING_SOURCE := testing/main.cpp

PERF_LIB := perfgaslib/perfgas.cpp
REAL_LIB := realgaslib/realgas.cpp

# === Object Directories ===
OBJ_DIR := build

# === Targets ===
PERF_TARGET := $(BIN_DIR)/solve_perf_gas
REAL_TARGET := $(BIN_DIR)/solve_real_gas
TESTING_TARGET := $(BIN_DIR)/test_program
GIBBS_TARGET := $(BIN_DIR)/gibbs_minimizer

# === Object Files ===
PERF_OBJECTS := $(COMMON_LIB_SOURCES:%.cpp=$(OBJ_DIR)/%.o) $(OBJ_DIR)/perfgaslib/perfgas.o $(OBJ_DIR)/perfgassolver_main.o
REAL_OBJECTS := $(COMMON_LIB_SOURCES:%.cpp=$(OBJ_DIR)/%.o) $(OBJ_DIR)/realgaslib/realgas.o $(OBJ_DIR)/realgassolver_main.o
GIBBS_OBJECTS := $(COMMON_LIB_SOURCES:%.cpp=$(OBJ_DIR)/%.o) $(OBJ_DIR)/realgaslib/realgas.o $(OBJ_DIR)/gibbs_main.o
TESTING_OBJECTS := $(COMMON_LIB_SOURCES:%.cpp=$(OBJ_DIR)/%.o) $(OBJ_DIR)/testing_main.o

# === Default Rule ===
all: $(PERF_TARGET) $(REAL_TARGET) $(TESTING_TARGET) $(GIBBS_TARGET)

# === Compile Main Files ===
$(OBJ_DIR)/perfgassolver_main.o: $(PERF_SOURCE)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(COMMON_INCLUDES) $(PERF_INCLUDES) -c $< -o $@

$(OBJ_DIR)/realgassolver_main.o: $(REAL_SOURCE)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(COMMON_INCLUDES) $(REAL_INCLUDES) -c $< -o $@

$(OBJ_DIR)/gibbs_main.o: $(GIBBS_SOURCE)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(COMMON_INCLUDES) $(REAL_INCLUDES) -c $< -o $@

$(OBJ_DIR)/testing_main.o: $(TESTING_SOURCE)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(COMMON_INCLUDES) -c $< -o $@

# === Compile Library Sources ===
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(COMMON_INCLUDES) $(PERF_INCLUDES) $(REAL_INCLUDES) -c $< -o $@

# === Link Executables ===
$(PERF_TARGET): $(PERF_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $^ -o $@

$(REAL_TARGET): $(REAL_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $^ -o $@

$(GIBBS_TARGET): $(GIBBS_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $^ -o $@

$(TESTING_TARGET): $(TESTING_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $^ -o $@

# === Clean Rule ===
clean:
	rm -rf $(OBJ_DIR) $(PERF_TARGET) $(REAL_TARGET) $(TESTING_TARGET) $(GIBBS_TARGET) $(BIN_DIR)

.PHONY: all clean
