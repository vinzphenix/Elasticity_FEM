# Nom du programme
PROG := deformation

SRC_DIR := src
INC_DIR := include
MODEL_DIR := models
BUILD_DIR := build

# Choix du compilateur
CC := gcc

# Localisation de gmsh (sdk, source ou autre)
GMSH_LIB_DIR := $(HOME)/gmsh/lib
GMSH_INC_DIR := $(HOME)/gmsh/api

# Flags de compilation
CFLAGS := -Wall -O3 #-fsanitize=address

# Chemins vers les dossiers `include`
INC_DIRS := -I $(INC_DIR) -I $(GMSH_INC_DIR)

# Chemins vers les dossiers `lib`
LIB_DIR := -L $(GMSH_LIB_DIR)

# Spécification du runtime path
LDFLAGS := -Wl,-rpath,$(GMSH_LIB_DIR) #-fsanitize=address

# Librairies à linker
LDLIBS := -lgmsh -lopenblas -lm -llapack


SRC_FILES := $(wildcard $(SRC_DIR)/*.c)
MODEL_FILES := $(wildcard $(MODEL_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC_FILES)) \
       $(patsubst $(MODEL_DIR)/%.c, $(BUILD_DIR)/%.o, $(MODEL_FILES))

# Compilation principale
all: | $(BUILD_DIR) $(PROG)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(MODEL_DIR)/%.c | $(BUILD_DIR)
	@echo "Compiling model $<..."
	@$(CC) -g -c $(CFLAGS) $(INC_DIRS) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	@echo "Compiling source $<..."
	@$(CC) -g -c $(CFLAGS) $(INC_DIRS) $< -o $@

# Link final
$(PROG): $(OBJS)
	@echo "Linking $(PROG)..."
	@$(CC) -g -o $@ $(OBJS) $(LIB_DIR) $(LDLIBS) $(LDFLAGS)

# Nettoyage
clean:
	rm -f $(PROG) $(OBJS)
