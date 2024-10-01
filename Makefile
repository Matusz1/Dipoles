CC = gcc
CFLAGS = -std=c11 -O3 -pedantic -Wall
LDLIBS = -lm -llapacke -lgsl -lgslcblas

SRC_DIR = ./src
OBJ_DIR = ./obj

SRCS := util.c polynomial.c hermite.c test.c state.c fock.c \
	workspace.c
SRCS := $(addprefix $(SRC_DIR)/,$(SRCS))
OBJS := $(SRCS:$(SRC_DIR)/%=$(OBJ_DIR)/%.o)

all: $(OBJ_DIR) run

run: $(OBJS) $(OBJ_DIR)/main.c.o
	$(CC) $^ $(LDLIBS) -o run

$(OBJ_DIR)/main.c.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.c.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm obj/*.o
