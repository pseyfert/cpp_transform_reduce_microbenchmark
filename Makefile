LDFLAGS  = -lstdc++ -m64 -g -march=native
CPPFLAGS = -march=native -std=c++17 -m64 -O3 -g -Wextra -Wall -Wshadow

CPPFLAGS += -isystem /home/pseyfert/.local/include
CPPFLAGS += -fopenmp

LDLIBS   += -L/home/pseyfert/.local/lib -lbenchmark
LDLIBS   += -ltbb

all: reduce
