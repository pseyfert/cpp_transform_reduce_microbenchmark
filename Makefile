LDFLAGS  = -lstdc++ -m64 -g -march=native -mavx2
CPPFLAGS = -march=native -std=c++17 -m64 -O3 -g -Wextra -Wall -Wshadow

CPPFLAGS += -isystem /home/pseyfert/.local/include
CPPFLAGS += -I.
CPPFLAGS += -fopenmp

LDLIBS   += -L/home/pseyfert/.local/lib -lbenchmark
LDLIBS   += -ltbb

all: reduce
