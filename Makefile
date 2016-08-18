CC = g++
CFLAGS = -I

all: mc-rect-lattice-func-linux1 mc-rect-lattice-func-linux2 mc-rect-lattice-func-linux4

mc1: mc-rect-lattice-func-linux1.cpp
	g++ -o mc-rect-lattice-func-linux1 mc-rect-lattice-func-linux1.cpp -I .

mc2: mc-rect-lattice-func-linux2.cpp
	g++ -o mc-rect-lattice-func-linux2 mc-rect-lattice-func-linux2.cpp -I .

mc4: mc-rect-lattice-func-linux4.cpp
	g++ -o mc-rect-lattice-func-linux4 mc-rect-lattice-func-linux4.cpp -I .

