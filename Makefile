CC = g++
CFLAGS = -I

all: mc-rect-lattice-func-linux1 mc-rect-lattice-func-linux2 mc-rect-lattice-func-linux4

mc1: mc-rect-lattice-func-linux1.cpp
	g++ -o mc-rect-lattice-func-linux1 mc-rect-lattice-func-linux1.cpp -I .

mc2: mc-rect-lattice-func-linux2.cpp
	g++ -o mc-rect-lattice-func-linux2 mc-rect-lattice-func-linux2.cpp -I .

mc4: mc-rect-lattice-func-linux4.cpp
	g++ -o mc-rect-lattice-func-linux4 mc-rect-lattice-func-linux4.cpp -I .

mc5: mc-rect-lattice-func-linux5.cpp
	g++ -o mc-rect-lattice-func-linux5 mc-rect-lattice-func-linux5.cpp -I .

mc6: mc-rect-lattice-func-linux6.cpp
	g++ -std=c++11 -o mc-rect-lattice-func-linux6 mc-rect-lattice-func-linux6.cpp -I .

mc7: mc-rect-lattice-func-linux7.cpp
	g++ -std=c++11 -o mc-rect-lattice-func-linux7 mc-rect-lattice-func-linux7.cpp -I .

mc8: mc-rect-lattice-func-linux8.cpp
	g++ -std=c++11 -o mc-rect-lattice-func-linux8 mc-rect-lattice-func-linux8.cpp -I .

mc9: mc-rect-lattice-func-linux9.cpp
	g++ -std=c++11 -o mc-rect-lattice-func-linux9 mc-rect-lattice-func-linux9.cpp -I .

test: test_read_file.cpp
	g++ -o test_read_file test_read_file.cpp -I .
