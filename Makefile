CC = gcc -Wall -O2

all:
	$(CC) fluid.c -o fluid -lGL -lGLU `sdl-config --cflags --libs`

clean:
	@echo Cleaning up...
	@rm fluid
	@echo Done.
