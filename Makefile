CC = gcc -Wall -O3 -funroll-loops -fomit-frame-pointer `./gcccpuopt`
CC = gcc -g
all:
	$(CC) liquid.c -o liquid -lGL -lGLU `sdl-config --cflags --libs`

clean:
	@echo Cleaning up...
	@rm liquid
	@echo Done.
