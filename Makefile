CC = gcc -Wall -O3 -funroll-all-loops -fomit-frame-pointer `./gcccpuopt`
#CC = gcc -O3 -march=core2 -mfpmath=sse,387 -fomit-frame-pointer -funroll-all-loops
#CC = gcc -g
all:
	$(CC) liquid.c -o liquid -lGL -lGLU `sdl-config --cflags --libs`

clean:
	@echo Cleaning up...
	@rm liquid
	@echo Done.
