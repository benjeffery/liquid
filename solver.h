#ifndef SOLVER_H
#define SOLVER_H

#include <stdlib.h>
#include <stdio.h>

#define TRUE 1
#define FALSE 0

#define NUM_PARTICLES_PER_CELL 4
#define SIZE 200
#define DT 0.1f
#define TOPBOTTOM 1
#define LEFTRIGHT 2
#define BOTH 3
#define DIFFUSION 0.0001

#define swap(x,y) {ValueArray* temp = x; x=y; y=temp;}

typedef float Value;
typedef Value ValueArray[SIZE][SIZE];
typedef int IntArray[SIZE][SIZE];

ValueArray* gu;
ValueArray* gv;
ValueArray* gu_old;
ValueArray* gv_old;
IntArray*   has_fluid;

typedef struct {
  float x;
  float y;
  float u;
  float v;
} Particle;

typedef Particle ParticleArray[50000];
ParticleArray* particles;

int num_particles;

float random_float(float min, float max);
void InitFluid();
void UpdateFluid();
void ResetFluid();
#endif
