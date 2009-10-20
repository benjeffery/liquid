#ifndef SOLVER_H
#define SOLVER_H

#include <stdlib.h>
#include <stdio.h>

#define NUM_PARTICLES 200000
#define SIZE 250
#define DT 0.1f
#define TOPBOTTOM 1
#define LEFTRIGHT 2
#define DIFFUSION 0.0001

#define swap(x,y) {ValueArray* temp = x; x=y; y=temp;}

typedef float Value;
typedef Value ValueArray[SIZE][SIZE];

ValueArray* gu;
ValueArray* gv;
ValueArray* gu_old;
ValueArray* gv_old;

typedef struct {
  float x;
  float y;
} Particle;

typedef Particle ParticleArray[NUM_PARTICLES];
ParticleArray* particles;

float random_float(float min, float max);
void InitFluid();
void UpdateFluid(int squirt, int reverse);
#endif
