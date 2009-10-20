#ifndef SOLVER_H
#define SOLVER_H

#include <stdlib.h>
#include <stdio.h>

#define NUM_PARTICLES 500000
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
void ImputeCorners(ValueArray* v);
void ApplyNeumannBoundary(ValueArray* v);
void ApplyDirichletBoundary(ValueArray* v, unsigned dir);
float InterpolateScalar(ValueArray array, float x, float y);
void JacobiPoissonSolver(int iterations, float a, float c, ValueArray* v, ValueArray* b);
void MoveParticles(ValueArray* u, ValueArray* v);
void Advect(ValueArray* dest, ValueArray* source, ValueArray* u, ValueArray* v);
void Diffuse(ValueArray* dest, ValueArray* source);
void Divergence(ValueArray* dest, ValueArray* u, ValueArray* v);
void GradientSubtract(ValueArray* u, ValueArray* v, ValueArray* pressure);
void Project(ValueArray* u, ValueArray* v, ValueArray* pressure, ValueArray* div);
void Squirt();
void UpdateFluid(int squirt, int reverse);
#endif
