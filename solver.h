#ifndef SOLVER_H
#define SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define NUM_PARTICLES_PER_CELL 4
#define SIZE 100
#define DT 0.1f

#define U 1
#define V 2

#define SOLID 1
#define FLUID 2
#define AIR 3

#define BOTH 3
#define DIFFUSION 0.0001

#define swap(x,y) {ValueArray* temp = x; x=y; y=temp;}
#define sqr(x) ((x)*(x))

typedef float Value;
typedef Value EdgeArray[SIZE+1][SIZE+1];
typedef Value CentreArray[SIZE][SIZE];
typedef int IntArray[SIZE][SIZE];

EdgeArray* u;
EdgeArray* v;
CentreArray* divergance;
CentreArray* pressure;
IntArray*   cell_type;
EdgeArray* distance;
CentreArray* a_diag;
CentreArray* a_plusi;
CentreArray* a_plusj;
CentreArray* aux;
CentreArray* search;
CentreArray* res;
CentreArray* weight;
CentreArray* precon;
CentreArray* m;


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
void UpdateFluid(int vel_on);
void ResetFluid();
#endif
