#include "solver.h"

float random_float(float min, float max)
{
  return ((float)rand()/(float)RAND_MAX) * (max-min) + min;
}

void InitFluid()
{
  gu = (ValueArray*) malloc(sizeof(ValueArray));
  gv = (ValueArray*) malloc(sizeof(ValueArray));
  gu_old = (ValueArray*) malloc(sizeof(ValueArray));
  gv_old = (ValueArray*) malloc(sizeof(ValueArray));
  particles = (ParticleArray*) malloc(sizeof(ParticleArray));

  int i,j,l;
  l = 0;
  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) {
      (*gu)[i][j] = 0;
      (*gv)[i][j] = 0;
      (*gu_old)[i][j] = 0;
      (*gv_old)[i][j] = 0;
    }
  }
  
  for (i = 0 ; i < NUM_PARTICLES ; i++) {
    (*particles)[i].x = random_float(1,SIZE-1);
    (*particles)[i].y = random_float(1,SIZE-1);
  }
}

void ImputeCorners(ValueArray* v)
{
  (*v)[0][0] = 0.5f*((*v)[1][0]+(*v)[0][1]);
  (*v)[0][SIZE-1] = 0.5f*((*v)[1][SIZE-1]+(*v)[0][SIZE-2]);
  (*v)[SIZE-1][0] = 0.5f*((*v)[SIZE-2][0]+(*v)[SIZE-1][1]);
  (*v)[SIZE-1][SIZE-1] = 0.5f*((*v)[SIZE-2][SIZE-1]+(*v)[SIZE-1][SIZE-2]);
}

void ApplyNeumannBoundary(ValueArray* v)
{
  unsigned x,y;
  for (x = 1; x < SIZE-1; ++x) {
    (*v)[x][0] = (*v)[x][1];
    (*v)[x][SIZE-1] = (*v)[x][SIZE-2];
  }
  for (y = 1; y < SIZE-1; ++y) {
    (*v)[0][y] = (*v)[1][y];
    (*v)[SIZE-1][y]= (*v)[SIZE-2][y];
  }
  ImputeCorners(v);
}

void ApplyDirichletBoundary(ValueArray* v, unsigned dir)
{
  unsigned x,y;
  for (x = 1; x < SIZE-1; ++x) {
    (*v)[x][0] = dir==TOPBOTTOM ? -(*v)[x][1] :  (*v)[x][1];
    (*v)[x][SIZE-1] = dir==TOPBOTTOM ? -(*v)[x][SIZE-2] : (*v)[x][SIZE-2];
  }
  for (y = 1; y < SIZE-1; ++y) {
    (*v)[0][y] = dir == LEFTRIGHT ? -(*v)[1][y] : (*v)[1][y];
    (*v)[SIZE-1][y] = dir == LEFTRIGHT ? -(*v)[SIZE-2][y] :  (*v)[SIZE-2][y];
  }
  ImputeCorners(v);
}

float InterpolateScalar(ValueArray array, float x, float y)
{
  if (x < 0.0f) x = 0.5f;
  if (y < 0.0f) y = 0.5f;
  if (x > SIZE-1.0f) x = SIZE-1.5f;
  if (y > SIZE-1.0f) y = SIZE-1.5f;
  unsigned f_x = x;
  unsigned f_y = y;
  float p_x = x-f_x;
  float p_y = y-f_y;
  return (array[f_x][f_y]*(1.0-p_x)*(1.0-p_y)) + (array[f_x+1][f_y]*p_x*(1.0-p_y)) + (array[f_x][f_y+1]*(1.0-p_x)*p_y) + (array[f_x+1][f_y+1]*p_x*p_y);
}

void JacobiPoissonSolver(int iterations, float a, float c, ValueArray* v, ValueArray* b)
{
  /* Solve Av = b Where A = del^2 */
  unsigned iteration,x,y;
  for (iteration = 0; iteration < iterations; ++iteration) {
    ApplyNeumannBoundary(v);
    for (x = 1; x < SIZE-1; ++x) {
      for (y = 1; y < SIZE-1; ++y) {
        (*v)[x][y] = ((a*((*v)[x-1][y] + (*v)[x+1][y]+ (*v)[x][y-1]+ (*v)[x][y+1])) + (*b)[x][y]) / c;
      }
    }
  }
}

void MoveParticles(ValueArray* u, ValueArray* v)
{
  int i;
  for (i = 0; i < NUM_PARTICLES; i++) {
    (*particles)[i].x += DT*SIZE*InterpolateScalar(*u, (*particles)[i].x, (*particles)[i].y);
    (*particles)[i].y += DT*SIZE*InterpolateScalar(*v, (*particles)[i].x, (*particles)[i].y);
  }
}

void Advect(ValueArray* dest, ValueArray* source, ValueArray* u, ValueArray* v)
{
  ApplyNeumannBoundary(source);
  unsigned x,y;
  for (x = 1; x < SIZE-1; ++x) {
    for (y = 1; y < SIZE-1; ++y) {
      float old_x = (float)x - (DT * SIZE * (*u)[x][y]);
      float old_y = (float)y - (DT * SIZE * (*v)[x][y]);
      (*dest)[x][y] = InterpolateScalar(*source, old_x, old_y);
    }
  }
}

void Diffuse(ValueArray* dest, ValueArray* source)
{
  float a=DT*DIFFUSION*SIZE*SIZE;
  JacobiPoissonSolver(40, a, 1+(4.0f*a), dest, source);
}

void Divergence(ValueArray* dest, ValueArray* u, ValueArray* v)
{
  unsigned x,y;
  for (x = 1; x < SIZE-1; ++x) {
    for (y = 1; y < SIZE-1; ++y) {
      (*dest)[x][y] = 0.5f * (((*u)[x-1][y] - (*u)[x+1][y]) + ((*v)[x][y-1] - (*v)[x][y+1])) / SIZE;
    }
  }
}

void GradientSubtract(ValueArray* u, ValueArray* v, ValueArray* pressure)
{
  unsigned x,y;
  for (x = 1; x < SIZE-1; ++x) {
    for (y = 1; y < SIZE-1; ++y) {
      (*u)[x][y] -= 0.5f * SIZE * ((*pressure)[x+1][y] - (*pressure)[x-1][y]);
      (*v)[x][y] -= 0.5f * SIZE * ((*pressure)[x][y+1] - (*pressure)[x][y-1]);
    }
  }
}

void Project(ValueArray* u, ValueArray* v, ValueArray* pressure, ValueArray* div)
{
    ApplyDirichletBoundary(u, LEFTRIGHT);
    ApplyDirichletBoundary(v, TOPBOTTOM);
    Divergence(div, u, v);
    ApplyNeumannBoundary(div);
    unsigned x,y;
    for (x = 0; x < SIZE; ++x) {
        for (y = 0; y < SIZE; ++y) {
            (*pressure)[x][y] = 0;
        }
    }
    JacobiPoissonSolver(100, 1.0f, 4.0f, pressure, div);
    ApplyNeumannBoundary(pressure);
    GradientSubtract(u, v, pressure);
}

void Squirt(int reverse)
{
  int k;
  for (k=-7 ;k < 8 ;k++) {
      (*gu)[SIZE/4][SIZE/2+k] = reverse ? .2f:-.2f;
      (*gu)[(SIZE/4)*3][SIZE/2+k] = reverse ? -.2f:.2f;
      (*gv)[SIZE/2+k][SIZE/4] = reverse ? .2f:-.2f;
      (*gv)[SIZE/2+k][(SIZE/4)*3] = reverse ? -.2f:.2f;
  }
}

void UpdateFluid(int squirt, int reverse)
{
  //VELOCITY COMPUTATIONS
  swap(gu, gu_old);
  swap(gv, gv_old);
  Advect(gu, gu_old, gu_old, gv_old);
  Advect(gv, gv_old, gu_old, gv_old);

  if (squirt)
    Squirt(reverse);

  //swap(gu, gu_old);
  //swap(gv, gv_old);
  //Diffuse(gu, gu_old);
  //Diffuse(gv, gv_old);

  Project(gu, gv, gu_old, gv_old);
  MoveParticles(gu,gv);
}
