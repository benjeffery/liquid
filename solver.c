#include "solver.h"

float random_float(float min, float max)
{
  return ((float)rand()/(float)RAND_MAX) * (max-min) + min;
}

void InitFluid()
{
  gu = (EdgeArray*) malloc(sizeof(EdgeArray));
  gv = (EdgeArray*) malloc(sizeof(EdgeArray));
  divergance = (CentreArray*) malloc(sizeof(CentreArray));
  pressure = (CentreArray*) malloc(sizeof(CentreArray));
  particles = (ParticleArray*) malloc(sizeof(ParticleArray));
  has_fluid = (IntArray*) malloc(sizeof(IntArray));
  ResetFluid();
}

void ResetFluid()
{
  int i,j,k;
  for (i = 0; i < SIZE+1; ++i) 
    for (j = 0; j < SIZE+1; ++j) 
      (*gu)[i][j] = 0;
  for (i = 0; i < SIZE+1; ++i) 
    for (j = 0; j < SIZE+1; ++j) 
      (*gv)[i][j] = 0;

  int c = 0;
  for (i = 0 ; i < NUM_PARTICLES_PER_CELL ; i++) {
    //Box to drop
    for (j = SIZE/2; j < SIZE/2+SIZE/6; j++) {
      for (k = 3*SIZE/4; k < 3*SIZE/4+SIZE/6; k++) {
        (*particles)[c].x = random_float((float)j-0.5f,(float)j+0.5f);
        (*particles)[c].y = random_float((float)k-0.5f,(float)k+0.5f);
        (*particles)[c].u = 0.0f;
        (*particles)[c].v = 0.0f;
        ++c;
      }
    }
    //Pool below
    for (j = 2; j < SIZE-2; j++) {
      for (k = 2; k < SIZE/6; k++) {
        (*particles)[c].x = random_float((float)j-0.5f,(float)j+0.5f);
        (*particles)[c].y = random_float((float)k-0.5f,(float)k+0.5f);
        (*particles)[c].u = 0.0f;
        (*particles)[c].v = 0.0f;
        ++c;
      }
    }
  }
  printf("%i particles created\n", c);
  num_particles = c;
}


//void ImputeCorners(ValueArray* v)
//{
//  (*v)[0][0] = 0.5f*((*v)[1][0]+(*v)[0][1]);
//  (*v)[0][SIZE-1] = 0.5f*((*v)[1][SIZE-1]+(*v)[0][SIZE-2]);
//  (*v)[SIZE-1][0] = 0.5f*((*v)[SIZE-2][0]+(*v)[SIZE-1][1]);
//  (*v)[SIZE-1][SIZE-1] = 0.5f*((*v)[SIZE-2][SIZE-1]+(*v)[SIZE-1][SIZE-2]);
//}


/* void ApplyNeumannBoundary(ValueArray* v) */
/* { */
/*   unsigned x,y; */
/*   for (x = 1; x < SIZE-1; ++x) { */
/*     (*v)[x][0] = (*v)[x][1]; */
/*     (*v)[x][SIZE-1] = (*v)[x][SIZE-2]; */
/*   } */
/*   for (y = 1; y < SIZE-1; ++y) { */
/*     (*v)[0][y] = (*v)[1][y]; */
/*     (*v)[SIZE-1][y]= (*v)[SIZE-2][y]; */
/*   } */
/*   ImputeCorners(v); */
/* } */


float InterpolateEdge(EdgeArray a, float x, float y, int type)
{
  if (x < 0.5f) x = 0.5f;
  if (y < 0.5f) y = 0.5f;
  if (x > SIZE-0.5f) x = SIZE-0.5f;
  if (y > SIZE-0.5f) y = SIZE-0.5f;
  if (type == U)
    y -= 0.5f;
  if (type == V)
    x -= 0.5f;
  unsigned idx_x = x;
  unsigned idx_y = y;
  float cell_x = x-idx_x;
  float cell_y = y-idx_y;
  return (a[idx_x][idx_y]*(1.0-cell_x)*(1.0-cell_y)) + 
         (a[idx_x+1][idx_y]*cell_x*(1.0-cell_y)) + 
         (a[idx_x][idx_y+1]*(1.0-cell_x)*cell_y) + 
         (a[idx_x+1][idx_y+1]*cell_x*cell_y);
}

void JacobiPoissonSolver(int iterations, float a, float c, CentreArray v, CentreArray b)
{
  //Should be replaced by a much better solver!
  /* Solve Av = b Where A = del^2 */
  unsigned iteration,x,y;
  for (iteration = 0; iteration < iterations; ++iteration) {
    for (x = 1; x < SIZE-1; ++x) {
      for (y = 1; y < SIZE-1; ++y) {
        // if ((*has_fluid)[x][y] == FLUID)
          v[x][y] = ((a*(v[x-1][y] + v[x+1][y]+ v[x][y-1]+ v[x][y+1])) + b[x][y]) / c;
      }
    }
  }
}

void MoveParticles(EdgeArray u, EdgeArray v)
{
  int i,k;
  //Should be replaced by a proper ODE integrator
  float p_u,p_v;
  for (i = 0; i < num_particles; i++) {
    for (k=0; k<4 ; k++) {
      p_u = InterpolateEdge(u, (*particles)[i].x, (*particles)[i].y, U);
      p_v = InterpolateEdge(v, (*particles)[i].x, (*particles)[i].y, V);
      (*particles)[i].x += (1.0f/4.0f)*DT*SIZE*p_u;
      (*particles)[i].y += (1.0f/4.0f)*DT*SIZE*p_v;
    }
    (*particles)[i].u = p_u;
    (*particles)[i].v = p_v;
  }
}

void Divergence(CentreArray dest, EdgeArray u, EdgeArray v)
{
  unsigned x,y;
  for (x = 0; x < SIZE; ++x) 
    for (y = 0; y < SIZE; ++y) 
      if (TRUE) //((*has_fluid)[x][y] == FLUID)
        dest[x][y] = ((u[x][y] - u[x+1][y]) + (v[x][y] - v[x][y+1])) / SIZE;
      else
        dest[x][y] = 0;
}

void GradientSubtract(EdgeArray u, EdgeArray v, CentreArray pressure)
{
  unsigned x,y;
  for (x = 0; x < SIZE-1; ++x) {
    for (y = 0; y < SIZE-1; ++y) {
      u[x][y] -= DT * SIZE * (pressure[x+1][y] - pressure[x][y]);
      v[x][y] -= DT * SIZE * (pressure[x][y+1] - pressure[x][y]);
    }
  }
}

void ApplyBoundary(EdgeArray u, EdgeArray v)
{
  unsigned x,y;
  for (x = 1; x < SIZE; ++x) {
    v[x][0] = -v[x][1];
    v[x][SIZE] = -v[x][SIZE-1];
  }
  for (y = 1; y < SIZE; ++y) {
    u[0][y] = -u[1][y];
    u[SIZE][y] = -u[SIZE-1][y];
  }
}

void Project(EdgeArray u, EdgeArray v, CentreArray pressure, CentreArray div)
{
  ApplyBoundary(u, v);
  Divergence(div, u, v);
  unsigned x,y;
  for (x = 0; x < SIZE; ++x)
    for (y = 0; y < SIZE; ++y)
      pressure[x][y] = 0;
  JacobiPoissonSolver(100, 1.0f, 4.0f, pressure, div);
  GradientSubtract(u, v, pressure);
}

void TransferParticlesToGrid(EdgeArray u, EdgeArray v, CentreArray sum_of_weights)
{
  int i,j;
  for (i = 0; i < SIZE; ++i) 
    for (j = 0; j < SIZE; ++j) 
      (*has_fluid)[i][j] = FALSE;
  for (i = 0; i < SIZE; ++i) 
    for (j = 0; j < SIZE; ++j) 
      sum_of_weights[i][j] = 0.0f;
  for (i = 0; i < SIZE+1; ++i) 
    for (j = 0; j < SIZE+1; ++j) 
      u[i][j] = 0.0f;
  for (i = 0; i < SIZE+1; ++i) 
    for (j = 0; j < SIZE+1; ++j) 
      v[i][j] = 0.0f;

  //TODO Particles should effect neighbouring cells and have proper weighting
  for (i = 0; i < num_particles; ++i) {
    unsigned x = (*particles)[i].x;
    unsigned y = (*particles)[i].y;
    //Mark the cell containing the box as fluid
    (*has_fluid)[x][y] = FLUID;
    u[x][y] += (*particles)[i].u;
    v[x][y] += (*particles)[i].v;
    sum_of_weights[x][y] += 1.0f;
  }
  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) { 
      if (sum_of_weights[i][j] > 0) {
        u[i][j] /= sum_of_weights[i][j]; 
        v[i][j] /= sum_of_weights[i][j]; 
      }
    }
  }    
}

void GravityToGrid(EdgeArray v)
{
  int i,j;
  for (i = 30; i < 32; ++i) 
    for (j = 1; j < 4; ++j) 
      //      if ((*has_fluid)[i][j] == FLUID)
        v[i][j] -= 0.1*DT;
}

void UpdateFluid()
{
  TransferParticlesToGrid(*gu, *gv, *divergance);
  GravityToGrid(*gv);
  Project(*gu, *gv, *divergance, *pressure);
  MoveParticles(*gu,*gv);
}
