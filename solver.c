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
  has_fluid = (IntArray*) malloc(sizeof(IntArray));
  ResetFluid();
}

void ResetFluid()
{
  int i,j,k;
  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) {
      (*gu)[i][j] = 0;
      (*gv)[i][j] = 0;
      (*gu_old)[i][j] = 0;
      (*gv_old)[i][j] = 0;
    }
  }
  
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
    (*v)[x][0] = (dir==TOPBOTTOM || dir == BOTH) ? -(*v)[x][1] :  0;
    (*v)[x][SIZE-1] = (dir==TOPBOTTOM || dir == BOTH) ? -(*v)[x][SIZE-2] : 0;
  }
  for (y = 1; y < SIZE-1; ++y) {
    (*v)[0][y] = (dir == LEFTRIGHT || dir == BOTH) ? -(*v)[1][y] : 0;
    (*v)[SIZE-1][y] = (dir == LEFTRIGHT || dir == BOTH) ? -(*v)[SIZE-2][y] : 0;
  }
  ImputeCorners(v);
}

float InterpolateScalar(ValueArray array, float x, float y)
{
  if (x < 0.5f) x = 0.5f;
  if (y < 0.5f) y = 0.5f;
  if (x > SIZE-1.5f) x = SIZE-1.5f;
  if (y > SIZE-1.5f) y = SIZE-1.5f;
  unsigned f_x = x;
  unsigned f_y = y;
  float p_x = x-f_x;
  float p_y = y-f_y;
  return (array[f_x][f_y]*(1.0-p_x)*(1.0-p_y)) + (array[f_x+1][f_y]*p_x*(1.0-p_y)) + (array[f_x][f_y+1]*(1.0-p_x)*p_y) + (array[f_x+1][f_y+1]*p_x*p_y);
}

void JacobiPoissonSolver(int iterations, float a, float c, ValueArray* v, ValueArray* b)
{
  //Should be replaced by a much better solver!
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
  int i,k;
  //Should be replaced by a proper ODE integrator
  float p_u,p_v;
  for (i = 0; i < num_particles; i++) {
    for (k=0; k<4 ; k++) {
      p_u = InterpolateScalar(*u, (*particles)[i].x, (*particles)[i].y);
      p_v = InterpolateScalar(*v, (*particles)[i].x, (*particles)[i].y);
      (*particles)[i].x += (1.0f/4.0f)*DT*SIZE*p_u;
      (*particles)[i].y += (1.0f/4.0f)*DT*SIZE*p_v;
    }
    (*particles)[i].u = p_u;
    (*particles)[i].v = p_v;
  }
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
    ApplyDirichletBoundary(div, BOTH);
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

void TransferParticlesToGrid(ValueArray* u, ValueArray* v, ValueArray* sum_of_weights)
{
  int i,j;
  for (i = 0; i < SIZE; ++i) 
    for (j = 0; j < SIZE; ++j) 
      (*has_fluid)[i][j] = FALSE;
  for (i = 0; i < SIZE; ++i) 
    for (j = 0; j < SIZE; ++j) 
      (*sum_of_weights)[i][j] = 0.0f;
  for (i = 0; i < SIZE; ++i) 
    for (j = 0; j < SIZE; ++j) 
      (*u)[i][j] = 0.0f;
  for (i = 0; i < SIZE; ++i) 
    for (j = 0; j < SIZE; ++j) 
      (*v)[i][j] = 0.0f;

  //TODO Particles should effect neighbouring cells and have proper weighting
  for (i = 0; i < num_particles; ++i) {
    unsigned x = (*particles)[i].x;
    unsigned y = (*particles)[i].y;
    //Mark the cell containing the box as fluid
    (*has_fluid)[x][y] = TRUE;
    (*u)[x][y] += (*particles)[i].u;
    (*v)[x][y] += (*particles)[i].v;
    (*sum_of_weights)[x][y] += 1.0f;
  }
  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) { 
      if ((*sum_of_weights)[i][j] > 0) {
        //printf("%f %f   ", (*u)[i][j], (*v)[i][j]);
        (*u)[i][j] /= (*sum_of_weights)[i][j]; 
        (*v)[i][j] /= (*sum_of_weights)[i][j]; 
        // printf("%f %f   ", (*u)[i][j], (*v)[i][j]);
      }
    }
  }    
}

void GravityToGrid(ValueArray* v)
{
  int i,j;
  for (i = 1; i < SIZE-1; ++i) 
    for (j = 1; j < SIZE-1; ++j) 
      (*v)[i][j] -= 0.1*DT;
}

void UpdateFluid()
{
  TransferParticlesToGrid(gu, gv, gu_old);
  GravityToGrid(gv);
  Project(gu, gv, gu_old, gv_old);
  MoveParticles(gu,gv);
}
