#include "solver.h"

float random_float(float min, float max)
{
  return ((float)rand()/(float)RAND_MAX) * (max-min) + min;
}

void InitFluid()
{
  u = (EdgeArray*) malloc(sizeof(EdgeArray));
  v = (EdgeArray*) malloc(sizeof(EdgeArray));
  divergance = (CentreArray*) malloc(sizeof(CentreArray));
  pressure = (CentreArray*) malloc(sizeof(CentreArray));
  particles = (ParticleArray*) malloc(sizeof(ParticleArray));
  cell_type = (IntArray*) malloc(sizeof(IntArray));
  distance = (EdgeArray*) malloc(sizeof(EdgeArray));
  a_diag = (CentreArray*) malloc(sizeof(CentreArray));
  a_plusi = (CentreArray*) malloc(sizeof(CentreArray));
  a_plusj = (CentreArray*) malloc(sizeof(CentreArray));
  aux = (CentreArray*) malloc(sizeof(CentreArray));
  search = (CentreArray*) malloc(sizeof(CentreArray));
  res = (CentreArray*) malloc(sizeof(CentreArray));
  weight = (CentreArray*) malloc(sizeof(CentreArray));
  precon = (CentreArray*) malloc(sizeof(CentreArray));
  m = (CentreArray*) malloc(sizeof(CentreArray));
  ResetFluid();
}

void PrintArray(CentreArray a)
{
  int i,j;
  for (j = SIZE-1; j > -1; --j) {
    for (i = 0; i < SIZE; ++i) {
      printf("%f ", a[i][j]);
    }
    printf("\n");
  }
}

void PrintIntArray(IntArray a)
{
  int i,j;
  for (i = SIZE-1; i > -1; --i) {
    for (j = 0; j < SIZE; ++j) {
      printf("%i ", a[j][i]);
    }
    printf("\n");
  }
}
void FastSweep()
{
  unsigned sweep_num,sweep_dir,i,j;
  //Set distance to be bigger than any line on grid
  for (i = 0; i < SIZE+1 ; i++)
    for (j = 0; j < SIZE+1 ; j++)
      (*distance)[i][j] = SIZE*SIZE;

  int sweep_dirs[4][2] = {{0,-1},{-1,0},{0,1},{1,0}};
  int inner_sweep_dir[2];
  for (sweep_num = 0; sweep_num < 4; sweep_num++) {
    for (sweep_dir = 0; sweep_dir < 4; sweep_dir++) {
      inner_sweep_dir[0] = sweep_dirs[sweep_dir][1] != 0;
      inner_sweep_dir[1] = sweep_dirs[sweep_dir][0] != 0;
      i = sweep_dirs[sweep_dir][0] == -1 ? SIZE-1 : 1;
      j = sweep_dirs[sweep_dir][1] == -1 ? SIZE-1 : 1;
      //Only need to check the upper bound as unsigned will wrap
      for (; i<SIZE+1 && j<SIZE+1; i += sweep_dirs[sweep_dir][0], j+= sweep_dirs[sweep_dir][1]) {
        for (i = inner_sweep_dir[0] ? 1 : i, j = inner_sweep_dir[1] ? 1 : j;
             i<SIZE+1 && j<SIZE+1;
             i += inner_sweep_dir[0], j+= inner_sweep_dir[1]) {
          if ((*cell_type)[i][j] == AIR) {
              float closest, c_u, c_v;
              closest = SIZE*SIZE+1;
              //Check each of our neighbours
              if ((*cell_type)[i-1][j] == FLUID) {
                  (*u)[i][j] = (*u)[i-1][j]; (*v)[i][j] = (*v)[i-1][j];
                  (*distance)[i][j] = 1;
                  continue;
              }
              if ((*cell_type)[i-1][j] == AIR) {
                if ((*distance)[i-1][j] < (*distance)[i][j]) {
                  closest = (*distance)[i-1][j];
                  c_u = (*u)[i-1][j]; c_v = (*v)[i-1][j];
                }
              }
              if ((*cell_type)[i+1][j] == FLUID) {
                  (*u)[i][j] = (*u)[i+1][j]; (*v)[i][j] = (*v)[i+1][j];
                  (*distance)[i][j] = 1;
                  continue;
              }
              if ((*cell_type)[i+1][j] == AIR) {
                if ((*distance)[i+1][j] < (*distance)[i][j]) {
                  closest = (*distance)[i+1][j];
                  c_u = (*u)[i+1][j]; c_v = (*v)[i+1][j];
                }
              }
              if ((*cell_type)[i][j-1] == FLUID) {
                (*u)[i][j] = (*u)[i][j-1]; (*v)[i][j] = (*v)[i][j-1];
                (*distance)[i][j] = 1;
                continue;
              }
              if ((*cell_type)[i][j-1] == AIR) {
                if ((*distance)[i][j-1] < (*distance)[i][j]) {
                  closest = (*distance)[i][j-1];
                  c_u = (*u)[i][j-1]; c_v = (*v)[i][j-1];
                }
              }
              if ((*cell_type)[i][j+1] == FLUID) {
                (*u)[i][j] = (*u)[i][j+1]; (*v)[i][j] = (*v)[i][j+1];
                (*distance)[i][j] = 1;
                continue;
              }
              if ((*cell_type)[i][j+1] == AIR) {
                if ((*distance)[i][j+1] < (*distance)[i][j]) {
                  closest = (*distance)[i][j+1];
                  c_u = (*u)[i][j+1]; c_v = (*v)[i][j+1];
                }
              }
              if (closest < (*distance)[i][j]) {
                (*distance)[i][j] = closest+1;
                (*u)[i][j] = c_u; (*v)[i][j] = c_v;
              }
            }
        }
        i = inner_sweep_dir[0] ? 1 : i, j = inner_sweep_dir[1] ? 1 : j;
      }
    }
  }
}
void ResetFluid()
{
  int i,j,k;
  for (i = 0; i < SIZE+1; ++i)
    for (j = 0; j < SIZE+1; ++j)
      (*u)[i][j] = 0;
  for (i = 0; i < SIZE+1; ++i)
    for (j = 0; j < SIZE+1; ++j)
      (*v)[i][j] = 0;

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
    for (j = 2; j < SIZE-3; j++) {
      for (k = 2; k < SIZE/5; k++) {
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

void bary_x(float x, int* i, float* fx)
{
  float sx=x;
  (*i)=(int)sx;
  (*fx)=sx-floor(sx);
}

void bary_x_centre(float x, int* i, float* fx)
{
  float sx=x-0.5;
  (*i)=(int)sx;
  if((*i)<0){ (*i)=0; (*fx)=0.0; }
  else if((*i)>SIZE-2){ (*i)=SIZE-2; (*fx)=1.0; }
  else{ (*fx)=sx-floor(sx); }
}

void bary_y(float y, int* j, float* fy)
{
  float sy=y;
  (*j)=(int)sy;
  (*fy)=sy-floor(sy);
}

void bary_y_centre(float y, int* j, float* fy)
{
  float sy=y-0.5;
  (*j)=(int)sy;
  if((*j)<0){ (*j)=0; (*fy)=0.0; }
  else if((*j)>SIZE-2){ (*j)=SIZE-2; (*fy)=1.0; }
  else{ (*fy)=sy-floor(sy); }
}

float InterpolateEdge(EdgeArray a, float x, float y, int type)
{
  int i,j;
  float cell_x,cell_y;
  if (type == U) {
    bary_x(x, &i, &cell_x);
    bary_y_centre(y, &j, &cell_y);
  }
  else {
    bary_x_centre(x, &i, &cell_x);
    bary_y(y, &j, &cell_y);
  }
  return (a[i][j]*(1.0-cell_x)*(1.0-cell_y)) +
         (a[i+1][j]*cell_x*(1.0-cell_y)) +
         (a[i][j+1]*(1.0-cell_x)*cell_y) +
         (a[i+1][j+1]*cell_x*cell_y);
}

void MoveParticles()
{
  int i,k;
  //Should be replaced by a proper ODE integrator
  float p_u,p_v;
  for (i = 0; i < num_particles; i++) {
    for (k=0; k<4 ; k++) {
      p_u = InterpolateEdge((*u), (*particles)[i].x, (*particles)[i].y, U);
      p_v = InterpolateEdge((*v), (*particles)[i].x, (*particles)[i].y, V);
      (*particles)[i].x += (1.0f/4.0f)*DT*p_u;
      (*particles)[i].y += (1.0f/4.0f)*DT*p_v;
    }
    (*particles)[i].u = p_u;
    (*particles)[i].v = p_v;
  }
}

void Divergence()
{
  unsigned x,y;
  for (x = 0; x < SIZE+1; ++x) {
    (*u)[0][x] = (*u)[1][x] = (*u)[SIZE][x] = (*u)[SIZE-1][x]=0;
    (*v)[x][0] = (*v)[x][1] = (*v)[x][SIZE] = (*v)[x][SIZE-1]=0;
  }

  for (x = 0; x < SIZE; ++x) {
    for (y = 0; y < SIZE; ++y) {
      if ((*cell_type)[x][y] == FLUID) {
        //printf("1 %f %f  %f %f  %f\n", (*u)[x][y], (*u)[x+1][y], (*v)[x][y], (*u)[x][y+1], (*divergance)[x][y]);
        (*divergance)[x][y] = -(((*u)[x+1][y]- (*u)[x][y]) + ((*v)[x][y+1] - (*v)[x][y]));
        //printf("2 %f %f  %f %f  %f\n", (*u)[x][y], (*u)[x+1][y], (*v)[x][y], (*u)[x][y+1], (*divergance)[x][y]);
      }
      else {
        (*divergance)[x][y] = 0.0;
      }
    }
  }

  //Correct for boundaries
  for (x = 0; x < SIZE; ++x) {
    for (y = 0; y < SIZE; ++y) {
      if ((*cell_type)[x][y] == FLUID) {
          //Where 0 read usolid
        if ((*cell_type)[x-1][y] == SOLID)
          (*divergance)[x][y] -= ((*u)[x][y] - 0);
        if ((*cell_type)[x+1][y] == SOLID)
          (*divergance)[x][y] += ((*u)[x+1][y] - 0);

        if ((*cell_type)[x][y-1] == SOLID)
          (*divergance)[x][y] -= ((*u)[x][y] - 0);
        if ((*cell_type)[x][y+1] == SOLID)
          (*divergance)[x][y] += ((*u)[x][y+1] - 0);
      }
    }
  }

}

void GradientSubtract()
{
  unsigned x,y;
  for (x = 1; x < SIZE; ++x) {
    for (y = 1; y < SIZE; ++y) {
      if (((*cell_type)[x][y] | (*cell_type)[x-1][y]) == FLUID) {
        (*u)[x][y] -= ((*pressure)[x][y] - (*pressure)[x-1][y]);
      }
      if (((*cell_type)[x][y] | (*cell_type)[x][y-1]) == FLUID) {
        (*v)[x][y] -= ((*pressure)[x][y] - (*pressure)[x][y-1]);
      }
    }
  }
}

void FormPreconditioner()
{
   const double mic_parameter=0.97;
   double d;
   unsigned i,j;
   for (j = 0; j < SIZE; ++j)
    for (i = 0; i < SIZE; ++i)
      (*precon)[i][j] = 0;
   for (j = 0; j < SIZE; ++j)
     for (i = 0; i < SIZE; ++i)
       if((*cell_type)[i][j] == FLUID) {
         d = (*a_diag)[i][j] - sqr((*a_plusi)[i-1][j] * (*precon)[i-1][j])
                             - sqr((*a_plusj)[i][j-1] * (*precon)[i][j-1])
                             - mic_parameter * ( (*a_plusi)[i-1][j]*(*a_plusj)[i-1][j]*sqr((*precon)[i-1][j])
                                               + (*a_plusj)[i][j-1]*(*a_plusi)[i][j-1]*sqr((*precon)[i][j-1]) );
         (*precon)[i][j]=1/sqrt(d+1e-6);
      }
}

void ApplyPreconditioner()
{
   int i, j;
   float d;
   for (j = 0; j < SIZE; ++j)
    for (i = 0; i < SIZE; ++i)
      (*m)[i][j] = 0;

   // Solve L*m=x
   for(j=0; j<SIZE; ++j) 
     for(i=0; i<SIZE; ++i)
       if((*cell_type)[i][j] == FLUID){
         d=(*res)[i][j] - (*a_plusi)[i-1][j] * (*precon)[i-1][j] * (*m)[i-1][j]
                        - (*a_plusj)[i][j-1] * (*precon)[i][j-1] * (*m)[i][j-1];
         (*m)[i][j] = (*precon)[i][j] * d;
       }
   // Solve L'*y=m
   for (j = 0; j < SIZE; ++j)
    for (i = 0; i < SIZE; ++i)
      (*aux)[i][j] = 0;

   for(j = SIZE-1; j > -1; --j) 
     for(i = SIZE-1; i > -1; --i)
       if((*cell_type)[i][j] == FLUID){
         d = (*m)[i][j] - (*a_plusi)[i][j] * (*precon)[i][j] * (*aux)[i+1][j]
                        - (*a_plusj)[i][j] * (*precon)[i][j] * (*aux)[i][j+1];
         (*aux)[i][j] = (*precon)[i][j] * d;
       }
}

void Solve()
{
  int i,j;
  float scale = 1.0f;
  //Form matrix A
  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) {
      (*a_diag)[i][j] = 0.0f;
      (*a_plusi)[i][j] = 0.0f;
      (*a_plusj)[i][j] = 0.0f;
      (*pressure)[i][j] = 0.0f; //Initial guess
      (*res)[i][j] = 0.0f;
    }
  }

  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) {
      if ((*cell_type)[i][j] == FLUID && (*cell_type)[i+1][j] == FLUID) {
        (*a_diag)[i][j] +=  scale;
        (*a_diag)[i+1][j] += scale;
        (*a_plusi)[i][j] = -scale;
      }
      if ((*cell_type)[i][j] == FLUID && (*cell_type)[i+1][j] == AIR) {
          (*a_diag)[i][j] += scale;
        }

      if ((*cell_type)[i][j] == FLUID && (*cell_type)[i][j+1] == FLUID) {
        (*a_diag)[i][j] +=  scale;
        (*a_diag)[i][j+1] += scale;
        (*a_plusj)[i][j] = -scale;
      }
      if ((*cell_type)[i][j] == FLUID && (*cell_type)[i][j+1] == AIR) {
        (*a_diag)[i][j] += scale;
      }
    }
  }

  FormPreconditioner();
  //r = b
  for (i = 0; i < SIZE; ++i)
    for (j = 0; j < SIZE; ++j)
      if ((*cell_type)[i][j] == FLUID)
        (*res)[i][j] = (*divergance)[i][j];
  
  //Check r = 0
  float max_res = 0;
  for (i = 0; i < SIZE; ++i)
    for (j = 0; j < SIZE; ++j)
      if (fabs((*res)[i][j]) > max_res)
        if ((*cell_type)[i][j] == FLUID)
          max_res = fabs((*res)[i][j]);
  if (max_res == 0) return;

  ApplyPreconditioner();

  for (i = 0; i < SIZE; ++i)
    for (j = 0; j < SIZE; ++j)
      if ((*cell_type)[i][j] == FLUID)
        (*search)[i][j] = (*aux)[i][j];

  double sigma = 0;
  double sigma_new = 0;
  for (i = 0; i < SIZE; ++i)
    for (j = 0; j < SIZE; ++j)
      if ((*cell_type)[i][j] == FLUID)
        sigma += (*aux)[i][j] * (*res)[i][j];

  int iterations = 0;
  while (iterations < 100) {
    ++iterations;
    for (i = 0; i < SIZE; ++i)
      for (j = 0; j < SIZE; ++j)
        if ((*cell_type)[i][j] == FLUID)
          (*aux)[i][j] = (*a_diag)[i][j] * (*search)[i][j]
            + (*a_plusi)[i-1][j] * (*search)[i-1][j]
            + (*a_plusi)[i][j] * (*search)[i+1][j]
            + (*a_plusj)[i][j-1] * (*search)[i][j-1]
            + (*a_plusj)[i][j] * (*search)[i][j+1];

    double temp = 0;
    for (i = 0; i < SIZE; ++i)
      for (j = 0; j < SIZE; ++j)
        if ((*cell_type)[i][j] == FLUID)
          temp += (*aux)[i][j] * (*search)[i][j];

    double alpha = sigma/temp;
    for (i = 0; i < SIZE; ++i)
      for (j = 0; j < SIZE; ++j)
        if ((*cell_type)[i][j] == FLUID)
          (*pressure)[i][j] += alpha * (*search)[i][j];
    for (i = 0; i < SIZE; ++i)
      for (j = 0; j < SIZE; ++j)
        if ((*cell_type)[i][j] == FLUID)
          (*res)[i][j] -= alpha * (*aux)[i][j];

    float max_res = 0;
    for (i = 0; i < SIZE; ++i)
      for (j = 0; j < SIZE; ++j)
        if ((*cell_type)[i][j] == FLUID)
          if (fabs((*res)[i][j]) > max_res)
            max_res = fabs((*res)[i][j]);

/*    printf("type\n");
    PrintIntArray(*cell_type);
    printf("divergance\n");
    PrintArray(*divergance);
    printf("pressure\n");
    PrintArray(*pressure);
    printf("a_diag\n");
    PrintArray(*a_diag);
    printf("a_plusi\n");
    PrintArray(*a_plusi);
    printf("a_plusj\n");
    PrintArray(*a_plusj);
    printf("aux\n");
    PrintArray(*aux);
    printf("res\n");
    PrintArray(*res);
    printf("search\n");
    PrintArray(*search);*/
    //printf("%f\n", max_res);
    //printf("%i\n", iterations);

    if (max_res < 0.000005) {
      printf("%i iterations %f\n", iterations, max_res);
      return;
    }

    ApplyPreconditioner();

    sigma_new = 0;
    for (i = 0; i < SIZE; ++i)
      for (j = 0; j < SIZE; ++j)
        if ((*cell_type)[i][j] == FLUID)
          sigma_new += (*aux)[i][j] * (*res)[i][j];

    double beta = sigma_new/sigma;

    for (i = 0; i < SIZE; ++i)
      for (j = 0; j < SIZE; ++j)
        if ((*cell_type)[i][j] == FLUID)
          (*search)[i][j] = (*aux)[i][j] + beta*(*search)[i][j];

    sigma = sigma_new;
  }
}

void Project()
{
  Divergence();
  Solve();
  GradientSubtract();
}

void TransferParticlesToGrid()
{
  int i,j;
  for (i = 0; i < SIZE; ++i)
    for (j = 0; j < SIZE; ++j)
      (*cell_type)[i][j] = AIR;
  for (i = 0; i < SIZE; ++i)
    for (j = 0; j < SIZE; ++j)
      (*weight)[i][j] = 0.0f;
  for (i = 0; i < SIZE+1; ++i)
    for (j = 0; j < SIZE+1; ++j)
      (*u)[i][j] = 0.0f;
  for (i = 0; i < SIZE+1; ++i)
    for (j = 0; j < SIZE+1; ++j)
      (*v)[i][j] = 0.0f;


  int x,y;
  float cell_x, cell_y, w;
  for (i = 0; i < num_particles; ++i) {
    bary_x((*particles)[i].x, &x, &cell_x);
    bary_y_centre((*particles)[i].y, &y, &cell_y);

    w=(1-cell_x)*(1-cell_y);
    (*u)[x][y] += w * (*particles)[i].u;
    (*weight)[x][y] += w;

    w=cell_x*(1-cell_y);
    (*u)[x+1][y] += w * (*particles)[i].u;
    (*weight)[x+1][y] += w;

    w=(1-cell_x)*cell_y;
    (*u)[x][y+1] += w * (*particles)[i].u;
    (*weight)[x][y+1] += w;

    w=cell_x*cell_y;
    (*u)[x+1][y+1] += w * (*particles)[i].u;
    (*weight)[x+1][y+1] += w;
  }
  for (j = 0; j < SIZE; ++j) 
    for (i = 0; i < SIZE; ++i) 
      if ((*weight)[i][j] > 0) {
        (*u)[i][j] /= (*weight)[i][j];
        (*weight)[i][j] = 0.0f;
      }

  for (i = 0; i < num_particles; ++i) {
    bary_x_centre((*particles)[i].x, &x, &cell_x);
    bary_y((*particles)[i].y, &y, &cell_y);

    w=(1-cell_x)*(1-cell_y);
    (*v)[x][y] += w * (*particles)[i].v;
    (*weight)[x][y] += w;

    w=cell_x*(1-cell_y);
    (*v)[x+1][y] += w * (*particles)[i].v;
    (*weight)[x+1][y] += w;

    w=(1-cell_x)*cell_y;
    (*v)[x][y+1] += w * (*particles)[i].v;
    (*weight)[x][y+1] += w;

    w=cell_x*cell_y;
    (*v)[x+1][y+1] += w * (*particles)[i].v;
    (*weight)[x+1][y+1] += w;
  }
  for (j = 0; j < SIZE; ++j) 
    for (i = 0; i < SIZE; ++i) 
      if ((*weight)[i][j] > 0) 
        (*v)[i][j] /= (*weight)[i][j];

  for (i = 0; i < num_particles; ++i) {
    bary_x((*particles)[i].x, &x, &cell_x);
    bary_y((*particles)[i].y, &y, &cell_y);
    (*cell_type)[x][y] = FLUID;
  }

  for (i = 0; i < SIZE; ++i) {
    (*cell_type)[0][i] = SOLID;
    (*cell_type)[SIZE-1][i] = SOLID;
  }
  for (j = 0; j < SIZE; ++j) {
    (*cell_type)[j][0] = SOLID;
    (*cell_type)[j][SIZE-1] = SOLID;
  }
}

void GravityToGrid()
{
  int i,j;
  for (i = 0; i < SIZE+1; ++i)
    for (j = 0; j < SIZE+1; ++j)
      (*v)[i][j] -= 0.1;
}

void UpdateFluid(int vel_on)
{
  TransferParticlesToGrid();
  GravityToGrid();
  FastSweep();
  Project();
  FastSweep();
  MoveParticles();
}
