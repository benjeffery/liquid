
#include <stdio.h>
#include <stdlib.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <unistd.h>
#include "SDL.h"

#define swap(x,y) {ValueArray* temp = x; x=y; y=temp;}

/* screen width, height, and bit depth */
#define SCREEN_WIDTH  500
#define SCREEN_HEIGHT 500
#define SCREEN_BPP     16
#define SIZE 250

#define DT 0.1f
#define DIFFUSION 0.0001f
#define VISCOSITY 0.0000001f 

/* Set up some booleans */
#define TRUE  1
#define FALSE 0

#define NEITHER 0
#define TOPBOTTOM 1
#define LEFTRIGHT 2

/* This is our SDL surface */
SDL_Surface *surface;

GLfloat xpos; 
GLfloat ypos; 
GLfloat zoom; 

typedef float Value;
typedef Value ValueArray[SIZE][SIZE];
typedef Value ValueCol[SIZE];
ValueArray gu_;
ValueArray* gu;
ValueArray gv_;
ValueArray* gv;
ValueArray gu_old_;
ValueArray* gu_old;
ValueArray gv_old_;
ValueArray* gv_old;

int vel_on;
int show_vel;
int mat_on;
int reverse;
int mat_dir;
int run;

typedef struct {
  float x;
  float y;
} Particle;

#define NUM_PARTICLES 500000

Particle particles[NUM_PARTICLES];

float random_float(float limit)
{
  return ((float)rand()/(float)RAND_MAX) * limit;
}

void InitFluid()
{
  
  xpos = 0;
  ypos = 0;
  zoom = -3.5;
  vel_on=TRUE;
  mat_on=FALSE;
  show_vel=FALSE;
  reverse=FALSE;
  mat_dir = 1;
  run = FALSE;

  gu = &gu_;
  gv = &gv_;
  gu_old = &gu_old_;
  gv_old = &gv_old_;


  int i,j,l;
  l = 0;
  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) {

      gu_[i][j] = 0;
      gv_[i][j] = 0;
      gu_old_[i][j] = 0;
      gv_old_[i][j] = 0;
    }
  }
  
  for (i = 0 ; i < NUM_PARTICLES ; i++)
    {
      particles[i].x = random_float(SIZE-2)+1;
      particles[i].y = random_float(SIZE-2)+1;
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
    particles[i].x += DT*SIZE*InterpolateScalar(*u, particles[i].x, particles[i].y);
    particles[i].y += DT*SIZE*InterpolateScalar(*v, particles[i].x, particles[i].y);
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

void UpdateFluid()
{
  //VELOCITY COMPUTATIONS
  swap(gu, gu_old);
  swap(gv, gv_old);
  Advect(gu, gu_old, gu_old, gv_old);
  Advect(gv, gv_old, gu_old, gv_old);

  unsigned k;
  for (k=0 ;k < 5 ;k++) {
    if (vel_on) {
      (*gu)[SIZE/4+k][SIZE/2] = reverse ? .5f:-.5f;
      (*gu)[(SIZE/4)*3+k][SIZE/2] = reverse ? -.5f:.5f;
      (*gv)[SIZE/2][SIZE/4+k] = reverse ? .5f:-.5f;
      (*gv)[SIZE/2][(SIZE/4)*3+k] = reverse ? -.5f:.5f;
    }    
  }
  // swap(gu, gu_old);
  //swap(gv, gv_old);
  //Diffuse(gu, gu_old);
  //Diffuse(gv, gv_old);

  Project(gu, gv, gu_old, gv_old);

  //MATERIAL CALCULATIONS
  
  //swap(gmaterial, gmaterial_old);
  //Diffuse(gmaterial, gmaterial_old);
  //  swap(gmaterial, gmaterial_old);
  //Advect(gmaterial, gmaterial_old, gu, gv);

  //swap(bmaterial, bmaterial_old);
  //Diffuse(bmaterial, bmaterial_old);
  //  swap(bmaterial, bmaterial_old);
  //Advect(bmaterial, bmaterial_old, gu, gv);

  // swap(rmaterial, rmaterial_old);
  //Diffuse(rmaterial, rmaterial_old);
  //swap(rmaterial, rmaterial_old);
  //Advect(rmaterial, rmaterial_old, gu, gv);
  
  MoveParticles(gu,gv);
  
  //      run = FALSE;
}

/* function to release/destroy our resources and restoring the old desktop */
void Quit( int returnCode )
{
    /* clean up the window */
    SDL_Quit( );

    /* and exit appropriately */
    exit( returnCode );
}


/* function to reset our viewport after a window resize */
int resizeWindow( int width, int height )
{
    /* Height / width ration */
    GLfloat ratio;

    /* Protect against a divide by zero */
    if ( height == 0 )
    height = 1;

    ratio = ( GLfloat )width / ( GLfloat )height;

    /* Setup our viewport. */
    glViewport( 0, 0, ( GLint )width, ( GLint )height );

    /*
     * change to the projection matrix and set
     * our viewing volume.
     */
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity( );

    /* Set our perspective */
    gluPerspective( 45.0f, ratio, 0.1f, 100.0f );

    /* Make sure we're chaning the model view and not the projection */
    glMatrixMode( GL_MODELVIEW );

    /* Reset The View */
    glLoadIdentity( );

    return( TRUE );
}

/* function to handle key press events */
void handleKeyPress( SDL_keysym *keysym )
{
  switch ( keysym->sym ) {
  case SDLK_ESCAPE:
    /* ESC key was pressed */
    Quit( 0 );
    break;
  case SDLK_F1:
    /* F1 key was pressed
     * this toggles fullscreen mode
     */
    SDL_WM_ToggleFullScreen( surface );
      break;
  case SDLK_v:
    vel_on = !vel_on;
    break;
  case SDLK_b:
    show_vel = !show_vel;
    break;
  case SDLK_m:
    mat_on = !mat_on;
    break;
  case SDLK_r:
    reverse = !reverse;
    break;
  case SDLK_n:
    mat_dir = -mat_dir;
    break;
  case SDLK_g:
    run = !run;
    break;
  case SDLK_w:
    ypos += 0.11;
    break;
  case SDLK_s:
    ypos -= 0.11;
    break;
  case SDLK_a:
    xpos -= 0.11;
    break;
  case SDLK_d:
    xpos += 0.11;
    break;
  case SDLK_q:
    zoom += 0.11;
    break;
  case SDLK_e:
    zoom -= 0.11;
    break;
  default:
    break;
  }
  return;
}

/* general OpenGL initialization function */
int initGL( GLvoid )
{

    InitFluid();

    /* Enable smooth shading */
    glShadeModel(GL_SMOOTH);

    /* Set the background black */
    glClearColor( 0.0f, 0.0f, 0.0f, 0.5f );

    /* Depth buffer setup */
    glClearDepth( 1.0f );

    /* Enables Depth Testing */
    glEnable( GL_DEPTH_TEST );

    /* The Type Of Depth Test To Do */
    glDepthFunc( GL_LEQUAL );

    glBlendFunc(GL_SRC_ALPHA,GL_ONE);
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);

    /* Really Nice Perspective Calculations */
    glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

    return( TRUE );
}

float toGLCoords(float x)
{
  return (x/(float)(SIZE))*2.0f - 1.0 ;
}

/* Here goes our drawing code */
int drawGLScene( GLvoid )
{
  //    usleep(1000000);

    /* These are to calculate our fps */
    static GLint T0     = 0;
    static GLint Frames = 0;

    /* Clear The Screen And The Depth Buffer */
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    /* Move Into The Screen 5 Units */
    glLoadIdentity( );
    glTranslatef( xpos, ypos, zoom );

    /* Select Our Texture */
    if (show_vel){
      glColor4f(1.0f,1.0f,1.0f,0.5f);
      glBegin(GL_LINES);{
        unsigned x,y;
        for (x = 0; x < SIZE; x+=1) {
          for (y = 0; y < SIZE; y+=1) {
            float x_off =  x + 0.5f;
            float y_off =  y + 0.5f;
            glVertex3f(toGLCoords(x_off), toGLCoords(y_off), 1.001f);
            glColor4f(1.0f,1.0f,1.0f,0.5f);
            glVertex3f(toGLCoords(x_off + (3.0f*DT*SIZE*(*gu)[x][y])), 
                       toGLCoords(y_off + (3.0f*DT*SIZE*(*gv)[x][y])),
                       1.001f);
            glColor4f(1.0f,1.0f,1.0f,1.0f);
          }
        }
      }
      glVertex3f( -1.0f,  1.0f, 1.0f );glVertex3f( -1.0f, -1.0f, 1.0f );
      glVertex3f( -1.0f, -1.0f, 1.0f );glVertex3f(  1.0f, -1.0f, 1.0f );
      glVertex3f(  1.0f, -1.0f, 1.0f );glVertex3f(  1.0f,  1.0f, 1.0f );
      glVertex3f(  1.0f,  1.0f, 1.0f );glVertex3f( -1.0f,  1.0f, 1.0f );
      glEnd();
    }
    glColor4f(1.0f,1.0f,1.0f,1.0f);
    glBegin(GL_POINTS);
    int i;
    for (i = 0; i < NUM_PARTICLES; i++)
      glVertex3f(toGLCoords(particles[i].x), toGLCoords(particles[i].y), 1.001f);
    glEnd();
    
    /* Draw it to the screen */
    SDL_GL_SwapBuffers( );

    /* Gather our frames per second */
    Frames++;
    {
    GLint t = SDL_GetTicks();
    if (t - T0 >= 5000) {
        GLfloat seconds = (t - T0) / 1000.0;
        GLfloat fps = Frames / seconds;
        printf("%d frames in %g seconds = %g FPS\n", Frames, seconds, fps);
        T0 = t;
        Frames = 0;
    }
    }

    if (run)
      UpdateFluid();
    return( TRUE );
}

int main( int argc, char **argv )
{
    /* Flags to pass to SDL_SetVideoMode */
    int videoFlags;
    /* main loop variable */
    int done = FALSE;
    /* used to collect events */
    SDL_Event event;
    /* this holds some info about our display */
    const SDL_VideoInfo *videoInfo;
    /* whether or not the window is active */
    int isActive = TRUE;

    /* initialize SDL */
    if ( SDL_Init( SDL_INIT_VIDEO ) < 0 )
    {
        fprintf( stderr, "Video initialization failed: %s\n",
             SDL_GetError( ) );
        Quit( 1 );
    }

    /* Fetch the video info */
    videoInfo = SDL_GetVideoInfo( );

    if ( !videoInfo )
    {
        fprintf( stderr, "Video query failed: %s\n",
             SDL_GetError( ) );
        Quit( 1 );
    }

    /* the flags to pass to SDL_SetVideoMode */
    videoFlags  = SDL_OPENGL;          /* Enable OpenGL in SDL */
    videoFlags |= SDL_GL_DOUBLEBUFFER; /* Enable double buffering */
    videoFlags |= SDL_HWPALETTE;       /* Store the palette in hardware */
    videoFlags |= SDL_RESIZABLE;       /* Enable window resizing */

    /* This checks to see if surfaces can be stored in memory */
    if ( videoInfo->hw_available )
    videoFlags |= SDL_HWSURFACE;
    else
    videoFlags |= SDL_SWSURFACE;

    /* This checks if hardware blits can be done */
    if ( videoInfo->blit_hw )
    videoFlags |= SDL_HWACCEL;

    /* Sets up OpenGL double buffering */
    SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );

    /* get a SDL surface */
    surface = SDL_SetVideoMode( SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP,
                videoFlags );

    /* Verify there is a surface */
    if ( !surface )
    {
        fprintf( stderr,  "Video mode set failed: %s\n", SDL_GetError( ) );
        Quit( 1 );
    }

    /* initialize OpenGL */
    initGL( );

    /* resize the initial window */
    resizeWindow( SCREEN_WIDTH, SCREEN_HEIGHT );

    /* wait for events */
    while ( !done )
    {
        /* handle the events in the queue */

        while ( SDL_PollEvent( &event ) )
        {
            switch( event.type )
            {
            case SDL_ACTIVEEVENT:
                /* Something's happend with our focus
                 * If we lost focus or we are iconified, we
                 * shouldn't draw the screen
                 */
                if ( event.active.gain == 0 )
                isActive = TRUE;
                else
                isActive = TRUE;
                break;
            case SDL_VIDEORESIZE:
                /* handle resize event */
                surface = SDL_SetVideoMode( event.resize.w,
                            event.resize.h,
                            16, videoFlags );
                if ( !surface )
                {
                    fprintf( stderr, "Could not get a surface after resize: %s\n", SDL_GetError( ) );
                    Quit( 1 );
                }
                resizeWindow( event.resize.w, event.resize.h );
                break;
            case SDL_KEYDOWN:
                /* handle key presses */
                handleKeyPress( &event.key.keysym );
                break;
            case SDL_QUIT:
                /* handle quit requests */
                done = TRUE;
                break;
            default:
                break;
            }
        }

        /* draw the scene */
        if ( isActive )
        drawGLScene( );
    }

    /* clean ourselves up and exit */
    Quit( 0 );

    /* Should never get here */
    return( 0 );
}
