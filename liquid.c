#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <unistd.h>
#include "SDL.h"
#include "solver.h"

/* screen width, height, and bit depth */
#define SCREEN_WIDTH  500
#define SCREEN_HEIGHT 500
#define FULL_WIDTH 1240
#define FULL_HEIGHT 1028
#define SCREEN_BPP     16

/* Set up some booleans */
#define TRUE  1
#define FALSE 0

/* This is our SDL surface */
SDL_Surface *surface;
/* Flags to pass to SDL_SetVideoMode */
int videoFlags;

GLfloat xpos; 
GLfloat ypos; 
GLfloat zoom; 

int vel_on;
int show_vel;
int reverse;
int run;

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
    /* Make sure we're changing the model view and not the projection */
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
    surface = SDL_SetVideoMode( FULL_WIDTH,
                                FULL_HEIGHT,
                                16, videoFlags );
    if ( !surface )
      {
        fprintf( stderr, "Could not get a surface after resize: %s\n", SDL_GetError( ) );
        Quit( 1 );
      }
    resizeWindow(FULL_WIDTH,  FULL_HEIGHT);
    SDL_WM_ToggleFullScreen( surface );
    break;
  case SDLK_v:
    vel_on = !vel_on;
    break;
  case SDLK_b:
    show_vel = !show_vel;
    break;
  case SDLK_g:
    run = !run;
    break;
  case SDLK_r:
    ResetFluid();
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
  xpos = 0;
  ypos = 0;
  zoom = -3.5;
  vel_on=TRUE;
  show_vel=FALSE;
  reverse=FALSE;
  run = FALSE;
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
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  /* Really Nice Perspective Calculations */
  glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );
  glEnable(GL_POINT_SMOOTH);
  return( TRUE );
}

float toGLCoords(float x)
{
  return (x/(float)(SIZE))*2.0f - 1.0 ;
}

/* Here goes our drawing code */
int drawGLScene( GLvoid )
{
  /* These are to calculate our fps */
  static GLint T0     = 0;
  static GLint Frames = 0;
  /* Clear The Screen And The Depth Buffer */
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  /* Move Into The Screen 5 Units */
  glLoadIdentity( );
  glTranslatef( xpos, ypos, zoom );

  if (show_vel){
    glColor4f(1.0f,1.0f,1.0f,0.5f);
    glBegin(GL_LINES);{
      unsigned x,y;
      for (x = 0; x < SIZE; x+=1) {
        for (y = 0; y < SIZE; y+=1) {
          float x_off =  x;
          float y_off =  y;
          glVertex3f(toGLCoords(x_off), toGLCoords(y_off), 1.001f);
          glColor4f(1.0f,1.0f,1.0f,0.5f);
          glVertex3f(toGLCoords(x_off + (1.0f*DT*SIZE*(*gu)[x][y])), 
                     toGLCoords(y_off + (1.0f*DT*SIZE*(*gv)[x][y])),
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
  glPointSize(1.0);
  glBegin(GL_POINTS);{
    int i;
    for (i = 0; i < num_particles; i++)
      glVertex3f(toGLCoords((*particles)[i].x), toGLCoords((*particles)[i].y), 1.001f);
    glEnd();
  }
  unsigned x,y;
  for (x = 0; x < SIZE-1; x+=1) {
    for (y = 0; y < SIZE-1; y+=1) {
      float x_off =  x + 0.5f;
      float y_off =  y + 0.5f;
      if ((*divergance)[x][y]*5000.0f > 0)
        glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
      else
        glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
      glPointSize(fabs((*divergance)[x][y]*5000.0f));
      glBegin(GL_POINTS);{
        glVertex3f(toGLCoords(x_off), toGLCoords(y_off), 1.001f);
        glEnd();
      }
    }
  }
  /* Draw it to the screen */
  SDL_GL_SwapBuffers( );
  /* Gather our frames per second */
  Frames++;
  
  GLint t = SDL_GetTicks();
  if (t - T0 >= 5000) {
    GLfloat seconds = (t - T0) / 1000.0;
    GLfloat fps = Frames / seconds;
    printf("%d frames in %g seconds = %g FPS\n", Frames, seconds, fps);
    T0 = t;
    Frames = 0;
  }
  if (run)
    UpdateFluid();
      run = FALSE;
  return(TRUE);
}

int main( int argc, char **argv )
{
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
