/*
 * This code was created by Jeff Molofee '99
 * (ported to Linux/SDL by Ti Leggett '01)
 *
 * If you've found this code useful, please let me know.
 *
 * Visit Jeff at http://nehe.gamedev.net/
 *
 * or for port-specific comments, questions, bugreports etc.
 * email to leggett@eecs.tulane.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <unistd.h>
#include "SDL.h"

/* screen width, height, and bit depth */
#define SCREEN_WIDTH  1000
#define SCREEN_HEIGHT 1000
#define SCREEN_BPP     16
#define HEIGHT 300
#define WIDTH 300
#define X 0
#define Y 1

#define DT 1.0
#define DX 1.0
#define RDX 1.0

/* Set up some booleans */
#define TRUE  1
#define FALSE 0

/* This is our SDL surface */
SDL_Surface *surface;

GLfloat xrot; /* X Rotation ( NEW ) */
GLfloat yrot; /* Y Rotation ( NEW ) */
GLfloat zrot; /* Z Rotation ( NEW ) */

GLuint material_tex;
GLuint velocity_tex;
GLuint pressure_tex;

typedef float Value;
typedef Value ValueArray[WIDTH][HEIGHT];
typedef Value ValueCol[HEIGHT];
ValueArray gu_;
ValueArray* gu;
ValueArray gv_;
ValueArray* gv;
ValueArray gu_old_;
ValueArray* gu_old;
ValueArray gv_old_;
ValueArray* gv_old;
ValueArray gmaterial_;
ValueArray* gmaterial;
ValueArray gmaterial_old_;
ValueArray* gmaterial_old;

int vel_on;

void InitFluid()
{
  vel_on=FALSE;
  gu = &gu_;
  gv = &gv_;
  gu_old = &gu_old_;
  gv_old = &gv_old_;
  gmaterial = &gmaterial_;
  gmaterial_old = &gmaterial_old_;

  int i,j,l;
  l = 0;
  for (i = 0; i < WIDTH; ++i) {
    for (j = 0; j < HEIGHT; ++j) {
      gmaterial_[i][j] = 0;
      gmaterial_old_[i][j] = 0;
      gu_[i][j] = 0;
      gv_[i][j] = 0;
      gu_old_[i][j] = 0;
      gv_old_[i][j] = 0;
    }
  }
  for (i = 1; i < WIDTH-1; ++i) {
    for (j = 1; j < HEIGHT-1; ++j) {
      gmaterial_[i][j] = l++%13 ? 0:1;
    }
  }
  gv_[WIDTH/4][HEIGHT/2] = .01;
  gv_[(WIDTH/4)*3][HEIGHT/2] = .01;
  gu_[WIDTH/2][HEIGHT/4] = .01;
  gu_[WIDTH/2][(HEIGHT/4)*3] = .01;
}

float InterpolateScalar(ValueArray array, float x, float y)
{
  if (x < 0) x = 0;
  if (y < 0) y = 0;
  if (x > WIDTH-1) x = WIDTH-1;
  if (y > HEIGHT-1) y = HEIGHT-1;
  unsigned f_x = x;
  unsigned f_y = y;
  float p_x = x-f_x;
  float p_y = y-f_y;
  return (array[f_x][f_y]*(1.0-p_x)*(1.0-p_y)) + (array[f_x+1][f_y]*p_x*(1.0-p_y)) + (array[f_x][f_y+1]*(1.0-p_x)*p_y) + (array[f_x+1][f_y+1]*p_x*p_y);
}

void JacobiPoissonSolver(int iterations, float alpha, float r_beta, ValueArray* k, ValueArray* b)
{
  /* Solve Ak = b Where A = del^2 */
  unsigned iteration,x,y;
  for (iteration = 0; iteration < iterations; ++iteration) {
    for (x = 1; x < WIDTH-1; ++x) {
      for (y = 1; y < HEIGHT-1; ++y) {
           //printf("%d %d %d\n", iteration, x,y);
           //printf("%d %d %d\n", iteration, x-1,y-1);
          (*k)[x][y] = ((*k)[x-1][y] + (*k)[x+1][y]+ (*k)[x][y-1]+ (*k)[x][y+1] + alpha * (*b)[x][y]) * r_beta;
      }
    }
  }
}

void Advect(ValueArray* dest, ValueArray* source, ValueArray* u, ValueArray* v, float dt)
{
  unsigned x,y;
  for (x = 1; x < WIDTH-1; ++x) {
    for (y = 1; y < HEIGHT-1; ++y) {
      float old_x = (float)x - (dt * RDX * (*u)[x][y]);
      float old_y = (float)y - (dt * RDX * (*v)[x][y]);
      (*dest)[x][y] = InterpolateScalar(*source, old_x, old_y);
    }
  }
}

void Diffuse(ValueArray* dest, ValueArray* source)
{
  JacobiPoissonSolver(20, (DX*DX)/DT, 1.0/(4.0 + (DX*DX)/DT), dest, source);
}

void Divergence(ValueArray* dest, ValueArray* u, ValueArray* v)
{
  unsigned x,y;
  for (x = 1; x < WIDTH-1; ++x) {
    for (y = 1; y < HEIGHT-1; ++y) {
      (*dest)[x][y] = (0.5 + RDX) * (((*u)[x-1][y] - (*u)[x+1][y]) + ((*v)[x][y-1] - (*v)[x][y+1]));
    }
  }
}

void GradientSubtract(ValueArray* u, ValueArray* v, ValueArray* pressure)
{
  unsigned x,y;
  for (x = 1; x < WIDTH-1; ++x) {
    for (y = 1; y < HEIGHT-1; ++y) {
      (*u)[x][y] -= 0.5 * RDX * ((*pressure)[x+1][y] - (*pressure)[x-1][y]);
      (*v)[x][y] -= 0.5 * RDX * ((*pressure)[x][y+1] - (*pressure)[x][y-1]);
    }
  }
}

void Project(ValueArray* u, ValueArray* v, ValueArray* pressure, ValueArray* div)
{
    Divergence(div, u, v);
    unsigned x,y;
    for (x = 0; x < WIDTH; ++x) {
        for (y = 0; y < HEIGHT; ++y) {
            (*pressure)[x][y] = 0;
        }
    }
    JacobiPoissonSolver(40, -(DX*DX), 1.0/4.0, pressure, div);
    GradientSubtract(u, v, pressure);
}

void ApplyPressureBoundary(ValueArray* pressure)
{
  unsigned x,y;
  for (x = 1; x < WIDTH-1; ++x) {
    (*pressure)[x][0] = (*pressure)[x][1];
    (*pressure)[x][HEIGHT-1] = (*pressure)[x][HEIGHT-2];
  }
  for (y = 1; y < HEIGHT-1; ++y) {
    (*pressure)[0][y] = (*pressure)[1][y];
    (*pressure)[WIDTH-1][y]= (*pressure)[WIDTH-2][y];
  }
}

void ApplyVelocityBoundary(ValueArray* u, ValueArray* v)
{
  unsigned x,y;
  for (x = 1; x < WIDTH-1; ++x) {
    (*u)[x][0] = -(*u)[x][1];
    (*v)[x][0] = -(*v)[x][1];
    (*u)[x][HEIGHT-1] = -(*u)[x][HEIGHT-2];
    (*v)[x][HEIGHT-1] = -(*v)[x][HEIGHT-2];
  }
  for (y = 1; y < HEIGHT-1; ++y) {
    (*u)[0][y] = -(*u)[1][y];
    (*v)[0][y] = -(*v)[1][y];
    (*u)[WIDTH-1][y] = -(*u)[WIDTH-2][y];
    (*v)[WIDTH-1][y] = -(*v)[WIDTH-2][y];
  }
}

void PrintValue(ValueArray* array)
{
  unsigned x,y;
  for (x = 0; x < WIDTH; ++x) {
    for (y = 0; y < HEIGHT; ++y) {
      printf("%2.1f %2.1f  ", (*array)[x][y], (*array)[x][y]);
    }
    printf("\n");
  }
}

/*void PrintDebug()
{
  printf("Velocity\n");
  PrintPair(velocity);
  printf("Pressure\n");
  PrintPair(velocity);
  printf("TempPair\n");
  PrintPair(temppair);
  printf("Material\n");
  PrintValue(material);
  printf("TempVal\n");
  PrintValue(tempval);
  printf("\n");
}*/

inline void Swap(ValueArray* x, ValueArray* y)
{
    ValueArray* temp = x;
    x=y;
    y=temp;
}

void UpdateFluid()
{
  //VELO
  //  PrintDebug();
  //Add mat and force;
  Swap(gu, gu_old);
  Swap(gv, gv_old);
  Diffuse(gu, gu_old);
  Diffuse(gv, gv_old);
  Project(gu, gv, gu_old, gv_old);
  Swap(gu, gu_old);
  Swap(gv, gv_old);
  Advect(gu, gu_old, gu_old, gv_old, DT);
  Advect(gv, gv_old, gu_old, gv_old, DT);
  Project(gu, gv, gu_old, gv_old);

  //MAT
  //add source
  Swap(gmaterial, gmaterial_old);
  Diffuse(gmaterial, gmaterial_old);
  Swap(gmaterial, gmaterial_old);
  Advect(gmaterial, gmaterial_old, gu, gv, DT);

  if (vel_on) {
    (*gu)[WIDTH/2][HEIGHT/2] = 0.3;
    (*gv)[WIDTH/2][HEIGHT/2] = 0.3;
  }


}

/* function to release/destroy our resources and restoring the old desktop */
void Quit( int returnCode )
{
    /* clean up the window */
    SDL_Quit( );

    /* and exit appropriately */
    exit( returnCode );
}

/* function to load in bitmap as a GL texture */
int LoadGLTextures()
{
    /* Create storage space for the texture */
    GLubyte TextureImage[WIDTH][HEIGHT][4];
    int i,j;
    for (i = 0; i < WIDTH; ++i) {
      for (j = 0; j < HEIGHT; ++j) {
        TextureImage[i][j][0] = (*gmaterial)[i][j] > 0 ? (*gmaterial)[i][j]*255 : 0;
        TextureImage[i][j][1] = (*gmaterial)[i][j] > 0 ? (*gmaterial)[i][j]*255 : 0;
        TextureImage[i][j][2] = (*gmaterial)[i][j] > 0 ? (*gmaterial)[i][j]*255 : 0;
        TextureImage[i][j][3] = 127;
      }
    }

    /* Create The Texture */
    glGenTextures( 1, &material_tex);

    /* Typical Texture Generation Using Data From The Bitmap */
    glBindTexture( GL_TEXTURE_2D, material_tex);

    /* Generate The Texture */
    glTexImage2D( GL_TEXTURE_2D, 0, 4, WIDTH,
                  HEIGHT, 0, GL_RGBA,
                  GL_UNSIGNED_BYTE, TextureImage );

    /* Linear Filtering */
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    for (i = 0; i < WIDTH; ++i) {
      for (j = 0; j < HEIGHT; ++j) {
        TextureImage[i][j][0] = ((*gu)[i][j]-0.5)*255;
        TextureImage[i][j][1] = ((*gv)[i][j]-0.5)*255;
        TextureImage[i][j][2] = 0;
        TextureImage[i][j][3] = 127;
      }
    }

    /* Create The Texture */
    glGenTextures( 1, &velocity_tex);

    /* Typical Texture Generation Using Data From The Bitmap */
    glBindTexture( GL_TEXTURE_2D, velocity_tex);

    /* Generate The Texture */
    glTexImage2D( GL_TEXTURE_2D, 0, 4, WIDTH,
                  HEIGHT, 0, GL_RGBA,
                  GL_UNSIGNED_BYTE, TextureImage );

    /* Linear Filtering */
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    for (i = 0; i < WIDTH; ++i) {
      for (j = 0; j < HEIGHT; ++j) {
        TextureImage[i][j][0] = 0;
        TextureImage[i][j][1] = ((*gu_old)[i][j]-0.5)*255;
        TextureImage[i][j][2] = 0;
        TextureImage[i][j][3] = 127;
      }
    }

    /* Create The Texture */
    glGenTextures( 1, &pressure_tex);

    /* Typical Texture Generation Using Data From The Bitmap */
    glBindTexture( GL_TEXTURE_2D, pressure_tex);

    /* Generate The Texture */
    glTexImage2D( GL_TEXTURE_2D, 0, 4, WIDTH,
                  HEIGHT, 0, GL_RGBA,
                  GL_UNSIGNED_BYTE, TextureImage );

    /* Linear Filtering */
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    return 1;
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
  case SDLK_F2:
    vel_on = !vel_on;
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
    /* Enable Texture Mapping ( NEW ) */
    glEnable( GL_TEXTURE_2D );

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

    /*glBlendFunc(GL_SRC_ALPHA,GL_ONE);
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);*/

    /* Really Nice Perspective Calculations */
    glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

    return( TRUE );
}

/* Here goes our drawing code */
int drawGLScene( GLvoid )
{
  //    usleep(1000000);

    /* These are to calculate our fps */
    LoadGLTextures();
    static GLint T0     = 0;
    static GLint Frames = 0;

    /* Clear The Screen And The Depth Buffer */
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    /* Move Into The Screen 5 Units */
    glLoadIdentity( );
    glTranslatef( 0.0f, 0.0f, -4.0f );

    glRotatef( xrot, 1.0f, 0.0f, 0.0f); /* Rotate On The X Axis */
    glRotatef( yrot, 0.0f, 1.0f, 0.0f); /* Rotate On The Y Axis */
    glRotatef( zrot, 0.0f, 0.0f, 1.0f); /* Rotate On The Z Axis */

    /* Select Our Texture */
    glBindTexture( GL_TEXTURE_2D, material_tex);
    glBegin(GL_QUADS);
      glTexCoord2f( 0.0f, 0.0f ); glVertex3f( 0.0f, 0.0f, 1.0f );
      glTexCoord2f( 1.0f, 0.0f ); glVertex3f(  1.0f, 0.0f, 1.0f );
      glTexCoord2f( 1.0f, 1.0f ); glVertex3f(  1.0f,  1.0f, 1.0f );
      glTexCoord2f( 0.0f, 1.0f ); glVertex3f( 0.0f,  1.0f, 1.0f );
    glEnd( );
    glBindTexture( GL_TEXTURE_2D, velocity_tex);
    glBegin(GL_QUADS);
      glTexCoord2f( 0.0f, 0.0f ); glVertex3f( 0.0f, -1.0f, 1.0f );
      glTexCoord2f( 1.0f, 0.0f ); glVertex3f(  1.0f, -1.0f, 1.0f );
      glTexCoord2f( 1.0f, 1.0f ); glVertex3f(  1.0f,  0.0f, 1.0f );
      glTexCoord2f( 0.0f, 1.0f ); glVertex3f( 0.0f,  0.0f, 1.0f );
    glEnd( );
    glBindTexture( GL_TEXTURE_2D, pressure_tex);
    glBegin(GL_QUADS);
      glTexCoord2f( 0.0f, 0.0f ); glVertex3f( -1.0f, -1.0f, 1.0f );
      glTexCoord2f( 1.0f, 0.0f ); glVertex3f(  0.0f, -1.0f, 1.0f );
      glTexCoord2f( 1.0f, 1.0f ); glVertex3f(  0.0f,  0.0f, 1.0f );
      glTexCoord2f( 0.0f, 1.0f ); glVertex3f( -1.0f,  0.0f, 1.0f );
    glEnd( );

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

    xrot += 0.00f; /* X Axis Rotation */
    yrot += 0.00f; /* Y Axis Rotation */
    zrot += 0.00f; /* Z Axis Rotation */

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
                isActive = FALSE;
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
