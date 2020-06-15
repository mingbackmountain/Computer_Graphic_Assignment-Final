//Nuengruethai Wutthisak 6088138
//Kanlayakorn Kesorn 6088057
//Thanakorn Pasangthien 6088109
/*____________________________________________________________________
|
| File: main.cpp
|
| Description: Hermite curve drawing.
|
| Functions:
|
| (C) Copyright 2007 Mores Prachyabrued
|___________________________________________________________________*/

// Enable older CRT functions (such as strcpy) without warnings from vc8 (vc 2005 .net)
#if _MSC_VER >= 1400
#define _CRT_SECURE_NO_DEPRECATE                // shut up the vc8 compiler
#define _CRT_NONSTDC_NO_DEPRECATE
#endif

/*___________________
|
| Include Files
|__________________*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gmtl/gmtl.h>

#include <GL/glut.h>

/*___________________
|
| Constants
|__________________*/
#define DEF_D 5

const int WIN_WIDTH_INIT = 800;
const int WIN_HEIGHT_INIT = 600;

const float C_TSTEP = 0.0005f;            // Curve time step (for rendering curve)

// plane
const float P_WIDTH = 3;
const float P_LENGTH = 3;
const float P_HEIGHT = 1.5f;

const gmtl::Vec3f PLANE_FORWARD(0, 0, 1.0f);            // Plane's forward translation vector (w.r.t. local frame)
const float PLANE_ROTATION = 5.0f;                      // Plane rotated by 5 degs per input


// Propeller dimensions (subpart)
const float PP_WIDTH = 1.5f;
const float PP_LENGTH = 1.5f;
const float PP_HEIGHT = 1.5f;

//Light transforms
const gmtl::Point3f LIGHT_POS(2.5, 1, 1.5);     // Propeller position on the plane (w.r.t. plane's frame)
const float LIGHT_ROTATION = 5.0f;                  // Propeller rotated by 5 degs per input

// Propeller transforms
const gmtl::Point3f PROPELLER_POS(2.5, 1, 1.5);     // Propeller position on the plane (w.r.t. plane's frame)
const float PROPELLER_ROTATION = 5.0f;                  // Propeller rotated by 5 degs per input


// Propeller rotation (subpart)
float pp_angle = 0;         // Rotation angle
float pp_angle_2 = 0;         // Rotation angle for subpart wing2
float pp_angle_3 = 0;         // Rotation angle for subpart wing3
float gun_angle = 0;          // Rotation angle for sub subpart gun

// Plane pose (position-quaternion pair)
gmtl::Point4f plane_p;      // Position (using explicit homogeneous form; see Quaternion example code)
gmtl::Quatf plane_q;        // Quaternion

gmtl::Point4f plane_p1;      // Position (using explicit homogeneous form; see Quaternion example code)
gmtl::Quatf plane_q1;        // Quaternion

// Quaternions to rotate plane
gmtl::Quatf zrotp_q;        // Positive and negative Z rotations
gmtl::Quatf zrotn_q;

gmtl::Quatf xrotp_q;
gmtl::Quatf xrotn_q;

gmtl::Quatf yrotp_q;
gmtl::Quatf yrotn_q;

const gmtl::Vec3f WORLD_UP(0, 1, 0);                  // World up axis (Y)

const int PT_NB = 6;                                  // Number of input points (equals #tangents and #segments)
const gmtl::Point3f input_pts[PT_NB] =
{
    gmtl::Point3f(100, 50, 10),
    gmtl::Point3f(0, 45, 0),
    gmtl::Point3f(-200, 0, -350),
    gmtl::Point3f(0, 35, -150),
    gmtl::Point3f(-190, -20, 0),
    gmtl::Point3f(0, -30, 160)
};

gmtl::Matrix44f MMAT;                                 // Basis matrix for Hermite curve form (const)

const float R_MAX = 2.5f * 1200.0f;                     // Expected maximum magnitude of local rightward component of plane's acceleration

const int PARTICLE_NB = 40;                              // Number of particles

const float VMAG_MEAN = 100.0f;                       // Velocity
const float VMAG_STD = 25.0f;
const float VDIR_STD = 0.25f;

const int TTL_BASE = 200;                          // Time to live
const int TTL_OFFSET = 50;

const float SMOKE_SIZE = 15;                         // Smoke size

const float S_TSTEP = 0.001f;                          // Simulation time step (for particle update)

const gmtl::Vec3f GRAVITY(0, -100.0f, 0);                 // w.r.t. world, can be upward for smoke

const gmtl::Vec3f V_WIND(100, 0, 0);                  // Wind velocity and drag coefficient (K)
const float K_COEF = 1.0f;

// Keyboard modifiers
enum KeyModifier { KM_SHIFT = 0, KM_CTRL, KM_ALT };

/*___________________
|
| Type Definitions
|__________________*/

// Particle structure storing position, velocity, and other properties
typedef struct _MyParticle {
    gmtl::Point3f p;      // Position
    gmtl::Vec3f v;        // Velocity
    float m;              // Mass
    int ttl;              // Time to live, decremented each iteration
    int ttl_init;                  // Initial ttl, used to fade color as the age increases
} MyParticle;

/*___________________
|
| Global variables
|__________________*/

int plane_id = 0;
// Track window dimensions, initialized to 800x600
int w_width = 800;
int w_height = 600;


// camera w.r.t. plane
float distance = 600.0f;
float elevation = -15.0f;                 // In degs
float azimuth = 180.0f;                 // In degs

// Mouse & keyboard
int mx_prev = 0, my_prev = 0;
bool mbuttons[3] = { false, false, false };
bool kmodifiers[3] = { false, false, false };



gmtl::Vec3f tangents[PT_NB];                          // Tangent at each input point
float s_tan = 1.0f;                                   // A scaling factor for tangents

gmtl::Matrix44f Cmats[PT_NB];                         // Coefficient matrix (C) for each Hermite curve segment (It's actually 4x3, ignore last column)

gmtl::Matrix44f ppose;                                // The plane's pose
gmtl::Matrix44f pposeadj;                             // Adjusted plane coordinate system that the (plane) camera will be attached to
gmtl::Matrix44f pposeadj_inv;                         // Adjusted plane coordinate system (plane's camera is attached to this frame), inverted
int ps = 0;                                        // The segment in which the plane currently belongs
float pt = 0;                                        // The current t-value for the plane
float pdt = 0.02f;                                    // delta_t for the plane

MyParticle particles[PARTICLE_NB];                          // Array of particles

enum TextureID {
    TID_SKYBACK = 0, TID_SKYLEFT, TID_SKYBOTTOM, TID_SKYFRONT, TID_SKYRIGHT, TID_SKYTOP,TEXTURN_BUILD1,TEXTURE_BUILD2,TEXTURE_TOPBUILD,
    TEXTURE_NB, PARTICLE1, PARTICLE2
};  // Texture IDs, with the last ID indicating the total number of textures

//skybox
const float SB_SIZE = 2000.0f;

//texture
GLuint textures[TEXTURE_NB];
GLuint texture;

// Lighting
gmtl::Point4f light_pos(-90.0, 100.0, 0.0, 1.0);
bool is_diffuse_on = true;
bool is_ambient_on = true;
bool is_specular_on = true;

// Rendering option
bool render_curve = true;
bool render_constraint = false;

// Lighting
const GLfloat NO_LIGHT[] = { 0.0, 0.0, 0.0, 1.0 };
const GLfloat AMBIENT_LIGHT[] = { 5, 5, 5, 1.0 };
const GLfloat DIFFUSE_LIGHT[] = { 7, 7, 7, 1.0 };
const GLfloat SPECULAR_LIGHT[] = { 10, 10, 10, 1.0 };

// Materials
const GLfloat DARKRED_COL[] = { 0.1, 0.0, 0.0, 1.0 };
const GLfloat DARKGREEN_COL[] = { 0.0, 0.0, 1.0, 1.0 };
const GLfloat BRIGHTGREEN_COL[] = { 0.0, 0.0, 2.0, 2.0 };
const GLfloat BRIGHTRED2_COL[] = { 2.0, 0.0, 0.0, 2.0 };
const GLfloat BRIGHTRED_COL[] = { 1.5, 0.0, 0.0, 1.0 };
const GLfloat DARKBLUE_COL[] = { 0.0, 0.0, 0.1, 1.0 };
const GLfloat BRIGHTBLUE_COL[] = { 0.0, 0.0, 0.7, 1.0 };
const GLfloat BRIGHTBLUE2_COL[] = { 0.0, 0.0, 2.0, 1.0 };
const GLfloat DARK_COL[] = { 0.1, 0.1, 0.1, 1.0 };
const GLfloat MEDIUMWHITE_COL[] = { 0.7, 0.7, 0.7, 1.0 };
const GLfloat SPECULAR_COL[] = { 0.7, 0.7, 0.7, 1.0 };

/*___________________
|
| Function Prototypes
|__________________*/

void Init();
void Init_AParticle(MyParticle& par);
void Compute_Tangents();
void Compute_Coefficiences();
void Display_Func();
void Idle_Func();
void Update_Plane();
void Update_PSystem();
void Keyboard_Func(unsigned char key, int x, int y);
void Reshape_Func(int width, int height);
void Mouse_Func(int button, int state, int x, int y);
void Motion_Func(int x, int y);
void Draw_World_Axes();
void Draw_Path();
void Draw_Path2();
void InitTransforms();
void Draw_Particles2();
void Draw_Particles();
void drawLight(const float radius);
void DrawPlaneBody(const float width, const float length, const float height);
float FastGauss(float mean, float std);
void DrawSkybox(const float s);
void LoadPPM(const char* fname, unsigned int* w, unsigned int* h, unsigned char** data, const int mallocflag);
void SetLight(const gmtl::Point4f& pos, const bool is_ambient, const bool is_diffuse, const bool is_specular);


/*____________________________________________________________________
|
| Function: main
|
| Input:
| Output: Program entry point.
|___________________________________________________________________*/

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(WIN_WIDTH_INIT, WIN_HEIGHT_INIT);
    glutCreateWindow("Final Assignment");

    glutDisplayFunc(Display_Func);
    glutIdleFunc(Idle_Func);
    glutKeyboardFunc(Keyboard_Func);
    glutReshapeFunc(Reshape_Func);
    glutMouseFunc(Mouse_Func);
    glutMotionFunc(Motion_Func);

    Init();

    glutMainLoop();

    return 0;
}

//|____________________________________________________________________
//|
//| Function: InitTransforms
//|
//! \param None.
//! \return None.
//!
//! Initializes all the transforms
//|____________________________________________________________________

void InitTransforms()
{
    const float COSTHETA_D2 = cos(gmtl::Math::deg2Rad(PLANE_ROTATION / 2));  // cos() and sin() expect radians
    const float SINTHETA_D2 = sin(gmtl::Math::deg2Rad(PLANE_ROTATION / 2));

    // Inits plane pose
    plane_p.set(1.0f, 0.0f, 4.0f, 1.0f);
    plane_q.set(0, 0, 0, 1);
    plane_p1.set(-4.0f, 0.0f, 4.0f, 1.0f);
    plane_q1.set(0, 0, 0, 1);

    // Z rotations (roll)
    zrotp_q.set(0, 0, SINTHETA_D2, COSTHETA_D2);      // +Z
    zrotn_q = gmtl::makeConj(zrotp_q);                // -Z

    // X rotation (pitch)
    xrotp_q.set(SINTHETA_D2, 0, 0, COSTHETA_D2);      // +X
    xrotn_q = gmtl::makeConj(zrotp_q);                 //-X

    // Y rotation (yaw)
    yrotp_q.set(0, SINTHETA_D2, 0, COSTHETA_D2);      // +Y
    yrotn_q = gmtl::makeConj(zrotp_q);                 //-Y


    // TODO: Initializes the remaining transforms
}

/*____________________________________________________________________
|
| Function: Init
|
| Input:
| Output: Initialization routine.
|___________________________________________________________________*/

void Init()
{
    int i;
    unsigned int texwidth, texheight;
    unsigned char* imagedata;

    // OpenGL
    glClearColor(0, 0, 0, 0);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);

    // Inits Hermite basis matrix
    MMAT.set(
        2, -2, 1, 1,
        -3, 3, -2, -1,
        0, 0, 1, 0,
        1, 0, 0, 0
    );

    // Inits tangents and coefficient matrices
    Compute_Tangents();
    Compute_Coefficiences();

    // Init adjusted plane-coordinate-system that the camera is attached to (HACK! by calling Update_Plane here)
    Update_Plane();

    // Inits particle system
    for (i = 0; i < PARTICLE_NB; i++) {
        Init_AParticle(particles[i]);
    }

    /*___________________________________________________________________
    |
    | Setup lighting
    |___________________________________________________________________*/

    // Disable global ambient
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, NO_LIGHT);

    // NOTE: for specular reflections, the "local viewer" model produces better
    // results than the default, but is slower. The default would not use the correct
    // vertex-to-eyepoint vector, treating it as always parallel to Z.
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);

    // Enable two sided lighting
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    // Enable lighting
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    /*____________________________________________________________________
    |
    | Load texture
    |___________________________________________________________________*/

    //particle 1

   // describe how data will be stored in memory
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // generate a new "texture object" and select it for setup
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, PARTICLE1);

    // load an image into memory
    LoadPPM("smoketex.ppm", &texwidth, &texheight, &imagedata, 1);

    // describe the image to the graphics system as a texture map
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);

    // select methods for "scaling" a texture region to a pixel
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // select the method for combining texture color with the lighting equation
    // (look up the third parameter)



    //particle 2



    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, PARTICLE2);
    LoadPPM("fire.PPM", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);


    // texture for skybox and building

    // Generate and setup texture objects
    glGenTextures(TEXTURE_NB, textures);

    // Skybox back wall
    glBindTexture(GL_TEXTURE_2D, textures[TID_SKYBACK]);
    LoadPPM("sky.ppm", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Skybox left wall
    glBindTexture(GL_TEXTURE_2D, textures[TID_SKYLEFT]);
    LoadPPM("sky.ppm", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Skybox bottom wall
    glBindTexture(GL_TEXTURE_2D, textures[TID_SKYBOTTOM]);
    LoadPPM("sky.ppm", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Skybox front wall
    glBindTexture(GL_TEXTURE_2D, textures[TID_SKYFRONT]);
    LoadPPM("sky.ppm", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Skybox right wall
    glBindTexture(GL_TEXTURE_2D, textures[TID_SKYRIGHT]);
    LoadPPM("sky.ppm", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Skybox top wall
    glBindTexture(GL_TEXTURE_2D, textures[TID_SKYTOP]);
    LoadPPM("sky.ppm", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Building1
    glBindTexture(GL_TEXTURE_2D, textures[TEXTURN_BUILD1]);
    LoadPPM("one.ppm", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Building2
    glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_BUILD2]);
    LoadPPM("two.ppm", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Building top
    glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_TOPBUILD]);
    LoadPPM("three.ppm", &texwidth, &texheight, &imagedata, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texwidth, texheight, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata);
    free(imagedata);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);


}

/*____________________________________________________________________
|
| Function: Init_AParticle
|
| Input:
| Output: Init a single particle.
|___________________________________________________________________*/

void Init_AParticle(MyParticle& par)
{
    gmtl::Vec3f v_dir;

    // position (= plane position)
    par.p.set(ppose[0][3], ppose[1][3], ppose[2][3]);

    // velocity (consider plane's -Z axis as mean velocity direction)
    v_dir.set(-ppose[0][2], -ppose[1][2], -ppose[2][2]);
    par.v.set(
        FastGauss(v_dir[0], VDIR_STD),
        FastGauss(v_dir[1], VDIR_STD),
        FastGauss(v_dir[2], VDIR_STD)
    );
    par.v = FastGauss(VMAG_MEAN, VMAG_STD) * par.v;

    // mass
    par.m = 1;

    // ttl
    par.ttl = par.ttl_init = (rand() % TTL_BASE) + TTL_OFFSET;
}

/*____________________________________________________________________
|
| Function: Compute_Tangents
|
| Input:
| Output: Computes tangents from input points and the scale factor
|___________________________________________________________________*/

void Compute_Tangents()
{
    int i;
    int prev, next;       // Adjacent points

    for (i = 0; i < PT_NB; i++) {
        // Compute adjacent points
        if (i == 0) { // First point
            prev = PT_NB - 1;
            next = i + 1;
        }
        else
            if (i == PT_NB - 1) { // Last point
                prev = i - 1;
                next = 0;
            }
            else { // Interior point
                prev = i - 1;
                next = i + 1;
            }

        tangents[i] = s_tan * (input_pts[next] - input_pts[prev]);
    }
}

/*____________________________________________________________________
|
| Function: Compute_Coefficiences
|
| Input:
| Output: Computes coefficiences matrices for curve segments
|___________________________________________________________________*/

void Compute_Coefficiences()
{
    int i;
    int n;                               // Next input point (index)
    gmtl::Matrix44f Gmat;                // Geometry matrix (It's actually 4x3, ignore last column)

    for (i = 0; i < PT_NB; i++) {
        if (i == PT_NB - 1) {
            n = 0;
        }
        else {
            n = i + 1;
        }

        Gmat.set(
            input_pts[i][0], input_pts[i][1], input_pts[i][2], 0,
            input_pts[n][0], input_pts[n][1], input_pts[n][2], 0,
            tangents[i][0], tangents[i][1], tangents[i][2], 0,
            tangents[n][0], tangents[n][1], tangents[n][2], 0
        );

        Cmats[i] = MMAT * Gmat;
    }
}

void Compute_Coefficiences2()
{
    int i;
    int n;                               // Next input point (index)
    gmtl::Matrix44f Gmat;                // Geometry matrix (It's actually 4x3, ignore last column)

    for (i = 0; i < PT_NB; i++) {
        if (i == PT_NB - 1) {
            n = 0;
        }
        else {
            n = i + 1;
        }

        Gmat.set(
            input_pts[i][0], input_pts[i][1], input_pts[i][2], 100,
            input_pts[n][0], input_pts[n][1], input_pts[n][2], 100,
            tangents[i][0], tangents[i][1], tangents[i][2], 20,
            tangents[n][0], tangents[n][1], tangents[n][2], 20
        );

        Cmats[i] = MMAT * Gmat;
    }
}

/*____________________________________________________________________
|
| Function: Display_Func
|
| Input:
| Output: GLUT display callback function.
|___________________________________________________________________*/

void Display_Func(void)
{
    gmtl::AxisAnglef aa;    // Converts plane's quaternion to axis-angle form to be used by glRotatef()
    gmtl::Vec3f axis;       // Axis component of axis-angle representation
    float angle;            // Angle component of axis-angle representation

    // Clear screen & depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

   

    // Add view matrix
    glTranslatef(0, 0, -distance);
    glRotatef(-elevation, 1, 0, 0);
    glRotatef(-azimuth, 0, 1, 0);
    glMultMatrixf(pposeadj_inv.mData);

   DrawSkybox(SB_SIZE);

    // Set light position wrt world

    SetLight(light_pos, is_ambient_on, is_diffuse_on, is_specular_on);

    glPushMatrix();
    glTranslatef(light_pos[0], light_pos[1], light_pos[2]);
    drawLight(9);
    glPopMatrix();

    Draw_Path();

    
    DrawPlaneBody(P_WIDTH, P_LENGTH, P_HEIGHT);
    Draw_Particles();
    Draw_Particles2();
    
    
    
 
    

    // Display!
    glutSwapBuffers();
}

/*____________________________________________________________________
|
| Function: Idle_Func
|
| Input:
| Output: GLUT idle callback function. Called for each simulation step.
|___________________________________________________________________*/

void Idle_Func()
{
    Update_Plane();
    Update_PSystem();

    glutPostRedisplay();
}

/*____________________________________________________________________
|
| Function: Update_Plane
|
| Input:
| Output: Called every simulation step to update the plane's pose.
|___________________________________________________________________*/

void Update_Plane()
{
    float x, y, z;                    // Plane's coordinates
    gmtl::Vec3f xa, ya, za;           // Plane's axes
    gmtl::Vec3f ac;                   // Acceleration vector
    float r;                          // Local rightward component of acceleration
    float ra;                         // Roll angle
    gmtl::AxisAnglef roll_aa;         // Local roll in axis angle and matrix forms
    gmtl::Matrix44f roll_mat;
    gmtl::Vec3f zadj;

    /*____________________________________________________________________
    |
    | Update t and possibly current segment
    |___________________________________________________________________*/

    pt += pdt;
    if (pt > 1) {
        pt -= 1;
        ps++;
        if (ps == PT_NB) {
            ps = 0;
        }
    }

    /*____________________________________________________________________
    |
    | Compute plane's new position by evaluating the polynomials (use Horner's rule for slight speedup)
    |___________________________________________________________________*/

    x = ((Cmats[ps][0][0] * pt + Cmats[ps][1][0]) * pt + Cmats[ps][2][0]) * pt + Cmats[ps][3][0];
    y = ((Cmats[ps][0][1] * pt + Cmats[ps][1][1]) * pt + Cmats[ps][2][1]) * pt + Cmats[ps][3][1];
    z = ((Cmats[ps][0][2] * pt + Cmats[ps][1][2]) * pt + Cmats[ps][2][2]) * pt + Cmats[ps][3][2];

    /*____________________________________________________________________
    |
    | Compute plane's orientation
    |___________________________________________________________________*/

    // Compute direction of motion (z) by evaluating the polynomial derivatives (use Horner's rule for slight speedup)
    za[0] = (3 * Cmats[ps][0][0] * pt + 2 * Cmats[ps][1][0]) * pt + Cmats[ps][2][0];
    za[1] = (3 * Cmats[ps][0][1] * pt + 2 * Cmats[ps][1][1]) * pt + Cmats[ps][2][1];
    za[2] = (3 * Cmats[ps][0][2] * pt + 2 * Cmats[ps][1][2]) * pt + Cmats[ps][2][2];
    gmtl::normalize(za);

    // Compute lateral axis (x)
    gmtl::cross(xa, WORLD_UP, za);
    gmtl::normalize(xa);

    // Compute up axis (x)
    gmtl::cross(ya, za, xa);
    gmtl::normalize(ya);

    /*____________________________________________________________________
    |
    | Compute banked turn (local roll)
    |___________________________________________________________________*/

    // Compute acceleration vector
    ac[0] = 6 * Cmats[ps][0][0] * pt + 2 * Cmats[ps][1][0];
    ac[1] = 6 * Cmats[ps][0][1] * pt + 2 * Cmats[ps][1][1];
    ac[2] = 6 * Cmats[ps][0][2] * pt + 2 * Cmats[ps][1][2];

    // Compute local rightward component of acceleration
    r = gmtl::dot(ac, xa);

    // Clamp r to R_MAX
    if (r > R_MAX) {
        r = R_MAX;
    }
    else
        if (r < -R_MAX) {
            r = -R_MAX;
        }

    // Compute roll angle
    ra = asin(r / R_MAX);

    //printf("%.2f\n", r);

  /*____________________________________________________________________
  |
  | Update plane's pose
  |___________________________________________________________________*/

    ppose.set(
        xa[0], ya[0], za[0], x,
        xa[1], ya[1], za[1], y,
        xa[2], ya[2], za[2], z,
        0, 0, 0, 1
    );
    ppose.setState(gmtl::Matrix44f::AFFINE);

    // Compute local roll (rotation about longitudinal axis (z))
    roll_aa.set(ra, za[0], za[1], za[2]);
    gmtl::set(roll_mat, roll_aa);

    // Apply local roll
    ppose = ppose * roll_mat;

    /*____________________________________________________________________
    |
    | Compute adjusted plane coordinate system and its inverse
    |___________________________________________________________________*/

    gmtl::cross(zadj, xa, WORLD_UP);
    gmtl::normalize(zadj);

    pposeadj.set(
        xa[0], WORLD_UP[0], zadj[0], x,
        xa[1], WORLD_UP[1], zadj[1], y,
        xa[2], WORLD_UP[2], zadj[2], z,
        0, 0, 0, 1
    );
    pposeadj.setState(gmtl::Matrix44f::AFFINE);

    gmtl::invert(pposeadj_inv, pposeadj);
}

/*____________________________________________________________________
|
| Function: Update_PSystem
|
| Input:
| Output: Called every simulation step to update the particle system.
|___________________________________________________________________*/

void Update_PSystem()
{
    int i;
    gmtl::Vec3f F;             // Net force
    gmtl::Vec3f a;             // Acceleration

    for (i = 0; i < PARTICLE_NB; i++) {
        // Life is shorten by one time unit
        particles[i].ttl--;

        if (particles[i].ttl > 0) { // Still active
            // Update position
            particles[i].p += S_TSTEP * particles[i].v;

            // Compute net force from gravity and vicous drag (from wind)
            F = (particles[i].m * GRAVITY) - (K_COEF * (particles[i].v - V_WIND));

            // Calculate acceleration from froce
            a = F / particles[i].m;

            // Update velocity
            particles[i].v += S_TSTEP * a;
        }
        else {  // Died, make it reborn
            Init_AParticle(particles[i]);
        }
    }
}

/*____________________________________________________________________
|
| Function: Keyboard_Func
|
| Input:
| Output: GLUT keyboard callback function.
|___________________________________________________________________*/

void Keyboard_Func(unsigned char key, int x, int y)
{
    switch ((char)key) {

    case 'c':       // Toggle curve drawing
        render_curve = !render_curve;
        break;
    case 'x':       // Toggle constraint drawing
        render_constraint = !render_constraint;
        break;
        /*____________________________________________________________________
        |
        |   Lighting Control
        |___________________________________________________________________*/

    case '1': // Light up (+Y translation)
        light_pos[0]++;
        printf("Light-Y = %.2f\n", light_pos[1]);
        break;
    case '2': // Light down (-Y translation)
        light_pos[0]--;
        printf("Light-Y = %.2f\n", light_pos[1]);
        break;

    case '3': // Light up (+Y translation)
        light_pos[1]++;
        printf("Light-Y = %.2f\n", light_pos[1]);
        break;
    case '4': // Light down (-Y translation)
        light_pos[1]--;
        printf("Light-Y = %.2f\n", light_pos[1]);
        break;

    case '5': // Light up (+Y translation)
        light_pos[2]++;
        printf("Light-Y = %.2f\n", light_pos[2]);
        break;
    case '6': // Light down (-Y translation)
        light_pos[2]--;
        printf("Light-Y = %.2f\n", light_pos[2]);
        break;

    case '7': // Toggles diffuse light ON/OFF
        is_diffuse_on = !is_diffuse_on;
        printf("Light-diffuse = %s\n", is_diffuse_on ? "ON" : "OFF");
        break;
    case '8': // Toggles diffuse light ON/OFF
        is_ambient_on = !is_ambient_on;
        printf("Light-diffuse = %s\n", is_ambient_on ? "ON" : "OFF");
        break;

    case '9': // Toggles diffuse light ON/OFF
        is_specular_on = !is_specular_on;
        printf("Light-diffuse = %s\n", is_ambient_on ? "ON" : "OFF");
        break;


    default:
        break;
    }

    glutPostRedisplay();
}

/*____________________________________________________________________
|
| Function: Reshape_Func
|
| Input:
| Output: GLUT reshape callback function.
|___________________________________________________________________*/

void Reshape_Func(int width, int height)
{
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, ((float)width) / height, 1.0f, 10000.0f);
}

//|____________________________________________________________________
//|
//| Function: Mouse_Func
//|
//! \param button     [in] one of GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON.
//! \param state      [in] one of GLUT_UP (event is due to release) or GLUT_DOWN (press).
//! \param x          [in] X-coordinate of mouse when an event occured.
//! \param y          [in] Y-coordinate of mouse when an event occured.
//! \return None.
//!
//! GLUT mouse-callback function: called for each mouse click.
//|____________________________________________________________________

void Mouse_Func(int button, int state, int x, int y)
{
    int km_state;

    // Updates button's sate and mouse coordinates
    if (state == GLUT_DOWN) {
        mbuttons[button] = true;
        mx_prev = x;
        my_prev = y;
    }
    else {
        mbuttons[button] = false;
    }

    // Updates keyboard modifiers
    km_state = glutGetModifiers();
    kmodifiers[KM_SHIFT] = km_state & GLUT_ACTIVE_SHIFT ? true : false;
    kmodifiers[KM_CTRL] = km_state & GLUT_ACTIVE_CTRL ? true : false;
    kmodifiers[KM_ALT] = km_state & GLUT_ACTIVE_ALT ? true : false;

    //glutPostRedisplay();      // Asks GLUT to redraw the screen
}

//|____________________________________________________________________
//|
//| Function: Motion_Func
//|
//! \param x      [in] X-coordinate of mouse when an event occured.
//! \param y      [in] Y-coordinate of mouse when an event occured.
//! \return None.
//!
//! GLUT motion-callback function: called for each mouse motion.
//|____________________________________________________________________

void Motion_Func(int x, int y)
{
    int dx, dy, d;

    if (mbuttons[GLUT_LEFT_BUTTON] || mbuttons[GLUT_RIGHT_BUTTON]) {
        // Computes distances the mouse has moved
        dx = x - mx_prev;
        dy = y - my_prev;

        // Updates mouse coordinates
        mx_prev = x;
        my_prev = y;

        // Hold left button to rotate camera
        if (mbuttons[GLUT_LEFT_BUTTON]) {
            if (!kmodifiers[KM_CTRL]) {
                elevation += dy;            // Elevation update
            }
            if (!kmodifiers[KM_SHIFT]) {
                azimuth += dx;              // Azimuth update
            }
        }

        // Hold right button to zoom
        if (mbuttons[GLUT_RIGHT_BUTTON]) {
            if (abs(dx) >= abs(dy)) {
                d = dx;
            }
            else {
                d = -dy;
            }
            distance += d;
        }

        glutPostRedisplay();      // Asks GLUT to redraw the screen
    }
}

/*____________________________________________________________________
|
| Function: Draw_World_Axes
|
| Input:
| Output: Draw world axes.
|___________________________________________________________________*/

void Draw_World_Axes()
{
    glBegin(GL_LINES);
    // x-axis is red
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(100.0, 0.0, 0.0);
    // y-axis is green
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 100.0, 0.0);
    // z-axis is blue
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 100.0);
    glEnd();
}

/*____________________________________________________________________
|
| Function: Draw_Path
|
| Input:
| Output: Draws a path.
|___________________________________________________________________*/

void Draw_Path()
{
    int i;
    float x, y, z, t;

    /*____________________________________________________________________
    |
    | Draw input points and tangents
    |___________________________________________________________________*/

    if (render_constraint) {
        for (i = 0; i < PT_NB; i++) {
            glPushMatrix();
            glTranslatef(input_pts[i][0], input_pts[i][1], input_pts[i][2]);
            glColor3f(1, 0, 0);
            glutSolidSphere(5, 16, 16);
            glColor3f(1, 1, 0);
            glBegin(GL_LINES);
            glVertex3f(0, 0, 0);
            glVertex3f(tangents[i][0], tangents[i][1], tangents[i][2]);
            glEnd();
            glPopMatrix();
        }
    }

    /*____________________________________________________________________
    |
    | Draw Hermite curve segments using line strips
    |___________________________________________________________________*/

    if (render_curve) {
        glColor3f(0, 1, 0);
        for (i = 0; i < PT_NB; i++) {
            // Draw each segment
            glBegin(GL_LINE_STRIP);
            for (t = 0; t <= 1; t += C_TSTEP) {
                // Simple polynomial evaluation
                /*float t2 = t*t;
                float t3 = t2*t;
                x = Cmats[i][0][0]*t3 + Cmats[i][1][0]*t2 + Cmats[i][2][0]*t + Cmats[i][3][0];
                y = Cmats[i][0][1]*t3 + Cmats[i][1][1]*t2 + Cmats[i][2][1]*t + Cmats[i][3][1];
                z = Cmats[i][0][2]*t3 + Cmats[i][1][2]*t2 + Cmats[i][2][2]*t + Cmats[i][3][2];*/

                // Use Horner's rule for slight speedup
                x = ((Cmats[i][0][0] * t + Cmats[i][1][0]) * t + Cmats[i][2][0]) * t + Cmats[i][3][0];
                y = ((Cmats[i][0][1] * t + Cmats[i][1][1]) * t + Cmats[i][2][1]) * t + Cmats[i][3][1];
                z = ((Cmats[i][0][2] * t + Cmats[i][1][2]) * t + Cmats[i][2][2]) * t + Cmats[i][3][2];
                glVertex3f(x, y, z);
            }
            glEnd();
        }
    }
}

void Draw_Path2()
{
    int i;
    float x, y, z, t;

    /*____________________________________________________________________
    |
    | Draw input points and tangents
    |___________________________________________________________________*/

    if (render_constraint) {
        for (i = 0; i < PT_NB; i++) {
            glPushMatrix();
            glTranslatef(input_pts[i][0], input_pts[i][1], input_pts[i][2]);
            glColor3f(1, 0, 0);
            glutSolidSphere(5, 16, 16);
            glColor3f(1, 1, 0);
            glBegin(GL_LINES);
            glVertex3f(0, 0, 0);
            glVertex3f(tangents[i][0], tangents[i][1], tangents[i][2]);
            glEnd();
            glPopMatrix();
        }
    }

    /*____________________________________________________________________
    |
    | Draw Hermite curve segments using line strips
    |___________________________________________________________________*/

    if (render_curve) {
        glColor3f(0, 1, 0);
        for (i = 0; i < PT_NB; i++) {
            // Draw each segment
            glBegin(GL_LINE_STRIP);
            for (t = 0; t <= 1; t += C_TSTEP) {
                // Simple polynomial evaluation
                float t2 = t * t;
                float t3 = t2 * t;
                x = Cmats[i][0][0] * t3 + Cmats[i][1][0] * t2 + Cmats[i][2][0] * t + Cmats[i][3][0];
                y = Cmats[i][0][1] * t3 + Cmats[i][1][1] * t2 + Cmats[i][2][1] * t + Cmats[i][3][1];
                z = Cmats[i][0][2] * t3 + Cmats[i][1][2] * t2 + Cmats[i][2][2] * t + Cmats[i][3][2];

                // Use Horner's rule for slight speedup
                /*x = ((Cmats[i][0][0] * t + Cmats[i][1][0]) * t + Cmats[i][2][0]) * t + Cmats[i][3][0];
                y = ((Cmats[i][0][1] * t + Cmats[i][1][1]) * t + Cmats[i][2][1]) * t + Cmats[i][3][1];
                z = ((Cmats[i][0][2] * t + Cmats[i][1][2]) * t + Cmats[i][2][2]) * t + Cmats[i][3][2];
                glVertex3f(x, y, z);*/
            }
            glEnd();
        }
    }
}


void drawLight(const float radius) {
    glDisable(GL_LIGHTING);
    glColor3f(1.0f, 1.0f, 1.0f);
    glutSolidSphere(radius, 20, 20);
    glEnable(GL_LIGHTING);
}



void DrawPlaneBody(const float width, const float length, const float height)
{
    
    glPushMatrix();
    glMultMatrixf(ppose.mData);
    glScalef(8, 8, 8); // If lighting enabled here, be sure to use GL_RESCALE_NORMAL
    
        float w = width/2;
        float l = length/2;
        float h = height / 2;
        
        // Sets materials
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 20.0);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SPECULAR_COL);
        
        glBegin(GL_TRIANGLES);
           // Wings are red
         glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, BRIGHTRED_COL);
         glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, BRIGHTRED_COL);
          
           //left wing
          glNormal3f(0.0, 1.0, 0.2);
           glVertex3f(w-1 , h, l - 2);
           glVertex3f(w-1, h, l-4);
           glVertex3f(w, h, l-4);
          
          // right wing
          glNormal3f(0.0, 1.0, 0.2);
          glVertex3f(-w+1, h, l - 2);
          glVertex3f(-w+1, h, l - 4);
          glVertex3f(-w, h, l - 4);
          
          //fin
          glNormal3f(0.0, 1.0, 0.2);
          glVertex3f(0, h, l - 3);
          glVertex3f(0, h, l - 4.5);
          glVertex3f(0, h+1, l - 4.5);
        glEnd();
          
             GLfloat x              = 0.0;
             GLfloat y              = 0.0;
             GLfloat angle          = 0.0;
             GLfloat angle_stepsize = 0.1;
             GLfloat radius = 1;

              // Sets materials
              glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 20.0);
              glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SPECULAR_COL);
        
          //body
              glBegin(GL_QUAD_STRIP);
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, BRIGHTBLUE_COL);
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, BRIGHTBLUE2_COL);
              angle = 0.0;
          while( angle < 2*gmtl::Math::PI ) {
                      x = radius * cos(angle);
                      y = radius * sin(angle);
                      glNormal3f(0.0, 1.0, 0.2);
                      glVertex3f(x, y , h-4);
                      glVertex3f(x, y , 0.0);
                      angle = angle + angle_stepsize;
                  }
                   glNormal3f(0.0, 1.0, 0.2);
                  glVertex3f(radius, 0.0, h);
                  glVertex3f(radius, 0.0, 0.0);
              glEnd();

      // Sets materials
              glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 20.0);
              glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SPECULAR_COL);
        
           glBegin(GL_POLYGON);
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, BRIGHTBLUE_COL);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, BRIGHTBLUE2_COL);
           angle = 0.0;
          while( angle < 2*gmtl::Math::PI ) {
                   x = radius * cos(angle);
                   y = radius * sin(angle);
                   glNormal3f(0.0, 1.0, 0.2);
                   glVertex3f(x, y , h-4);
                   angle = angle + angle_stepsize;
               }
                glNormal3f(0.0, 1.0, 0.2);
               glVertex3f(radius, 0.0, h);
           glEnd();
          
          // Sets materials
                  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 20.0);
                  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SPECULAR_COL);
          //front
          int k;
          glBegin(GL_TRIANGLES);
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, BRIGHTGREEN_COL);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, BRIGHTGREEN_COL);
          for (k=0;k<=360;k+=DEF_D){
            glNormal3f(0.0, 1.0, 0.2);
            glVertex3f(0,0,1);
            glVertex3f(cos(k),sin(k),0);
            glVertex3f(cos(k+DEF_D),sin(k+DEF_D),0);
          }
          glEnd();
}

/*____________________________________________________________________
|
| Function: Draw_Particles
|
| Input:
| Output: Draw particles as view-oriented billboards.
|___________________________________________________________________*/

void Draw_Particles()
{
    int i;
    gmtl::Matrix44f cam_mat;                    // Camera matrix
    gmtl::Matrix44f dt_mat, el_mat, az_mat;     // distance, elevation, and azimuth matrices, initialized to IDENTITY
    float age_scale;                            // Age factor

  /*____________________________________________________________________
  |
  | Enable texturing and blending
  |___________________________________________________________________*/

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, PARTICLE1);
    glDepthMask(GL_FALSE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);


        // Set distance matrix
        dt_mat[2][3] = distance;

    // Set elevation matrix
    gmtl::set(el_mat,
        gmtl::AxisAnglef(gmtl::Math::deg2Rad(elevation), 1, 0, 0)
    );

    // Set azimuth matrix
    gmtl::set(az_mat,
        gmtl::AxisAnglef(gmtl::Math::deg2Rad(azimuth), 0, 1, 0)
    );

    // Compute camera w.r.t. world and discard translation
    cam_mat = pposeadj * az_mat * el_mat * dt_mat;
    cam_mat[0][3] = cam_mat[1][3] = cam_mat[2][3] = 0;

    /*____________________________________________________________________
    |
    | Render particles as billboards
    |___________________________________________________________________*/

    for (i = 0; i < PARTICLE_NB; i++) {
        glPushMatrix();
        glTranslatef(particles[i].p[0], particles[i].p[1], particles[i].p[2]);
        glMultMatrixf(cam_mat.mData);       // Orient billboards to face camera
        glBegin(GL_QUADS);
        age_scale = ((float)particles[i].ttl) / particles[i].ttl_init;
        glColor3f(age_scale, age_scale, age_scale);
        glTexCoord2f(0.0, 0.0);
        glVertex3f(-SMOKE_SIZE, -SMOKE_SIZE, 0.0);
        glTexCoord2f(1.0, 0.0);
        glVertex3f(SMOKE_SIZE, -SMOKE_SIZE, 0.0);
        glTexCoord2f(1.0, 1.0);
        glVertex3f(SMOKE_SIZE, SMOKE_SIZE, 0.0);
        glTexCoord2f(0.0, 1.0);
        glVertex3f(-SMOKE_SIZE, SMOKE_SIZE, 0.0);
        glEnd();
        glPopMatrix();
    }

    /*____________________________________________________________________
    |
    | Restore rendering states
    |___________________________________________________________________*/

    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    glDisable(GL_TEXTURE_2D);
}

void Draw_Particles2()
{
    int i;
    gmtl::Matrix44f cam_mat;                    // Camera matrix
    gmtl::Matrix44f dt_mat, el_mat, az_mat;     // distance, elevation, and azimuth matrices, initialized to IDENTITY
    float age_scale;                            // Age factor

  /*____________________________________________________________________
  |
  | Enable texturing and blending
  |___________________________________________________________________*/

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, PARTICLE2);
    glDepthMask(GL_FALSE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    /*____________________________________________________________________
    |
    | Orient billboards to face the camera
    |___________________________________________________________________*/

    // Set distance matrix
    dt_mat[2][3] = distance;

    // Set elevation matrix
    gmtl::set(el_mat,
        gmtl::AxisAnglef(gmtl::Math::deg2Rad(elevation), 1, 0, 0)
    );

    // Set azimuth matrix
    gmtl::set(az_mat,
        gmtl::AxisAnglef(gmtl::Math::deg2Rad(azimuth), 0, 1, 0)
    );

    // Compute camera w.r.t. world and discard translation
    cam_mat = pposeadj * az_mat * el_mat * dt_mat;
    cam_mat[0][3] = cam_mat[1][3] = cam_mat[2][3] = 0;

    /*____________________________________________________________________
    |
    | Render particles as billboards
    |___________________________________________________________________*/

    for (i = 0; i < PARTICLE_NB; i++) {
        glPushMatrix();
        glTranslatef(particles[i].p[0], particles[i].p[1] + 45, particles[i].p[2]);
        glMultMatrixf(cam_mat.mData);       // Orient billboards to face camera
        glBegin(GL_QUADS);
        age_scale = ((float)particles[i].ttl) / particles[i].ttl_init;
        glColor3f(age_scale, age_scale, age_scale);
        glTexCoord2f(0.0, 0.0);
        glVertex3f(-SMOKE_SIZE, -SMOKE_SIZE, 0.0);
        glTexCoord2f(1.0, 0.0);
        glVertex3f(SMOKE_SIZE, -SMOKE_SIZE, 0.0);
        glTexCoord2f(1.0, 1.0);
        glVertex3f(SMOKE_SIZE, SMOKE_SIZE, 0.0);
        glTexCoord2f(0.0, 1.0);
        glVertex3f(-SMOKE_SIZE, SMOKE_SIZE, 0.0);
        glEnd();
        glPopMatrix();
    }

    /*____________________________________________________________________
    |
    | Restore rendering states
    |___________________________________________________________________*/

    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    glDisable(GL_TEXTURE_2D);
}

/*
    Random number generator by Donald H. House.
    Modified to compile on Visual Studio.Net 2003
*/
float FastGauss(float mean, float std)
{
#define RESOLUTION 2500
    static float lookup[RESOLUTION + 1];

#define itblmax    20
    /* length - 1 of table describing F inverse */
#define didu    40.0
/* delta table position / delta ind. variable           itblmax / 0.5 */

    static float tbl[] =
    { 0.00000E+00, 6.27500E-02, 1.25641E-01, 1.89000E-01,
     2.53333E-01, 3.18684E-01, 3.85405E-01, 4.53889E-01,
     5.24412E-01, 5.97647E-01, 6.74375E-01, 7.55333E-01,
     8.41482E-01, 9.34615E-01, 1.03652E+00, 1.15048E+00,
     1.28167E+00, 1.43933E+00, 1.64500E+00, 1.96000E+00,
     3.87000E+00 };

    static int hereb4;
    //static struct timeval tv;

    float u, di, delta/*, result*/;
    int i, index, minus = 1;

    if (!hereb4) {
        for (i = 0; i <= RESOLUTION; i++) {
            if ((u = i / (float)RESOLUTION) > 0.5) {
                minus = 0;
                u -= 0.5;
            }

            /* interpolate gaussian random number using table */

            index = (int)(di = (didu * u));
            di -= (float)index;
            delta = tbl[index] + (tbl[index + 1] - tbl[index]) * di;
            lookup[i] = (minus ? -delta : delta);
        }

        /*gettimeofday(&tv, NULL);
        srand((unsigned int)tv.tv_usec);*/
        srand((unsigned)time(NULL));
        hereb4 = 1;
    }

    i = rand() / (RAND_MAX / RESOLUTION) + 1;
    if (i > RESOLUTION) {
        i = RESOLUTION;
    }

    return(mean + std * lookup[i]);
}

/*____________________________________________________________________
|
| Function: LoadPPM
|
| Input:
| Output: Draw world axes.
|  LoadPPM - a minimal Portable Pixelformat image file loader
|  fname: name of file to load (input)
|  w: width of loaded image in pixels (output)
|  h: height of loaded image in pixels (output)
|  data: image data address (input or output depending on mallocflag)
|  mallocflag: 1 if memory not pre-allocated, 0 if data already points
|              to allocated memory that can hold the image. Note that
|              if new memory is allocated, free() should be used to
|              deallocate when it is no longer needed.
|___________________________________________________________________*/

void LoadPPM(const char* fname, unsigned int* w, unsigned int* h, unsigned char** data, const int mallocflag)
{
    FILE* fp;
    char P, num;
    int max;
    char s[1000];

    if (!(fp = fopen(fname, "rb")))
    {
        perror("cannot open image file\n");
        exit(0);
    }

    fscanf(fp, "%c%c\n", &P, &num);
    if ((P != 'P') || (num != '6'))
    {
        perror("unknown file format for image\n");
        exit(0);
    }

    do
    {
        fgets(s, 999, fp);
    } while (s[0] == '#');


    sscanf(s, "%d%d", w, h);
    fgets(s, 999, fp);
    sscanf(s, "%d", &max);

    if (mallocflag)
        if (!(*data = (unsigned char*)malloc(*w * *h * 3)))
        {
            perror("cannot allocate memory for image data\n");
            exit(0);
        }

    fread(*data, 3, *w * *h, fp);

    fclose(fp);
}

/*____________________________________________________________________
|
| Function: DrawSkybox
|
| ! \param s      [in] Skybox size.
| ! \return None.
| !
| ! Draws a skybox.
|____________________________________________________________________*/

void DrawSkybox(const float s)
{
    float s2 = s / 2;

       // Turn on texture mapping and disable lighting
       glEnable(GL_TEXTURE_2D);
       glDisable(GL_LIGHTING);

       // Back wall
       glBindTexture(GL_TEXTURE_2D, textures[TID_SKYBACK]);  // Specify which texture will be used
       glBegin(GL_QUADS);
       glColor3f(1.0f, 1.0f, 1.0f);
       glTexCoord2f(0.0, 1.0);
       glVertex3f(-s2, -s2, -s2);
       glTexCoord2f(1.0, 1.0);
       glVertex3f(s2, -s2, -s2);
       glTexCoord2f(1.0, 0.0);
       glVertex3f(s2, s2, -s2);
       glTexCoord2f(0.0, 0.0);
       glVertex3f(-s2, s2, -s2);
       glEnd();

       // Left wall
       glBindTexture(GL_TEXTURE_2D, textures[TID_SKYLEFT]);
       glBegin(GL_QUADS);
       glColor3f(1.0f, 1.0f, 1.0f);
       glTexCoord2f(0.0, 1.0);
       glVertex3f(-s2, -s2, s2);
       glTexCoord2f(1.0, 1.0);
       glVertex3f(-s2, -s2, -s2);
       glTexCoord2f(1.0, 0.0);
       glVertex3f(-s2, s2, -s2);
       glTexCoord2f(0.0, 0.0);
       glVertex3f(-s2, s2, s2);
       glEnd();

       // Bottom wall
       glBindTexture(GL_TEXTURE_2D, textures[TID_SKYBOTTOM]);
       glBegin(GL_QUADS);
       glColor3f(1.0f, 1.0f, 1.0f);
       glTexCoord2f(0.0, 1.0);
       glVertex3f(-s2, -s2, s2);
       glTexCoord2f(1.0, 1.0);
       glVertex3f(s2, -s2, s2);
       glTexCoord2f(1.0, 0.0);
       glVertex3f(s2, -s2, -s2);
       glTexCoord2f(0.0, 0.0);
       glVertex3f(-s2, -s2, -s2);
       glEnd();

       // Front wall
       glBindTexture(GL_TEXTURE_2D, textures[TID_SKYFRONT]);  // Specify which texture will be used
       glBegin(GL_QUADS);
       glColor3f(1.0f, 1.0f, 1.0f);
       glTexCoord2f(0.0, 1.0);
       glVertex3f(-s2, -s2, s2);
       glTexCoord2f(1.0, 1.0);
       glVertex3f(s2, -s2, s2);
       glTexCoord2f(1.0, 0.0);
       glVertex3f(s2, s2, s2);
       glTexCoord2f(0.0, 0.0);
       glVertex3f(-s2, s2, s2);
       glEnd();

       // Right wall
       glBindTexture(GL_TEXTURE_2D, textures[TID_SKYRIGHT]);
       glBegin(GL_QUADS);
       glColor3f(1.0f, 1.0f, 1.0f);
       glTexCoord2f(0.0, 1.0);
       glVertex3f(s2, -s2, s2);
       glTexCoord2f(1.0, 1.0);
       glVertex3f(s2, -s2, -s2);
       glTexCoord2f(1.0, 0.0);
       glVertex3f(s2, s2, -s2);
       glTexCoord2f(0.0, 0.0);
       glVertex3f(s2, s2, s2);
       glEnd();

       // Top wall
       glBindTexture(GL_TEXTURE_2D, textures[TID_SKYTOP]);
       glBegin(GL_QUADS);
       glColor3f(1.0f, 1.0f, 1.0f);
       glTexCoord2f(0.0, 1.0);
       glVertex3f(-s2, s2, s2);
       glTexCoord2f(1.0, 1.0);
       glVertex3f(s2, s2, s2);
       glTexCoord2f(1.0, 0.0);
       glVertex3f(s2, s2, -s2);
       glTexCoord2f(0.0, 0.0);
       glVertex3f(-s2, s2, -s2);
       glEnd();

      // The building Top
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_TOPBUILD]);
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(-s2 / 50, s2 / 50 - 11, s2 / 50);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(s2 / 50, s2 / 50 - 11, s2 / 50);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(s2 / 50, s2 / 50 - 11, -s2 / 50);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(-s2 / 50, s2 / 50 - 11, -s2 / 50);
           glEnd();
       
           // The building side1
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURN_BUILD1]);  // Specify which texture will be used
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(-s2 / 50, -s2, -s2 / 50);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(s2 / 50, -s2, -s2 / 50);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(s2 / 50, s2 / 50 - 11, -s2 / 50);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(-s2 / 50, s2 / 50 - 11, -s2 / 50);
           glEnd();
       
           // The building side2
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURN_BUILD1]);  // Specify which texture will be used
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(s2 / 50, -s2, -s2 / 50);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(s2 / 50, -s2, s2 / 50);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(s2 / 50, s2 / 50 - 11, s2 / 50);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(s2 / 50, s2 / 50 - 11, -s2 / 50);
           glEnd();
       
           // The building side3
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURN_BUILD1]);  // Specify which texture will be used
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(s2 / 50, -s2, s2 / 50);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(-s2 / 50, -s2, s2 / 50);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(-s2 / 50, s2 / 50 - 11, s2 / 50);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(s2 / 50, s2 / 50 - 11, s2 / 50);
           glEnd();
       
           // The building side4
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURN_BUILD1]);  // Specify which texture will be used
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(-s2 / 50, -s2, s2 / 50);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(-s2 / 50, -s2, -s2 / 50);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(-s2 / 50, s2 / 50 - 11, -s2 / 50);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(-s2 / 50, s2 / 50 - 11, s2 / 50);
           glEnd();
       
           // The building2 Top
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_TOPBUILD]);
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(-s2 / 30, s2 / 70 - 20, s2 / 30 + 50);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(s2 / 30, s2 / 70 - 20, s2 / 30 + 50);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(s2 / 30, s2 / 70 - 20, s2 / 30 + 25);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(-s2 / 30, s2 / 70 - 20, s2 / 30 + 25);
           glEnd();
           // The buildin2 side1
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_BUILD2]);  // Specify which texture will be used
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(-s2 / 30, -s2, s2 / 30 +25);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(s2 / 30, -s2, s2 / 30 +25);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(s2 / 30, s2 / 70-20, s2 / 30 + 25);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(-s2 / 30, s2 / 70-20, s2 / 30 + 25);
           glEnd();
       
           // The building2 side2
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_BUILD2]);  // Specify which texture will be used
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(-s2 / 30, -s2, s2 / 30 + 50);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(s2 / 30, -s2, s2 / 30 + 50);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(s2 / 30, s2 / 70 - 20, s2 / 30 + 50);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(-s2 / 30, s2 / 70 - 20, s2 / 30 + 50);
           glEnd();
       
           // The building2 side3
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_BUILD2]);  // Specify which texture will be used
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(-s2 / 30, -s2, s2 / 30 + 50);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(-s2 / 30, -s2, s2 / 30 + 25);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(-s2 / 30, s2 / 70 - 20, s2 / 30 + 25);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(-s2 / 30, s2 / 70 - 20, s2 / 30 + 50);
           glEnd();
       
           // The building2 side4
           glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_BUILD2]);  // Specify which texture will be used
           glBegin(GL_QUADS);
           glColor3f(1.0f, 1.0f, 1.0f);
           glTexCoord2f(0.0, 1.0);
           glVertex3f(s2 / 30, -s2, s2 / 30 + 50);
           glTexCoord2f(1.0, 1.0);
           glVertex3f(s2 / 30, -s2, s2 / 30 + 25);
           glTexCoord2f(1.0, 0.0);
           glVertex3f(s2 / 30, s2 / 70 - 20, s2 / 30 + 25);
           glTexCoord2f(0.0, 0.0);
           glVertex3f(s2 / 30, s2 / 70 - 20, s2 / 30 + 50);
           glEnd();

       // Turn off texture mapping and enable lighting
       glEnable(GL_LIGHTING);
       glDisable(GL_TEXTURE_2D);
}

/*____________________________________________________________________
|
| Function: SetLight
|
| ! \param pos          [in] Light position.
| ! \param is_ambient   [in] Is ambient enabled?
| ! \param is_diffuse   [in] Is diffuse enabled?
| ! \param is_specular  [in] Is specular enabled?
| ! \return None.
| !
| ! Set light properties.
|____________________________________________________________________*/

void SetLight(const gmtl::Point4f& pos, const bool is_ambient, const bool is_diffuse, const bool is_specular)
{
    glLightfv(GL_LIGHT0, GL_POSITION, pos.mData);
    if (is_ambient) {
        glLightfv(GL_LIGHT0, GL_AMBIENT, AMBIENT_LIGHT);
    }
    else {
        glLightfv(GL_LIGHT0, GL_AMBIENT, NO_LIGHT);
    }
    if (is_diffuse) {
        glLightfv(GL_LIGHT0, GL_DIFFUSE, DIFFUSE_LIGHT);
    }
    else {
        glLightfv(GL_LIGHT0, GL_DIFFUSE, NO_LIGHT);
    }
    if (is_specular) {
        glLightfv(GL_LIGHT0, GL_SPECULAR, SPECULAR_LIGHT);
    }
    else {
        glLightfv(GL_LIGHT0, GL_SPECULAR, NO_LIGHT);
    }
}


