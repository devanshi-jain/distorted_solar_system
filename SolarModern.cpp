/*
 * SolarModern.cpp - Version 1.3 - October 10, 2020
 *
 * Example program illustrating a simple solar system with a sun,
 *    earth, and moon rendered as wireframe spheres.
 * 
 * USAGE:
 *    Press "r" or "R" key to toggle (off and on) running the animation
 *    Press "s" or "S" key to single-step animation
 *    Press "t" or "T" key to toggle between animating at a fixed time step
 *                                     or at real elapsed time.
 *    Press "c" or "C" to toggle culling backfaces
 *    The up arrow key and down array key control the
 *			time step used in the animation rate.  Each key
 *			press multiplies or divides the times by a factor 
 *			of two (2).
 *	  Press ESCAPE to exit.
 *
 * Author: Sam Buss
 *
 * Software accompanying POSSIBLE SECOND EDITION TO the book
 *		3D Computer Graphics: A Mathematical Introduction with OpenGL,
 *		by S. Buss, Cambridge University Press, 2003.
 *
 * Software is "as-is" and carries no warranty.  It may be used without
 *   restriction, but if you modify it, please change the filenames to
 *   prevent confusion between different versions.
 * Bug reports: Sam Buss, sbuss@ucsd.edu.
 * Web page: http://math.ucsd.edu/~sbuss/MathCG2
 */

// These libraries are needed to link the program.
#pragma comment(lib,"opengl32.lib")
#pragma comment(lib,"glu32.lib")
#pragma comment(lib,"glfw3.lib")
#pragma comment(lib,"glew32s.lib")
#pragma comment(lib,"glew32.lib")

// Use the static library (so glew32.dll is not needed):
#define GLEW_STATIC

#include "GL/glew.h" 
#include "GLFW/glfw3.h"

#include "MathMisc.h"
#include "LinearR3.h"
#include "LinearR4.h"		
#include "GlGeomSphere.h"
#include "ShaderMgrSLR.h"
bool check_for_opengl_errors();     // Function prototype (should really go in a header file)

// Enable standard input and output via printf(), etc.
// Put this include *after* the includes for glew and GLFW!
#include <stdio.h>

// ********************
// Animation controls and state infornation
// ********************

// These two variables control whether running or paused.
bool spinMode = true;
bool singleStep = false;
bool cullBackFaces = false; // Equals true to cull backfaces. Equals false to not cull backfaces.   
bool UseRealTime = false;   // Initially use a fixed animation increment step.
double PreviousTime = 0.0;

// These three variables control the animation's state and speed.
double HourOfDay = 0.0;
double DayOfYear = 0.0;
double DayOfYearX = 0.0;
double AnimateIncrement = 24.0;  // Time step for animation (in units of hours)

double viewAzimuth = 0.25;	// Angle of view up/down (in radians)
LinearMapR4 viewMatrix;		// The current view matrix, based on viewAzimuth and viewDirection.

// ************************
// General data helping with setting up VAO (Vertex Array Objects)
//    and Vertex Buffer Objects.
// ***********************

// These objects take care of generating and loading VAO's, VBO's and EBO's,
//    rendering spheres for the moon, earch and sun.
// They render as radius 1 spheres (but will be scaled by the Model matrix)
GlGeomSphere Moon1(6, 6);    // A sphere with 6 slices and 6 stacks
GlGeomSphere Earth(8, 12);    // A sphere with 8 slices and 12 stacks
GlGeomSphere PlanetX(9, 11);    // A sphere with 9 slices and 11 stacks
GlGeomSphere Sun1(10, 10);    // A sphere with 10 slices and 10 stacks
GlGeomSphere Sun2(10, 10);    // A sphere with 10 slices and 10 stacks
GlGeomSphere SD(4,8);    // A sphere with 10 slices and 10 stacks
GlGeomSphere moonmoon(4, 6);    // A sphere with 10 slices and 10 stacks


// We create one shader program: consisting of a vertex shader and a fragment shader
unsigned int shaderProgram1;
const unsigned int vertPos_loc = 0;   // Corresponds to "location = 0" in the verter shader definitions
const unsigned int vertColor_loc = 1; // Corresponds to "location = 1" in the verter shader definitions
const char* projMatName = "projectionMatrix";		// Name of the uniform variable projectionMatrix
int projMatLocation;						// Location of the projectionMatrix in the "smooth" shader program.
const char* modelviewMatName = "modelviewMatrix";	// Name of the uniform variable modelviewMatrix
int modelviewMatLocation;					// Location of the modelviewMatrix in the "smooth" shader program.

//  The Projection matrix: Controls the "camera view/field-of-view" transformation
//     Generally is the same for all objects in the scene.
LinearMapR4 theProjectionMatrix;		//  The Projection matrix: Controls the "camera/view" transformation

// A ModelView matrix controls the placement of a particular object in 3-space.
//     It is generally different for each object.
// The array matEntries holds the matrix values as floats to be loaded into the shader program. 
float matEntries[16];		// Holds 16 floats (since cannot load doubles into a shader that uses floats)

// *****************************
// These variables set the dimensions of the perspective region we wish to view.
// They are used to help form the projection matrix and the view matrix
// All rendered objects lie in the rectangular prism centered on the z-axis
//     equal to (-Xmax,Xmax)x(-Ymax,Ymax)x(Zmin,Zmax)
// The model/view matrix can be used to move objects to this position
// Increase SceneXZradius if the solar system becomes bigger
const double SceneXZradius = 7.5;       
const double Xmax = SceneXZradius;                
const double Ymax = 2.0;
const double Zmax = SceneXZradius;
const double Zmin = -SceneXZradius;
// CameraDistance gives distance from the camera to the origin (usually center of the scene)
//      This controls the perpective of the camera, so that that at the origin,
//          all of (-Xmax,Xmax)x(-Ymax,Ymax)x{-CameraDistance} is visible.
//      MakeCameraDistance larger or smaller to affect field of view.
// ZnearMin = minimum value for distance to near clipping plane
const double CameraDistance = 30.0;   
const double ZnearMin = 1.0;

// *************************
// mySetupGeometries defines the scene data, especially vertex  positions and colors.
//    - It also loads all the data into the VAO's (Vertex Array Objects) and
//      into the VBO's (Vertex Buffer Objects).
// This routine is only called once to initialize the data.
// *************************
void mySetupGeometries() {

	// These routines take care of loading info into their VAO's, VBO's and EBO's.
	Sun1.InitializeAttribLocations(vertPos_loc);
    Sun2.InitializeAttribLocations(vertPos_loc);
	Earth.InitializeAttribLocations(vertPos_loc);
    PlanetX.InitializeAttribLocations(vertPos_loc);
	Moon1.InitializeAttribLocations(vertPos_loc);
    SD.InitializeAttribLocations(vertPos_loc);
    moonmoon.InitializeAttribLocations(vertPos_loc);

	// Set the initial view matrix
	viewMatrix.Set_glTranslate(0.0, 0.0, -CameraDistance);        // Translate to be in front of the camera
    viewMatrix.Mult_glRotate(viewAzimuth, 1.0, 0.0, 0.0);         // Rotate to view from slightly above   

	check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}
// *************************************
// Main routine for rendering the scene
// myRenderScene() is called every time the scene needs to be redrawn.
// mySetupGeometries() has already created the vertex and buffer objects
//    and the model view matrices.
// setup_shaders() has already created the shader programs.
// *************************************
void myRenderScene() {

	// Clear the rendering window
	static const float black[] = { 0.0f, 0.0f, 0.0f, 0.0f };
	const float clearDepth = 1.0f;
	glClearBufferfv(GL_COLOR, 0, black);
	glClearBufferfv(GL_DEPTH, 0, &clearDepth);	// Must pass in a pointer to the depth value!

    if (spinMode) {
        // Calculate the time step, either as AnimateIncrement or based on actual elapsed time
        double thisAnimateIncrement = AnimateIncrement;
        if (UseRealTime && !singleStep) {
            double curTime = glfwGetTime();
            thisAnimateIncrement *= (curTime - PreviousTime)*60.0;
            PreviousTime = curTime;
        }
        // Update the animation state
        HourOfDay += thisAnimateIncrement;
        DayOfYear += thisAnimateIncrement / 24.0;
        DayOfYearX += thisAnimateIncrement / 24.0;

        HourOfDay = HourOfDay - ((int)(HourOfDay / 24)) * 24;       // Wrap back to be in range [0,24)
        DayOfYear = DayOfYear - ((int)(DayOfYear / 365)) * 366;     // Wrap back to be in range [0,365)
        DayOfYearX = DayOfYearX - ((int)(DayOfYearX / 365)) * 500;     // Wrap back to be in range [0,500)

        if (singleStep) {
            spinMode = false;       // If in single step mode, turn off future animation
        }
    }


	glUseProgram(shaderProgram1);

	LinearMapR4 Sun1PosMatrix = viewMatrix;		                    // Place Sun1 at center of the scene
    double revolveAnglesun1 = (DayOfYear / 365.0) * PI2;
    Sun1PosMatrix.Mult_glRotate(revolveAnglesun1, 0.0, 1.0, 0.0);   // Revolve the sun1 around sun2
    Sun1PosMatrix.Mult_glTranslate(0.0, 0.0, 2.0);		            // Place the earth five units away from the sun
    Sun1PosMatrix.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    glVertexAttrib3f(vertColor_loc, 1.0f, 1.0f, 0.0f);     // Make the sun yellow
    Sun1.Render();

    LinearMapR4 Sun2PosMatrix = viewMatrix;			// Place Sun2 at center of the scene
    // SunPosMatrix.Mult_glScale(1.0);				// Scaling by (1, 1, 1) has no effect
    double revolveAnglesun2 = (DayOfYear / 365.0) * PI2;
    Sun2PosMatrix.Mult_glRotate(revolveAnglesun2, 0.0, 1.0, 0.0);   // Revolve the sun1 around sun2
    Sun2PosMatrix.Mult_glTranslate(0.0, 0.0, -2.0);		        // Place the earth five units away from the sun
    // These two lines load the matrix into the shader
    Sun2PosMatrix.DumpByColumns(matEntries);           // These two lines load the matrix into the shader
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    glVertexAttrib3f(vertColor_loc, 1.0f, 1.0f, 0.0f);     // Make the sun yellow
    Sun2.Render();

    
    // EarthPosMatrix - specifies position of the earth
    // EarthMatrix - specifies the size of the earth and its rotation on its axis
    LinearMapR4 EarthPosMatrix = viewMatrix;//SunPosMatrix;
    double revolveAngle = -1 * (DayOfYear / 365.0) * PI2;//0.4 //
    EarthPosMatrix.Mult_glTranslate(5.0*cos(revolveAngle), 0.0, 5.0*sin(revolveAngle));		// Place the earth five units away from the sun
    double tiltAngle = 0.41; //radians
    EarthPosMatrix.Mult_glRotate(-1*tiltAngle, 0.0, 0.0, 1.0); //IMP (Z)  // Revolve the earth around the sun

	LinearMapR4 EarthMatrix = EarthPosMatrix;
    double earthRotationAngle =  (HourOfDay / 24.0) * PI2;
    EarthMatrix.Mult_glRotate(earthRotationAngle, 0.0, 1.0, 0.0);   // Rotate earth on y-axis
	EarthMatrix.Mult_glScale(0.5);                                  // Make radius 0.5.
	EarthMatrix.DumpByColumns(matEntries);
	glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
	glVertexAttrib3f(vertColor_loc, 0.2f, 0.4f, 1.0f);     // Make the earth bright cyan-blue
	Earth.Render();

    // MoonMatrix - control placement, and size of the moon.
 	LinearMapR4 MoonMatrix = EarthPosMatrix;        // Base the moon's matrix off the earth's *POS* matrix (EarthPosMatrix)
    double moonRotationAngle = (DayOfYear*12.0 / 365.0)*PI2;
    MoonMatrix.Mult_glRotate(moonRotationAngle, 0.0, 1.0, 0.0);  // Revolving around the earth twelve times per year
	MoonMatrix.Mult_glTranslate(0.0, 0.0, 1.0);	    // Place the Moon one unit away from the earth
    MoonMatrix.Mult_glScale(0.2);                   // Moon has radius 0.2
    MoonMatrix.DumpByColumns(matEntries);
	glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
	glVertexAttrib3f(vertColor_loc, 0.9f, 0.9f, 0.9f);     // Make the moon bright gray
	Moon1.Render();

    LinearMapR4 PlanetXPosMatrix = viewMatrix;			// Place Sun2 at center of the scene
    // SunPosMatrix.Mult_glScale(1.0);				// Scaling by (1, 1, 1) has no effect
    double revolveAnglePlanetX =  -1 * (DayOfYearX  / 500.0) * PI2;
    PlanetXPosMatrix.Mult_glRotate(revolveAnglePlanetX, 0.0, 1.0, 0.0);   // Revolve the sun1 around sun2
    PlanetXPosMatrix.Mult_glTranslate(0.0, 0.0, 7.0);		        // Place the earth five units away from the sun
    PlanetXPosMatrix.Mult_glScale(0.6);
    PlanetXPosMatrix.DumpByColumns(matEntries);           // These two lines load the matrix into the shader
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    glVertexAttrib3f(vertColor_loc, 1.0f, 0.0f, 0.0f);     // Make the sun yellow
    Earth.Render();
    
    LinearMapR4 SDMatrix = EarthMatrix;//bcoz earthmatrix contains rotation
    //double SDRotationAngle = (HourOfDay / 24.0) * PI2;
    //SDMatrix.Mult_glRotate(SDRotationAngle, 0.0, 1.0, 0.0);   // Rotate earth on y-axis
                                    // Make radius 0.5.
    SDMatrix.Mult_glTranslate(0.84, 0.55, 0.0);		// Place the earth five units away from the sun
    
    SDMatrix.Mult_glScale(0.1); 
    SDMatrix.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    glVertexAttrib3f(vertColor_loc, 1.0f, 0.0f, 1.0f);     // Make the earth bright cyan-blue
    SD.Render();

    // MoonMatrix - control placement, and size of the moon.
    LinearMapR4 moonmoonMatrix = MoonMatrix;        // Base the moon's matrix off the earth's *POS* matrix (EarthPosMatrix)
    double moonmoonRotationAngle = (DayOfYear * 12.0 / 365.0)* PI2;
    moonmoonMatrix.Mult_glRotate(moonmoonRotationAngle, 0.0, 1.0, 0.0);  // Revolving around the earth twelve times per year
    moonmoonMatrix.Mult_glTranslate(0.0, 0.0, 3.0);	    // Place the Moon one unit away from the earth
    moonmoonMatrix.Mult_glScale(0.6);                   // Moon has radius 0.2
    moonmoonMatrix.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    glVertexAttrib3f(vertColor_loc, 0.0f, 1.0f, 1.0f);     // Make the moon bright gray
    moonmoon.Render();


	check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}

void my_setup_SceneData() {
	mySetupGeometries();
	setup_shaders();

	// Get the locations of the projection and model view matrices in the shader programs.
	projMatLocation = glGetUniformLocation(shaderProgram1, projMatName);
	modelviewMatLocation = glGetUniformLocation(shaderProgram1, modelviewMatName);

    // Initialize for animation (not really necessary)
    glfwSetTime(PreviousTime);     // PreviousTime is equal to 0.0 when this line reached.

	check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}

// *******************************************************
// Process all key press events.
// This routine is called each time a key is pressed or released.
// *******************************************************
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (action == GLFW_RELEASE) {
		return;			// Ignore key up (key release) events
	}
	switch (key) {
	case GLFW_KEY_ESCAPE:
		glfwSetWindowShouldClose(window, true);
		break;
    case GLFW_KEY_R:
        if (singleStep) {			// If ending single step mode
            singleStep = false;
            spinMode = true;		// Restart animation
        }
        else {
            spinMode = !spinMode;	// Toggle animation on and off.
        }
        if (spinMode && UseRealTime) {
            PreviousTime = glfwGetTime();   // Reset base time in case of "real time" animation
        }
        break;
    case GLFW_KEY_S:
        singleStep = true;
        spinMode = true;
        break;
    case GLFW_KEY_C:		// Toggle backface culling
        cullBackFaces = !cullBackFaces;     // Negate truth value of cullBackFaces
        if (cullBackFaces) {
            glEnable(GL_CULL_FACE);
        }
        else {
            glDisable(GL_CULL_FACE);
        }
        break;
    case GLFW_KEY_T:		// Toggle using real time versus fixed time step
        UseRealTime = !UseRealTime;
        if (UseRealTime) {
            glfwSetTime(0.0);
            PreviousTime = 0.0;
        }
        break;
    // Next two cases: should check for underflow or overflow
    case GLFW_KEY_UP:
        AnimateIncrement *= 2.0;			// Double the animation time step
		break;
	case GLFW_KEY_DOWN:
        AnimateIncrement *= 0.5;			// Halve the animation time step
		break;
	}
}

// *************************************************
// This function is called with the graphics window is first created,
//    and again whenever it is resized.
// The Projection Matrix is set here.
// Conventions for how the projection matrix is set
// The "viewscreen" is a rectangle occupying  (-Xmax,Xmax)x(-Ymax,Ymax)x{-CameraDistance}
//    Any object with z<Zmin-CameraDistAnce or z>Zmax-CameraDistance is clipped.
//    Unless Zmax-CameraDistance is >-ZnearMin, in which case z>-ZnearMin is used as near clipping plane. 
// *************************************************
void window_size_callback(GLFWwindow* window, int width, int height) {
	// Define the portion of the window used for OpenGL rendering.
	glViewport(0, 0, width, height);

	// Setup the projection matrix as a perspective view.
	// The complication is that the aspect ratio of the window may not match the
	//		aspect ratio of the scene we want to view.
	double w = (width == 0) ? 1.0 : (double)width;
	double h = (height == 0) ? 1.0 : (double)height;
	double windowXmax, windowYmax;
    double aspectFactor = w * Ymax / (h * Xmax);   // == (w/h)/(Xmax/Ymax), ratio of aspect ratios
	if (aspectFactor>1) {
		windowXmax = Xmax * aspectFactor;
		windowYmax = Ymax;
	}
	else {
		windowYmax = Ymax / aspectFactor;
		windowXmax = Xmax;
	}

	// Using the max & min values for x & y & z that should be visible in the window,
	//		we set up the perspective projection.
    double zFar = -Zmin + CameraDistance;
    double zNear = Max(CameraDistance - Zmax, ZnearMin);
	theProjectionMatrix.Set_glFrustum(-windowXmax, windowXmax, -windowYmax, windowYmax, zNear, zFar);

	if (glIsProgram(shaderProgram1)) {
		glUseProgram(shaderProgram1);
		theProjectionMatrix.DumpByColumns(matEntries);
		glUniformMatrix4fv(projMatLocation, 1, false, matEntries);
	}
	check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}

void my_setup_OpenGL() {
	
	glEnable(GL_DEPTH_TEST);	// Enable depth buffering
	glDepthFunc(GL_LEQUAL);		// Useful for multipass shaders
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glCullFace(GL_BACK);        // GL_BACK is the default anyway

	check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}

void error_callback(int error, const char* description)
{
	// Print error
	fputs(description, stderr);
}

void setup_callbacks(GLFWwindow* window) {
	// Set callback function for resizing the window
	glfwSetFramebufferSizeCallback(window, window_size_callback);

	// Set callback for key up/down/repeat events
	glfwSetKeyCallback(window, key_callback);

	// Set callbacks for mouse movement (cursor position) and mouse botton up/down events.
	// glfwSetCursorPosCallback(window, cursor_pos_callback);
	// glfwSetMouseButtonCallback(window, mouse_button_callback);
}

int main() {
	glfwSetErrorCallback(error_callback);	// Supposed to be called in event of errors. (doesn't work?)
	glfwInit();
#if defined(__APPLE__) || defined(__linux__)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    //glfwWindowHint(GLFW_SAMPLES, 4);        // Invoke Multisample Antialiasing (MSAA)
    //glEnable(GL_MULTISAMPLE);               // Usually is enabled by default, but call just in case.

	const int initWidth = 800;
	const int initHeight = 600;
	GLFWwindow* window = glfwCreateWindow(initWidth, initHeight, "SolarModern", NULL, NULL);
	if (window == NULL) {
		printf("Failed to create GLFW window!\n");
		return -1;
	}
	glfwMakeContextCurrent(window);

	if (GLEW_OK != glewInit()) {
		printf("Failed to initialize GLEW!.\n");
		return -1;
	}

	// Print info of GPU and supported OpenGL version
	printf("Renderer: %s\n", glGetString(GL_RENDERER));
	printf("OpenGL version supported %s\n", glGetString(GL_VERSION));
#ifdef GL_SHADING_LANGUAGE_VERSION
	printf("Supported GLSL version is %s.\n", (char *)glGetString(GL_SHADING_LANGUAGE_VERSION));
#endif
    printf("Using GLEW version %s.\n", glewGetString(GLEW_VERSION));

	printf("------------------------------\n");
	printf("Press 'r' or 'R' (Run) to toggle(off and on) running the animation.\n");
    printf("Press 's' or 'S' to single-step animation.\n");
    printf("Press 't' or 'T' (Time) to toggle animation between fixed time step, and real elapsed time.\n");
    printf("Press up and down arrow keys to increase and decrease animation rate.\n   ");
    printf("    - animation step size is doubled or halved with each press.\n");
	printf("Press 'c' or 'C' (Cull) to toggle whether back faces are culled.\n");
	printf("Press ESCAPE to exit.\n");
	
    setup_callbacks(window);
   
	// Initialize OpenGL, the scene and the shaders
    my_setup_OpenGL();
	my_setup_SceneData();
 	window_size_callback(window, initWidth, initHeight);

    // Loop while program is not terminated.
	while (!glfwWindowShouldClose(window)) {
	
		myRenderScene();				// Render into the current buffer
		glfwSwapBuffers(window);		// Displays what was just rendered (using double buffering).

		// Poll events (key presses, mouse events)
		glfwWaitEventsTimeout(1.0/60.0);	    // Use this to animate at 60 frames/sec (timing is NOT reliable)
		// glfwWaitEvents();					// Or, Use this instead if no animation.
		// glfwPollEvents();					// Use this version when animating as fast as possible
	}

	glfwTerminate();
	return 0;
}

// If an error is found, it could have been caused by any command since the
//   previous call to check_for_opengl_errors()
// To find what generated the error, you can try adding more calls to
//   check_for_opengl_errors().
char errNames[9][36] = {
    "Unknown OpenGL error",
    "GL_INVALID_ENUM", "GL_INVALID_VALUE", "GL_INVALID_OPERATION",
    "GL_INVALID_FRAMEBUFFER_OPERATION", "GL_OUT_OF_MEMORY",
    "GL_STACK_UNDERFLOW", "GL_STACK_OVERFLOW", "GL_CONTEXT_LOST" };
bool check_for_opengl_errors() {
    int numErrors = 0;
    GLenum err;
    while ((err = glGetError()) != GL_NO_ERROR) {
        numErrors++;
        int errNum = 0;
        switch (err) {
        case GL_INVALID_ENUM:
            errNum = 1;
            break;
        case GL_INVALID_VALUE:
            errNum = 2;
            break;
        case GL_INVALID_OPERATION:
            errNum = 3;
            break;
        case GL_INVALID_FRAMEBUFFER_OPERATION:
            errNum = 4;
            break;
        case GL_OUT_OF_MEMORY:
            errNum = 5;
            break;
        case GL_STACK_UNDERFLOW:
            errNum = 6;
            break;
        case GL_STACK_OVERFLOW:
            errNum = 7;
            break;
        case GL_CONTEXT_LOST:
            errNum = 8;
            break;
        }
        printf("OpenGL ERROR: %s.\n", errNames[errNum]);
    }
    return (numErrors != 0);
}
