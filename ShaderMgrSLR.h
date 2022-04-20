// *******************************
// ShaderMgrSLR.h - Version 1.2 - March 30, 2020
//
// ShaderMgrSLRcpp code defines, compiles and manages shaders
//           for the SolarModern.cpp program
//
// For use with Math 155A - Winter 2017.
//
// Software is "as-is" and carries no warranty. It may be used without
//  restriction, but if you modify it, please change the filenames to
//  prevent confusion between different versions.
// Bug reports: Sam Buss, sbuss@ucsd.edu
// *******************************


#pragma once

void setup_shaders();
unsigned int setup_shader_vertfrag(const char* vertexShaderSource, const char* fragmentShaderSource);

GLuint check_compilation_shader(GLuint shader);
GLuint check_link_status(GLuint program);



