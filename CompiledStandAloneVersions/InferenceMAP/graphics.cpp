/*
 *
 *  InferenceMAP v1.0
 *  4/3/2015
 *
 *  Author: Mohamed El Beheiry, Physico-Chimie Curie, Institut Curie
 *  		Jean-Baptiste Masson, Physics of Biological Systems, Institut Pasteur
 *  Contact e-mail: mohamed.elbeheiry@gmail.com
 *  Copyright (c) 2015, Mohamed El Beheiry, Jean-Baptiste Masson, Institut Curie, Institut Pasteur
 *  All rights Reserved.
 *
 *  InferenceMAP is released under an "academic use only" licence.
 *  Details of which are provided in the InferenceMAP_License.doc file.
 *  Usage of InferenceMAP requires acceptance of this license.
 *
 *  User instructions for using InferenceMAP are provided in the InferenceMAP User Manual.
 *
 */

#include "stdafx.h"

#include "globals.h"
#include "file.h"
#include "graphics.h"
#include "draw.h"
#include "annotations.h"

#define PI 3.1415926535897932384626433832795028841971693993751

// global variable class
extern Globals *iMAP;
extern File *file;
extern void stereoVisionCallback(Fl_Widget*w,void*v);

class iMAPGlWindow;

// OpenGL Extensions
#ifdef __APPLE__
	PFNGLACTIVETEXTUREPROC glActiveTexture;
	PFNGLATTACHSHADERPROC glAttachShader;
	PFNGLBINDBUFFERPROC glBindBuffer;
	PFNGLBINDFRAMEBUFFERPROC glBindFramebuffer;
	PFNGLBINDRENDERBUFFERPROC glBindRenderbuffer;
	PFNGLBINDVERTEXARRAYPROC glBindVertexArray;
	PFNGLBUFFERDATAPROC glBufferData;
	PFNGLCHECKFRAMEBUFFERSTATUSPROC glCheckFramebufferStatus;
	PFNGLCOLORTABLEPROC glColorTable;
	PFNGLCOMPILESHADERPROC glCompileShader;
	PFNGLCREATEPROGRAMPROC glCreateProgram;
	PFNGLCREATESHADERPROC glCreateShader;
	PFNGLDELETEBUFFERSPROC glDeleteBuffers;
	PFNGLDELETEPROGRAMPROC glDeleteProgram;
	PFNGLDELETESHADERPROC glDeleteShader;
	PFNGLDETACHSHADERPROC glDetachShader;
	PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray;
	PFNGLFRAMEBUFFERRENDERBUFFERPROC glFramebufferRenderbuffer;
	PFNGLFRAMEBUFFERTEXTURE2DPROC glFramebufferTexture2D;
	PFNGLGENBUFFERSPROC glGenBuffers;
	PFNGLGENFRAMEBUFFERSPROC glGenFramebuffers;
	PFNGLGENRENDERBUFFERSPROC glGenRenderbuffers;
	PFNGLGENVERTEXARRAYSPROC glGenVertexArrays;
	PFNGLGETATTRIBLOCATIONPROC glGetAttribLocation;
	PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog;
	PFNGLGETPROGRAMIVPROC glGetProgramiv;
	PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
	PFNGLGETSHADERIVPROC glGetShaderiv;
	PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
	PFNGLISSHADERPROC glIsShader;
	PFNGLLINKPROGRAMPROC glLinkProgram;
	PFNGLMULTITEXCOORD3FPROC glMultiTexCoord3f;
	PFNGLPROGRAMPARAMETERIEXTPROC glProgramParameteriEXT;
	PFNGLRENDERBUFFERSTORAGEPROC glRenderbufferStorage;
	PFNGLSHADERSOURCEPROC glShaderSource;
	PFNGLTEXIMAGE3DPROC glTexImage3D;
	PFNGLUNIFORM1IPROC glUniform1i;
	PFNGLUNIFORM1FPROC glUniform1f;
	PFNGLUNIFORM2FPROC glUniform2f;
	PFNGLUNIFORM3FPROC glUniform3f;
	PFNGLUNIFORMMATRIX3FVPROC glUniformMatrix3fv;
	PFNGLUNIFORMMATRIX4FVPROC glUniformMatrix4fv;
	PFNGLUSEPROGRAMPROC glUseProgram;
	PFNGLVERTEXATTRIB1FPROC glVertexAttrib1f;
	PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer;
#endif

void setupGLExtensions() {

//	glewExperimental = GL_TRUE;
//	glewInit();

	#ifdef __APPLE__
		glActiveTexture = (PFNGLACTIVETEXTUREPROC) glutGetProcAddress("glActiveTexture");
		glAttachShader = (PFNGLATTACHSHADERPROC) glutGetProcAddress("glAttachShader");
		glBindBuffer = (PFNGLBINDBUFFERPROC) glutGetProcAddress("glBindBuffer");
		glBindFramebuffer = (PFNGLBINDFRAMEBUFFERPROC) glutGetProcAddress("glBindFramebuffer");
		glBindRenderbuffer = (PFNGLBINDRENDERBUFFERPROC) glutGetProcAddress("glBindRenderbuffer");
		glBindVertexArray = (PFNGLBINDVERTEXARRAYPROC) glutGetProcAddress("glBindVertexArray");
		glBufferData = (PFNGLBUFFERDATAPROC) glutGetProcAddress("glBufferData");
		glCheckFramebufferStatus = (PFNGLCHECKFRAMEBUFFERSTATUSPROC) glutGetProcAddress("glCheckFramebufferStatus");
		glColorTable = (PFNGLCOLORTABLEPROC) glutGetProcAddress("glColorTable");
		glCompileShader = (PFNGLCOMPILESHADERPROC) glutGetProcAddress("glCompileShader");
		glCreateProgram = (PFNGLCREATEPROGRAMPROC) glutGetProcAddress("glCreateProgram");
		glCreateShader = (PFNGLCREATESHADERPROC) glutGetProcAddress("glCreateShader");
		glDeleteBuffers = (PFNGLDELETEBUFFERSPROC) glutGetProcAddress("glDeleteBuffers");
		glDeleteProgram = (PFNGLDELETEPROGRAMPROC) glutGetProcAddress("glDeleteProgram");
		glDeleteShader = (PFNGLDELETESHADERPROC) glutGetProcAddress("glDeleteShader");
		glDetachShader = (PFNGLDETACHSHADERPROC) glutGetProcAddress("glDetachShader");
		glEnableVertexAttribArray = (PFNGLENABLEVERTEXATTRIBARRAYPROC) glutGetProcAddress("glEnableVertexAttribArray");
		glFramebufferRenderbuffer = (PFNGLFRAMEBUFFERRENDERBUFFERPROC) glutGetProcAddress("glFramebufferRenderbuffer");
		glFramebufferTexture2D = (PFNGLFRAMEBUFFERTEXTURE2DPROC) glutGetProcAddress("glFramebufferTexture2D");
		glGetAttribLocation = (PFNGLGETATTRIBLOCATIONPROC) glutGetProcAddress("glGetAttribLocation");
		glGenBuffers = (PFNGLGENBUFFERSPROC) glutGetProcAddress("glGenBuffers");
		glGenFramebuffers = (PFNGLGENFRAMEBUFFERSPROC) glutGetProcAddress("glGenFramebuffers");
		glGenRenderbuffers = (PFNGLGENRENDERBUFFERSPROC) glutGetProcAddress("glGenRenderbuffers");
		glGenVertexArrays = (PFNGLGENVERTEXARRAYSPROC) glutGetProcAddress("glGenVertexArrays");
		glGetProgramInfoLog = (PFNGLGETPROGRAMINFOLOGPROC) glutGetProcAddress("glGetProgramInfoLog");
		glGetProgramiv = (PFNGLGETPROGRAMIVPROC) glutGetProcAddress("glGetProgramiv");
		glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) glutGetProcAddress("glGetShaderInfoLog");
		glGetShaderiv = (PFNGLGETSHADERIVPROC) glutGetProcAddress("glGetShaderiv");
		glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) glutGetProcAddress("glGetUniformLocation");
		glIsShader = (PFNGLISSHADERPROC) glutGetProcAddress("glIsShader");
		glLinkProgram = (PFNGLLINKPROGRAMPROC) glutGetProcAddress("glLinkProgram");
		glMultiTexCoord3f = (PFNGLMULTITEXCOORD3FPROC) glutGetProcAddress("glMultiTexCoord3f");
		glProgramParameteriEXT = (PFNGLPROGRAMPARAMETERIEXTPROC) glutGetProcAddress("glProgramParameteriEXT");
		glRenderbufferStorage = (PFNGLRENDERBUFFERSTORAGEPROC) glutGetProcAddress("glRenderbufferStorage");
		glShaderSource = (PFNGLSHADERSOURCEPROC) glutGetProcAddress("glShaderSource");
		glTexImage3D = (PFNGLTEXIMAGE3DPROC) glutGetProcAddress("glTexImage3D");
		glUniform1i = (PFNGLUNIFORM1IPROC) glutGetProcAddress("glUniform1i");
		glUniform1f = (PFNGLUNIFORM1FPROC) glutGetProcAddress("glUniform1f");
		glUniform2f = (PFNGLUNIFORM2FPROC) glutGetProcAddress("glUniform2f");
		glUniform3f = (PFNGLUNIFORM3FPROC) glutGetProcAddress("glUniform3f");
		glUniformMatrix3fv = (PFNGLUNIFORMMATRIX3FVPROC) glutGetProcAddress("glUniformMatrix3fv");
		glUniformMatrix4fv = (PFNGLUNIFORMMATRIX4FVPROC) glutGetProcAddress("glUniformMatrix4fv");
		glUseProgram = (PFNGLUSEPROGRAMPROC) glutGetProcAddress("glUseProgram");
		glVertexAttrib1f = (PFNGLVERTEXATTRIB1FPROC) glutGetProcAddress("glVertexAttrib1f");
		glVertexAttribPointer = (PFNGLVERTEXATTRIBPOINTERPROC) glutGetProcAddress("glVertexAttribPointer");
	#elif _WIN32
		glActiveTexture = (PFNGLACTIVETEXTUREPROC) wglGetProcAddress("glActiveTexture");
		glAttachShader = (PFNGLATTACHSHADERPROC) wglGetProcAddress("glAttachShader");
		glBindBuffer = (PFNGLBINDBUFFERPROC) wglGetProcAddress("glBindBuffer");
		glBindFramebuffer = (PFNGLBINDFRAMEBUFFERPROC) wglGetProcAddress("glBindFramebuffer");
		glBindRenderbuffer = (PFNGLBINDRENDERBUFFERPROC) wglGetProcAddress("glBindRenderbuffer");
		glBindVertexArray = (PFNGLBINDVERTEXARRAYPROC) wglGetProcAddress("glBindVertexArray");
		glBufferData = (PFNGLBUFFERDATAPROC) wglGetProcAddress("glBufferData");
		glCheckFramebufferStatus = (PFNGLCHECKFRAMEBUFFERSTATUSPROC) wglGetProcAddress("glCheckFramebufferStatus");
		glColorTable = (PFNGLCOLORTABLEPROC) wglGetProcAddress("glColorTable");
		glCompileShader = (PFNGLCOMPILESHADERPROC) wglGetProcAddress("glCompileShader");
		glCreateProgram = (PFNGLCREATEPROGRAMPROC) wglGetProcAddress("glCreateProgram");
		glCreateShader = (PFNGLCREATESHADERPROC) wglGetProcAddress("glCreateShader");
		glDeleteBuffers = (PFNGLDELETEBUFFERSPROC) wglGetProcAddress("glDeleteBuffers");
		glDeleteProgram = (PFNGLDELETEPROGRAMPROC) wglGetProcAddress("glDeleteProgram");
		glDeleteShader = (PFNGLDELETESHADERPROC) wglGetProcAddress("glDeleteShader");
		glDetachShader = (PFNGLDETACHSHADERPROC) wglGetProcAddress("glDetachShader");
		glEnableVertexAttribArray = (PFNGLENABLEVERTEXATTRIBARRAYPROC) wglGetProcAddress("glEnableVertexAttribArray");
		glFramebufferRenderbuffer = (PFNGLFRAMEBUFFERRENDERBUFFERPROC) wglGetProcAddress("glFramebufferRenderbuffer");
		glFramebufferTexture2D = (PFNGLFRAMEBUFFERTEXTURE2DPROC) wglGetProcAddress("glFramebufferTexture2D");
		glGetAttribLocation = (PFNGLGETATTRIBLOCATIONPROC) wglGetProcAddress("glGetAttribLocation");
		glGenBuffers = (PFNGLGENBUFFERSPROC) wglGetProcAddress("glGenBuffers");
		glGenFramebuffers = (PFNGLGENFRAMEBUFFERSPROC) wglGetProcAddress("glGenFramebuffers");
		glGenRenderbuffers = (PFNGLGENRENDERBUFFERSPROC) wglGetProcAddress("glGenRenderbuffers");
		glGenVertexArrays = (PFNGLGENVERTEXARRAYSPROC) wglGetProcAddress("glGenVertexArrays");
		glGetProgramInfoLog = (PFNGLGETPROGRAMINFOLOGPROC) wglGetProcAddress("glGetProgramInfoLog");
		glGetProgramiv = (PFNGLGETPROGRAMIVPROC) wglGetProcAddress("glGetProgramiv");
		glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) wglGetProcAddress("glGetShaderInfoLog");
		glGetShaderiv = (PFNGLGETSHADERIVPROC) wglGetProcAddress("glGetShaderiv");
		glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) wglGetProcAddress("glGetUniformLocation");
		glIsShader = (PFNGLISSHADERPROC) wglGetProcAddress("glIsShader");
		glLinkProgram = (PFNGLLINKPROGRAMPROC) wglGetProcAddress("glLinkProgram");
		glMultiTexCoord3f = (PFNGLMULTITEXCOORD3FPROC) wglGetProcAddress("glMultiTexCoord3f");
		glProgramParameteriEXT = (PFNGLPROGRAMPARAMETERIEXTPROC) wglGetProcAddress("glProgramParameteriEXT");
		glRenderbufferStorage = (PFNGLRENDERBUFFERSTORAGEPROC) wglGetProcAddress("glRenderbufferStorage");
		glShaderSource = (PFNGLSHADERSOURCEPROC) wglGetProcAddress("glShaderSource");
		glTexImage3D = (PFNGLTEXIMAGE3DPROC) wglGetProcAddress("glTexImage3D");
		glUniform1i = (PFNGLUNIFORM1IPROC) wglGetProcAddress("glUniform1i");
		glUniform1f = (PFNGLUNIFORM1FPROC) wglGetProcAddress("glUniform1f");
		glUniform2f = (PFNGLUNIFORM2FPROC) wglGetProcAddress("glUniform2f");
		glUniform3f = (PFNGLUNIFORM3FPROC) wglGetProcAddress("glUniform3f");
		glUniformMatrix3fv = (PFNGLUNIFORMMATRIX3FVPROC) wglGetProcAddress("glUniformMatrix3fv");
		glUniformMatrix4fv = (PFNGLUNIFORMMATRIX4FVPROC) wglGetProcAddress("glUniformMatrix4fv");
		glUseProgram = (PFNGLUSEPROGRAMPROC) wglGetProcAddress("glUseProgram");
		glVertexAttrib1f = (PFNGLVERTEXATTRIB1FPROC) wglGetProcAddress("glVertexAttrib1f");
		glVertexAttribPointer = (PFNGLVERTEXATTRIBPOINTERPROC) wglGetProcAddress("glVertexAttribPointer");
	#endif

}
// create gaussian detection texture
float* detectionImage(float sigma, int dimension) {
	float *image = new float[dimension*dimension];
	for (int y = 0; y < dimension; y++) {
		for (int x = 0; x < dimension; x++) {
			image[y*dimension+x] = exp( -( (float)(y-dimension/2)*(y-dimension/2) + (float)(x-dimension/2)*(x-dimension/2) )/(2*sigma) );
		}
	}
	return image;
}

void colormap(File *fileObject, int index, float value) {

	if (fileObject->colormapFlip) {
		value = 1-value;
	}

	float max2,max3,max6;
	float rgb[] = { 0, 0, 0 };

	const float cMin = fileObject->cMin;
	float cMax = fileObject->cMax;

	switch (fileObject->colormapType) {
	// autumn
	case 1:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=0;rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		break;
	// blue
	case 2:
		value-=cMin;
		if(value>1)
			{rgb[0]=rgb[1]=rgb[2]=1;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=value/cMax;rgb[2]=1;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// blue red
	case 3:
		max3=(cMax-cMin)/2;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=0;rgb[2]=0;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<max3)
			{rgb[0]=(value/max3);rgb[1]=(value/max3);rgb[2]=1;}
		else if(value<2*max3)
			{rgb[0]=1;rgb[1]=(1-(value-max3)/max3);rgb[2]=(1-(value-max3)/max3);}
		else {rgb[0]=1;rgb[1]=0;rgb[2]=0;}
		break;
	// cool
	case 4:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=1;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=(1-value/cMax);rgb[2]=1;}
		else {rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		break;
	// cyan hot
	case 5:
		max2=(cMax-cMin)/2;
			value-=cMin;
		if(value>1)
			{rgb[0]=rgb[1]=rgb[2]=1;}
		else if(value<0)
			{rgb[0]=rgb[1]=0;rgb[2]=1;}
		else if(value<max2)
		{rgb[0]=0;rgb[1]=(value/max2);rgb[2]=1;}
		else if(value<2*max2)
			{rgb[0]=(value-max2)/max2;rgb[1]=1;rgb[2]=1;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// fire
	case 6:
		max3=(cMax-cMin)/3;
			value-=cMin;
		if(value>1)
			{rgb[0]=rgb[1]=rgb[2]=1;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		else if(value<max3)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1-(value/max3);}
		else if(value<2*max3)
			{rgb[0]=1;rgb[1]=(value-max3)/max3;rgb[2]=0;}
		else if(value<3*max3)
			{rgb[0]=1;rgb[1]=1;rgb[2]=(value-2*max3)/max3;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// gray
	case 7:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=rgb[1]=rgb[2]=1;}
		else if(value<0)
			{rgb[0]=rgb[1]=rgb[2]=0;}
		else if(value<cMax) {
			rgb[0]=rgb[1]=rgb[2]=value/cMax;
		}
		break;
	// green
	case 8:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=1;rgb[2]=1;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=value/cMax;}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=1;}
		break;
	// hsv
	case 9:
		max6=(cMax-cMin)/6;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<max6)
			{rgb[0]=1;rgb[1]=0;rgb[2]=value/max6;}
		else if(value<2*max6)
			{rgb[0]=(1-(value-max6)/max6);rgb[1]=0;rgb[2]=1;}
		else if(value<3*max6)
			{rgb[0]=0;rgb[1]=(value-2*max6)/max6;rgb[2]=1;}
		else if(value<4*max6)
			{rgb[0]=0;rgb[1]=1;rgb[2]=(1-(value-3*max6)/max6);}
		else if(value<5*max6)
			{rgb[0]=(value-4*max6)/max6;rgb[1]=1;rgb[2]=0;}
		else if(value<6*max6)
			{rgb[0]=1;rgb[1]=(1-(value-5*max6)/max6);rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=rgb[2]=0;}
		break;
	// jet
	case 10:
		max3=(cMax-cMin)/3;
		value-=cMin;
		if(value>1)
			{rgb[0]=rgb[1]=rgb[2]=1;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<max3)
			{rgb[0]=0;rgb[1]=(value/max3);rgb[2]=1;}
		else if(value<2*max3)
			{rgb[0]=(value-max3)/max3;rgb[1]=1;rgb[2]=1-rgb[0];}
		else if(value<3*max3)
			{rgb[0]=1;rgb[1]=(1-(value-2*max3)/max3);rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=rgb[2]=0;}
		break;
	// red
	case 11:
		value-=cMin;
		if(value>1)
			{rgb[0]=rgb[1]=rgb[2]=1;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=value/cMax;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// redhot
	case 12:
		max2=(cMax-cMin)/2;
			value-=cMin;
		if(value>1)
			{rgb[0]=rgb[1]=rgb[2]=1;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<max2)
		{rgb[0]=1;rgb[1]=(value/max2);rgb[2]=0;}
		else if(value<2*max2)
			{rgb[0]=1;rgb[1]=1;rgb[2]=(value-max2)/max2;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// spring
	case 13:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		break;
	// summer
	case 14:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		break;
	// winter
	case 15:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=0;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else {rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		break;
	}

	// store r
	fileObject->colorArray[index*4] = rgb[0];
	// store g
	fileObject->colorArray[index*4+1] = rgb[1];
	// store b
	fileObject->colorArray[index*4+2] = rgb[2];
}

float* colormap(int colormapType, float value, float cMin, float cMax, bool flip) {

	if (flip) {
		value = 1.0f-value;
	}

	float max2,max3,max6;
	float rgb[] = { 0.0f, 0.0f, 0.0f };

	switch (colormapType) {
	// autumn
	case 0:
		cMax-=cMin;
		value-=cMin;
		if(value>1)	{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=0;}
		else if(value<0) {rgb[0]=1;rgb[1]=0;rgb[2]=0;}
		else if(value<cMax)	{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		break;
	// blue
	case 1:
		value-=cMin;
		if(value>1)	{rgb[0]=value/cMax;rgb[1]=value/cMax;rgb[2]=1;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=value/cMax;rgb[2]=1;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// blue red
	case 2:
		max3=(cMax-cMin)/2;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=(1-(value-max3)/max3);rgb[2]=(1-(value-max3)/max3);}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<max3)
			{rgb[0]=(value/max3);rgb[1]=(value/max3);rgb[2]=1;}
		else if(value<2*max3)
			{rgb[0]=1;rgb[1]=(1-(value-max3)/max3);rgb[2]=(1-(value-max3)/max3);}
		else {rgb[0]=1;rgb[1]=0;rgb[2]=0;}
		break;
	// cool
	case 3:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=value/cMax;rgb[1]=(1-value/cMax);rgb[2]=1;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=1;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=(1-value/cMax);rgb[2]=1;}
		else {rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		break;
	// cyan hot
	case 4:
		max2=(cMax-cMin)/2;
			value-=cMin;
		if(value>1)
			{rgb[0]=(value-max2)/max2;rgb[1]=1;rgb[2]=1;}
		else if(value<0)
			{rgb[0]=rgb[1]=0;rgb[2]=1;}
		else if(value<max2)
		{rgb[0]=0;rgb[1]=(value/max2);rgb[2]=1;}
		else if(value<2*max2)
			{rgb[0]=(value-max2)/max2;rgb[1]=1;rgb[2]=1;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// fire
	case 5:
		max3=(cMax-cMin)/3;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=1;rgb[2]=(value-2*max3)/max3;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		else if(value<max3)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1-(value/max3);}
		else if(value<2*max3)
			{rgb[0]=1;rgb[1]=(value-max3)/max3;rgb[2]=0;}
		else if(value<3*max3)
			{rgb[0]=1;rgb[1]=1;rgb[2]=(value-2*max3)/max3;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// gray
	case 6:
		cMax-=cMin;
		value-=cMin;
		if(value>1)	{ rgb[0]=rgb[1]=rgb[2]=value/cMax;}
		else if(value<0) {rgb[0]=rgb[1]=rgb[2]=0;}
		else if(value<cMax) { rgb[0]=rgb[1]=rgb[2]=value/cMax;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// green
	case 7:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=value/cMax;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=value/cMax;}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=1;}
		break;
	// hsv
	case 8:
		max6=(cMax-cMin)/6;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=(1-(value-5*max6)/max6);rgb[2]=0;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<max6)
			{rgb[0]=1;rgb[1]=0;rgb[2]=value/max6;}
		else if(value<2*max6)
			{rgb[0]=(1-(value-max6)/max6);rgb[1]=0;rgb[2]=1;}
		else if(value<3*max6)
			{rgb[0]=0;rgb[1]=(value-2*max6)/max6;rgb[2]=1;}
		else if(value<4*max6)
			{rgb[0]=0;rgb[1]=1;rgb[2]=(1-(value-3*max6)/max6);}
		else if(value<5*max6)
			{rgb[0]=(value-4*max6)/max6;rgb[1]=1;rgb[2]=0;}
		else if(value<6*max6)
			{rgb[0]=1;rgb[1]=(1-(value-5*max6)/max6);rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=rgb[2]=0;}
		break;
	// jet
	case 9:
		max3=(cMax-cMin)/3;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=(1-(value-2*max3)/max3);rgb[2]=0;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<max3)
			{rgb[0]=0;rgb[1]=(value/max3);rgb[2]=1;}
		else if(value<2*max3)
			{rgb[0]=(value-max3)/max3;rgb[1]=1;rgb[2]=1-rgb[0];}
		else if(value<3*max3)
			{rgb[0]=1;rgb[1]=(1-(value-2*max3)/max3);rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=rgb[2]=0;}
		break;
	// red
	case 10:
		value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=value/cMax;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=value/cMax;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// redhot
	case 11:
		max2=(cMax-cMin)/2;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=1;rgb[2]=(value-max2)/max2;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<max2)
		{rgb[0]=1;rgb[1]=(value/max2);rgb[2]=0;}
		else if(value<2*max2)
			{rgb[0]=1;rgb[1]=1;rgb[2]=(value-max2)/max2;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// spring
	case 12:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else if(value<0)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		break;
	// summer
	case 13:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=0;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		break;
	// winter
	case 14:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=0;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=0;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else {rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		break;
	}

	return &rgb[0];
}

void colormap(float *rgb, int colormapType, float value, float cMin, float cMax, bool flip) {

	if (flip) {
		value = 1.0-value;
	}

	float max2,max3,max6;

	switch (colormapType) {
	// autumn
	case 0:
		cMax-=cMin;
		value-=cMin;
		if(value>1)	{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=0;}
		else if(value<0) {rgb[0]=1;rgb[1]=0;rgb[2]=0;}
		else if(value<cMax)	{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		break;
	// blue
	case 1:
		value-=cMin;
		if(value>1)	{rgb[0]=value/cMax;rgb[1]=value/cMax;rgb[2]=1;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=value/cMax;rgb[2]=1;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// blue red
	case 2:
		max3=(cMax-cMin)/2;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=(1-(value-max3)/max3);rgb[2]=(1-(value-max3)/max3);}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<max3)
			{rgb[0]=(value/max3);rgb[1]=(value/max3);rgb[2]=1;}
		else if(value<2*max3)
			{rgb[0]=1;rgb[1]=(1-(value-max3)/max3);rgb[2]=(1-(value-max3)/max3);}
		else {rgb[0]=1;rgb[1]=0;rgb[2]=0;}
		break;
	// cool
	case 3:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=value/cMax;rgb[1]=(1-value/cMax);rgb[2]=1;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=1;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=(1-value/cMax);rgb[2]=1;}
		else {rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		break;
	// cyan hot
	case 4:
		max2=(cMax-cMin)/2;
			value-=cMin;
		if(value>1)
			{rgb[0]=(value-max2)/max2;rgb[1]=1;rgb[2]=1;}
		else if(value<0)
			{rgb[0]=rgb[1]=0;rgb[2]=1;}
		else if(value<max2)
		{rgb[0]=0;rgb[1]=(value/max2);rgb[2]=1;}
		else if(value<2*max2)
			{rgb[0]=(value-max2)/max2;rgb[1]=1;rgb[2]=1;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// fire
	case 5:
		max3=(cMax-cMin)/3;
			value-=cMin;
		if(value>1)
			{rgb[0]=rgb[1]=rgb[2]=1;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		else if(value<max3)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1-(value/max3);}
		else if(value<2*max3)
			{rgb[0]=1;rgb[1]=(value-max3)/max3;rgb[2]=0;}
		else if(value<3*max3)
			{rgb[0]=1;rgb[1]=1;rgb[2]=(value-2*max3)/max3;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// gray
	case 6:
		cMax-=cMin;
		value-=cMin;
		if(value>1)	{ rgb[0]=rgb[1]=rgb[2]=value/cMax;}
		else if(value<0) {rgb[0]=rgb[1]=rgb[2]=0;}
		else if(value<cMax) { rgb[0]=rgb[1]=rgb[2]=value/cMax;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// green
	case 7:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=value/cMax;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=value/cMax;}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=1;}
		break;
	// hsv
	case 8:
		max6=(cMax-cMin)/6;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=(1-(value-5*max6)/max6);rgb[2]=0;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<max6)
			{rgb[0]=1;rgb[1]=0;rgb[2]=value/max6;}
		else if(value<2*max6)
			{rgb[0]=(1-(value-max6)/max6);rgb[1]=0;rgb[2]=1;}
		else if(value<3*max6)
			{rgb[0]=0;rgb[1]=(value-2*max6)/max6;rgb[2]=1;}
		else if(value<4*max6)
			{rgb[0]=0;rgb[1]=1;rgb[2]=(1-(value-3*max6)/max6);}
		else if(value<5*max6)
			{rgb[0]=(value-4*max6)/max6;rgb[1]=1;rgb[2]=0;}
		else if(value<6*max6)
			{rgb[0]=1;rgb[1]=(1-(value-5*max6)/max6);rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=rgb[2]=0;}
		break;
	// jet
	case 9:
		max3=(cMax-cMin)/3;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=(1-(value-2*max3)/max3);rgb[2]=0;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<max3)
			{rgb[0]=0;rgb[1]=(value/max3);rgb[2]=1;}
		else if(value<2*max3)
			{rgb[0]=(value-max3)/max3;rgb[1]=1;rgb[2]=1-rgb[0];}
		else if(value<3*max3)
			{rgb[0]=1;rgb[1]=(1-(value-2*max3)/max3);rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=rgb[2]=0;}
		break;
	// red
	case 10:
		value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=value/cMax;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=value/cMax;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// redhot
	case 11:
		max2=(cMax-cMin)/2;
			value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=1;rgb[2]=(value-max2)/max2;}
		else if(value<0)
			{rgb[0]=1;rgb[1]=rgb[2]=0;}
		else if(value<max2)
		{rgb[0]=1;rgb[1]=(value/max2);rgb[2]=0;}
		else if(value<2*max2)
			{rgb[0]=1;rgb[1]=1;rgb[2]=(value-max2)/max2;}
		else {rgb[0]=rgb[1]=rgb[2]=1;}
		break;
	// spring
	case 12:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else if(value<0)
			{rgb[0]=1;rgb[1]=0;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=1;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		break;
	// summer
	case 13:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=0;}
		else if(value<0)
			{rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		else if(value<cMax)
			{rgb[0]=value/cMax;rgb[1]=1;rgb[2]=0;}
		else {rgb[0]=1;rgb[1]=1;rgb[2]=0;}
		break;
	// winter
	case 14:
		cMax-=cMin;
		value-=cMin;
		if(value>1)
			{rgb[0]=0;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else if(value<0)
			{rgb[0]=0;rgb[1]=0;rgb[2]=1;}
		else if(value<cMax)
			{rgb[0]=0;rgb[1]=value/cMax;rgb[2]=(1-value/cMax);}
		else {rgb[0]=0;rgb[1]=1;rgb[2]=0;}
		break;
	}
}

void defaultScreen() {

	if (iMAP->showLogo) { 
		glClearColor(iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		GLUquadricObj *quadObj = gluNewQuadric();
		gluQuadricNormals(quadObj, GLU_SMOOTH);

		// create gaussian bell
		if (iMAP->gaussianBell == false) { drawGaussianBell(); }

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		
		iMAP->defaultZoom = 25.0;

		gluPerspective(iMAP->defaultZoom, 1.0, 1, 50);

		glMatrixMode(GL_MODELVIEW);

		glPushMatrix();

			glTranslatef(2.5,-0.5,0.0);
			glRotatef(15.0,1.0,0.0,0.0);

			const GLfloat mat_ambient[]    = { 0.20f,0.20f,0.20f,1.0f };  // RGBA
			const GLfloat mat_diffuse[]    = { iMAP->defaultRotate/100.0f,iMAP->defaultRotate/100.0f,iMAP->defaultRotate/100.0f,1.0f };  // RGBA
			const GLfloat mat_specular[]   = { iMAP->defaultRotate/100.0f,iMAP->defaultRotate/100.0f,iMAP->defaultRotate/100.0f,1.0f };  // RGBA
			const GLfloat light_position0[] = { -iMAP->defaultRotate, 50.0, iMAP->defaultRotate, 0.0 };  // XYZ

			glShadeModel(GL_SMOOTH);
			glEnable(GL_COLOR_MATERIAL);

			glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
			glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
			glMaterialf(GL_FRONT,  GL_SHININESS, 100.0);

			glEnable(GL_LIGHT0);
			glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
			glLightfv(GL_LIGHT0, GL_AMBIENT, mat_ambient);
			glLightfv(GL_LIGHT0, GL_DIFFUSE,   mat_diffuse);
			glLightfv(GL_LIGHT0, GL_SPECULAR, mat_ambient);
			glEnable(GL_LIGHTING);
			glEnable(GL_DEPTH_TEST);

			glColor4f(1.0f,1.0f,0.0f,1.0f);

			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GEQUAL, 0.35f);
			glEnable(GL_NORMALIZE);

			const int slices = 30;

			// sphere at origin
			glRotatef(160,0.0,1.0,0.0);
			glTranslatef(6.0,-1.0,0.0);

			glPushMatrix();

				gluSphere(quadObj,0.4f,slices,slices);

				glRotatef(-90,0,1,0);
				gluCylinder(quadObj, 0.2f, 0.2f, 4.0f, slices, slices);
				glPushMatrix();
					gluSphere(quadObj,0.2f,slices,slices);
					glTranslatef(0,0,4);
					glRotatef(60,0,1,0);
					gluCylinder(quadObj, 0.2f, 0.2f, 3.0f, slices, slices);
					glPushMatrix();
						gluSphere(quadObj,0.2f,slices,slices);
						glTranslatef(0,0,3);
						glRotatef(70,-0.3,1,0.6);
						gluCylinder(quadObj, 0.2f, 0.2f, 2.0f, slices, slices);
						glPushMatrix();
						gluSphere(quadObj,0.2f,slices,slices);
							glTranslatef(0,0,2);
							glRotatef(45,0.0,1,-0.3);
							gluCylinder(quadObj, 0.2f, 0.2f, 2.0f, slices, slices);
							glPushMatrix();
								gluSphere(quadObj,0.2f,slices,slices);
								glTranslatef(0,0,2);
								glRotatef(50,0.0,1,0.3);
								gluCylinder(quadObj, 0.2f, 0.2f, 3.0f, slices, slices);
								glPushMatrix();
									gluSphere(quadObj,0.2f,slices,slices);
									glTranslatef(0,0,3);
									glRotatef(85,0.0,1,0.2);
									gluCylinder(quadObj, 0.2f, 0.2f, 3.0f, slices, slices);
									glPushMatrix();
										gluSphere(quadObj,0.2f,slices,slices);
										glTranslatef(0,0,3);
										glRotatef(100,0.0,1,0.4);
										gluCylinder(quadObj, 0.2f, 0.2f, 1.5f, slices, slices);
										glPushMatrix();
											gluSphere(quadObj,0.2f,slices,slices);
											glTranslatef(0,0,1.5);
											gluCylinder(quadObj, 0.5, 0, 1, slices, slices);
										glPopMatrix();
									glPopMatrix();
								glPopMatrix();
							glPopMatrix();
						glPopMatrix();
					glPopMatrix();
				glPopMatrix();
			glPopMatrix ();

			glPushMatrix();
				glScalef(1.0,1.0f,1.0);
				glRotatef(-90,1,0,0);
				glTranslatef(-6.7,-6.7,-1.0);
				glEnableClientState(GL_VERTEX_ARRAY);
				glEnableClientState(GL_COLOR_ARRAY);
				glEnableClientState(GL_NORMAL_ARRAY);

				glVertexPointer(3,GL_FLOAT,0,iMAP->bellQuads);
				glColorPointer(4,GL_FLOAT,0,iMAP->bellColors);
				glNormalPointer(GL_FLOAT,0,iMAP->bellNormals);

				glDrawArrays(GL_QUADS,0,iMAP->bellDimensions);

				glDisableClientState(GL_VERTEX_ARRAY);
				glDisableClientState(GL_COLOR_ARRAY);
				glDisableClientState(GL_NORMAL_ARRAY);
			glPopMatrix();

		glPopMatrix();

		glDisable(GL_ALPHA_TEST);


		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
		glDisable(GL_COLOR_MATERIAL);
		glDisable(GL_DEPTH_TEST);

		gluDeleteQuadric(quadObj);

		// draw text logo
		int width,ascent,descent;

		char label[20];
		glEnable(GL_TEXTURE_2D);
		glColor4f(1.0f,1.0f,1.0f,1.0f);
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER,0.5f);
		sprintf(label,"InferenceMAP");
		txfEstablishTexture(iMAP->fontTex, 0, GL_TRUE);
		txfGetStringMetrics(iMAP->fontTex,label,strlen(label),&width,&ascent,&descent);
		glPushMatrix();
			glTranslatef(-2.0f,-3.0f,0.0f);
			glScalef(60.0f/(float)width/15.0,60.0f/(float)width/15.0,60.0f/(float)width/15.0);
			txfRenderString(iMAP->fontTex, label, strlen(label));
		glPopMatrix();
		glDisable(GL_ALPHA_TEST);
		glDisable(GL_TEXTURE_2D);

		if (iMAP->defaultRotate <= 100) { iMAP->defaultRotate += 1.0; }

	}
}

void drawGaussianBell() {
	const int dimension = 60;
	const float sigma = 70.0;
	const float step = 0.14;
	float **bell = new float*[dimension];
	for (int k = 0; k < dimension; k++) {
		bell[k] = new float[dimension];
	}

	for (int x = 0; x < dimension; x++) {
		for (int y = 0; y < dimension; y++) {
			bell[x][y] = exp( -( (float)(y-dimension/2)*(y-dimension/2) + (float)(x-dimension/2)*(x-dimension/2) )/(2*sigma) );
		}
	}

	iMAP->bellDimensions = 4*(dimension-1)*(dimension-1);
	iMAP->bellQuads = new float[3*iMAP->bellDimensions];
	iMAP->bellColors = new float[4*iMAP->bellDimensions];
	iMAP->bellNormals = new float[3*iMAP->bellDimensions];

	int i = 0;
	float v1[3];
	float v2[3];
	for (int a = 0; a < dimension-1; a++) {
		for (int b = 0; b < dimension-1; b++) {
			// bottom left corner
			iMAP->bellQuads[12*i+0] = a*step;
			iMAP->bellQuads[12*i+1] = b*step;
			iMAP->bellQuads[12*i+2] = bell[a][b];
			colormap(&iMAP->bellColors[16*i+0],4,iMAP->bellQuads[12*i+2],0.0,1.0,false);
			iMAP->bellColors[16*i+3] = 8.0*iMAP->bellQuads[12*i+2];
			iMAP->bellQuads[12*i+2] *= 6.0;
			// bottom right corner
			iMAP->bellQuads[12*i+3] = (a+1)*step;
			iMAP->bellQuads[12*i+4] = b*step;
			iMAP->bellQuads[12*i+5] = bell[a+1][b];
			colormap(&iMAP->bellColors[16*i+4],4,iMAP->bellQuads[12*i+5],0.0,1.0,false);
			iMAP->bellColors[16*i+7] = 8.0*iMAP->bellQuads[12*i+5];
			iMAP->bellQuads[12*i+5] *= 6.0;
			// top right corner
			iMAP->bellQuads[12*i+6] = (a+1)*step;
			iMAP->bellQuads[12*i+7] = (b+1)*step;
			iMAP->bellQuads[12*i+8] = bell[a+1][b+1];
			colormap(&iMAP->bellColors[16*i+8],4,iMAP->bellQuads[12*i+8],0.0,1.0,false);
			iMAP->bellColors[16*i+11] = 8.0*iMAP->bellQuads[12*i+8];
			iMAP->bellQuads[12*i+8] *= 6.0;
			// top left corner
			iMAP->bellQuads[12*i+9] = (a)*step;
			iMAP->bellQuads[12*i+10] = (b+1)*step;
			iMAP->bellQuads[12*i+11] = bell[a][b+1];
			colormap(&iMAP->bellColors[16*i+12],4,iMAP->bellQuads[12*i+11],0.0,1.0,false);
			iMAP->bellColors[16*i+15] = 8.0*iMAP->bellQuads[12*i+11];
			iMAP->bellQuads[12*i+11] *= 6.0;

			// normals
			v1[0] = iMAP->bellQuads[12*i+0]-iMAP->bellQuads[12*i+9];
			v1[1] = iMAP->bellQuads[12*i+1]-iMAP->bellQuads[12*i+10];
			v1[2] = iMAP->bellQuads[12*i+2]-iMAP->bellQuads[12*i+11];
			v2[0] = iMAP->bellQuads[12*i+3]-iMAP->bellQuads[12*i+0];
			v2[1] = iMAP->bellQuads[12*i+4]-iMAP->bellQuads[12*i+1];
			v2[2] = iMAP->bellQuads[12*i+5]-iMAP->bellQuads[12*i+2];
			crossProduct(&iMAP->bellNormals[12*i+0],v1,v2);
			v1[0] = iMAP->bellQuads[12*i+3]-iMAP->bellQuads[12*i+0];
			v1[1] = iMAP->bellQuads[12*i+4]-iMAP->bellQuads[12*i+1];
			v1[2] = iMAP->bellQuads[12*i+5]-iMAP->bellQuads[12*i+2];
			v2[0] = iMAP->bellQuads[12*i+6]-iMAP->bellQuads[12*i+3];
			v2[1] = iMAP->bellQuads[12*i+7]-iMAP->bellQuads[12*i+4];
			v2[2] = iMAP->bellQuads[12*i+8]-iMAP->bellQuads[12*i+5];
			crossProduct(&iMAP->bellNormals[12*i+3],v1,v2);
			v1[0] = iMAP->bellQuads[12*i+6]-iMAP->bellQuads[12*i+3];
			v1[1] = iMAP->bellQuads[12*i+7]-iMAP->bellQuads[12*i+4];
			v1[2] = iMAP->bellQuads[12*i+8]-iMAP->bellQuads[12*i+5];
			v2[0] = iMAP->bellQuads[12*i+9]-iMAP->bellQuads[12*i+6];
			v2[1] = iMAP->bellQuads[12*i+10]-iMAP->bellQuads[12*i+7];
			v2[2] = iMAP->bellQuads[12*i+11]-iMAP->bellQuads[12*i+8];
			crossProduct(&iMAP->bellNormals[12*i+6],v1,v2);
			v1[0] = iMAP->bellQuads[12*i+9]-iMAP->bellQuads[12*i+6];
			v1[1] = iMAP->bellQuads[12*i+10]-iMAP->bellQuads[12*i+7];
			v1[2] = iMAP->bellQuads[12*i+11]-iMAP->bellQuads[12*i+8];
			v2[0] = iMAP->bellQuads[12*i+0]-iMAP->bellQuads[12*i+9];
			v2[1] = iMAP->bellQuads[12*i+1]-iMAP->bellQuads[12*i+10];
			v2[2] = iMAP->bellQuads[12*i+2]-iMAP->bellQuads[12*i+11];
			crossProduct(&iMAP->bellNormals[12*i+9],v1,v2);

			i++;
		}
	}

	float *normalsTemp = new float[3*iMAP->bellDimensions];
	for (int w = 0; w < 3*iMAP->bellDimensions; w++) { normalsTemp[w] = 0.0; }
	int num = 0;
	for (int r = 0; r < 4*(dimension-1)*(dimension-1); r++) {
		const float x = iMAP->bellQuads[3*r+0];
		const float y = iMAP->bellQuads[3*r+1];
		num = 0;
		for (int q = 0; q < 4*(dimension-1)*(dimension-1); q++) {
			if (fabs(iMAP->bellQuads[3*q+0]-x) < 0.001 && fabs(iMAP->bellQuads[3*q+1]-y) < 0.001) {
				normalsTemp[3*r+0] += iMAP->bellNormals[3*q+0];
				normalsTemp[3*r+1] += iMAP->bellNormals[3*q+1];
				normalsTemp[3*r+2] += iMAP->bellNormals[3*q+2];
				num++;
			}
		}
		normalsTemp[3*r+0] /= (float)num;
		normalsTemp[3*r+1] /= (float)num;
		normalsTemp[3*r+2] /= (float)num;
	}

	for (int e = 0; e < 3*iMAP->bellDimensions; e++) {
		iMAP->bellNormals[e] = normalsTemp[e];
	}

	for (int h = 0; h < dimension; h++) {
		delete [] bell[h];
	}

	delete [] bell;
	delete [] normalsTemp;

	iMAP->gaussianBell = true;
}

// initialize OpenGL display parameters
void iMAPGlWindow::GlInit() {

	// Load png file for drawParticle
	glEnable (GL_DEPTH_TEST);

	glShadeModel (GL_SMOOTH);

	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glGenTextures(1, &iMAP->gaussianTex);
	glBindTexture (GL_TEXTURE_2D, iMAP->gaussianTex);
	glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glPixelStorei (GL_UNPACK_ALIGNMENT, 1);

	float sigma = 300;
	float dim = 100;

	float* gaussian = detectionImage(sigma,dim);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, dim, dim, 0, GL_ALPHA, GL_FLOAT, gaussian);
	delete [] gaussian;

	glutInitDisplayMode(GLUT_ALPHA | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO);

	// Clear The Screen And The Depth Buffer
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	glClearColor(iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0.0f);

	gluLookAt(	0.0, 0.0, 20,	/* eye is at (0, 0, 5) */
			0.0, 0.0, 0.0,	/* center is at (0, 0, 0) */
			0.0, 1.0, 0.0);	/* up is in +Y direction */

}

void iMAPGlWindow::TimerCallback(void *userdata) {
	iMAPGlWindow *mygl = (iMAPGlWindow*)userdata;

	if (iMAP->fileNumber > 0) {
		if (file->animationIndex < (int)((file->tMax-file->tMin)/file->exposureTime)) {
			file->animationIndex++;
			if (iMAP->synchronizeTiff) {
				iMAP->tiffSequence->setSlice(file->animationIndex);
				iMAP->overlayTiffGui->sequenceSlider->value(file->animationIndex);
			}
		}
		else {
			file->animationIndex = 0;
		}
	}

	// file needs to be loaded for these functions
	if (iMAP->fileNumber > 0) {
		// movie maker routine
		if ( (iMAP->movieMakerPlayEnable == true) && (file->fileLoaded == true) ) {
			if (iMAP->currentRow <= iMAP->movieMakerGui->table->getRows()) {
				if (iMAP->framesArray[iMAP->currentRow] > iMAP->framesPlayedArray[iMAP->currentRow]) {
					movieMakerOperation(iMAP->operationArray[iMAP->currentRow],iMAP->axisArray[iMAP->currentRow],iMAP->incrementArray[iMAP->currentRow]);
					iMAP->framesPlayedArray[iMAP->currentRow]++;
					iMAP->movieMakerGui->table->highlightRow(iMAP->currentRow);
					if (iMAP->movieMakerRecordEnable == true) {
						saveScreenSequence();
					}
				}
				else {
					iMAP->framesPlayedArray[iMAP->currentRow] = 0;
					iMAP->currentRow++;
				}
				// for overlap animation (makes for smoother transition between steps)
				if ( ((iMAP->currentRow+1) <= iMAP->movieMakerGui->table->getRows()) && (iMAP->framesArray[iMAP->currentRow]-iMAP->framesPlayedArray[iMAP->currentRow] < (int)( (iMAP->overlapArray[iMAP->currentRow+1])*(float)iMAP->framesArray[iMAP->currentRow+1]) ) ) {
					if (iMAP->framesPlayedArray[iMAP->currentRow+1] < iMAP->framesArray[iMAP->currentRow+1]) {
						movieMakerOperation(iMAP->operationArray[iMAP->currentRow+1],iMAP->axisArray[iMAP->currentRow+1],iMAP->incrementArray[iMAP->currentRow+1]);
						iMAP->framesPlayedArray[iMAP->currentRow+1]++;
					}
				}
			}
			// handle end of movie
			else {
				iMAP->movieMakerGui->table->highlightRow(0);
				// have to run callback twice for whatever reason
				iMAP->movieMakerGui->stopButton->do_callback();
				iMAP->movieMakerGui->stopButton->do_callback();
			}
		}
	}

	mygl->redraw();
	Fl::repeat_timeout(1.0f/iMAP->fpsSlider->value(),TimerCallback,mygl);
}

void iMAPGlWindow::draw() {
	// for stereo vision
	if (iMAP->stereoVision == true) {
		// left screen
		drawMaster();
		// right screen
		drawSlave();
	}
	// for normal vision
	else {
		drawMaster();
	}
}

void iMAPGlWindow::drawMaster() {

	// Initialize/handle reshaped viewport
	if ( !valid() ) {
		valid(1);
		GlInit();
		glViewport(0,0,w(),h());
	}

	if (iMAP->fileNumber > 0) {
		glViewport(0,0,w(),h());

		/********** 2D Trajectories **********/
		if (file->landscapePlot == false && file->file3D == false) {
			glLoadIdentity();
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glOrtho(-file->orthoLimit, file->orthoLimit, -file->orthoLimit, file->orthoLimit, 1.0f, 50.0f);
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			gluLookAt(0.0f, 0.0f, 20.0f,
					  0.0f, 0.0f, 0.0f,
					  0.0f, 1.0f, 0.0f);
			glClearColor(iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0.0f);
			glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

			// overlay TIFF
			if (iMAP->tiffOverlay) { iMAP->tiffSequence->draw(0); }
			// draw bounding box
			if (file->boundingBoxEnable) { drawBoundingBox(); }
			// draw dimensions
			if (file->dimensionsEnable) { drawDimensions(); }
			// draw grid
			if (file->gridEnable) {	drawGrid(); }
			// draw ticks
			if (file->ticksEnable) { drawTicks(); }
			// draw trajectory
			if (file->drawTrajectories) { drawTrajectories(); }
			// Select plot type
			if (file->drawLocalizations) {
				switch(file->plotType) {
					case 0:
						drawFrame();
						break;
					case 1:
						drawIntensity();
						break;
				}
			}
			// animate trajectory
			if (file->animateTrajectories) { animateTrajectories(); }

			// draw single trajectory overlay
			if (iMAP->singleTrajectoryInferenceMode && !file->animateTrajectories) {
				iMAP->singleTrajectoryInferenceGui->draw();
			}

			// mesh overlay graphics
			if (file->squareMeshOverlay) {
				drawSquareMesh();
				// draw posteriori highlight
				drawSquarePosteriorHighlight();
				if (file->drawRandomizedOptimizationZones) {
					drawSquareRandomizedOptimizationZones();
				}
			}

			if (file->voronoiMeshOverlay) {
				drawVoronoiMesh();
				drawVoronoiPosteriorHighlight();
				if (file->drawRandomizedOptimizationZones) {
					drawVoronoiRandomizedOptimizationZones();
				}
			}

			if (file->treeMeshOverlay) {
				drawQuadTreeMesh();
				drawQuadtreePosteriorHighlight();

				if (file->drawRandomizedOptimizationZones) {
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					drawQuadTreeRandomizedOptimizationZones();
					glDisable(GL_BLEND);
				}

				// spot visualization
				if (file->treeMeshGui->spotVisualizationButton->value()) {
					glEnable(GL_POINT_SPRITE);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

					glEnable(GL_TEXTURE_2D);
					glBindTexture(GL_TEXTURE_2D, iMAP->gaussianTex);
					glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);

					glPushMatrix();
						updateView();
						drawQuadTreeSpotVisualization(file->treeMesh->quadTree);
					glPopMatrix();

					glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_FALSE);

					glBindTexture(GL_TEXTURE_2D,0);
					glDisable(GL_TEXTURE_2D);
					glDisable(GL_POINT_SPRITE);
					glDisable(GL_BLEND);
				}
			}

			glEnable(GL_BLEND);
			glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			if (iMAP->selectionButtonPressed) {
				if (iMAP->selectionMade) {
					Node *temp = file->selection->head->next;
					if (temp != NULL) {
						glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f);
						glPushMatrix();
//							updateView();
							glLineWidth(1);
							glTranslatef(file->xTranslate,file->yTranslate,0.0);
							glBegin(GL_LINE_STRIP);
								glVertex2d(file->selection->head->xCoord,file->selection->head->yCoord);
								while (temp->next != NULL) {
									glVertex2d(temp->xCoord,temp->yCoord);
									temp = temp->next;
								}
							glEnd();
						glPopMatrix();
					}
				}
				else if (file->selection != NULL) {
					Node *temp = file->selection->head->next;
					if (temp != NULL) {

						glEnableClientState(GL_VERTEX_ARRAY);

						if (iMAP->customSelectionInferenceGui->DDrinferred) { glColor4f(0.0,1.0,0.0,0.8f); }
						else if (iMAP->customSelectionInferenceGui->DFinferred) { glColor4f(1.0,1.0,0.0,0.8f); }
						else { glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f); }

						glPushMatrix();
							updateView();
							glLineWidth(2);
							glVertexPointer(2,GL_DOUBLE,0,file->selection->vertices);
							glDrawArrays(GL_LINE_LOOP,0,file->selection->n-1);
						glPopMatrix();

						glDisableClientState(GL_VERTEX_ARRAY);

						if (iMAP->customSelectionInferenceGui->DFinferred) {
							glColor4f(1.0,1.0,1.0,0.4);
							glPushMatrix();
								updateView();
								file->selection->draw();
								file->selection->drawForce();
							glPopMatrix();
						}
						else if (iMAP->customSelectionInferenceGui->DDrinferred) {
							glColor4f(1.0,1.0,1.0,0.4);
							glPushMatrix();
								updateView();
								file->selection->draw();
								file->selection->drawForce();
							glPopMatrix();
						}

					}
				}
			}
			else {
				x1 = (float) pushPosition[0];
				y1 = (float) pushPosition[1];
				x2 = (float) releasePosition[0];
				y2 = (float) releasePosition[1];

				iMAP->zoomBox[0] = x1;
				iMAP->zoomBox[1] = x2;
				iMAP->zoomBox[2] = y1;
				iMAP->zoomBox[3] = y2;

				glPushMatrix();
					glLineWidth(1);
					glBegin(GL_LINE_LOOP);
						glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f);
						glVertex2f(x1,-y1);
						glVertex2f(x1,-y2);
						glVertex2f(x2,-y2);
						glVertex2f(x2,-y1);
//						glVertex2f(x1-file->xTranslate,-y1-file->yTranslate);
//						glVertex2f(x1-file->xTranslate,-y2-file->yTranslate);
//						glVertex2f(x2-file->xTranslate,-y2-file->yTranslate);
//						glVertex2f(x2-file->xTranslate,-y1-file->yTranslate);
					glEnd();
				glPopMatrix();
				// annotation
				drawZoomBox(fabsf(x1-x2),fabsf(-y1+y2));
			}
			glDisable(GL_BLEND);
		}
		/********** SURFACE PLOT VIEWING FOR 2D MODE **********/
		else if (file->landscapePlot == true && file->file3D == false) {

			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(iMAP->fovy, 1.0, 1.0, 100.0);

			if (iMAP->stereoVision == true) {
				xyz right;
				xyz vd = {0.0,0.0,1.0};
				xyz vu = {0.0,1.0,0.0};

				CROSSPROD(vd,vu,right);
				Normalise(&right);

				const float eyeSeparation = iMAP->stereoVisionGui->getEyeSeparation();

				right.x *= eyeSeparation / 2.0;
				right.y *= eyeSeparation / 2.0;
				right.z *= eyeSeparation / 2.0;

				glMatrixMode(GL_MODELVIEW);
				glDrawBuffer(GL_BACK_LEFT);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				glLoadIdentity();

				gluLookAt(right.x,
						  right.y,
						  iMAP->initialLookAt+right.z,
						  0.0,0.0,0.0,
						  0.0,1.0,0);

				const int stereoCentreOffset = (h()-w()/2)/2;
				glViewport(0,stereoCentreOffset,w()/2,w()/2);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			} else {
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				gluLookAt(0,0, iMAP->initialLookAt,
						  0,0,0,
						  0.0, 1.0, 0.0);
				glViewport(0,0,w(),h());
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			}

			glClearColor(iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0.0f);

			// overlay TIFF
			if (iMAP->tiffOverlay) { iMAP->tiffSequence->draw(0); }
			// draw bounding box
			if (file->boundingBoxEnable) { drawBoundingBox(); }
			// draw dimensions
			if (file->dimensionsEnable) { drawDimensions(); }
			// draw grid
			if (file->gridEnable) {	drawGrid(); }
			// draw ticks
			if (file->ticksEnable) { drawTicks(); }
			// draw trajectory
			if (file->drawTrajectories) { drawTrajectories(); }
			// Select plot type
			if (file->drawLocalizations) {
				switch(file->plotType) {
					case 0:
						drawFrame();
						break;
					case 1:
						drawIntensity();
						break;
				}
			}
			// animate trajectory
			if (file->animateTrajectories) { animateTrajectories(); }

			// mesh overlay graphics
			if (file->squareMeshOverlay) {
				drawSquareMesh();
				if (file->landscapePlot) {
					drawSquareLandscape();
					drawLandscapeAxes(0);
				}
			}

			if (file->voronoiMeshOverlay) {
				drawVoronoiMesh();
				drawVoronoiPosteriorHighlight();
				if (file->landscapePlot) {
					drawVoronoiLandscape();
					drawLandscapeAxes(0);
				}
			}

			if (file->treeMeshOverlay) {
				drawQuadTreeMesh();
				if (file->landscapePlot) {
					drawQuadTreeLandscape();
					drawLandscapeAxes(0);
				}
				drawQuadtreePosteriorHighlight();

				// spot visualization
				if (file->treeMeshGui->spotVisualizationButton->value()) {
					glEnable(GL_POINT_SPRITE);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

					glEnable(GL_TEXTURE_2D);
					glBindTexture(GL_TEXTURE_2D, iMAP->gaussianTex);
					glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);

					glPushMatrix();
						updateView();
						drawQuadTreeSpotVisualization(file->treeMesh->quadTree);
					glPopMatrix();

					glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_FALSE);

					glBindTexture(GL_TEXTURE_2D,0);
					glDisable(GL_TEXTURE_2D);
					glDisable(GL_POINT_SPRITE);
					glDisable(GL_BLEND);
				}
			}
		}
		else if (file->file3D == true) {
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(iMAP->fovy, 1.0, 1.0, 100.0);

			if (iMAP->stereoVision == true) {
				xyz right;
				xyz vd = {0.0,0.0,1.0};
				xyz vu = {0.0,1.0,0.0};

				CROSSPROD(vd,vu,right);
				Normalise(&right);

				const float eyeSeparation = iMAP->stereoVisionGui->getEyeSeparation();

				right.x *= eyeSeparation / 2.0;
				right.y *= eyeSeparation / 2.0;
				right.z *= eyeSeparation / 2.0;

				glMatrixMode(GL_MODELVIEW);
				glDrawBuffer(GL_BACK_LEFT);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				glLoadIdentity();

				gluLookAt(right.x,
						  right.y,
						  iMAP->initialLookAt+right.z,
						  0.0,0.0,0.0,
						  0.0,1.0,0);

				const int stereoCentreOffset = (h()-w()/2)/2;
				glViewport(0,stereoCentreOffset,w()/2,w()/2);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			} else {
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				gluLookAt(0,0, iMAP->initialLookAt,
						  0,0,0,
						  0.0, 1.0, 0.0);
				glViewport(0,0,w(),h());
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			}

			glClearColor(iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0.0f);


			// draw bounding box
			if (file->boundingBoxEnable) { drawBoundingBox3D(); }
			// draw trajectory
			if (file->drawTrajectories) { drawTrajectories3D(); }
			// draw localizations
			if (file->drawLocalizations) { drawFrame(); }
			// animate trajectory
			if (file->animateTrajectories) { animateTrajectories3D(); }
			// draw axes
			drawAxes(0);
		}
	}
	else {
		defaultScreen();
	}
}

void iMAPGlWindow::drawSlave() {

	// Initialize/handle reshaped viewport
	if ( !valid() ) {
		valid(1);
		GlInit();
		glViewport(0,0,w(),h());
	}

	if (iMAP->fileNumber > 0) {
		glViewport(0,0,w(),h());

		/********** 2D Trajectories **********/
		if (file->landscapePlot == false && file->file3D == false) {
			glLoadIdentity();
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glOrtho(-file->orthoLimit, file->orthoLimit, -file->orthoLimit, file->orthoLimit, 1.0f, 50.0f);
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			gluLookAt(0.0f, 0.0f, 20.0f,
					  0.0f, 0.0f, 0.0f,
					  0.0f, 1.0f, 0.0f);
			glClearColor(iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0.0f);
			glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

			// overlay TIFF
			if (iMAP->tiffOverlay) { iMAP->tiffSequence->draw(0); }
			// draw bounding box
			if (file->boundingBoxEnable) { drawBoundingBox(); }
			// draw dimensions
			if (file->dimensionsEnable) { drawDimensions(); }
			// draw grid
			if (file->gridEnable) {	drawGrid(); }
			// draw ticks
			if (file->ticksEnable) { drawTicks(); }
			// draw trajectory
			if (file->drawTrajectories) { drawTrajectories(); }
			// Select plot type
			if (file->drawLocalizations) {
				switch(file->plotType) {
					case 0:
						drawFrame();
						break;
					case 1:
						drawIntensity();
						break;
				}
			}
			// animate trajectory
			if (file->animateTrajectories) { animateTrajectories(); }

			// draw single trajectory overlay
			if (iMAP->singleTrajectoryInferenceMode && !file->animateTrajectories) {
				iMAP->singleTrajectoryInferenceGui->draw();
			}

			// mesh overlay graphics
			if (file->squareMeshOverlay) {
				drawSquareMesh();
				// draw posteriori highlight
				drawSquarePosteriorHighlight();
				if (file->drawRandomizedOptimizationZones) {
					drawSquareRandomizedOptimizationZones();
				}
			}

			if (file->voronoiMeshOverlay) {
				drawVoronoiMesh();
				drawVoronoiPosteriorHighlight();
				if (file->drawRandomizedOptimizationZones) {
					drawVoronoiRandomizedOptimizationZones();
				}
			}

			if (file->treeMeshOverlay) {
				drawQuadTreeMesh();
				drawQuadtreePosteriorHighlight();

				if (file->drawRandomizedOptimizationZones) {
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					drawQuadTreeRandomizedOptimizationZones();
					glDisable(GL_BLEND);
				}

				// spot visualization
				if (file->treeMeshGui->spotVisualizationButton->value()) {
					glEnable(GL_POINT_SPRITE);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

					glEnable(GL_TEXTURE_2D);
					glBindTexture(GL_TEXTURE_2D, iMAP->gaussianTex);
					glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);

					glPushMatrix();
						updateView();
						drawQuadTreeSpotVisualization(file->treeMesh->quadTree);
					glPopMatrix();

					glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_FALSE);

					glBindTexture(GL_TEXTURE_2D,0);
					glDisable(GL_TEXTURE_2D);
					glDisable(GL_POINT_SPRITE);
					glDisable(GL_BLEND);
				}
			}

			glEnable(GL_BLEND);
			glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			if (iMAP->selectionButtonPressed) {
				if (iMAP->selectionMade) {
					Node *temp = file->selection->head->next;
					if (temp != NULL) {
						glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f);
						glPushMatrix();
//							updateView();
							glLineWidth(1);
							glTranslatef(file->xTranslate,file->yTranslate,0.0);
							glBegin(GL_LINE_STRIP);
								glVertex2d(file->selection->head->xCoord,file->selection->head->yCoord);
								while (temp->next != NULL) {
									glVertex2d(temp->xCoord,temp->yCoord);
									temp = temp->next;
								}
							glEnd();
						glPopMatrix();
					}
				}
				else if (file->selection != NULL) {
					Node *temp = file->selection->head->next;
					if (temp != NULL) {

						glEnableClientState(GL_VERTEX_ARRAY);
						glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f);
						glPushMatrix();
							updateView();
							glLineWidth(1);
							glVertexPointer(2,GL_DOUBLE,0,file->selection->vertices);
							glDrawArrays(GL_LINE_LOOP,0,file->selection->n-1);
						glPopMatrix();

						glDisableClientState(GL_VERTEX_ARRAY);

						if (iMAP->customSelectionInferenceGui->DFinferred) {
							glColor4f(0.2,1.0,0.2,0.3);
							glPushMatrix();
								updateView();
								file->selection->draw();
								glColor4f(1.0,1.0,0.0,0.9);
								file->selection->drawForce();
							glPopMatrix();
						}
					}
				}
			}
			else {
				x1 = (float) pushPosition[0];
				y1 = (float) pushPosition[1];
				x2 = (float) releasePosition[0];
				y2 = (float) releasePosition[1];

				iMAP->zoomBox[0] = x1;
				iMAP->zoomBox[1] = x2;
				iMAP->zoomBox[2] = y1;
				iMAP->zoomBox[3] = y2;

				glPushMatrix();
					glLineWidth(1);
					glBegin(GL_LINE_LOOP);
						glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f);
						glVertex2f(x1,-y1);
						glVertex2f(x1,-y2);
						glVertex2f(x2,-y2);
						glVertex2f(x2,-y1);
//						glVertex2f(x1-file->xTranslate,-y1-file->yTranslate);
//						glVertex2f(x1-file->xTranslate,-y2-file->yTranslate);
//						glVertex2f(x2-file->xTranslate,-y2-file->yTranslate);
//						glVertex2f(x2-file->xTranslate,-y1-file->yTranslate);
					glEnd();
				glPopMatrix();
				// annotation
				drawZoomBox(fabsf(x1-x2),fabsf(-y1+y2));
			}
			glDisable(GL_BLEND);
		}
		/********** SURFACE PLOT VIEWING FOR 2D MODE **********/
		else if (file->landscapePlot == true && file->file3D == false) {

			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(iMAP->fovy, 1.0, 1, 100);

			xyz right;
			xyz vd = {0.0,0.0,1.0};
			xyz vu = {0.0,1.0,0.0};

			CROSSPROD(vd,vu,right);
			Normalise(&right);

			const float eyeSeparation = iMAP->stereoVisionGui->getEyeSeparation();

			right.x *= eyeSeparation / 2.0;
			right.y *= eyeSeparation / 2.0;
			right.z *= eyeSeparation / 2.0;

			glMatrixMode(GL_MODELVIEW);
			glDrawBuffer(GL_BACK_RIGHT);
//				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glLoadIdentity();
			gluLookAt(-right.x,
					  -right.y,
					  iMAP->initialLookAt-right.z,
					  0.0,0.0,0.0,
					  0.0,1.0,0);

			glClearColor(iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0.0f);

			const int stereoCentreOffset = (h()-w()/2)/2;
	//		glViewport(0,stereoCentreOffset,ViSP->screenDimensions[2]/2,ViSP->screenDimensions[2]/2);
			glViewport(w()/2,stereoCentreOffset,w()/2,w()/2);

			// overlay TIFF
			if (iMAP->tiffOverlay) { iMAP->tiffSequence->draw(0); }
			// draw bounding box
			if (file->boundingBoxEnable) { drawBoundingBox(); }
			// draw dimensions
			if (file->dimensionsEnable) { drawDimensions(); }
			// draw grid
			if (file->gridEnable) {	drawGrid(); }
			// draw ticks
			if (file->ticksEnable) { drawTicks(); }
			// draw trajectory
			if (file->drawTrajectories) { drawTrajectories(); }
			// Select plot type
			if (file->drawLocalizations) {
				switch(file->plotType) {
					case 0:
						drawFrame();
						break;
					case 1:
						drawIntensity();
						break;
				}
			}
			// animate trajectory
			if (file->animateTrajectories) { animateTrajectories(); }

			// mesh overlay graphics
			if (file->squareMeshOverlay) {
				drawSquareMesh();
				if (file->landscapePlot) {
					drawSquareLandscape();
					drawLandscapeAxes(1);
				}
			}

			if (file->voronoiMeshOverlay) {
				drawVoronoiMesh();
				drawVoronoiPosteriorHighlight();
				if (file->landscapePlot) {
					drawVoronoiLandscape();
					drawLandscapeAxes(1);
				}
			}

			if (file->treeMeshOverlay) {
				drawQuadTreeMesh();
				if (file->landscapePlot) {
					drawQuadTreeLandscape();
					drawLandscapeAxes(1);
				}
				drawQuadtreePosteriorHighlight();

				// spot visualization
				if (file->treeMeshGui->spotVisualizationButton->value()) {
					glEnable(GL_POINT_SPRITE);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

					glEnable(GL_TEXTURE_2D);
					glBindTexture(GL_TEXTURE_2D, iMAP->gaussianTex);
					glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);

					glPushMatrix();
						updateView();
						drawQuadTreeSpotVisualization(file->treeMesh->quadTree);
					glPopMatrix();

					glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_FALSE);

					glBindTexture(GL_TEXTURE_2D,0);
					glDisable(GL_TEXTURE_2D);
					glDisable(GL_POINT_SPRITE);
					glDisable(GL_BLEND);
				}

			}
		}
		else if (file->file3D == true) {
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(iMAP->fovy, 1.0, 1, 100);
			if (iMAP->stereoVision == true) {

				xyz right;
				xyz vd = {0.0,0.0,1.0};
				xyz vu = {0.0,1.0,0.0};

				CROSSPROD(vd,vu,right);
				Normalise(&right);

				const float eyeSeparation = iMAP->stereoVisionGui->getEyeSeparation();

				right.x *= eyeSeparation / 2.0;
				right.y *= eyeSeparation / 2.0;
				right.z *= eyeSeparation / 2.0;

				glMatrixMode(GL_MODELVIEW);
				glDrawBuffer(GL_BACK_RIGHT);
	//				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				glLoadIdentity();
				gluLookAt(-right.x,
						  -right.y,
						  iMAP->initialLookAt-right.z,
						  0.0,0.0,0.0,
						  0.0,1.0,0);

				glClearColor(iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0.0f);

				const int stereoCentreOffset = (h()-w()/2)/2;
		//		glViewport(0,stereoCentreOffset,ViSP->screenDimensions[2]/2,ViSP->screenDimensions[2]/2);
				glViewport(w()/2,stereoCentreOffset,w()/2,w()/2);

			} else {
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				gluLookAt(0,0, iMAP->initialLookAt,
						  0,0,0,
						  0.0, 1.0, 0.0);
				glViewport(0,0,w(),h());
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			}

			glClearColor(iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0.0f);

			// draw bounding box
			if (file->boundingBoxEnable) { drawBoundingBox3D(); }
			// draw trajectory
			if (file->drawTrajectories) { drawTrajectories3D(); }
			// draw localizations
			if (file->drawLocalizations) { drawFrame(); }
			// animate trajectory
			if (file->animateTrajectories) { animateTrajectories3D(); }
			// draw axes
			drawAxes(1);
		}
	} else {
		defaultScreen();
	}
}

int iMAPGlWindow::handle(int e) {

	if (iMAP->fileNumber > 0) {
		int ret = Fl_Gl_Window::handle(e);

		if (file->landscapePlot == false && file->file3D == false) {
			// for removing crop/delete box after region has been cropped
			if (boxChanged == true) {
				boxChanged = false;
				pushPosition[0] = pushPosition[1] = 0;
				releasePosition[0] = releasePosition[1] = 0;
			}

			switch (e) {
				case FL_MOUSEWHEEL:
					return(1);
				case FL_PUSH:
					// create selection box
					glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
					glGetIntegerv(GL_VIEWPORT,viewport);
					// right click to select/deselect cells
					if (file->squareMeshApplied()) {
						float xMin,xMax,yMin,yMax;
						if (file->squareMesh->selectionMode()) {
							xMin = (float)file->squareMesh->selection.xMin;//-file->xTranslate;
							xMax = (float)file->squareMesh->selection.xMax;//-file->xTranslate;
							yMin = (float)file->squareMesh->selection.yMin;//-file->yTranslate;
							yMax = (float)file->squareMesh->selection.yMax;//-file->yTranslate;
						} else {
							xMin = (float)file->xMin;//-file->xTranslate;
							xMax = (float)file->xMax;//-file->xTranslate;
							yMin = (float)file->yMin;//-file->yTranslate;
							yMax = (float)file->yMax;//-file->yTranslate;
						}
						if (Fl::event_button() == FL_RIGHT_MOUSE) {
							gluUnProject(
								Fl::event_x(),
								viewport[3]-Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							const int x = (int)( (cellSelect[0] + fabsf(xMin))/file->squareMesh->getDx() );
							const int y = (int)( (cellSelect[1] + fabsf(yMin))/file->squareMesh->getDy() );

							// activate or deactivate cell selected
							if (x <= file->squareMesh->getXCells() && y <= file->squareMesh->getYCells() && x >= 0 && y >= 0) {
								if (file->squareMesh->active(x,y)) { file->squareMesh->deactivateCell(x,y); }
								else { file->squareMesh->activateCell(x,y); }
								file->squareMesh->updateNeighbours();
								file->squareMeshGui->updateVariables(iMAP->activeZones);
							}
							iMAP->overlayAdjustment = true;
						}
						if (file->squareMeshGui->posterioriGroup->active() && file->squareMeshGui->posterioriGroup->visible()) {
							gluUnProject(
								Fl::event_x(),
								viewport[3]-Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							const int x = (int)( (cellSelect[0] + fabsf(xMin))/file->squareMesh->getDx() );
							const int y = (int)( (cellSelect[1] + fabsf(yMin))/file->squareMesh->getDy() );

							if (file->squareMeshGui->potentialReferenceButton->value()) {
								if (file->squareMesh->active(x,y)) {
									file->squareMesh->setCurrentZone2(x,y);

									const double min = 0.0;
									const double max = 20.0;
									file->squareMeshGui->minPosteriorSlider->value(0.0);
									file->squareMeshGui->minPosteriorSlider->bounds(min,max);
									file->squareMeshGui->minPosteriorSlider->deactivate();
									file->squareMeshGui->maxPosteriorSlider->value(max/2.0);
									file->squareMeshGui->maxPosteriorSlider->bounds(min,max);

								}
							} else {
								if (file->squareMesh->active(x,y)) {
									file->squareMesh->setCurrentZone(x,y);
								}
								if (file->squareMeshGui->dPosterioriButton->value()) {
									file->squareMeshGui->dPosterioriButton->do_callback();
								}
								else if (file->squareMeshGui->fPosteriorButton->value()) {
									file->squareMeshGui->fPosteriorButton->do_callback();
								}
								else if (file->squareMeshGui->vPosteriorButton->value()) {
									file->squareMeshGui->vPosteriorButton->do_callback();
								}
							}
						}
					}

					// right click to select/deselect cells
					if (file->voronoiMeshApplied()) {
						if (Fl::event_button() == FL_RIGHT_MOUSE) {
							gluUnProject(
								Fl::event_x(),
								viewport[3]-Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							// find which cell was clicked
							const int selection = file->voronoiMesh->getCellSelection(cellSelect[0],cellSelect[1]);

							// activate or deactivate cell selected
							if (file->voronoiMesh->active(selection)) { file->voronoiMesh->deactivateCellClick(selection); }
							else { file->voronoiMesh->activateCellClick(selection); }
							file->voronoiMesh->updateNeighbours();

						}
						if (file->voronoiMeshGui->posterioriGroup->active() && file->voronoiMeshGui->posterioriGroup->visible()) {
							gluUnProject(
								Fl::event_x(),
								viewport[3]-Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							const int selection = file->voronoiMesh->getCellSelection(cellSelect[0],cellSelect[1]);

							if (file->voronoiMeshGui->potentialReferenceButton->value()) {
								if (file->voronoiMesh->active(selection)) {
									file->voronoiMesh->setCurrentZone2(selection);
								}
							} else {
								file->voronoiMesh->setCurrentZone(selection);
								if (file->voronoiMeshGui->dPosterioriButton->value()) {
									file->voronoiMeshGui->dPosterioriButton->do_callback();
								}
								else if (file->voronoiMeshGui->fPosteriorButton->value()) {
									file->voronoiMeshGui->fPosteriorButton->do_callback();
								}
								else if (file->voronoiMeshGui->vPosteriorButton->value()) {
									file->voronoiMeshGui->vPosteriorButton->do_callback();
								}
							}
						}
					}

					// right click to select/deselect cells
					if (file->treeMeshApplied()) {
						float xMin,xMax,yMin,yMax;
						if (file->treeMesh->selectionMode()) {
							xMin = (float)file->treeMesh->selection.xMin;
							xMax = (float)file->treeMesh->selection.xMax;
							yMin = (float)file->treeMesh->selection.yMin;
							yMax = (float)file->treeMesh->selection.yMax;
						} else {
							xMin = (float)file->xMin;
							xMax = (float)file->xMax;
							yMin = (float)file->yMin;
							yMax = (float)file->yMax;
						}
						if (Fl::event_button() == FL_RIGHT_MOUSE) {
							gluUnProject(
								Fl::event_x(),
								viewport[3]-Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							// find which cell was clicked
							file->treeMesh->selectedQuadLeaf = 0;
							bool found = false;
							file->treeMesh->getLeafSelection(cellSelect[0],cellSelect[1],file->treeMesh->quadTree,&found);

							if (file->treeMesh->selectedQuadLeaf != 0) {
								file->treeMesh->leafSelected = true;
								// activate or deactivate cell selected
								if (file->treeMesh->selectedQuadLeaf->active()) {
									file->treeMesh->selectedQuadLeaf->deactivate();
									file->treeMesh->totalVariables--;
								}
								else {
									file->treeMesh->selectedQuadLeaf->activate();
									file->treeMesh->totalVariables++;
								}
								file->treeMeshGui->updateVariables(file->treeMesh->totalVariables);
								file->treeMesh->updateNeighbours(file->treeMesh->quadTree);
							} else {
								file->treeMesh->leafSelected = false;
							}
							iMAP->overlayAdjustment = true;
						}
						if (file->treeMeshGui->posterioriGroup->visible() && file->treeMeshGui->posterioriGroup->active()) {
							gluUnProject(
								Fl::event_x(),
								viewport[3]-Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							// find which cell was clicked

							if (file->treeMeshGui->potentialReferenceButton->value()) {
//								if (file->treeMesh->selectedQuadLeaf2 != 0) {
									bool found = false;
									file->treeMesh->getLeafSelection2(cellSelect[0],cellSelect[1],file->treeMesh->quadTree,&found);
									file->treeMesh->leaf2Selected = true;
//								}
							} else {
								if (file->treeMesh->selectedQuadLeaf != 0) {
									bool found = false;
									file->treeMesh->getLeafSelection(cellSelect[0],cellSelect[1],file->treeMesh->quadTree,&found);
									file->treeMesh->leafSelected = true;
									if (file->treeMeshGui->dPosteriorButton->value()) {
										file->treeMeshGui->dPosteriorButton->do_callback();
									}
									else if (file->treeMeshGui->fPosteriorButton->value()) {
										file->treeMeshGui->fPosteriorButton->do_callback();
									}
									else if (file->treeMeshGui->vPosteriorButton->value()) {
										file->treeMeshGui->vPosteriorButton->do_callback();
									}
								}
							}
						}
					}

					if (Fl::event_button() == FL_LEFT_MOUSE) {
						if (iMAP->selectionButtonPressed) {
							iMAP->selectionMade = true;
							gluUnProject(
								Fl::event_x(),
								Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							cellSelect[0] += (double)file->xTranslate;
							cellSelect[1] += (double)file->yTranslate;

							if (file->selection != NULL) {
								delete file->selection;
								file->selection = NULL;
							}

							file->selection = new Selection(cellSelect[0],-cellSelect[1]);
							iMAP->customSelectionInferenceGui->deactivate();

//	//						printf("%f\t%f\n",file->selection->head->xCoord,file->selection->head->yCoord);
						}
						else {
							// assign x1,x2 coords
							gluUnProject(
								Fl::event_x(),
								Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&pushPosition[0],
								&pushPosition[1],
								&pushPosition[2]);

							pushPosition[0] += file->xTranslate;
							pushPosition[1] += file->yTranslate;

							releasePosition[0] = pushPosition[0];
							releasePosition[1] = pushPosition[1];
							releasePosition[2] = pushPosition[2];
						}
					}

					// reset view on double click
					if (Fl::event_clicks() == 1 && Fl::event_button() == FL_LEFT_MOUSE) {
						doubleClick = true;
						file->xTranslate = 0.0f;
						file->yTranslate = 0.0f;
						file->orthoLimit = 0.6*file->maxRange;
					}

					return(1);
				case FL_DRAG:				// handle mouse drag
					if (Fl::event_button() == FL_LEFT_MOUSE) {
						if (iMAP->selectionButtonPressed) {
							gluUnProject(
								Fl::event_x(),
								Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&releasePosition[0],
								&releasePosition[1],
								&releasePosition[2]);

							releasePosition[0] += (double)file->xTranslate;
							releasePosition[1] += (double)file->yTranslate;

							// only add node if greater than 3 nm away
							if ( file->selection->farEnough(releasePosition[0],-releasePosition[1]) ) {
								file->selection->add(releasePosition[0],-releasePosition[1]);
							}
							dragging = true;
						}
						else {// assign x2,y2 in realtime
							gluUnProject(
								Fl::event_x(),
								Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&releasePosition[0],
								&releasePosition[1],
								&releasePosition[2]);

							releasePosition[0] += file->xTranslate;
							releasePosition[1] += file->yTranslate;

							lastX = currentX;
							lastY = currentY;

							dragging = true;
						}
						iMAP->zoomBoxEnable = true;
						return(1);
					}
					dragging = false;
					return(1);
				case FL_RELEASE:			// handle mouse release

					if (iMAP->selectionButtonPressed) {
						iMAP->selectionMade = false;
						iMAP->customSelectionInferenceGui->DFinferred = false;
						iMAP->customSelectionInferenceGui->DDrinferred = false;
						file->selection->assignPoints();
						if (file->selection->cell.count < 2) {
							fl_alert("Not enough points for selection.");
							iMAP->customSelectionInferenceGui->selectionButton->value(0);
							iMAP->customSelectionInferenceGui->selectionButton->do_callback();
							delete file->selection;
							file->selection = NULL;
							iMAP->selectionButtonPressed = false;
							iMAP->customSelectionInferenceGui->inferDFButton->deactivate();
							iMAP->customSelectionInferenceGui->inferDDrButton->deactivate();
						} 
						else {
							iMAP->customSelectionInferenceGui->inferDFButton->activate();
							iMAP->customSelectionInferenceGui->inferDDrButton->activate();
						}
						iMAP->customSelectionInferenceGui->activate();
					}

					if (doubleClick) { doubleClick = false; }
					else if (dragging) {
						const float w = (float)fabs(iMAP->zoomBox[1] - iMAP->zoomBox[0]);
						const float h = (float)fabs(iMAP->zoomBox[3] - iMAP->zoomBox[2]);

						const float xMid = (float)(iMAP->zoomBox[1] + iMAP->zoomBox[0])/2.0f;
						const float yMid = (float)(iMAP->zoomBox[3] + iMAP->zoomBox[2])/2.0f;

						if (w != 0.0f || h != 0.0f) {
							file->xTranslate += -xMid;
							file->yTranslate += yMid;
							if (w >= h) { file->orthoLimit = 0.6*w; }
							else { file->orthoLimit = 0.6*h; }
						}

						iMAP->zoomBoxEnable = false;

						boxChanged = true;
						dragging = false;

						//iMAP->overlayAdjustment = true;
					}
					// square meshing
					else if (file->meshType == 0) {
						if (file->inferred && file->squareMeshGui->infoOverlayButton->value() && Fl::event_button() == FL_LEFT_MOUSE) {
							gluUnProject(
								Fl::event_x(),
								viewport[3]-Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							float xMin,xMax,yMin,yMax;
							if (file->squareMesh->selectionMode()) {
								xMin = (float)file->squareMesh->selection.xMin;
								xMax = (float)file->squareMesh->selection.xMax;
								yMin = (float)file->squareMesh->selection.yMin;
								yMax = (float)file->squareMesh->selection.yMax;
							} else {
								xMin = (float)file->xMin;
								xMax = (float)file->xMax;
								yMin = (float)file->yMin;
								yMax = (float)file->yMax;
							}

							const int x = (int)( (cellSelect[0] + fabsf(xMin))/file->squareMesh->getDx() );
							const int y = (int)( (cellSelect[1] + fabsf(yMin))/file->squareMesh->getDy() );

							if (x < file->squareMesh->getXCells() && y < file->squareMesh->getYCells() && x >= 0 && y >= 0) {
								if (file->squareMesh->getCell(x,y)->active()) {
									char clusterString [50];
									char localizationsString [50];
									char diffusionString [50];
									char potentialString [50];
									char forceXString [50];
									char forceYString [50];
									char forceMagString [50];

									if (file->squareMesh->active(x,y)) {
										sprintf(clusterString,"Zone (%i,%i)",x,y);
										sprintf(localizationsString,"%i Localizations",file->squareMesh->getCell(x,y)->getCount());
										sprintf(diffusionString,"D = %f [um\262/s]",file->squareMesh->getCell(x,y)->getDiffusion());
										sprintf(potentialString,"V = %f [kT]",file->squareMesh->getCell(x,y)->getPotential());
										sprintf(forceXString,"Fx = %f [pN]",file->squareMesh->getCell(x,y)->getForceX());
										sprintf(forceYString,"Fy = %f [pN]",file->squareMesh->getCell(x,y)->getForceY());
										sprintf(forceMagString,"||F|| = %f [pN]",file->squareMesh->getCell(x,y)->getForceMagnitude());

										Fl_Menu_Item lclick_menu[] = {
											{clusterString,0,nullCallback,0,FL_MENU_DIVIDER},
											{localizationsString, 0, nullCallback},
											{diffusionString, 0, nullCallback},
											{potentialString, 0,  nullCallback},
											{forceXString, 0, nullCallback},
											{forceYString, 0, nullCallback},
											{forceMagString, 0, nullCallback},
											{0}
										};

										const Fl_Menu_Item *m = lclick_menu->popup(Fl::event_x(), Fl::event_y(), 0, 0, 0);
										if ( m ) m->do_callback(0, m->user_data());
									} else {
										sprintf(clusterString,"Zone (%i,%i)",x,y);
										sprintf(localizationsString,"%i Localizations",file->squareMesh->getCell(x,y)->getCount());

										Fl_Menu_Item lclick_menu[] = {
											{clusterString,0,nullCallback,0,FL_MENU_DIVIDER},
											{localizationsString, 0, nullCallback},
											{0}
										};

										const Fl_Menu_Item *m = lclick_menu->popup(Fl::event_x(), Fl::event_y(), 0, 0, 0);
										if ( m ) m->do_callback(0, m->user_data());
									}
								}
							}
						}
					}
					// voronoi tessellation
					else if (file->meshType == 1) {
						if (file->inferred && file->voronoiMeshGui->infoOverlayButton->value() && Fl::event_button() == FL_LEFT_MOUSE) {
							gluUnProject(
								Fl::event_x(),
								viewport[3]-Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							const int selection = file->voronoiMesh->getCellSelection(cellSelect[0],cellSelect[1]);

							if (file->voronoiMesh->getCell(selection)->active()) {
								char clusterString [50];
								char localizationsString [50];
								char diffusionString [50];
								char potentialString [50];
								char forceXString [50];
								char forceYString [50];
								char forceMagString [50];

								sprintf(clusterString,"Cluster %i",selection);
								sprintf(localizationsString,"%i Localizations",file->voronoiMesh->getCell(selection)->getCount());
								sprintf(diffusionString,"D = %f [um\262/s]",file->voronoiMesh->getCell(selection)->getDiffusion());
								sprintf(potentialString,"V = %f [kT]",file->voronoiMesh->getCell(selection)->getPotential());
								sprintf(forceXString,"Fx = %f [a.u.]",file->voronoiMesh->getCell(selection)->getForceX());
								sprintf(forceYString,"Fy = %f [a.u.]",file->voronoiMesh->getCell(selection)->getForceY());
								sprintf(forceMagString,"||F|| = %f [a.u]",file->voronoiMesh->getCell(selection)->getForceMagnitude());

								Fl_Menu_Item lclick_menu[] = {
									{clusterString,0,nullCallback,0,FL_MENU_DIVIDER},
									{localizationsString, 0, nullCallback},
									{diffusionString, 0, nullCallback},
									{potentialString, 0,  nullCallback},
									{forceXString, 0, nullCallback},
									{forceYString, 0, nullCallback},
									{forceMagString, 0, nullCallback},
									{0}
								};

								const Fl_Menu_Item *m = lclick_menu->popup(Fl::event_x(), Fl::event_y(), 0, 0, 0);
								if ( m ) m->do_callback(0, m->user_data());
							}
						}
					}
					// quad-tree meshing
					else if (file->meshType == 2) {
						if (file->inferred && file->treeMeshGui->infoOverlayButton->value() && Fl::event_button() == FL_LEFT_MOUSE) {
							gluUnProject(
								Fl::event_x(),
								viewport[3]-Fl::event_y(), // important!!
								1.0f, // magic z value
								iMAP->mv,
								projMatrix,
								viewport,
								&cellSelect[0],
								&cellSelect[1],
								&cellSelect[2]);

							file->treeMesh->selectedQuadLeaf = 0;
							bool found = false;
							file->treeMesh->getLeafSelection(cellSelect[0],cellSelect[1],file->treeMesh->quadTree,&found);

							char localizationsString [50];
							char diffusionString [50];
							char potentialString [50];
							char forceXString [50];
							char forceYString [50];
							char forceMagString [50];

							if (file->treeMesh->selectedQuadLeaf != 0 && file->treeMesh->selectedQuadLeaf->active()) {
								sprintf(localizationsString,"%i Localizations",file->treeMesh->selectedQuadLeaf->getCount());
								sprintf(diffusionString,"D = %f [um\262/s]",file->treeMesh->selectedQuadLeaf->getDiffusion());
								sprintf(potentialString,"V = %f [kT]",file->treeMesh->selectedQuadLeaf->getPotential());
								sprintf(forceXString,"Fx = %f [a.u.]",file->treeMesh->selectedQuadLeaf->getForceX());
								sprintf(forceYString,"Fy = %f [a.u.]",file->treeMesh->selectedQuadLeaf->getForceY());
								sprintf(forceMagString,"||F|| = %f [a.u]",file->treeMesh->selectedQuadLeaf->getForce());

								Fl_Menu_Item lclick_menu[] = {
									{localizationsString, 0, nullCallback},
									{diffusionString, 0, nullCallback},
									{potentialString, 0,  nullCallback},
									{forceXString, 0, nullCallback},
									{forceYString, 0, nullCallback},
									{forceMagString, 0, nullCallback},
									{0}
								};

								const Fl_Menu_Item *m = lclick_menu->popup(Fl::event_x(), Fl::event_y(), 0, 0, 0);
								if ( m ) m->do_callback(0, m->user_data());
							}
						}
					}

					return(1);
				case FL_ENTER:
					window()->cursor(FL_CURSOR_CROSS);

					ret = 1;
					break;
				case FL_LEAVE:
					window()->cursor(FL_CURSOR_DEFAULT);
					ret = 1;
					break;
			}
		} else if (file->landscapePlot == true || file->file3D == true) {
			switch (e) {
				case FL_MOUSEWHEEL:
					#ifdef _WIN32
						if (Fl::event_dy() < 0 && iMAP->fovy < 180) { iMAP->fovy += iMAP->fovy/25; }
						else if (Fl::event_dy() > 0 && iMAP->fovy > 0) { iMAP->fovy -= iMAP->fovy/25; }
					#elif __APPLE__
					if (Fl::event_dy() < 0 && iMAP->fovy < 180) { iMAP->fovy += iMAP->fovy/100; }
					else if (Fl::event_dy() > 0 && iMAP->fovy > 0) { iMAP->fovy -= iMAP->fovy/100; }
					#endif
					glMatrixMode(GL_PROJECTION);
					glLoadIdentity();
					gluPerspective(iMAP->fovy, 1.0, 1, 100);
					glMatrixMode(GL_MODELVIEW);
					glLoadIdentity();
					gluLookAt(0,0, iMAP->initialLookAt,
							  0,0,0,
							  0.0, 1.0, 0.0);
					return(1);
				case FL_PUSH:
					lastX = Fl::event_x();
					lastY = Fl::event_y();

		        	if (Fl::event_clicks() == 1) {
		        		file->xTranslate = file->yTranslate = file->zTranslate = 0.0;
		        		file->xRotate = file->yRotate = file->zRotate = 0;
		        		iMAP->fovy = 2*atan(file->maxRange/2/15)*180/PI;
		        	}
					return(1);
				case FL_DRAG:
					currentX = Fl::event_x();
					currentY = Fl::event_y();

					// rotation
					if (Fl::event_button() == FL_LEFT_MOUSE) {
						// y increment
						if ((float) (currentX-lastX) < 0) {	file->yRotate += (360 + (float) (currentX-lastX)); }
						else { file->yRotate += (float) (currentX-lastX); }
						// x increment
						if ((float) (currentY-lastY) < 0) { file->xRotate += (360 + (float) (currentY-lastY)); }
						else { file->xRotate += (float) (currentY-lastY); }
					}

					// panning/translation
					else if (Fl::event_button() == FL_RIGHT_MOUSE) {
						// create selection box
		    			glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
		    			glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
		    			glGetIntegerv(GL_VIEWPORT,viewport);

		    			// assign last position
		    			gluUnProject(
							(float)lastX,
							(float)lastY,
							0.5, // magic z value
							modelMatrix,
							projMatrix,
							viewport,
							&pushPosition[0],
							&pushPosition[1],
							&pushPosition[2]);

		    			// assign current position
		    			gluUnProject(
							(float)currentX,
							(float)currentY,
							0.5, // magic z value
							modelMatrix,
							projMatrix,
							viewport,
							&releasePosition[0],
							&releasePosition[1],
							&releasePosition[2]);

		    			file->xTranslate += releasePosition[0]-pushPosition[0];
		    			file->yTranslate -= releasePosition[1]-pushPosition[1];
					}
					lastX = currentX;
					lastY = currentY;
					return(1);
				case FL_RELEASE:			// handle mouse release
					return(1);
				case FL_ENTER:
					break;
				case FL_LEAVE:
					break;
				case FL_SHORTCUT:
					if (Fl::event_key() == 65307) {
						if (iMAP->stereoVision == true) {
//							iMAP->stereoVision = false;
							stereoVisionCallback((Fl_Widget*)NULL,(void*)0);
						}
					}
					return(1);
			}
		}
		return(ret);
	}
	return(Fl_Gl_Window::handle(e));
}

void drawLandscapeAxes(int screen) {
    if (iMAP->landscapeAxesEnable == true) {

            GLUquadricObj *quadObj = gluNewQuadric();
            gluQuadricNormals(quadObj, GLU_SMOOTH);

    		switch(screen) {
				case 0: // left screen
		            glViewport(0,0,iMAP->screenDimensions[2]/10,iMAP->screenDimensions[2]/10);
					break;
				case 1: // right screen (for stereo viewing)
					glViewport(iMAP->screenDimensions[2]/2.0,0,iMAP->screenDimensions[2]/10,iMAP->screenDimensions[2]/10);
					break;
    		}

            // prevent overlapping viewports
            glClear(GL_DEPTH_BUFFER_BIT);

            glPushMatrix();

                    updateViewAxis();

                    glMatrixMode(GL_PROJECTION);
                    glLoadIdentity();
                    gluPerspective(45, 1.0, 1, 50);

                    glMatrixMode(GL_MODELVIEW);

                    const GLfloat mat_ambient[]    = { 0.20f,0.20f,0.20f,1.0f };  // RGBA
                    const GLfloat mat_diffuse[]    = { 1.0f/2.0f,1.0f/2.0f,1.0f/2.0f,1.0f };  // RGBA
                    const GLfloat mat_specular[]   = { 1.0f/1.5f,1.0f/1.5f,1.0f/1.5f,1.0f };  // RGBA
                    const GLfloat light_position0[] = { 50.0, 50.0, 50.0, 0.0 };  // XYZ

                    glShadeModel(GL_SMOOTH);
                    glEnable(GL_COLOR_MATERIAL);

                    glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
//                      glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
                    glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
                    glMaterialf(GL_FRONT,  GL_SHININESS, 100.0);

                    glEnable(GL_LIGHT0);
                    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
                    glLightfv(GL_LIGHT0, GL_AMBIENT, mat_ambient);
                    glLightfv(GL_LIGHT0, GL_DIFFUSE,   mat_diffuse);
                    glLightfv(GL_LIGHT0, GL_SPECULAR, mat_ambient);
                    glEnable(GL_LIGHTING);
                    glEnable(GL_DEPTH_TEST);

                    glColor4f(1.0,1.0,0.0,0.5);
                    const float scale = 0.05f;

                    glEnable(GL_ALPHA_TEST);
                    glAlphaFunc(GL_GEQUAL, 0.35f);

                    txfEstablishTexture(iMAP->fontTex, 0, GL_TRUE);

                    const int slices = 20;

                    // sphere at origin
                    gluSphere(quadObj,0.5f,slices,slices);


                    // x axis
                    glRotatef(90,0,1,0);
                    gluCylinder(quadObj, 0.2f, 0.2f, 5.0f, slices, slices);
                    glPushMatrix();
                            glTranslatef(0,0,5);
                            gluCylinder(quadObj, 0.5, 0, 1, slices, slices);
                            // draw x
                            glTranslatef(0.0f,0.0f,1.5f);
                            glEnable(GL_TEXTURE_2D);
                            glDisable(GL_LIGHTING);
                            billboardBegin();
                            glTranslatef(-0.75f,-1.0f,0.0f);
                            glScalef(scale,scale,scale);
                            txfRenderString(iMAP->fontTex, "X", strlen("X"));
                            billboardEnd();
                            glEnable(GL_LIGHTING);
                            glDisable(GL_TEXTURE_2D);
                    glPopMatrix();

                    // y axis
                    glRotatef(270,1,0,0);
                    gluCylinder(quadObj, 0.2f, 0.2f, 5.0f, slices, slices);
                    glPushMatrix();
                            glTranslatef(0,0,5);
                            gluCylinder(quadObj, 0.5, 0, 1, slices, slices);
                            // draw y
                            glTranslatef(0.0f,0.0f,1.5f);
                            glEnable(GL_TEXTURE_2D);
                            glDisable(GL_LIGHTING);
                            billboardBegin();
                            glTranslatef(-0.75f,-0.5f,0.0f);
                            glScalef(scale,scale,scale);
                            txfRenderString(iMAP->fontTex, "Y", strlen("Y"));
                            billboardEnd();
                            glEnable(GL_LIGHTING);
                            glDisable(GL_TEXTURE_2D);
                    glPopMatrix();

                    // z axis
                    glRotatef(-90,0,1,0);
                    gluCylinder(quadObj, 0.2f, 0.2f, 5.0f, slices, slices);
                    glPushMatrix();
                            glTranslatef(0,0,5);
                            gluCylinder(quadObj, 0.5, 0, 1, slices, slices);
                            // draw z
                            glTranslatef(0.0f,0.0f,1.5f);
                            glEnable(GL_TEXTURE_2D);
                            glDisable(GL_LIGHTING);
                            billboardBegin();
                            glTranslatef(-0.75f,-1.0f,0.0f);
                            glScalef(scale,scale,scale);

                            switch(file->meshType) {
                            case 0: // square mesh
                            	if (file->squareMeshGui->overlayPointNumberButton->value()) { txfRenderString(iMAP->fontTex, "#", strlen("#")); }
                            	else if (file->squareMeshGui->overlayDiffusionButton->value()) { txfRenderString(iMAP->fontTex, "D", strlen("D")); }
                            	else if (file->squareMeshGui->overlayPotentialButton->value()) { txfRenderString(iMAP->fontTex, "V", strlen("V")); }
                            	else if (file->squareMeshGui->overlayForceMagnitudeButton->value()) {
                            		if (file->optimizationMode == 2) { txfRenderString(iMAP->fontTex, "||Dr||", strlen("||Dr||")); }
                            		else { txfRenderString(iMAP->fontTex, "||F||", strlen("||F||")); }
                            	}
                            	else {  txfRenderString(iMAP->fontTex, "?", strlen("?")); }
                            	break;
                            case 1: // voronoi mesh
                            	if (file->voronoiMeshGui->overlayPointNumberButton->value()) { txfRenderString(iMAP->fontTex, "#", strlen("#")); }
                            	else if (file->voronoiMeshGui->overlayDiffusionButton->value()) { txfRenderString(iMAP->fontTex, "D", strlen("D")); }
                            	else if (file->voronoiMeshGui->overlayPotentialButton->value()) { txfRenderString(iMAP->fontTex, "V", strlen("V")); }
                            	else if (file->voronoiMeshGui->overlayForceMagnitudeButton->value()) {
                            		if (file->optimizationMode == 2) { txfRenderString(iMAP->fontTex, "||Dr||", strlen("||Dr||")); }
                            		else { txfRenderString(iMAP->fontTex, "||F||", strlen("||F||")); }
                            	}
                            	else {  txfRenderString(iMAP->fontTex, "?", strlen("?")); }
                            	break;
                            case 2: // quadtree mesh
                            	if (file->treeMeshGui->overlayPointNumberButton->value()) { txfRenderString(iMAP->fontTex, "#", strlen("#")); }
                            	else if (file->treeMeshGui->overlayDiffusionButton->value()) { txfRenderString(iMAP->fontTex, "D", strlen("D")); }
                            	else if (file->treeMeshGui->overlayPotentialButton->value()) { txfRenderString(iMAP->fontTex, "V", strlen("V")); }
                            	else if (file->treeMeshGui->overlayForceMagnitudeButton->value()) {
                            		if (file->optimizationMode == 2) { txfRenderString(iMAP->fontTex, "||Dr||", strlen("||Dr||")); }
                            		else { txfRenderString(iMAP->fontTex, "||F||", strlen("||F||")); }
                            	}
                            	else {  txfRenderString(iMAP->fontTex, "?", strlen("?")); }
                            	break;
                            }

                            billboardEnd();
                            glEnable(GL_LIGHTING);
                            glDisable(GL_TEXTURE_2D);
                    glPopMatrix();
            glPopMatrix ();

            glDisable(GL_ALPHA_TEST);

            glDisable(GL_LIGHTING);
            glDisable(GL_LIGHT0);
            glDisable(GL_COLOR_MATERIAL);
            glDisable(GL_DEPTH_TEST);

            gluDeleteQuadric(quadObj);
    }
}

void drawAxes(int screen) {
    if (file->file3D == true) {

            GLUquadricObj *quadObj = gluNewQuadric();
            gluQuadricNormals(quadObj, GLU_SMOOTH);

			int dims[4];
			Fl::screen_xywh(dims[0],dims[1],dims[2],dims[3],0);

    		switch(screen) {
				case 0: // left screen
		            glViewport(0,0,dims[2]/10,dims[2]/10);
					break;
				case 1: // right screen (for stereo viewing)
					glViewport(dims[2]/2.0,0,dims[2]/10,dims[2]/10);
					break;
    		}

            // prevent overlapping viewports
            glClear(GL_DEPTH_BUFFER_BIT);

            glPushMatrix();

                    updateViewAxis();

                    glMatrixMode(GL_PROJECTION);
                    glLoadIdentity();
                    gluPerspective(45, 1.0, 1, 50);

                    glMatrixMode(GL_MODELVIEW);

                    const GLfloat mat_ambient[]    = { 0.20f,0.20f,0.20f,1.0f };  // RGBA
                    const GLfloat mat_diffuse[]    = { 1.0f/2.0f,1.0f/2.0f,1.0f/2.0f,1.0f };  // RGBA
                    const GLfloat mat_specular[]   = { 1.0f/1.5f,1.0f/1.5f,1.0f/1.5f,1.0f };  // RGBA
                    const GLfloat light_position0[] = { 50.0, 50.0, 50.0, 0.0 };  // XYZ

                    glShadeModel(GL_SMOOTH);
                    glEnable(GL_COLOR_MATERIAL);

                    glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
//                      glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
                    glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
                    glMaterialf(GL_FRONT,  GL_SHININESS, 100.0);

                    glEnable(GL_LIGHT0);
                    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
                    glLightfv(GL_LIGHT0, GL_AMBIENT, mat_ambient);
                    glLightfv(GL_LIGHT0, GL_DIFFUSE,   mat_diffuse);
                    glLightfv(GL_LIGHT0, GL_SPECULAR, mat_ambient);
                    glEnable(GL_LIGHTING);
                    glEnable(GL_DEPTH_TEST);

                    glColor4f(1.0,1.0,0.0,0.5);
                    const float scale = 0.05f;

                    glEnable(GL_ALPHA_TEST);
                    glAlphaFunc(GL_GEQUAL, 0.35f);

                    txfEstablishTexture(iMAP->fontTex, 0, GL_TRUE);

                    const int slices = 20;

                    // sphere at origin
                    gluSphere(quadObj,0.5f,slices,slices);


                    // x axis
                    glRotatef(90,0,1,0);
                    gluCylinder(quadObj, 0.2f, 0.2f, 5.0f, slices, slices);
                    glPushMatrix();
                            glTranslatef(0,0,5);
                            gluCylinder(quadObj, 0.5, 0, 1, slices, slices);
                            // draw x
                            glTranslatef(0.0f,0.0f,1.5f);
                            glEnable(GL_TEXTURE_2D);
                            glDisable(GL_LIGHTING);
                            billboardBegin();
                            glTranslatef(-0.75f,-1.0f,0.0f);
                            glScalef(scale,scale,scale);
                            txfRenderString(iMAP->fontTex, "X", strlen("X"));
                            billboardEnd();
                            glEnable(GL_LIGHTING);
                            glDisable(GL_TEXTURE_2D);
                    glPopMatrix();

                    // y axis
                    glRotatef(270,1,0,0);
                    gluCylinder(quadObj, 0.2f, 0.2f, 5.0f, slices, slices);
                    glPushMatrix();
                            glTranslatef(0,0,5);
                            gluCylinder(quadObj, 0.5, 0, 1, slices, slices);
                            // draw y
                            glTranslatef(0.0f,0.0f,1.5f);
                            glEnable(GL_TEXTURE_2D);
                            glDisable(GL_LIGHTING);
                            billboardBegin();
                            glTranslatef(-0.75f,-0.5f,0.0f);
                            glScalef(scale,scale,scale);
                            txfRenderString(iMAP->fontTex, "Y", strlen("Y"));
                            billboardEnd();
                            glEnable(GL_LIGHTING);
                            glDisable(GL_TEXTURE_2D);
                    glPopMatrix();

                    // z axis
                    glRotatef(-90,0,1,0);
                    gluCylinder(quadObj, 0.2f, 0.2f, 5.0f, slices, slices);
                    glPushMatrix();
                            glTranslatef(0,0,5);
                            gluCylinder(quadObj, 0.5, 0, 1, slices, slices);
                            // draw z
                            glTranslatef(0.0f,0.0f,1.5f);
                            glEnable(GL_TEXTURE_2D);
                            glDisable(GL_LIGHTING);
                            billboardBegin();
                            glTranslatef(-0.75f,-1.0f,0.0f);
                            glScalef(scale,scale,scale);
                            txfRenderString(iMAP->fontTex, "Z", strlen("Z"));

                            billboardEnd();
                            glEnable(GL_LIGHTING);
                            glDisable(GL_TEXTURE_2D);
                    glPopMatrix();
            glPopMatrix ();

            glDisable(GL_ALPHA_TEST);

            glDisable(GL_LIGHTING);
            glDisable(GL_LIGHT0);
            glDisable(GL_COLOR_MATERIAL);
            glDisable(GL_DEPTH_TEST);

            gluDeleteQuadric(quadObj);

    }

}

void drawZoomBox(float xRange, float yRange) {
	if (iMAP->zoomBoxEnable) {
		glViewport(iMAP->glWindow->w()/2.35,iMAP->glWindow->h()/1.04,iMAP->screenDimensions[2]/10,iMAP->screenDimensions[2]/30);

		glLoadIdentity();
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-1.0, 1.0, -1.0, 1.0, 1.0f, 50.0f);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(0.0f, 0.0f, 20.0f,
				  0.0f, 0.0f, 0.0f,
				  0.0f, 1.0f, 0.0f);

		// prevent overlapping viewports
		glClear(GL_DEPTH_BUFFER_BIT);

		glColor4f(1.0,1.0,1.0,1.0);
		const float scale = 0.015f;
		char label[260];
		sprintf(label,"%.3f x %.3f",xRange,yRange);

		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GEQUAL, 0.35f);

		txfEstablishTexture(iMAP->fontTex, 0, GL_TRUE);

		glPushMatrix();
			glEnable(GL_TEXTURE_2D);
			glTranslatef(-0.8f,-0.6f,0.0f);
			glScalef(scale/3.0,scale,scale);
			txfRenderString(iMAP->fontTex, label, strlen(label));
			glDisable(GL_TEXTURE_2D);
		glPopMatrix();

		glDisable(GL_ALPHA_TEST);
	}
}

void Normalise(xyz *p) {
   double length;

   length = p->x * p->x + p->y * p->y + p->z * p->z;
   if (length > 0) {
      length = sqrt(length);
      p->x /= length;
      p->y /= length;
      p->z /= length;
   } else {
      p->x = 0;
      p->y = 0;
      p->z = 0;
   }
}
