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

#ifndef TEXFONT_H_
#define TEXFONT_H_

// STANDARD LIBRARIES
#include <stdio.h>

// OpenGL Libraries
#ifdef __APPLE__
	#include <GL/glew.h>
	#include <OpenGL/OpenGL.h>
//	#include <OpenGL/gl.h>
//	#include <OpenGL/glu.h>
	#include <GLUT/glut.H>
	#include <FL/gl.h>
#elif _WIN32
	#include <GL/glew.h>
	#include <GL/glext.h>
	#include <glut.h>
#endif


#define TXF_FORMAT_BYTE		0
#define TXF_FORMAT_BITMAP	1

typedef struct {
  unsigned short c;       /* Potentially support 16-bit glyphs. */
  unsigned char width;
  unsigned char height;
  signed char xoffset;
  signed char yoffset;
  signed char advance;
  char dummy;           /* Space holder for alignment reasons. */
  short x;
  short y;
} TexGlyphInfo;

typedef struct {
  GLfloat t0[2];
  GLshort v0[2];
  GLfloat t1[2];
  GLshort v1[2];
  GLfloat t2[2];
  GLshort v2[2];
  GLfloat t3[2];
  GLshort v3[2];
  GLfloat advance;
} TexGlyphVertexInfo;

typedef struct {
  GLuint texobj;
  int tex_width;
  int tex_height;
  int max_ascent;
  int max_descent;
  int num_glyphs;
  int min_glyph;
  int range;
  unsigned char *teximage;
  TexGlyphInfo *tgi;
  TexGlyphVertexInfo *tgvi;
  TexGlyphVertexInfo **lut;
} TexFont;

extern char *txfErrorString(void);

extern TexFont *txfLoadFont(
  char *filename);

extern void txfUnloadFont(
  TexFont * txf);

extern GLuint txfEstablishTexture(
  TexFont * txf,
  GLuint texobj,
  GLboolean setupMipmaps);

extern void txfBindFontTexture(
  TexFont * txf);

extern void txfGetStringMetrics(
  TexFont * txf,
  char *string,
  int len,
  int *width,
  int *max_ascent,
  int *max_descent);

extern void txfRenderGlyph(
  TexFont * txf,
  int c);

extern void txfRenderString(
  TexFont * txf,
  char *string,
  int len);

extern void txfRenderFancyString(
  TexFont * txf,
  char *string,
  int len);

#endif /* __TEXFONT_H__ */
