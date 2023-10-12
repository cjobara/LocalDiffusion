///*
// *
// *  InferenceMAP v1.0
// *  22/05/2014
// *
// *  Author: Mohamed El Beheiry, Physico-Chimie Curie, Institut Curie
// *  		Jean-Baptiste Masson, Physics of Biological Systems, Institut Pasteur
// *  Contact e-mail: mohamed.elbeheiry@gmail.com
// *  Copyright (c) 2014, Mohamed El Beheiry, Jean-Baptiste Masson, Institut Curie, Institut Pasteur
// *  All rights Reserved.
// *
// *  InferenceMAP is released under an "academic use only" licence.
// *  Details of which are provided in the xxx.doc file.
// *  Usage of InferenceMAP requires acceptance of this license.
// *
// *  User instructions for using InferenceMAP are provided in the InferenceMAP User Manual.
// *
// */
//
//#ifndef UTILITIES_H_
//#define UTILITIES_H_
//
//// OpenGL Libraries
//#ifdef __APPLE__
//	#include <GL/glew.h>
//	#include <OpenGL/OpenGL.h>
////	#include <OpenGL/gl.h>
////	#include <OpenGL/glu.h>
//	#include <GLUT/glut.H>
//	#include <FL/gl.h>
//#elif _WIN32
//	#include <GL/glew.h>
//	#include <GL/glext.h>
//	#include <glut.h>
//#endif
//
//#include <stdio.h>
//#include <math.h>
//#include <stddef.h>
//#include <stdlib.h>
//
//#include <iostream>
//#include <fstream>
//#include <sstream>
//
//#define PI 3.1415926535897932384626433832795
//#define PIdiv180 (PI/180.0)
//
////#  if HAVE_PTHREAD_H
//// Use POSIX threading...
//
//#include <pthread.h>
//
//typedef pthread_t Fl_Thread;
//
//static int fl_create_thread(Fl_Thread& t, void *(*f) (void *), void* p) {
//  return pthread_create((pthread_t*)&t, 0, f, p);
//}
//
//static void fl_exit_thread(void *value) {
//	pthread_exit(NULL);
//}
//
//static int fl_join_thread(Fl_Thread th, void **thread_return) {
//	return pthread_join(th,thread_return);
//}
//
////static int fl_destroy_thread(Fl_Thread& t);
////#  elif defined(WIN32) && !defined(__WATCOMC__) // Use Windows threading...
////
////#    include <windows.h>
////#    include <process.h>
////
////typedef unsigned long Fl_Thread;
////
////static int fl_create_thread(Fl_Thread& t, void *(*f) (void *), void* p) {
////  return t = (Fl_Thread)_beginthread((void( __cdecl * )( void * ))f, 0, p);
////}
////
////#  elif defined(__WATCOMC__)
////#    include <process.h>
////
////typedef unsigned long Fl_Thread;
////
////static int fl_create_thread(Fl_Thread& t, void *(*f) (void *), void* p) {
////  return t = (Fl_Thread)_beginthread((void(* )( void * ))f, 32000, p);
////}
////#  endif // !HAVE_PTHREAD_H
//
//void ind2sub(int *siz, int N, int idx, int *sub);
//int sub2ind(int *siz, int N, int *sub);
//
//struct Vector3D {
//	GLfloat x,y,z;
//	Vector3D(GLfloat xx, GLfloat yy, GLfloat zz) : x(xx), y(yy), z(zz) {}
//	Vector3D() : x(0.0), y(0.0), z(0.0) {}
//};
//
//float DotProduct (Vector3D * u, Vector3D * v);
//Vector3D Normalize3dVector(Vector3D v);
//
//class Camera {
//	private:
//
//		Vector3D ViewDir;
//		Vector3D RightVector;
//		Vector3D UpVector;
//		Vector3D Position;
//
//		GLfloat RotatedX, RotatedY, RotatedZ;
//
//	public:
//		Camera();				//inits the values (Position: (0|0|0) Target: (0|0|-1) )
//		void Render ( void );	//executes some glRotates and a glTranslate command
//								//Note: You should call glLoadIdentity before using Render
//
//		void Move ( Vector3D Direction );
//		void RotateX ( GLfloat Angle );
//		void RotateY ( GLfloat Angle );
//		void RotateZ ( GLfloat Angle );
//
//		void MoveForward ( GLfloat Distance );
//		void MoveUpward ( GLfloat Distance );
//		void StrafeRight ( GLfloat Distance );
//
//		void Reset( void );
//};
//
//#endif /* UTILITIES_H_ */
//
////// "$Id: threads.cxx 8864 2011-07-19 04:49:30Z greg.ercolano $"
//////
////// Threading example program for the Fast Light Tool Kit (FLTK).
//////
////// Copyright 1998-2010 by Bill Spitzak and others.
//////
////// This library is free software. Distribution and use rights are outlined in
////// the file "COPYING" which should have been included with this file.  If this
////// file is missing or damaged, see the license at:
//////
//////     http://www.fltk.org/COPYING.php
//////
////// Please report all bugs and problems on the following page:
//////
//////     http://www.fltk.org/str.php
//////
////
////#include <config.h>
////
////#if HAVE_PTHREAD || defined(WIN32)
////#  include <FL/Fl.H>
////#  include <FL/Fl_Double_Window.H>
////#  include <FL/Fl_Browser.H>
////#  include <FL/Fl_Value_Output.H>
////#  include <FL/fl_ask.H>
////#  include "threads.h"
////#  include <stdio.h>
////#  include <math.h>
////
////Fl_Thread prime_thread;
////
////Fl_Browser *browser1, *browser2;
////Fl_Value_Output *value1, *value2;
////int start2 = 3;
////
////void magic_number_cb(void *p)
////{
////  Fl_Value_Output *w = (Fl_Value_Output*)p;
////  w->labelcolor(FL_RED);
////  w->redraw_label();
////}
////
////void* prime_func(void* p)
////{
////  Fl_Browser* browser = (Fl_Browser*) p;
////  Fl_Value_Output *value;
////  int n;
////  int step;
////  char proud = 0;
////
////  if (browser == browser2) {
////    n      = start2;
////    start2 += 2;
////    step   = 12;
////    value  = value2;
////  } else {
////    n     = 3;
////    step  = 2;
////    value = value1;
////  }
////
////  // very simple prime number calculator !
////  //
////  // The return at the end of this function can never be reached and thus
////  // will generate a warning with some compilers, however we need to have
////  // a return statement or other compilers will complain there is no return
////  // statement. To avoid warnings on all compilers, we fool the smart ones
////  // into beleiving that there is a chance that we reach the end by testing
////  // n>=0, knowing that logically, n will never be negative in this context.
////  if (n>=0) for (;;) {
////    int pp;
////    int hn = (int)sqrt((double)n);
////
////    for (pp=3; pp<=hn; pp+=2) if ( n%pp == 0 ) break;
////    if (pp >= hn) {
////      char s[128];
////      sprintf(s, "%d", n);
////
////      // Obtain a lock before we access the browser widget...
////      Fl::lock();
////
////      browser->add(s);
////      browser->bottomline(browser->size());
////      if (n > value->value()) value->value(n);
////      n += step;
////
////      // Release the lock...
////      Fl::unlock();
////
////      // Send a message to the main thread, at which point it will
////      // process any pending redraws for our browser widget.  The
////      // message we pass here isn't used for anything, so we could also
////      // just pass NULL.
////      Fl::awake(p);
////      if (n>10000 && !proud) {
////        proud = 1;
////        Fl::awake(magic_number_cb, value);
////      }
////    } else {
////      // This should not be necessary since "n" and "step" are local variables,
////      // however it appears that at least MacOS X has some threading issues
////      // that cause semi-random corruption of the (stack) variables.
////      Fl::lock();
////      n += step;
////      Fl::unlock();
////    }
////  }
////  return 0L;
////}
////
////int main(int argc, char **argv)
////{
////  Fl_Double_Window* w = new Fl_Double_Window(200, 200, "Single Thread");
////  browser1 = new Fl_Browser(0, 0, 200, 175);
////  w->resizable(browser1);
////  value1 = new Fl_Value_Output(100, 175, 200, 25, "Max Prime:");
////  w->end();
////  w->show(argc, argv);
////  w = new Fl_Double_Window(200, 200, "Six Threads");
////  browser2 = new Fl_Browser(0, 0, 200, 175);
////  w->resizable(browser2);
////  value2 = new Fl_Value_Output(100, 175, 200, 25, "Max Prime:");
////  w->end();
////  w->show();
////
////  browser1->add("Prime numbers:");
////  browser2->add("Prime numbers:");
////
////  // Enable multi-thread support by locking from the main
////  // thread.  Fl::wait() and Fl::run() call Fl::unlock() and
////  // Fl::lock() as needed to release control to the child threads
////  // when it is safe to do so...
////  Fl::lock();
////
////  // Start threads...
////
////  // One thread displaying in one browser
////  fl_create_thread(prime_thread, prime_func, browser1);
////
////  // Several threads displaying in another browser
////  fl_create_thread(prime_thread, prime_func, browser2);
////  fl_create_thread(prime_thread, prime_func, browser2);
////  fl_create_thread(prime_thread, prime_func, browser2);
////  fl_create_thread(prime_thread, prime_func, browser2);
////  fl_create_thread(prime_thread, prime_func, browser2);
////  fl_create_thread(prime_thread, prime_func, browser2);
////
////  Fl::run();
////
////  return 0;
////}
////#else
////#  include <FL/fl_ask.H>
////
////int main() {
////  fl_alert("Sorry, threading not supported on this platform!");
////}
////#endif // HAVE_PTHREAD || WIN32
////
////
//////
////// End of "$Id: threads.cxx 8864 2011-07-19 04:49:30Z greg.ercolano $".
//////
