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

#ifndef STRUCTS_H_
#define STRUCTS_H_

#define CROSSPROD(p1,p2,p3) \
   p3.x = p1.y*p2.z - p1.z*p2.y; \
   p3.y = p1.z*p2.x - p1.x*p2.z; \
   p3.z = p1.x*p2.y - p1.y*p2.x

#include <stddef.h>

struct Detection {
	double x;
	double y;
	double z;
	double i;
	double t;
//	bool visited;
	Detection() {
		x = 0.0;
		y = 0.0;
		z = 0.0;
		i = 0.0;
		t = 0.0;
//		visited = false;
	}
};

struct Trajectory {
	int identifier;
	int n; // number of detections in the trajectory
	Detection* d;
	double *xy;
	double *xyz;
	double averageDx;
	double averageDy;
	double averageDz;
	double averageDt;
	double xMax;
	double xMin;
	double yMax;
	double yMin;
	double zMax;
	double zMin;
	double length;
	float rgb[3];
	// constructor
	Trajectory() {
		identifier = -1;
		rgb[0] = 1.0;
		rgb[1] = 1.0;
		rgb[2] = 1.0;
		n = 0;
		d = NULL;
		length = 0.0;
		averageDx = 0.0;
		averageDy = 0.0;
		averageDz = 0.0;
		averageDt = 0.0;
		xMin = yMin = zMin = 1000000.0;
		xMax = yMax = zMax = -1000000.0;
		xy = NULL;
		xyz = NULL;
	}
};

struct Square {
	Detection p1,p2,p3,p4;
	bool inside;
	Square() {
		inside = false;
	}
	/*
	 *  p1 ------ p2
	 *    |      |
	 *    |      |
	 *  p3 ------ p4
	 */
};

struct Zone {
	float value;
	bool active;
	float xCentre;
	float yCentre;
	Zone() {
		value = 0.0;
		active = false;
		xCentre = 0.0;
		yCentre = 0.0;
	}
};

struct DataNode {
	double data;
	DataNode *next;
};

struct Coordinate {
	float x;
	float y;
//	float z;
	int index;
	int id;
	// constructor
	Coordinate() {
		x = 0.0f;
		y = 0.0f;
//		z = 0.0f;
		index = 0;
		id = -1;
	}
};

struct Partition {
	int count;
	int vCount;
	int progress;
	Coordinate *coords;
//	Vertex3d *vertices;
	// constructor
	Partition() {
		vCount = 0;
		count = 0;
		progress = 0;
		coords = NULL;
//		vertices = NULL;
	}
};

typedef struct {
	double x,y,z;
} xyz;

#endif /* STRUCTS_H_ */
