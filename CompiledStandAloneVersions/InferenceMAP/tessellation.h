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

#ifndef TESSELLATION_H_
#define TESSELLATION_H_

#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "cell.h"

#ifndef NULL
#define NULL 0
#endif
#define DELETED -2

#define le 0
#define re 1

using namespace std;

//namespace abelian
//{

struct PointVDG
{
	float x,y;
};

#define FLOAT_EPSILON	0.00001
#define THETA_EPSILON	0.0001
class Point
{
	public:

	Point(float _x, float _y) {x = _x; y = _y;}

	bool equals(const Point& a)
	{
		float diffx = a.x - x;
		diffx = (diffx < 0) ? -diffx : diffx;

		float diffy = a.y - y;
		diffy = (diffy < 0) ? -diffy : diffy;

		if(diffx < FLOAT_EPSILON && diffy < FLOAT_EPSILON) {
			return true;
		}
		else return false;
	}

	bool operator==(const Point& a)
	{
		return equals(a);
	}

	struct comparator
	{
		//this function acts as both the equals and less operator
		//the == and < of the class don't get used when std::set has a comparator
		bool operator() (const Point& lhs, const Point& rhs) const
		{
			float angle_lhs = atan2(lhs.y, lhs.x);
			float angle_rhs = atan2(rhs.y, rhs.x);

			float diff = angle_lhs - angle_rhs;
			diff = (diff < 0) ? -diff : diff;

			if(diff < THETA_EPSILON) return false;	//always the case (takes precedence)
			else if(angle_lhs < angle_rhs) return true;
			else return false;
		}
	};

	float x;
	float y;
};

class Cell
{
	public:
		Cell();
		~Cell();

		struct PointVDG		coord;
		set< Point, Point::comparator >		border;
		//set< Point >		border;

		float height;

		void copy(VoronoiCell *ic);
		void draw();
		void debugDraw();
		void cleanup();

	protected:
};

//}	//namespace


struct	Freenode
{
	struct	Freenode *nextfree;
};

struct FreeNodeArrayList
{
	struct	Freenode* memory;
	struct	FreeNodeArrayList* next;

};

struct	Freelist
{
	struct	Freenode	*head;
	int		nodesize;
};

struct Point3
{
	float x,y,z;
	int count;
};

struct VertexLink
{
	PointVDG coord;
	PointVDG v[3];
	int count;
};

// structure used both for sites and for vertices
struct Site
{
	struct	PointVDG	coord;
	int		sitenbr;
	int		refcnt;
	int		overallRefcnt;
};

struct Edge
{
	float   a,b,c;
	struct	Site 	*ep[2];
	struct	Site	*reg[2];
	int		edgenbr;

};

struct GraphEdge
{
	float x1,y1,x2,y2;
//	long v1,v2; //vertices that this was created from
	struct GraphEdge* next;
};

struct Halfedge
{
	struct	Halfedge	*ELleft, *ELright;
	struct	Edge	*ELedge;
	int		ELrefcnt;
	char	ELpm;
	struct	Site	*vertex;
	float	ystar;
	struct	Halfedge *PQnext;
};

class VoronoiDiagramGenerator
{
public:
	VoronoiDiagramGenerator();
	virtual ~VoronoiDiagramGenerator();

	bool generateVoronoi(float *xValues, float *yValues, int numPoints,
		float minX, float maxX, float minY, float maxY, float minDist, bool genVectorInfo=true);

	//By default, the delaunay triangulation is NOT generated
	void setGenerateDelaunay(bool genDel);

	//By default, the voronoi diagram IS generated
	void setGenerateVoronoi(bool genVor);

	void resetIterator()
	{
		iteratorEdges = allEdges;
	}

	bool getNext(float& x1, float& y1, float& x2, float& y2)
	{
		if(iteratorEdges == 0)
			return false;

		x1 = iteratorEdges->x1;
		x2 = iteratorEdges->x2;
		y1 = iteratorEdges->y1;
		y2 = iteratorEdges->y2;

//		//LOG<<"getNext returned the edge ("<<x1<<","<<y1<<") -> ("<<x2<<","<<y2<<")";

		iteratorEdges = iteratorEdges->next;

		return true;
	}

	void resetDelaunayEdgesIterator()
	{
		iteratorDelaunayEdges = delaunayEdges;
		//LOG<<"resetDelaunayEdgesIterator set iteratorDelaunayEdges = "<<iteratorDelaunayEdges;
	}

	bool getNextDelaunay(float& x1, float& y1, float& x2, float& y2)
	{
		if(iteratorDelaunayEdges == 0)
		{
			//LOG<<"iteratorDelaunayEdges = 0, returning false";
			return false;
		}
		x1 = iteratorDelaunayEdges->x1;
		x2 = iteratorDelaunayEdges->x2;
		y1 = iteratorDelaunayEdges->y1;
		y2 = iteratorDelaunayEdges->y2;

		iteratorDelaunayEdges = iteratorDelaunayEdges->next;

		//LOG<<"getNextDelaunay returned the edge ("<<x1<<","<<y1<<") -> ("<<x2<<","<<y2<<")";

		return true;
	}

	void resetVertexPairIterator()
	{
		currentVertexLink = 0;
	}

	bool getNextVertexPair(float& x1, float& y1, float& x2, float& y2);

	void resetVerticesIterator()
	{
		currentVertex = 0;
	}

	bool getNextVertex(float& x, float& y)
	{
		if(finalVertices == 0)
			return false;

		if(currentVertex >= sizeOfFinalVertices) return false;
		x = finalVertices[currentVertex].x;
		y = finalVertices[currentVertex].y;
		currentVertex++;
		return true;
	}

	void reset();


public:
	void cleanup();
	void cleanupEdges();
	char *getfree(struct Freelist *fl);
	struct	Halfedge *PQfind();
	int PQempty();

	struct	Halfedge **ELhash;
	struct	Halfedge *HEcreate(), *ELleft(), *ELright(), *ELleftbnd();
	struct	Halfedge *HEcreate(struct Edge *e,int pm);


	struct PointVDG PQ_min();
	struct Halfedge *PQextractmin();
	void freeinit(struct Freelist *fl,int size);
	void makefree(struct Freenode *curr,struct Freelist *fl);
	void geominit();
	void plotinit();
	bool voronoi(bool genVectorInfo);
	void ref(struct Site *v);
	void deref(struct Site *v);
	void endpoint(struct Edge *e,int lr,struct Site * s);

	void ELdelete(struct Halfedge *he);
	struct Halfedge *ELleftbnd(struct PointVDG *p);
	struct Halfedge *ELright(struct Halfedge *he);
	void makevertex(struct Site *v);
	void out_triple(struct Site *s1, struct Site *s2,struct Site * s3);

	void		PQinsert(struct Halfedge *he,struct Site * v, float offset);
	void		PQdelete(struct Halfedge *he);
	bool		ELinitialize();
	void		ELinsert(struct	Halfedge *lb, struct Halfedge *newHe);
	struct Halfedge * ELgethash(int b);
	struct Halfedge *ELleft(struct Halfedge *he);
	struct Site *leftreg(struct Halfedge *he);
	void		out_site(struct Site *s);
	bool		PQinitialize();
	int			PQbucket(struct Halfedge *he);
	void		clip_line(struct Edge *e);
	char		*myalloc(unsigned n);
	int			right_of(struct Halfedge *el,struct PointVDG *p);

	struct Site *rightreg(struct Halfedge *he);
	struct Edge *bisect(struct	Site *s1,struct	Site *s2);
	float dist(struct Site *s,struct Site *t);
	struct Site *intersect(struct Halfedge *el1, struct Halfedge *el2, struct PointVDG *p=0);

	void		out_bisector(struct Edge *e);
	void		out_ep(struct Edge *e);
	void		out_vertex(struct Site *v);
	struct Site *nextone();

	void		pushGraphEdge(float x1, float y1, float x2, float y2);
	void		pushDelaunayGraphEdge(float x1, float y1, float x2, float y2);


	void		openpl();
	void		line(float x1, float y1, float x2, float y2);
	void		circle(float x, float y, float radius);
	void		range(float minX, float minY, float maxX, float maxY);

	void  		insertVertexAddress(long vertexNum, struct Site* address);
	void		insertVertexLink(long vertexNum, long vertexLinkedTo);
	void		generateVertexLinks();

	bool		genDelaunay;
	bool		genVoronoi;

	struct		Freelist	hfl;
	struct		Halfedge *ELleftend, *ELrightend;
	int 		ELhashsize;

	int			triangulate, sorted, plot, debug;
	float		xmin, xmax, ymin, ymax, deltax, deltay;

	struct		Site	*sites;
	int			nsites;
	int			siteidx;
	int			sqrt_nsites;
	int			nvertices;
	struct 		Freelist sfl;
	struct		Site	*bottomsite;

	int			nedges;
	struct		Freelist efl;
	int			PQhashsize;
	struct		Halfedge *PQhash;
	int			PQcount;
	int			PQmin;

	int			ntry, totalsearch;
	float		pxmin, pxmax, pymin, pymax, cradius;
	int			total_alloc;

	float		borderMinX, borderMaxX, borderMinY, borderMaxY;

	FreeNodeArrayList* allMemoryList;
	FreeNodeArrayList* currentMemoryBlock;

	GraphEdge*	allEdges;
	GraphEdge*	iteratorEdges;

	GraphEdge*	delaunayEdges;
	GraphEdge*	iteratorDelaunayEdges;

	Point3*		vertexLinks; //lists all the vectors that each vector is directly connected to
	long		sizeOfVertexLinks;

	Site**		vertices;
	long		sizeOfVertices ;

	VertexLink* finalVertexLinks;
	long 		sizeOfFinalVertexLinks;
	long		currentVertexLink;

	PointVDG	*finalVertices;
	long		sizeOfFinalVertices ;
	long 		currentVertex;

	float		minDistanceBetweenSites;

//	abelian::MemoryPool< abelian::Cell > cellPool;
	Cell *cells;
};

int scomp(const void *p1,const void *p2);

#endif /* TESSELLATION_H_ */
