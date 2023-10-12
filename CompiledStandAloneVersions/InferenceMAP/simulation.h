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

#ifndef SIMULATION_H_
#define SIMULATION_H_

class Simulation;

#include "mesh.h"
#include "trees.h"
#include "structs.h"

class Simulation {
	public:
		Simulation(SquareMesh *mesh);
		Simulation(VoronoiMesh *mesh);
		Simulation(TreeMesh *mesh);

		int numberOfTrajectories;
		double dt;
		int steps;
		double duration;

		void createTrajectories(SquareMesh *mesh);
		void createTrajectories(VoronoiMesh *mesh);
		void createTrajectories(TreeMesh *mesh);

		/* Quad-Tree Functions */
		QuadTree *initialLeaf;
		void quadTreeInitialPosition(QuadTree *tree, int *increment, int *index);

		void saveTrajectories();

		Trajectory *trajectories;

		long seed;

		int n_traj_sortie;
		int n_sortie_regularise;

		double **traj_sortie;
		double **traj_sortie_regularise;

		friend class SquareMesh;
		friend class VoronoiMesh;
		friend class TreeMesh;
};

#endif /* SIMULATION_H_ */
