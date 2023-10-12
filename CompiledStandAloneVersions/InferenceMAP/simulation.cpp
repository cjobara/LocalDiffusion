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

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <cmath>
#include <stdlib.h>
#include <string.h>

#include "globals.h"
#include "optimization.h"
#include "file.h"

extern Globals *iMAP;
extern File *file;

Simulation::Simulation(SquareMesh *mesh) {
	numberOfTrajectories = (int)file->squareMeshGui->numberOfTrajectoriesSlider->value();
	dt = file->squareMeshGui->deltaTSlider->value()/1000.0;
	steps = (int)file->squareMeshGui->timeStepsSlider->value();

	duration = (double)steps*dt;

	trajectories = new Trajectory[numberOfTrajectories];

	seed = -time(NULL);

	n_traj_sortie =(int) floor(5.*floor(duration/dt));
	n_sortie_regularise = (int) floor(duration/dt);

	traj_sortie = new double*[steps+1];
	traj_sortie_regularise = new double*[steps+1];
	for (int i = 0; i < steps+1; i++) {
		traj_sortie[i] = new double[3];
		traj_sortie_regularise[i] = new double[3];
	}
	for (int i = 0; i < steps+1; i++) {
		traj_sortie[i][0] = 0.0;
		traj_sortie[i][1] = 0.0;
		traj_sortie[i][2] = 0.0;
		traj_sortie_regularise[i][0] = 0.0;
		traj_sortie_regularise[i][1] = 0.0;
		traj_sortie_regularise[i][2] = 0.0;
	}
}

Simulation::Simulation(VoronoiMesh *mesh) {
	numberOfTrajectories = (int)file->voronoiMeshGui->numberOfTrajectoriesSlider->value();
	dt = file->voronoiMeshGui->deltaTSlider->value()/1000.0;
	steps = (int)file->voronoiMeshGui->timeStepsSlider->value();

	duration = (double)steps*dt;

	trajectories = new Trajectory[numberOfTrajectories];

	seed = -time(NULL);
	n_traj_sortie =(int) floor(5.*floor(duration/dt));
	n_sortie_regularise = (int) floor(duration/dt);

	traj_sortie = new double*[steps+1];
	traj_sortie_regularise = new double*[steps+1];
	for (int i = 0; i < steps+1; i++) {
		traj_sortie[i] = new double[3];
		traj_sortie_regularise[i] = new double[3];
	}
	for (int i = 0; i < steps+1; i++) {
		traj_sortie[i][0] = 0.0;
		traj_sortie[i][1] = 0.0;
		traj_sortie[i][2] = 0.0;
		traj_sortie_regularise[i][0] = 0.0;
		traj_sortie_regularise[i][1] = 0.0;
		traj_sortie_regularise[i][2] = 0.0;
	}
}

Simulation::Simulation(TreeMesh *mesh) {
	numberOfTrajectories = (int)file->treeMeshGui->numberOfTrajectoriesSlider->value();
	dt = file->treeMeshGui->deltaTSlider->value()/1000.0;
	steps = (int)file->treeMeshGui->timeStepsSlider->value();

	duration = (double)steps*dt;

	trajectories = new Trajectory[numberOfTrajectories];

	seed = -time(NULL);
	n_traj_sortie =(int) floor(5.*floor(duration/dt));
	n_sortie_regularise = (int) floor(duration/dt);

	traj_sortie = new double*[steps+1];
	traj_sortie_regularise = new double*[steps+1];
	for (int i = 0; i < steps+1; i++) {
		traj_sortie[i] = new double[3];
		traj_sortie_regularise[i] = new double[3];
	}
	for (int i = 0; i < steps+1; i++) {
		traj_sortie[i][0] = 0.0;
		traj_sortie[i][1] = 0.0;
		traj_sortie[i][2] = 0.0;
		traj_sortie_regularise[i][0] = 0.0;
		traj_sortie_regularise[i][1] = 0.0;
		traj_sortie_regularise[i][2] = 0.0;
	}
}

void Simulation::createTrajectories(SquareMesh *mesh) {
	int i,j,indice, indice_loc, i_init, j_init;
	double r1, r2,t, tau;
	double a0, a1, a2, a3, a4;
	int iOld,jOld;

	for (int h = 0; h < numberOfTrajectories; h++) {

		// set initial values to zero
		for (int c = 0; c < steps; c++) {
			traj_sortie[c][0] = 0.0;
			traj_sortie[c][1] = 0.0;
			traj_sortie[c][2] = 0.0;
		}

		// determine starting zone
		indice =(int) floor(mesh->dimensions*ran1(&seed));

	//	fprintf(stderr,"cell %i of %i\n",indice,mesh->dimensions);

		int xPos,yPos;
		int c = 0;
		for (int a = 0; a < mesh->xCells; a++) {
			for (int b = 0; b < mesh->yCells; b++) {
				if (mesh->active(a,b)) {
					c++;
				}
				/* PROBLEM IS HERE */
				if (c == indice && c > 0) {
					xPos = a;
					yPos = b;
					indice = -1;
					break;
				}
				if (indice == 0 && c == 1) {
					xPos = a;
					yPos = b;
					indice = -1;
					break;
				}
			}
		}

		i = xPos; i_init = i;
		j = yPos; j_init = j;

		t = 0.0;
		indice_loc = 0;

		// only 4 transition probabilities for regular meshing
		a1 = -10.;
		a2 = -10.;
		a3 = -10.;
		a4 = -10.;

		while (t < duration) {
			r1 = ran1(&seed);
			r2 = ran1(&seed);
			iOld = i;
			jOld = j;
			if ((i ==0) && (j ==0)){
	//			fprintf(stderr, "condition 1 (%i,%i)\n",i,j);;
				a1 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-0.5*mesh->dy*(mesh->cells[i][j+1].potential-mesh->cells[i][j].potential)*4e-21/mesh->dy);
				a2 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-0.5*mesh->dx*(mesh->cells[i+1][j].potential-mesh->cells[i][j].potential)*4e-21/mesh->dx);
				a0 = a1+a2;
				tau = 1/a0*log(1/r1);

				if (tau>100.){ break; }
				t += tau;
				if (t > duration){ }
				else {
				if (r2*a0 <= a1) { j += 1; }
				else { i += 1; }
				}
			}
			else if ((i == mesh->xCells-1) && (j ==0)){
	//			fprintf(stderr, "condition 2 (%i,%i)\n",i,j);
				a1 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j+1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a4 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i-1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a0 = a1 + a4;
				tau = 1/a0*log(1/r1);

				if (tau>100.){break;}
				t += tau;
				if (t>duration){ } else{
				if (r2*a0 <= a1){ j += 1;}
				else {i -= 1;}
					}
			}
			else if ((i ==mesh->xCells-1)&&(j ==mesh->yCells-1)){
	//			fprintf(stderr, "condition 3 (%i,%i)\n",i,j);
				a3 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j-1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a4 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i-1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a0 = a3 +a4;
				tau = 1/a0*log(1/r1);
				if (tau>100.){break;}
				t += tau;
				if (t>duration){} else{
				if (r2*a0 <= a3){ j -= 1;}
				else {i -= 1;}
					}
			}
			else if ((i ==0)&&(j == mesh->yCells-1)){
	//			fprintf(stderr, "condition 4 (%i,%i)\n",i,j);
				a2 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i+1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a3 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j-1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a0 = a2 +a3;
				tau = 1/a0*log(1/r1);
				if (tau>100.){break;}
				t += tau;
				if (t>duration){} else{
				if (r2*a0 <= a2){ i += 1;}
				else {j -= 1;}
					}
			}
			else if (i==0){
	//		fprintf(stderr, "condition 5 (%i,%i)\n",i,j);
				a1 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j+1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a2 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i+1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a3 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j-1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));

				a0 = a1 +a2 +a3;
				tau = 1/a0*log(1/r1);
				if (tau>100.){break;}
				t += tau;
				if (t>duration){} else{
				if (r2*a0 <= a1){ j += 1;}
				else if ( (a1<r2*a0)&&(r2*a0<=a1+a2) ){i += 1;}
				else {j -= 1;}
					}
			}
			else if (i == mesh->xCells-1){
	//			fprintf(stderr, "condition 6 (%i,%i)\n",i,j);
				a1 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j+1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a3 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j-1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a4 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i-1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a0 = a1 +a3 +a4;
				tau = 1/a0*log(1/r1);
				if (tau>100.){break;}
				t += tau;
				if (t>duration){} else{
				if (r2*a0 <= a1){ j += 1;}
				else if ( (a1<r2*a0)&&(r2*a0<=a1+a3) ){j -= 1;}
				else {i -= 1;}
					}
			}
			else if (j == 0){
	//			fprintf(stderr, "condition 7 (%i,%i)\n",i,j);
				a1 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j+1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a2 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i+1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a4 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i-1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a0 = a1 +a2 +a4;
				tau = 1/a0*log(1/r1);
				if (tau>100.){break;}
				t += tau;
				if (t>duration){} else{
				if (r2*a0 <= a1){ j += 1;}
				else if ( (a1<r2*a0)&&(r2*a0<=a1+a2) ){i += 1;}
				else {i -= 1;}
					}
			}
			else if (j == mesh->yCells-1){
	//			fprintf(stderr, "condition 8 (%i,%i)\n",i,j);
				a2 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i+1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a3 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j-1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a4 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i-1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a0 = a2 +a3 +a4;
				tau = 1/a0*log(1/r1);
				if (tau>100.){break;}
				t += tau;
				if (t>duration){} else{
				if (r2*a0 <= a2){ i += 1;}
				else if ( (a2<r2*a0)&&(r2*a0<=a2+a3) ){j -= 1;}
				else {i -= 1;}
					}
			}
			else {
	//			fprintf(stderr, "condition 9 (%i,%i)\n",i,j);
				a1 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j+1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a2 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i+1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a3 = mesh->cells[i][j].diffusion/(mesh->dy*mesh->dy)*exp(-1./2*mesh->dy*((mesh->cells[i][j-1].potential - mesh->cells[i][j].potential)*4e-21/mesh->dy));
				a4 = mesh->cells[i][j].diffusion/(mesh->dx*mesh->dx)*exp(-1./2*mesh->dx*((mesh->cells[i-1][j].potential - mesh->cells[i][j].potential)*4e-21/mesh->dx));
				a0 = a1 +a2 +a3 +a4;

	//			fprintf(stderr,"%f + %f + %f + %f = %f\n",a1,a2,a3,a4,a0);

				tau = 1/a0*log(1/r1);
	//			fprintf(stderr,"tau = %f\n",tau);

				if (tau>100.){break;}
				t += tau;
				if (t > duration){} else{
				if (r2*a0 <= a1){ j += 1;}
				else if ( (a1<r2*a0)&&(r2*a0<=a1+a2) ){i += 1;}
				else if ((a1+a2<r2*a0)&&(r2*a0<=a1+a2+a3)){j -= 1;}
				else {i -= 1;}
				}
			}

			indice_loc++;
			const int pos = (int)floor(t/dt);

//			fprintf(stderr,"%f\t%f\t%f\n",traj_sortie[indice_loc][0],traj_sortie[indice_loc][1],traj_sortie[indice_loc][2]);

			if (pos < steps && t < duration) {
				if (indice_loc < steps && fabs(traj_sortie[pos][0]) < 0.000001) {
		//			fprintf(stderr,"%i/%i\n",indice_loc,n_traj_sortie);

	//				traj_sortie[indice_loc][0] = t;
	//				traj_sortie[indice_loc][1] = mesh->cells[i][j].getXCentre();
	//				traj_sortie[indice_loc][2] = mesh->cells[i][j].getYCentre();
	//
	//				traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][0] = floor(traj_sortie[indice_loc][0]/dt)*dt;
	//				traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][1] = mesh->cells[i][j].getXCentre();
	//				traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][2] = mesh->cells[i][j].getYCentre();

	//				traj_sortie[indice_loc][0] = t;
					traj_sortie[pos][0] = (double)pos*dt;
					traj_sortie[pos][1] = mesh->cells[i][j].getXCentre();
					traj_sortie[pos][2] = mesh->cells[i][j].getYCentre();

	//				if (traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][0] == 0.0) {
	//					traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][0] = floor(traj_sortie[indice_loc][0]/dt)*dt;
	//					traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][1] = mesh->cells[i][j].getXCentre();
	//					traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][2] = mesh->cells[i][j].getYCentre();
	//				}
	//				fprintf(stderr,"%i\t%f\t%f\n",indice_loc,t,floor(traj_sortie[indice_loc][0]/dt)*dt);

	//				fprintf(stderr,"%i\t%i\n",indice_loc,(int) floor(traj_sortie[indice_loc][0]/dt));

				} else {
					i = iOld;
					j = jOld;
//					fprintf(stderr,"trajectory %i\n",h);
				}
			}
		}

		trajectories[h].n = steps;
		trajectories[h].d = new Detection[steps];

		int u = 0, v = 0, w = 0;
		for (v = 0; v < steps; v++){

			if (fabs(traj_sortie[v][1]) < 0.000001 && fabs(traj_sortie[v][2]) < 0.000001) {
				if (v > 0) {
					traj_sortie[v][1] = traj_sortie[v-1][1];
					traj_sortie[v][2] = traj_sortie[v-1][2];
				}
				else {
					u = 0;
					while (fabs(traj_sortie[u][1]) < 0.000001 && fabs(traj_sortie[u][2]) < 0.000001 && u < steps) {
//						fprintf(stderr,"traj_sortie[%i][1] = %f\t traj_sortie[%i][2]%f\n",u,traj_sortie[u][1],u,traj_sortie[u][2]);
						u++;
					}
					traj_sortie[0][0] = traj_sortie[u][0];
					traj_sortie[0][1] = traj_sortie[u][1];
					traj_sortie[0][2] = traj_sortie[u][2];
				}
			}
			traj_sortie[v][0] = v*dt;

			trajectories[h].d[w].i = 1.0;
			trajectories[h].d[w].x = traj_sortie[v][1];
			trajectories[h].d[w].y = traj_sortie[v][2];
			trajectories[h].d[w].t = traj_sortie[v][0];

			w++;

//			fprintf(stderr,"[u,v,w] = [%i,%i,%i] (%i)\n",u,v,w,steps);
		}

		//n_sortie_regularise = (int) floor(t/dt);
//		fprintf(stderr,"Trajectory %i of %i\n",h+1,numberOfTrajectories);
	}
//	file->squareMeshGui->saveTrajectoriesButton->activate();
}

void Simulation::createTrajectories(VoronoiMesh *mesh) {
	int i, indice, indice_loc, i_init;
	double r1, r2,t, tau;
	double a0;
	int iOld;

	for (int h = 0; h < numberOfTrajectories; h++) {
		
	fprintf(stderr,"vor 0 : %i\n",h);
		// set initial values to zero
		for (int c = 0; c < steps; c++) {
			traj_sortie[c][0] = 0.0;
			traj_sortie[c][1] = 0.0;
			traj_sortie[c][2] = 0.0;
		}

		// determine starting zone
		indice = (int) floor(mesh->dimensions*ran1(&seed));

		int initPos;
		int c = 0;
		for (int b = 0; b < mesh->numberOfClusters; b++) {
			if (mesh->active(b)) {
				c++;
			}
			/* PROBLEM IS HERE */
			if (c == indice && c > 0) {
				initPos = b;
				indice = -1;
				break;
			}
			if (indice == 0 && c == 1) {
				initPos = b;
				indice = -1;
				break;
			}
		}

		i = initPos; i_init = i;

		t = 0.0;
		indice_loc = 0;

		while (t < duration) {
	fprintf(stderr,"vorr 0 : %i\n",h);
//			fprintf(stderr,"%f/%f\n",t,duration);
			int n = mesh->getNumberOfNeighbours(i);
			double* a = new double[n];
			int* aIndex = new int[n];
			
			r1 = ran1(&seed);
			r2 = ran1(&seed);
			iOld = i;
			a0 = 0.0;
			
			for (int r = 0; r < n; r++) {
				const double dx = fabs(mesh->cells[i].getXCentre()-mesh->cells[mesh->cells[i].getNeighbour(r)].getXCentre());
				const double dy = fabs(mesh->cells[i].getYCentre()-mesh->cells[mesh->cells[i].getNeighbour(r)].getYCentre());

				const double dr = sqrt(dx*dx+dy*dy);

				a[r] = mesh->cells[i].diffusion/(dr*dr)*exp(-1./2*dr*((mesh->cells[i].potential - mesh->cells[mesh->cells[i].getNeighbour(r)].potential)*4e-21/dr));
				aIndex[r] = mesh->cells[i].getNeighbour(r);
				a0 += a[r];
			}
			
			tau = 1.0/a0*log(1.0/r1);

			if (tau > 100.0){ break; }
			t += tau;
			if (t > duration) {	}
			else {
				double *aSum = new double[n+1];
				aSum[0] = 0.0;
				// fill transition sums array (used to determine next reaction site)
				for (int m = 0; m < n; m++) {
					aSum[m+1] = aSum[m] + a[m];
				}
				// find reaction site
				for (int p = 0; p < n; p++) {
					if ( aSum[p] <= r2*a0 && r2*a0 < aSum[p+1]) {
						i = aIndex[p];
						delete [] aSum;
						break;
					}
				}
				//delete [] aSum;
			}
			
			indice_loc++;
			const int pos = (int)floor(t/dt);

//			fprintf(stderr,"%f\t%f\t%f\n",traj_sortie[indice_loc][0],traj_sortie[indice_loc][1],traj_sortie[indice_loc][2]);

			if (pos < steps && t < duration) {
				if (/*indice_loc < steps && */fabs(traj_sortie[pos][0]) < 0.000001) {
		//			fprintf(stderr,"%i/%i\n",indice_loc,n_traj_sortie);

	//				traj_sortie[indice_loc][0] = t;
	//				traj_sortie[indice_loc][1] = mesh->cells[i][j].getXCentre();
	//				traj_sortie[indice_loc][2] = mesh->cells[i][j].getYCentre();
	//
	//				traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][0] = floor(traj_sortie[indice_loc][0]/dt)*dt;
	//				traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][1] = mesh->cells[i][j].getXCentre();
	//				traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][2] = mesh->cells[i][j].getYCentre();

	//				traj_sortie[indice_loc][0] = t;
					traj_sortie[pos][0] = (double)pos*dt;
					traj_sortie[pos][1] = mesh->cells[i].getXCentre();
					traj_sortie[pos][2] = mesh->cells[i].getYCentre();

	//				if (traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][0] == 0.0) {
	//					traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][0] = floor(traj_sortie[indice_loc][0]/dt)*dt;
	//					traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][1] = mesh->cells[i][j].getXCentre();
	//					traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][2] = mesh->cells[i][j].getYCentre();
	//				}
	//				fprintf(stderr,"%i\t%f\t%f\n",indice_loc,t,floor(traj_sortie[indice_loc][0]/dt)*dt);

	//				fprintf(stderr,"%i\t%i\n",indice_loc,(int) floor(traj_sortie[indice_loc][0]/dt));

				} else {
					i = iOld;
	//					fprintf(stderr,"trajectory %i\n",h);
				}
			}
			
			delete [] a;
			delete [] aIndex;
		}
		
		trajectories[h].n = steps;
		trajectories[h].d = new Detection[steps];

		int u = 0, v = 0, w = 0;
		for (v = 0; v < steps; v++){

			if (fabs(traj_sortie[v][1]) < 0.000001 && fabs(traj_sortie[v][2]) < 0.000001) {
				if (v > 0) {
					traj_sortie[v][1] = traj_sortie[v-1][1];
					traj_sortie[v][2] = traj_sortie[v-1][2];
				}
				else {
					u = 0;
					while (fabs(traj_sortie[u][1]) < 0.000001 && fabs(traj_sortie[u][2]) < 0.000001 && u < steps) {
//						fprintf(stderr,"traj_sortie[%i][1] = %f\t traj_sortie[%i][2]%f\n",u,traj_sortie[u][1],u,traj_sortie[u][2]);
						u++;
					}
					traj_sortie[0][0] = traj_sortie[u][0];
					traj_sortie[0][1] = traj_sortie[u][1];
					traj_sortie[0][2] = traj_sortie[u][2];
				}
			}
			traj_sortie[v][0] = v*dt;

			trajectories[h].d[w].i = 1.0;
			trajectories[h].d[w].x = traj_sortie[v][1];
			trajectories[h].d[w].y = traj_sortie[v][2];
			trajectories[h].d[w].t = traj_sortie[v][0];

			w++;

//			fprintf(stderr,"[u,v,w] = [%i,%i,%i] (%i)\n",u,v,w,steps);
		}
		
	fprintf(stderr,"vorr 3 : %i\n",h);
		//n_sortie_regularise = (int) floor(t/dt);
//		fprintf(stderr,"Trajectory %i of %i\n",h+1,numberOfTrajectories);
	}
//	file->voronoiMeshGui->saveTrajectoriesButton->activate();
}

void Simulation::createTrajectories(TreeMesh *mesh) {
	int indice, indice_loc;
	double r1, r2,t, tau;
	double a0;
	QuadTree *iOld;
	QuadTree *i;

	for (int h = 0; h < numberOfTrajectories; h++) {

		// set initial values to zero
		for (int c = 0; c < steps; c++) {
			traj_sortie[c][0] = 0.0;
			traj_sortie[c][1] = 0.0;
			traj_sortie[c][2] = 0.0;
		}

		// determine starting zone
		indice = (int) floor(mesh->dimensions*ran1(&seed));

		int c = 0;
		quadTreeInitialPosition(mesh->quadTree, &c, &indice);
//		fprintf(stderr,"%i %i %i\n",c,indice,initLeaf->identifier);

		i = initialLeaf;

		t = 0.0;
		indice_loc = 0;

		while (t < duration) {
//			fprintf(stderr,"%f/%f\n",t,duration);

			int n = i->getNumberOfNeighbours();

			// construct single Quad-Tree array to store neighbours
			QuadTree** neighbours = new QuadTree*[n];

			int e = 0;
			for (int g = 0; g < i->getNumberOfLeftNeighbours(); g++) {
				neighbours[e] = i->getLeftNeighbour(g);
				e++;
			}
			for (int g = 0; g < i->getNumberOfRightNeighbours(); g++) {
				neighbours[e] = i->getRightNeighbour(g);
				e++;
			}
			for (int g = 0; g < i->getNumberOfTopNeighbours(); g++) {
				neighbours[e] = i->getTopNeighbour(g);
				e++;
			}
			for (int g = 0; g < i->getNumberOfBottomNeighbours(); g++) {
				neighbours[e] = i->getBottomNeighbour(g);
				e++;
			}

			double* a = new double[n];

			r1 = ran1(&seed);
			r2 = ran1(&seed);
			iOld = i;
			a0 = 0.0;

			for (int r = 0; r < n; r++) {
				const double dx = fabs(i->getXCentre()-neighbours[r]->getXCentre());
				const double dy = fabs(i->getYCentre()-neighbours[r]->getYCentre());

				const double dr = sqrt(dx*dx+dy*dy);

				a[r] = i->getD()/(dr*dr)*exp(-1./2*dr*((i->getEp() - neighbours[r]->getEp())*4e-21/dr));
				a0 += a[r];
			}

			tau = 1/a0*log(1/r1);

			if (tau > 100.0){ break; }
			t += tau;
			if (t > duration) {	}
			else {
				double* aSum = new double[n+1];
				aSum[0] = 0.0;
				// fill transition sums array (used to determine next reaction site)
				for (int m = 0; m < n; m++) {
					aSum[m+1] = aSum[m] + a[m];
				}
				// find reaction site
				for (int p = 0; p < n; p++) {
					if ( aSum[p] <= r2*a0 && r2*a0 < aSum[p+1]) {
						i = neighbours[p];
						delete [] aSum;
						break;
					}
				}
				//delete [] aSum;
			}

			indice_loc++;
			const int pos = (int)floor(t/dt);

//			fprintf(stderr,"%f\t%f\t%f\n",traj_sortie[indice_loc][0],traj_sortie[indice_loc][1],traj_sortie[indice_loc][2]);

			if (pos < steps && t < duration) {
				if (/*indice_loc < steps && */fabs(traj_sortie[pos][0]) < 0.000001) {
		//			fprintf(stderr,"%i/%i\n",indice_loc,n_traj_sortie);

	//				traj_sortie[indice_loc][0] = t;
	//				traj_sortie[indice_loc][1] = mesh->cells[i][j].getXCentre();
	//				traj_sortie[indice_loc][2] = mesh->cells[i][j].getYCentre();
	//
	//				traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][0] = floor(traj_sortie[indice_loc][0]/dt)*dt;
	//				traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][1] = mesh->cells[i][j].getXCentre();
	//				traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][2] = mesh->cells[i][j].getYCentre();

	//				traj_sortie[indice_loc][0] = t;
					traj_sortie[pos][0] = (double)pos*dt;
					traj_sortie[pos][1] = i->getXCentre();
					traj_sortie[pos][2] = i->getYCentre();

	//				if (traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][0] == 0.0) {
	//					traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][0] = floor(traj_sortie[indice_loc][0]/dt)*dt;
	//					traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][1] = mesh->cells[i][j].getXCentre();
	//					traj_sortie_regularise[(int) floor(traj_sortie[indice_loc][0]/dt)][2] = mesh->cells[i][j].getYCentre();
	//				}
	//				fprintf(stderr,"%i\t%f\t%f\n",indice_loc,t,floor(traj_sortie[indice_loc][0]/dt)*dt);

	//				fprintf(stderr,"%i\t%i\n",indice_loc,(int) floor(traj_sortie[indice_loc][0]/dt));

				} else {
					i = iOld;
	//					fprintf(stderr,"trajectory %i\n",h);
				}
			}
			delete [] a;
			delete [] neighbours;
		}

		trajectories[h].n = steps;
		trajectories[h].d = new Detection[steps];

		int u = 0, v = 0, w = 0;
		for (v = 0; v < steps; v++){

			if (fabs(traj_sortie[v][1]) < 0.000001 && fabs(traj_sortie[v][2]) < 0.000001) {
				if (v > 0) {
					traj_sortie[v][1] = traj_sortie[v-1][1];
					traj_sortie[v][2] = traj_sortie[v-1][2];
				}
				else {
					u = 0;
					while (fabs(traj_sortie[u][1]) < 0.000001 && fabs(traj_sortie[u][2]) < 0.000001 && u < steps) {
//						fprintf(stderr,"traj_sortie[%i][1] = %f\t traj_sortie[%i][2]%f\n",u,traj_sortie[u][1],u,traj_sortie[u][2]);
						u++;
					}
					traj_sortie[0][0] = traj_sortie[u][0];
					traj_sortie[0][1] = traj_sortie[u][1];
					traj_sortie[0][2] = traj_sortie[u][2];
				}
			}
			traj_sortie[v][0] = v*dt;

			trajectories[h].d[w].i = 1.0;
			trajectories[h].d[w].x = traj_sortie[v][1];
			trajectories[h].d[w].y = traj_sortie[v][2];
			trajectories[h].d[w].t = traj_sortie[v][0];

			w++;

//			fprintf(stderr,"[u,v,w] = [%i,%i,%i] (%i)\n",u,v,w,steps);
		}

		n_sortie_regularise = (int) floor(t/dt);
	}
//	file->treeMeshGui->saveTrajectoriesButton->activate();
}

void Simulation::quadTreeInitialPosition(QuadTree *tree, int *increment, int *idx) {

//	int initPos;
//	int c = 0;
//	for (int b = 0; b < mesh->numberOfClusters; b++) {
//		if (mesh->activated(b)) {
//			c++;
//		}
//		/* PROBLEM IS HERE */
//		if (c == indice && c > 0) {
//			initPos = b;
//			indice = -1;
//			break;
//		}
//		if (indice == 0 && c == 1) {
//			initPos = b;
//			indice = -1;
//			break;
//		}
//	}

	// only count points in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) { *increment += 1; }
		if (*increment == *idx && *increment > 0) {
			initialLeaf = tree;
			*idx = -1;
			return;
		}
		if (*idx == 0 && *increment == 1) {
			initialLeaf = tree;
			*idx = -1;
			return;
		}
		return;
	}
	quadTreeInitialPosition(tree->nw,increment,idx);
	quadTreeInitialPosition(tree->ne,increment,idx);
	quadTreeInitialPosition(tree->sw,increment,idx);
	quadTreeInitialPosition(tree->se,increment,idx);

	return;
}

void Simulation::saveTrajectories() {
	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Simulated Trajectories Data");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Multi-Trajectory File\t*.trxyt\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".trxyt");
		FILE * writeFile;

		writeFile = fopen (nativeFilename,"wb");

		for (int b = 0; b < numberOfTrajectories; b++) {
			for (int a = 0; a < trajectories[b].n; a++) {
				fprintf(writeFile,"%i\t%lf\t%lf\t%lf\n",
						b,
						trajectories[b].d[a].x,
						trajectories[b].d[a].y,
						trajectories[b].d[a].t
				);
			}
		}

		fclose (writeFile);
	}
}

