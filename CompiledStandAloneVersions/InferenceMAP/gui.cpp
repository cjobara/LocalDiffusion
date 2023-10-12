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

#include "gui.h"
#include "file.h"
#include "globals.h"
#include "draw.h"
#include "inference.h"
#include "simulation.h"
#include "density.h"

extern Globals *iMAP;
extern File *file;

#define GTOL 1.e-16
#define K_R 4.e-7
#define PI 3.1415926535897932384626433832795028841971693993751

const Fl_Menu_Item optmizationMenu[] = {
	{"(D) Inference"},
	{"(D,F) Inference"},
	{"(D,Drift) Inference"},
	{"(D,V) Inference"},
	{"Polynomial Potential Inference"},
	{0}
};

const Fl_Menu_Item colormapItemsSquare[] = {
		{"autumn",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)0},
		{"blue",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)1},
		{"blue red",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)2},
		{"cool",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)3},
		{"cyan hot",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)4},
		{"fire",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)5},
		{"gray",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)6},
		{"green",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)7},
		{"hsv",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)8},
		{"jet",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)9},
		{"red",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)10},
		{"red hot",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)11},
		{"spring",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)12},
		{"summer",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)13},
		{"winter",0,(Fl_Callback*)colormapChoiceSquareCallback,(int*)14},
	{0},
};

const Fl_Menu_Item colormapItemsVoronoi[] = {
		{"autumn",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)0},
		{"blue",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)1},
		{"blue red",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)2},
		{"cool",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)3},
		{"cyan hot",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)4},
		{"fire",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)5},
		{"gray",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)6},
		{"green",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)7},
		{"hsv",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)8},
		{"jet",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)9},
		{"red",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)10},
		{"red hot",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)11},
		{"spring",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)12},
		{"summer",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)13},
		{"winter",0,(Fl_Callback*)colormapChoiceVoronoiCallback,(int*)14},
	{0},
};

const Fl_Menu_Item colormapItemsQuadTree[] = {
		{"autumn",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)0},
		{"blue",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)1},
		{"blue red",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)2},
		{"cool",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)3},
		{"cyan hot",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)4},
		{"fire",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)5},
		{"gray",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)6},
		{"green",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)7},
		{"hsv",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)8},
		{"jet",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)9},
		{"red",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)10},
		{"red hot",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)11},
		{"spring",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)12},
		{"summer",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)13},
		{"winter",0,(Fl_Callback*)colormapChoiceQuadTreeCallback,(int*)14},
	{0},
};

const Fl_Menu_Item clusteringMenu[] = {
	{"K-Means"},
	{"H-Means 1"},
	{"H-Means 2"}
};

void FileInfo::update() {

		char label [100];

		for (int c = 0; c < 100; c++) { label[c] = '\0'; }
		sprintf(label,"%i",file->numberOfFiles);
		filesBox->copy_label(label);

		for (int c = 0; c < 100; c++) { label[c] = '\0'; }
		sprintf(label,"%.1f",file->tMax - file->tMin);
		durationBox->copy_label(label);

		for (int c = 0; c < 100; c++) { label[c] = '\0'; }

		sprintf(label,"%.0f",1000.0f * (iMAP->exposureTime)); // fix
		acquisitionTimeBox->copy_label(label);

		for (int c = 0; c < 100; c++) { label[c] = '\0'; }
		sprintf(label,"%i / %i",file->localizationCountVoronoi,file->localizationCountDisplay);
		this->totalDetectionsBox->copy_label(label);

		for (int c = 0; c < 100; c++) { label[c] = '\0'; }
		sprintf(label,"%.1f x %.1f",file->xRange,file->yRange);
		this->dimensionsBox->copy_label(label);

		for (int c = 0; c < 100; c++) { label[c] = '\0'; }
		sprintf(label,"%.0f",1000.0*sqrt(file->averageDx*file->averageDx+file->averageDy*file->averageDy));
		this->averageStepBox->copy_label(label);

}

void FileInfo::updatePoints() {

		char label [100];
		for (int c = 0; c < 100; c++) { label[c] = '\0'; }
		sprintf(label,"%i / %i",file->localizationCountVoronoi,file->localizationCountDisplay);
		this->totalDetectionsBox->copy_label(label);

//		Fl::redraw();
}

void dPosteriorCallback(Fl_Button*w,void*v) {
	switch(file->meshType) {
		case 0: // square
		{
			if (w->value()) {
				file->squareMeshGui->fPosteriorButton->value(0);
				file->squareMeshGui->vPosteriorButton->value(0);
				file->squareMeshGui->potentialReferenceButton->value(0);
				file->squareMeshGui->potentialReferenceButton->deactivate();
				file->squareMeshGui->minPosteriorSlider->activate();

//				const double Dtemp = (file->squareMesh->getDiffusionMax()-file->squareMesh->getDiffusionMin())/2.0;
//				const int pts = file->squareMesh->getMinPointsPerCell();
//				const double min = Dtemp-0.8/sqrt((double)pts);
//				const double max = Dtemp+0.8/sqrt((double)pts);

				const double min = file->squareMesh->getDiffusionMin()*0.75;
				const double max = file->squareMesh->getDiffusionMax()*1.33;

				file->squareMeshGui->minPosteriorSlider->value(min);
				file->squareMeshGui->minPosteriorSlider->bounds(min,max);
				file->squareMeshGui->maxPosteriorSlider->value(max);
				file->squareMeshGui->maxPosteriorSlider->bounds(min,max);

				file->squareMeshGui->savePosteriorButton->deactivate();

			} else {
				w->value(1);
			}
		}
		break;
		case 1: // voronoi
		{
			if (w->value()) {
				file->voronoiMeshGui->fPosteriorButton->value(0);
				file->voronoiMeshGui->vPosteriorButton->value(0);
				file->voronoiMeshGui->potentialReferenceButton->value(0);
				file->voronoiMeshGui->potentialReferenceButton->deactivate();
				file->voronoiMeshGui->minPosteriorSlider->activate();

//				const double Dtemp = (file->voronoiMesh->getDiffusionMax()-file->voronoiMesh->getDiffusionMin())/2.0;
//				const int pts = file->voronoiMesh->getMinPointsPerCell();
//				const double min = Dtemp-0.8/sqrt((double)pts);
//				const double max = Dtemp+0.8/sqrt((double)pts);

				const double min = file->voronoiMesh->getDiffusionMin()*0.75;
				const double max = file->voronoiMesh->getDiffusionMax()*1.33;

				file->voronoiMeshGui->minPosteriorSlider->value(min);
				file->voronoiMeshGui->minPosteriorSlider->bounds(min,max);
				file->voronoiMeshGui->maxPosteriorSlider->value(max);
				file->voronoiMeshGui->maxPosteriorSlider->bounds(min,max);

				file->voronoiMeshGui->savePosteriorButton->deactivate();

//				fprintf(stderr,"current zone : %i\ndiffusion : %f\nforce : %f\n",file->voronoiMesh->getCurrentZone(),file->voronoiMesh->getDiffusion(file->voronoiMesh->getCurrentZone()),file->voronoiMesh->getForce(file->voronoiMesh->getCurrentZone()));
			} else {
				w->value(1);
			}
		}
		break;
		case 2: // quad tree
		{
			if (w->value()) {
				file->treeMeshGui->fPosteriorButton->value(0);
				file->treeMeshGui->vPosteriorButton->value(0);
				file->treeMeshGui->potentialReferenceButton->value(0);
				file->treeMeshGui->potentialReferenceButton->deactivate();
				file->treeMeshGui->minPosteriorSlider->activate();

//				const double Dtemp = (file->treeMesh->getDiffusionMax()-file->treeMesh->getDiffusionMin())/2.0;
//				const int pts = file->treeMesh->getMinPointsPerCell();
//				const double min = Dtemp-0.8/sqrt((double)pts);
//				const double max = Dtemp+0.8/sqrt((double)pts);

				const double min = file->treeMesh->getDiffusionMin()*0.75;
				const double max = file->treeMesh->getDiffusionMax()*1.33;

				file->treeMeshGui->minPosteriorSlider->value(min);
				file->treeMeshGui->minPosteriorSlider->bounds(min,max);
				file->treeMeshGui->maxPosteriorSlider->value(max);
				file->treeMeshGui->maxPosteriorSlider->bounds(min,max);

				file->treeMeshGui->savePosteriorButton->deactivate();

			} else {
				w->value(1);
			}
		}
		break;
	}
	Fl::redraw();
}

void fPosterioriCallback(Fl_Button*w,void*v) {
	switch(file->meshType) {
		case 0: // square
			if (w->value()) {
				file->squareMeshGui->dPosterioriButton->value(0);
				file->squareMeshGui->vPosteriorButton->value(0);
				file->squareMeshGui->potentialReferenceButton->value(0);
				file->squareMeshGui->potentialReferenceButton->deactivate();
				file->squareMeshGui->minPosteriorSlider->deactivate();

				const double min = 0.0;
				const double max = 20.0;
				file->squareMeshGui->minPosteriorSlider->value(min);
				file->squareMeshGui->minPosteriorSlider->bounds(min,max);
				file->squareMeshGui->maxPosteriorSlider->value(max);
				file->squareMeshGui->maxPosteriorSlider->bounds(min,max);

				file->squareMeshGui->savePosteriorButton->deactivate();

			} else {
				w->value(1);
			}
			break;
		case 1: // voronoi
			if (w->value()) {
				file->voronoiMeshGui->dPosterioriButton->value(0);
				file->voronoiMeshGui->vPosteriorButton->value(0);
				file->voronoiMeshGui->potentialReferenceButton->value(0);
				file->voronoiMeshGui->potentialReferenceButton->deactivate();
				file->voronoiMeshGui->minPosteriorSlider->deactivate();

				const double min = 0.0;
				const double max = 20.0;
				file->voronoiMeshGui->minPosteriorSlider->value(min);
				file->voronoiMeshGui->minPosteriorSlider->bounds(min,max);
				file->voronoiMeshGui->maxPosteriorSlider->value(max);
				file->voronoiMeshGui->maxPosteriorSlider->bounds(min,max);

				file->voronoiMeshGui->savePosteriorButton->deactivate();

			} else {
				w->value(1);
			}
			break;
		case 2: // quad tree
			if (w->value()) {
				file->treeMeshGui->dPosteriorButton->value(0);
				file->treeMeshGui->vPosteriorButton->value(0);
				file->treeMeshGui->potentialReferenceButton->value(0);
				file->treeMeshGui->potentialReferenceButton->deactivate();
				file->treeMeshGui->minPosteriorSlider->deactivate();

				const double min = 0.0;
				const double max = 20.0;
				file->treeMeshGui->minPosteriorSlider->value(min);
				file->treeMeshGui->minPosteriorSlider->bounds(min,max);
				file->treeMeshGui->maxPosteriorSlider->value(max);
				file->treeMeshGui->maxPosteriorSlider->bounds(min,max);

				file->treeMeshGui->savePosteriorButton->deactivate();

			} else {
				w->value(1);
			}
			break;
	}
	Fl::redraw();
}

void vPosterioriCallback(Fl_Button*w,void*v) {
	switch(file->meshType) {
		case 0: // square
			if (w->value()) {

				file->squareMeshGui->dPosterioriButton->value(0);
				file->squareMeshGui->fPosteriorButton->value(0);
				file->squareMeshGui->potentialReferenceButton->activate();

				file->squareMeshGui->minPosteriorSlider->activate();

//				const double Vtemp = (file->squareMesh->getPotentialMax()-file->squareMesh->getPotentialMin())/2.0;
//				const int pts = file->squareMesh->getMinPointsPerCell();
//				const double min = Vtemp-30.0/sqrt((double)pts);
//				const double max = Vtemp+30.0/sqrt((double)pts);

				const double min = file->squareMesh->getPotentialMin()-2.0;
				const double max = file->squareMesh->getPotentialMax()+2.0;

				file->squareMeshGui->minPosteriorSlider->value(min);
				file->squareMeshGui->minPosteriorSlider->bounds(min,max);
				file->squareMeshGui->maxPosteriorSlider->value(max);
				file->squareMeshGui->maxPosteriorSlider->bounds(min,max);

			} else {
				file->squareMeshGui->potentialReferenceButton->deactivate();
				file->squareMeshGui->potentialReferenceButton->value(0);
				w->value(1);
			}
			break;
		case 1: // voronoi
			if (w->value()) {
				file->voronoiMeshGui->dPosterioriButton->value(0);
				file->voronoiMeshGui->fPosteriorButton->value(0);
				file->voronoiMeshGui->potentialReferenceButton->activate();

				file->voronoiMeshGui->minPosteriorSlider->activate();

//				const double Vtemp = (file->voronoiMesh->getPotentialMax()-file->voronoiMesh->getPotentialMin())/2.0;
//				const int pts = file->voronoiMesh->getMinPointsPerCell();
//				const double min = Vtemp-30.0/sqrt((double)pts);
//				const double max = Vtemp+30.0/sqrt((double)pts);

				const double min = file->voronoiMesh->getPotentialMin()-2.0;
				const double max = file->voronoiMesh->getPotentialMax()+2.0;

				file->voronoiMeshGui->minPosteriorSlider->value(min);
				file->voronoiMeshGui->minPosteriorSlider->bounds(min,max);
				file->voronoiMeshGui->maxPosteriorSlider->value(max);
				file->voronoiMeshGui->maxPosteriorSlider->bounds(min,max);

			} else {
				file->voronoiMeshGui->potentialReferenceButton->deactivate();
				file->voronoiMeshGui->potentialReferenceButton->value(0);
				w->value(1);
			}
			break;
		case 2: // quad tree
			if (w->value()) {
				file->treeMeshGui->dPosteriorButton->value(0);
				file->treeMeshGui->fPosteriorButton->value(0);
				file->treeMeshGui->potentialReferenceButton->activate();

				file->treeMeshGui->minPosteriorSlider->activate();

//				const double Vtemp = (file->treeMesh->getPotentialMax()-file->treeMesh->getPotentialMin())/2.0;
//				const int pts = file->treeMesh->getMinPointsPerCell();
//				const double min = Vtemp-30.0/sqrt((double)pts);
//				const double max = Vtemp+30.0/sqrt((double)pts);

				const double min = file->treeMesh->getPotentialMin()-2.0;
				const double max = file->treeMesh->getPotentialMax()+2.0;

				file->treeMeshGui->minPosteriorSlider->value(min);
				file->treeMeshGui->minPosteriorSlider->bounds(min,max);
				file->treeMeshGui->maxPosteriorSlider->value(max);
				file->treeMeshGui->maxPosteriorSlider->bounds(min,max);

			} else {
				file->treeMeshGui->potentialReferenceButton->deactivate();
				file->treeMeshGui->potentialReferenceButton->value(0);
				w->value(1);
			}
			break;
	}
	Fl::redraw();
}

void potentialReferenceCallback(Fl_Button*w,void*v) {
	switch(file->meshType) {
		case 0: // square
			if (w->value()) {
				// get inferred delta based on two points
				const double min = 0.0;
				const double max = 10.0;
				file->squareMeshGui->minPosteriorSlider->value(0.0);
				file->squareMeshGui->minPosteriorSlider->bounds(min,max);
				file->squareMeshGui->maxPosteriorSlider->value(max/2.0);
				file->squareMeshGui->maxPosteriorSlider->bounds(min,max);

				file->squareMeshGui->samplePosteriorButton->activate();
				file->squareMeshGui->minPosteriorSlider->deactivate();
			}
			else {
				file->squareMeshGui->minPosteriorSlider->activate();
				const double min = file->squareMesh->getPotentialMin()-2.0;
				const double max = file->squareMesh->getPotentialMax()+2.0;

				file->squareMeshGui->minPosteriorSlider->value(min);
				file->squareMeshGui->minPosteriorSlider->bounds(min,max);
				file->squareMeshGui->maxPosteriorSlider->value(max);
				file->squareMeshGui->maxPosteriorSlider->bounds(min,max);
			}
			break;
		case 1: // voronoi
			if (w->value()) {
				const double min = 0.0;
				const double max = 10.0;
				file->voronoiMeshGui->minPosteriorSlider->value(0.0);
				file->voronoiMeshGui->minPosteriorSlider->bounds(min,max);
				file->voronoiMeshGui->maxPosteriorSlider->value(max/2.0);
				file->voronoiMeshGui->maxPosteriorSlider->bounds(min,max);

				file->voronoiMeshGui->samplePosteriorButton->activate();
				file->voronoiMeshGui->minPosteriorSlider->deactivate();
			}
			else {
				file->voronoiMeshGui->minPosteriorSlider->activate();
				const double min = file->voronoiMesh->getPotentialMin()-2.0;
				const double max = file->voronoiMesh->getPotentialMax()+2.0;

				file->voronoiMeshGui->minPosteriorSlider->value(min);
				file->voronoiMeshGui->minPosteriorSlider->bounds(min,max);
				file->voronoiMeshGui->maxPosteriorSlider->value(max);
				file->voronoiMeshGui->maxPosteriorSlider->bounds(min,max);
			}
			break;
		case 2: // quad tree
			if (w->value()) {
				const double min = 0.0;
				const double max = 10.0;
				file->treeMeshGui->minPosteriorSlider->value(0.0);
				file->treeMeshGui->minPosteriorSlider->bounds(min,max);
				file->treeMeshGui->maxPosteriorSlider->value(max/2.0);
				file->treeMeshGui->maxPosteriorSlider->bounds(min,max);

				file->treeMeshGui->samplePosteriorButton->activate();
				file->treeMeshGui->minPosteriorSlider->deactivate();
			}
			else {
				file->treeMeshGui->minPosteriorSlider->activate();
				const double min = file->treeMesh->getPotentialMin()-2.0;
				const double max = file->treeMesh->getPotentialMax()+2.0;

				file->treeMeshGui->minPosteriorSlider->value(min);
				file->treeMeshGui->minPosteriorSlider->bounds(min,max);
				file->treeMeshGui->maxPosteriorSlider->value(max);
				file->treeMeshGui->maxPosteriorSlider->bounds(min,max);
			}
			break;
	}
	Fl::redraw();
}

SquareMeshGui::SquareMeshGui() {

	if (iMAP->fileNumber > 0) {

		colormap = 9;
		
		const int w = 500;
		const int h = 300;

		window = new Fl_Window(w,h,"Square Meshing");
		window->callback(squareMeshWindowCallback);
		window->color(iMAP->bgColor);
		window->begin();
		window->set_non_modal();

		// Permanent interface buttons
		{
			overlayPointNumberButton = new Fl_Check_Button(5,245,65,20,"Points");
			overlayPointNumberButton->labelsize(12);
			overlayPointNumberButton->value(1);
			overlayPointNumberButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayPointNumberButton->labelfont(iMAP->normalFont);
			overlayPointNumberButton->labelcolor(FL_WHITE);
			overlayPointNumberButton->callback((Fl_Callback*)pointNumberOverlayCallback,(int*)0);
			overlayPointNumberButton->show();

			overlayDiffusionButton = new Fl_Check_Button(70,245,165,20,"Diffusion Coefficient");
			overlayDiffusionButton->labelsize(12);
			overlayDiffusionButton->callback((Fl_Callback*)diffusionOverlayCallback,(int*)0);
			overlayDiffusionButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayDiffusionButton->labelfont(iMAP->normalFont);
			overlayDiffusionButton->labelcolor(iMAP->bgColor);
			overlayDiffusionButton->deactivate();
			overlayDiffusionButton->show();

			overlayDiffusionLogButton = new Fl_Check_Button(70,260,165,20,"Log");
			overlayDiffusionLogButton->labelsize(12);
			overlayDiffusionLogButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayDiffusionLogButton->callback((Fl_Callback*)diffusionLogOverlayCallback,(int*)0);
			overlayDiffusionLogButton->labelfont(iMAP->normalFont);
			overlayDiffusionLogButton->labelcolor(iMAP->bgColor);
			overlayDiffusionLogButton->deactivate();
			overlayDiffusionLogButton->show();

			overlayForceMagnitudeButton = new Fl_Check_Button(240,245,135,20,"Force Magnitude");
			overlayForceMagnitudeButton->labelsize(12);
			overlayForceMagnitudeButton->callback((Fl_Callback*)forceMagnitudeOverlayCallback,(int*)0);
			overlayForceMagnitudeButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayForceMagnitudeButton->labelfont(iMAP->normalFont);
			overlayForceMagnitudeButton->labelcolor(iMAP->bgColor);
			overlayForceMagnitudeButton->deactivate();
			overlayForceMagnitudeButton->show();

			overlayForceArrowsButton = new Fl_Check_Button(240,260,135,20,"Arrows");
			overlayForceArrowsButton->labelsize(12);
			overlayForceArrowsButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayForceArrowsButton->labelfont(iMAP->normalFont);
			overlayForceArrowsButton->labelcolor(iMAP->bgColor);
			overlayForceArrowsButton->deactivate();
			overlayForceArrowsButton->show();

			overlayPotentialButton = new Fl_Check_Button(365,245,140,20,"Potential Energy");
			overlayPotentialButton->labelsize(12);
			overlayPotentialButton->callback((Fl_Callback*)potentialOverlayCallback,(int*)0);
			overlayPotentialButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayPotentialButton->labelfont(iMAP->normalFont);
			overlayPotentialButton->labelcolor(iMAP->bgColor);
			overlayPotentialButton->deactivate();
			overlayPotentialButton->show();

			variablesBox = new Fl_Box(0,275,w,25,"0 Zones");
			variablesBox->labelsize(14);
			variablesBox->align(FL_ALIGN_CENTER);
			variablesBox->labelfont(iMAP->boldFont);
			variablesBox->labelcolor(FL_WHITE);
			variablesBox->show();
		}

		tabs = new Fl_Tabs(0,0,w,h-60);
		tabs->labelsize(12);
//		tabs->color(FL_DARK1,FL_DARK1);
		tabs->box(FL_BORDER_FRAME);
		tabs->labelcolor(FL_WHITE);
		tabs->labelfont(iMAP->normalFont);
		tabs->begin();

		inferenceGroup = new Fl_Group(0,25,w,h-85,"Inference");
		inferenceGroup->labelsize(12);
//		inferenceGroup->color(FL_WHITE,iMAP->bgColor);
		inferenceGroup->labelcolor(FL_DARK1);
		inferenceGroup->labelfont(iMAP->boldFont);
		inferenceGroup->box(FL_BORDER_FRAME);
		inferenceGroup->begin();
		{
			dxdySlider = new Slider(10,45,w-20,20,"Side Length [nm]");
			dxdySlider->precision(0);
			dxdySlider->labelsize(11);
			dxdySlider->labelcolor(FL_WHITE);
			dxdySlider->value(sqrt(file->averageDx*file->averageDx+file->averageDy*file->averageDy)*1000);
			dxdySlider->labelfont(iMAP->boldFont);
			dxdySlider->bounds(sqrt(file->averageDx*file->averageDx+file->averageDy*file->averageDy)*1000,2*sqrt(file->averageDx*file->averageDx+file->averageDy*file->averageDy)*1000);
			dxdySlider->callback((Fl_Callback*)colorOverlayCallback,(int*)0);
			dxdySlider->show();

			minPointsSlider = new Slider(10,85,w-20,20,"Minimum Points / Zone");
			minPointsSlider->precision(0);
			minPointsSlider->labelsize(11);
			minPointsSlider->labelcolor(FL_WHITE);
			minPointsSlider->labelfont(iMAP->boldFont);
			minPointsSlider->value(20);
			minPointsSlider->callback((Fl_Callback*)minPointsSliderCallback);
			minPointsSlider->deactivate();
			minPointsSlider->show();

			inferenceModeChoice = new Fl_Choice(10,170,(w-24)/2,20,"Mode");
			inferenceModeChoice->menu(optmizationMenu);
			inferenceModeChoice->align(FL_ALIGN_TOP);
			inferenceModeChoice->textsize(11);
			inferenceModeChoice->textfont(iMAP->boldFont);
			inferenceModeChoice->labelcolor(FL_WHITE);
			inferenceModeChoice->labelfont(iMAP->boldFont);
//			inferenceModeChoice->box(FL_BORDER_FRAME);
			inferenceModeChoice->labelsize(11);
			inferenceModeChoice->color(FL_GRAY,FL_GRAY);
			inferenceModeChoice->callback((Fl_Callback*)optimizationChoiceCallback,(int*)0);
			inferenceModeChoice->deactivate();
			inferenceModeChoice->show();

			localizationPrecisionSlider = new Slider(15+(w-24)/2,170,(w-24)/2,20,"Localization Precision [nm]");
			localizationPrecisionSlider->labelsize(12);
			localizationPrecisionSlider->labelcolor(FL_WHITE);
			localizationPrecisionSlider->labelfont(iMAP->boldFont);
			localizationPrecisionSlider->precision(0);
			localizationPrecisionSlider->value(30);
			localizationPrecisionSlider->bounds(0,100);
			localizationPrecisionSlider->show();

			applyButton = new Fl_Button(10,200,(w-40)/6,30,"Apply");
			applyButton->labelsize(11);
			applyButton->labelfont(iMAP->boldFont);
			applyButton->callback((Fl_Callback*)applyMeshCallback,(int*)0);
			applyButton->show();

			const int buttonWidth = 76;
			resetButton = new Fl_Button(10+buttonWidth+5,200,buttonWidth,30,"Reset");
			resetButton->labelsize(11);
			resetButton->labelfont(iMAP->boldFont);
			resetButton->callback((Fl_Callback*)resetMeshCallback,(int*)0);
			resetButton->show();
			resetButton->deactivate();

			inferButton = new Fl_Button(10+2*buttonWidth+10,200,buttonWidth,30,"Infer");
			inferButton->labelsize(11);
			inferButton->labelfont(iMAP->boldFont);
			inferButton->callback((Fl_Callback*)inferMeshCallback,(int*)0);
			inferButton->show();
			inferButton->deactivate();

			pauseButton = new Fl_Button(10+3*buttonWidth+15,200,buttonWidth,30,"Pause");
			pauseButton->labelsize(14);
			pauseButton->label("@||");
			pauseButton->callback((Fl_Callback*)pauseCalculationButtonCallback,(int*)0);
			pauseButton->deactivate();
			pauseButton->show();

			stopButton = new Fl_Button(10+4*buttonWidth+20,200,buttonWidth,30,"Stop");
			stopButton->labelsize(12);
			stopButton->label("@square");
			stopButton->labelcolor(FL_RED);
			stopButton->callback((Fl_Callback*)stopCalculationButtonCallback,(int*)0);
			stopButton->deactivate();
			stopButton->show();

			saveButton = new Fl_Button(10+5*buttonWidth+25,200,buttonWidth,30,"Save");
			saveButton->labelsize(12);
			saveButton->labelfont(iMAP->boldFont);
			saveButton->callback((Fl_Callback*)saveSquareMeshCallback);
			saveButton->deactivate();
			saveButton->show();
		}
		inferenceGroup->end();

		annotationsGroup = new Fl_Group(0,25,w,h-85,"Overlay");
		annotationsGroup->labelsize(12);
		annotationsGroup->labelcolor(FL_DARK1);
		annotationsGroup->box(FL_BORDER_FRAME);
		annotationsGroup->labelfont(iMAP->normalFont);
		annotationsGroup->begin();
		{
			// font options
			localizationNumberLabelButton = new Fl_Check_Button(10,45,150,20,"Localization Number");
			localizationNumberLabelButton->labelsize(12);
			localizationNumberLabelButton->labelcolor(FL_WHITE);
			localizationNumberLabelButton->labelfont(iMAP->boldFont);
			localizationNumberLabelButton->value(1);
			localizationNumberLabelButton->show();

			fontScaleSlider = new Slider(175,45,w-185,20,"Font Scale");
			fontScaleSlider->labelsize(12);
			fontScaleSlider->precision(2);
			fontScaleSlider->bounds(0.5f,1.5f);
			fontScaleSlider->value(1.0f);
			fontScaleSlider->labelfont(iMAP->boldFont);
			fontScaleSlider->labelcolor(FL_WHITE);
			fontScaleSlider->show();

			// grid options
			displayGridButton = new Fl_Check_Button(10,85,105,20,"Display Grid");
			displayGridButton->labelsize(12);
			displayGridButton->labelcolor(FL_WHITE);
			displayGridButton->labelfont(iMAP->boldFont);
			displayGridButton->value(1);
			displayGridButton->show();

			gridColorButton = new Fl_Button(120,85,85,20,"Grid Color");
			gridColorButton->labelsize(12);
			gridColorButton->labelfont(iMAP->boldFont);
			gridColorButton->callback((Fl_Callback*)gridColorButtonCallback,(int*)0);
			gridColorButton->show();

			gridRGB[0] = gridRGB[1] = gridRGB[2] = 1.0;

			gridAlphaSlider = new Slider(215,85,w-225,20,"Grid Alpha");
			gridAlphaSlider->labelsize(12);
			gridAlphaSlider->precision(2);
			gridAlphaSlider->labelcolor(FL_WHITE);
			gridAlphaSlider->labelfont(iMAP->boldFont);
			gridAlphaSlider->value(1.0f);
			gridAlphaSlider->bounds(0.0,1.0);
			gridAlphaSlider->show();

			// color overlay options
			overlayButton = new Fl_Check_Button(10,125,115,20,"Color Overlay");
			overlayButton->labelsize(12);
			overlayButton->labelcolor(FL_WHITE);
			overlayButton->labelfont(iMAP->boldFont);
			overlayButton->value(1);
			overlayButton->show();

			overlayColorChoice = new Fl_Choice(155,125,70,20,"Map");
			overlayColorChoice->labelsize(12);
			overlayColorChoice->menu(colormapItemsSquare);
			overlayColorChoice->textsize(12);
			overlayColorChoice->labelcolor(FL_WHITE);
			overlayColorChoice->textfont(iMAP->boldFont);
			overlayColorChoice->labelfont(iMAP->boldFont);
			overlayColorChoice->callback((Fl_Callback*)colorOverlayCallback,(int*)0);
			overlayColorChoice->value(colormap); // jet by default
			overlayColorChoice->show();

			flipColormapButton = new Fl_Check_Button(225,125,40,20,"Flip");
			flipColormapButton->callback((Fl_Callback*)flipColormapCallback,(int*)0);
			flipColormapButton->labelsize(12);
			flipColormapButton->labelcolor(FL_WHITE);
			flipColormapButton->labelfont(iMAP->boldFont);
			flipColormapButton->value(0);
			flipColormapButton->show();

			overlayAlphaSlider = new Slider(280,125,w-290,20,"Overlay Alpha");
			overlayAlphaSlider->labelsize(12);
			overlayAlphaSlider->precision(2);
			overlayAlphaSlider->bounds(0.0f,1.0f);
			overlayAlphaSlider->value(0.5f);
			overlayAlphaSlider->labelfont(iMAP->boldFont);
			overlayAlphaSlider->labelcolor(FL_WHITE);
			overlayAlphaSlider->callback((Fl_Callback*)colorOverlayCallback,(int*)0);
			overlayAlphaSlider->show();

			// spot visualization options
			spotVisualizationButton = new Fl_Check_Button(10,165,135,20,"Spot Visualization");
			spotVisualizationButton->labelsize(12);
			spotVisualizationButton->labelcolor(FL_WHITE);
			spotVisualizationButton->labelfont(iMAP->boldFont);
			spotVisualizationButton->value(0);
			spotVisualizationButton->show();

			spotVisualizationScaleSlider = new Slider(160,165,w-170,20,"Scale");
			spotVisualizationScaleSlider->labelsize(12);
			spotVisualizationScaleSlider->precision(2);
			spotVisualizationScaleSlider->bounds(0.0f,5.0f);
			spotVisualizationScaleSlider->labelcolor(FL_WHITE);
			spotVisualizationScaleSlider->labelfont(iMAP->boldFont);
			spotVisualizationScaleSlider->value(1.0f);
			spotVisualizationScaleSlider->show();

			// hover information
			infoOverlayButton = new Fl_Check_Button(10,190,140,20,"Hover Information");
			infoOverlayButton->labelsize(12);
			infoOverlayButton->value(0);
			infoOverlayButton->labelfont(iMAP->boldFont);
			infoOverlayButton->labelcolor(FL_WHITE);
			infoOverlayButton->deactivate();
			infoOverlayButton->show();

			neighbourDistanceViewButton = new Fl_Check_Button(10,210,195,20,"View Neighbor Connections");
			neighbourDistanceViewButton->labelcolor(FL_WHITE);
			neighbourDistanceViewButton->labelfont(iMAP->boldFont);
			neighbourDistanceViewButton->labelsize(12);
			neighbourDistanceViewButton->value(1);
			neighbourDistanceViewButton->show();

		}
		annotationsGroup->deactivate();
		annotationsGroup->end();

		advancedGroup = new Fl_Group(0,25,w,h-85,"Advanced");
		advancedGroup->labelsize(12);
		advancedGroup->labelcolor(FL_DARK1);
		advancedGroup->box(FL_BORDER_FRAME);
		advancedGroup->labelfont(iMAP->normalFont);
		advancedGroup->begin();
		{

			roButton = new Fl_Check_Button(10,30,180,20,"Randomized Optimization");
			roButton->labelcolor(FL_WHITE);
			roButton->labelfont(iMAP->boldFont);
			roButton->labelsize(12);
			roButton->value(0);
			roButton->callback((Fl_Callback*)roButtonCallback,(int*)0);
			roButton->deactivate();
			roButton->show();

			roRadiusSlider = new Slider(10,65,(w-20)/3-5,20,"Selection Radius [nm]");
			roRadiusSlider->labelsize(12);
			roRadiusSlider->labelcolor(FL_WHITE);
			roRadiusSlider->labelfont(iMAP->boldFont);
			roRadiusSlider->precision(0);
			roRadiusSlider->value(500);
			roRadiusSlider->bounds(100,2000);
			roRadiusSlider->deactivate();
			roRadiusSlider->show();

			roToleranceSlider = new Slider(12+(w-20)/3,65,w/3-10,20,"Cost Tolerance [%]");
			roToleranceSlider->labelsize(12);
			roToleranceSlider->labelcolor(FL_WHITE);
			roToleranceSlider->labelfont(iMAP->boldFont);
			roToleranceSlider->precision(4);
			roToleranceSlider->value(0.001);
			roToleranceSlider->bounds(0.0,1.0);
			roToleranceSlider->callback((Fl_Callback*)randomizedOptimizationToleranceCallback,(int*)0);
			roToleranceSlider->deactivate();
			roToleranceSlider->show();

			roMaximumIterationsSlider = new Slider(14+2*(w-20)/3,65,w/3-10,20,"Maximum Iterations");
			roMaximumIterationsSlider->labelsize(12);
			roMaximumIterationsSlider->labelcolor(FL_WHITE);
			roMaximumIterationsSlider->labelfont(iMAP->boldFont);
			roMaximumIterationsSlider->precision(0);
			roMaximumIterationsSlider->value(5);
			roMaximumIterationsSlider->bounds(3,10);
			roMaximumIterationsSlider->deactivate();
			roMaximumIterationsSlider->show();

			neighbourDistanceSlider = new Slider(10,105,w-20,20,"Maximum Neighbor Distance [nm]");
			neighbourDistanceSlider->labelsize(12);
			neighbourDistanceSlider->labelcolor(FL_WHITE);
			neighbourDistanceSlider->labelfont(iMAP->boldFont);
			neighbourDistanceSlider->precision(1);
			neighbourDistanceSlider->value(1000.0);
			neighbourDistanceSlider->bounds(50,2000);
			neighbourDistanceSlider->callback((Fl_Callback*)maximumNeighbourDistanceCallback,(int*)0);
			neighbourDistanceSlider->show();

			Fl_Box *dfLabel = new Fl_Box(15,145,100,20,"(D,F) Inference:");
			dfLabel->labelsize(12);
			dfLabel->labelcolor(FL_WHITE);
			dfLabel->labelfont(iMAP->boldFont);
			dfLabel->show();

			betaSlider = new Slider(130,145,w-140,20,"Potential Energy Penalization (Beta)");
			betaSlider->precision(1);
			betaSlider->labelcolor(FL_WHITE);
			betaSlider->labelfont(iMAP->boldFont);
			betaSlider->bounds(0.0f,5.0f);
			betaSlider->value(2.0f);
			betaSlider->deactivate();
			betaSlider->show();

			Fl_Box *polynomialLabel = new Fl_Box(10,185,220,20,"Polynomial Potential Inference:");
			polynomialLabel->labelsize(12);
			polynomialLabel->labelcolor(FL_WHITE);
			polynomialLabel->labelfont(iMAP->boldFont);
			polynomialLabel->show();

			polynomialOrderSlider = new Slider(235,185,w-245,20,"Polynomial Order");
			polynomialOrderSlider->precision(0);
			polynomialOrderSlider->labelsize(12);
			polynomialOrderSlider->labelcolor(FL_WHITE);
			polynomialOrderSlider->value(2);
			polynomialOrderSlider->labelfont(iMAP->boldFont);
			polynomialOrderSlider->bounds(2,8);
			polynomialOrderSlider->show();
			polynomialOrderSlider->deactivate();

		}
		advancedGroup->deactivate();
		advancedGroup->end();

		priorGroup = new Fl_Group(0,25,w,h-85,"Prior");
		priorGroup->labelsize(12);
		priorGroup->labelcolor(FL_DARK1);
		priorGroup->box(FL_BORDER_FRAME);
		priorGroup->labelfont(iMAP->normalFont);
		priorGroup->begin();
		{

			jeffreysPriorButton = new Fl_Check_Button(10,30,180,25,"Enable Jeffreys Prior");
			jeffreysPriorButton->labelsize(12);
			jeffreysPriorButton->labelcolor(FL_WHITE);
			jeffreysPriorButton->labelfont(iMAP->boldFont);
			jeffreysPriorButton->value(1);
			jeffreysPriorButton->callback((Fl_Callback*)jeffreysPriorButtonCallback,(int*)0);
			jeffreysPriorButton->show();

			smoothingPriorButton = new Fl_Check_Button(10,50,180,25,"Enable Smoothing Prior");
			smoothingPriorButton->labelsize(12);
			smoothingPriorButton->labelcolor(FL_WHITE);
			smoothingPriorButton->labelfont(iMAP->boldFont);
			smoothingPriorButton->value(0);
			smoothingPriorButton->callback((Fl_Callback*)smoothingPriorButtonCallback,(int*)0);
			smoothingPriorButton->show();
//			smoothingPriorButton->deactivate();

			lambdaSlider = new Slider(10,85,(w-30)/2,20,"V Penalisation (lambda)");
			lambdaSlider->precision(3);
			lambdaSlider->labelcolor(FL_WHITE);
			lambdaSlider->labelfont(iMAP->boldFont);
			lambdaSlider->labelsize(12);
			lambdaSlider->value(0.1);
			lambdaSlider->bounds(0.0,10);
			lambdaSlider->deactivate();
			lambdaSlider->show();

			muSlider = new Slider((w-30)/2+15,85,(w-30)/2,20,"D Penalization (mu)");
			muSlider->precision(3);
			muSlider->labelsize(12);
			muSlider->labelcolor(FL_WHITE);
			muSlider->labelfont(iMAP->boldFont);
			muSlider->value(0.1);
			muSlider->bounds(0.0,10.0);
			muSlider->deactivate();
			muSlider->show();

		}
//		priorGroup->deactivate();
		priorGroup->end();

		posterioriGroup = new Fl_Group(0,25,w,h-85,"Posterior");
		posterioriGroup->labelsize(12);
		posterioriGroup->labelcolor(FL_DARK1);
		posterioriGroup->box(FL_BORDER_FRAME);
		posterioriGroup->labelfont(iMAP->normalFont);
		posterioriGroup->begin();
		{
			dPosterioriButton = new Fl_Round_Button(10,30,155,20,"Diffusion");
			dPosterioriButton->labelfont(1);
			dPosterioriButton->labelsize(12);
			dPosterioriButton->labelfont(iMAP->boldFont);
			dPosterioriButton->type(FL_TOGGLE_BUTTON);
			dPosterioriButton->labelcolor(FL_WHITE);
			dPosterioriButton->callback((Fl_Callback*)dPosteriorCallback);
			dPosterioriButton->value(0);
			dPosterioriButton->deactivate();
			dPosterioriButton->show();

			fPosteriorButton = new Fl_Round_Button(10,50,155,20,"Force");
			fPosteriorButton->labelfont(iMAP->boldFont);
			fPosteriorButton->labelsize(12);
			fPosteriorButton->type(FL_TOGGLE_BUTTON);
			fPosteriorButton->labelcolor(FL_WHITE);
			fPosteriorButton->callback((Fl_Callback*)fPosterioriCallback);
			fPosteriorButton->value(0);
			fPosteriorButton->deactivate();
			fPosteriorButton->show();

			vPosteriorButton = new Fl_Round_Button(10,70,82,20,"Potential");
			vPosteriorButton->labelfont(iMAP->boldFont);
			vPosteriorButton->labelsize(12);
			vPosteriorButton->type(FL_TOGGLE_BUTTON);
			vPosteriorButton->labelcolor(FL_WHITE);
			vPosteriorButton->callback((Fl_Callback*)vPosterioriCallback);
			vPosteriorButton->value(0);
			vPosteriorButton->deactivate();
			vPosteriorButton->show();

			potentialReferenceButton = new Fl_Button(95,70,70,20,"Reference");
			potentialReferenceButton->labelfont(iMAP->boldFont);
			potentialReferenceButton->labelsize(12);
			potentialReferenceButton->type(FL_TOGGLE_BUTTON);
			potentialReferenceButton->callback((Fl_Callback*)potentialReferenceCallback);
			potentialReferenceButton->value(0);
			potentialReferenceButton->deactivate();
			potentialReferenceButton->show();

			minPosteriorSlider = new Slider(10,105,155,20,"Minimum");
			minPosteriorSlider->labelsize(12);
			minPosteriorSlider->labelfont(iMAP->boldFont);
			minPosteriorSlider->precision(3);
			minPosteriorSlider->labelcolor(FL_WHITE);
			minPosteriorSlider->value(0.001f);
			minPosteriorSlider->bounds(0.001,1.0);
			minPosteriorSlider->deactivate();
			minPosteriorSlider->show();

			maxPosteriorSlider = new Slider(10,140,155,20,"Maximum");
			maxPosteriorSlider->labelsize(12);
			maxPosteriorSlider->labelfont(iMAP->boldFont);
			maxPosteriorSlider->precision(3);
			maxPosteriorSlider->labelcolor(FL_WHITE);
			maxPosteriorSlider->value(1.0f);
			maxPosteriorSlider->bounds(0.001,1.0);
			maxPosteriorSlider->deactivate();
			maxPosteriorSlider->show();

			posterioriSampleNumberSlider = new Slider(10,175,155,20,"Samples");
			posterioriSampleNumberSlider->labelsize(12);
			posterioriSampleNumberSlider->labelfont(iMAP->boldFont);
			posterioriSampleNumberSlider->precision(0);
			posterioriSampleNumberSlider->labelcolor(FL_WHITE);
			posterioriSampleNumberSlider->bounds(100,1000);
			posterioriSampleNumberSlider->value(100);
			posterioriSampleNumberSlider->deactivate();
			posterioriSampleNumberSlider->show();

			samplePosteriorButton = new Fl_Button(10,205,75,25,"Sample");
			samplePosteriorButton->labelfont(iMAP->boldFont);
			samplePosteriorButton->labelsize(12);
			samplePosteriorButton->callback((Fl_Callback*)samplePosteriorCallback,(int*)0);
//			samplePosterioriButton->deactivate();
			samplePosteriorButton->show();

			savePosteriorButton = new Fl_Button(90,205,75,25,"Save");
			savePosteriorButton->labelfont(iMAP->boldFont);
			savePosteriorButton->labelsize(12);
			savePosteriorButton->callback((Fl_Callback*)savePosteriorCallback,(int*)0);
			savePosteriorButton->deactivate();
			savePosteriorButton->show();

			// draw box
			Fl_Box *posterioriBox = new Fl_Box(170,35,320,195);
			posterioriBox->color(FL_WHITE);

		}
		posterioriGroup->deactivate();
		posterioriGroup->end();

		simulationGroup = new Fl_Group(0,25,w,h-85,"Simulation");
		simulationGroup->labelsize(12);
//		simulationGroup->color(FL_WHITE,FL_WHITE);
		simulationGroup->labelcolor(FL_DARK1);
		simulationGroup->box(FL_BORDER_FRAME);
		simulationGroup->labelfont(iMAP->normalFont);
		simulationGroup->begin();
		{
			numberOfTrajectoriesSlider = new Slider(10,45,w-20,20,"Number of Trajectories");
			numberOfTrajectoriesSlider->precision(0);
			numberOfTrajectoriesSlider->labelfont(iMAP->boldFont);
			numberOfTrajectoriesSlider->labelcolor(FL_WHITE);
			numberOfTrajectoriesSlider->bounds(100,1000);
			numberOfTrajectoriesSlider->value(500);
//			numberOfTrajectoriesSlider->callback((Fl_Callback*)simulationWidgetsCallback,(int*)0);
			numberOfTrajectoriesSlider->activate();

			deltaTSlider = new Slider(10,85,(w-30)/2,20,"Delta [ms]");
			deltaTSlider->labelsize(12);
			deltaTSlider->precision(0);
			deltaTSlider->labelfont(iMAP->boldFont);
			deltaTSlider->labelcolor(FL_WHITE);
			deltaTSlider->value(25);
			deltaTSlider->bounds(0,100);
//			deltaTSlider->callback((Fl_Callback*)simulationWidgetsCallback,(int*)0);
			deltaTSlider->activate();

			timeStepsSlider = new Slider((w-30)/2+20,85,(w-30)/2,20,"Maximum Time Steps");
			timeStepsSlider->labelsize(12);
			timeStepsSlider->labelfont(iMAP->boldFont);
			timeStepsSlider->labelcolor(FL_WHITE);
			timeStepsSlider->precision(0);
			timeStepsSlider->value(25);
			timeStepsSlider->bounds(0,100);
//			timeStepsSlider->callback((Fl_Callback*)simulationWidgetsCallback,(int*)0);
			timeStepsSlider->activate();

			generateTrajectoriesButton = new Fl_Button(w/2-50,115,100,40,"Save\nTrajectories");
			generateTrajectoriesButton->labelsize(12);
			generateTrajectoriesButton->labelfont(iMAP->boldFont);
			generateTrajectoriesButton->callback((Fl_Callback*)generateTrajectoriesCallback,(int*)0);
			generateTrajectoriesButton->show();
			generateTrajectoriesButton->activate();

		}
		simulationGroup->deactivate();
		simulationGroup->end();

		landscapeGroup = new Fl_Group(0,25,w,h-85,"Landscape");
		landscapeGroup->labelsize(12);
//		landscapeGroup->color(FL_WHITE,FL_WHITE);
		landscapeGroup->labelcolor(FL_DARK1);
		landscapeGroup->box(FL_BORDER_FRAME);
		landscapeGroup->labelfont(iMAP->normalFont);
		landscapeGroup->begin();
		{
			landscapeButton = new Fl_Check_Button(10,30,120,25,"View Landscape");
			landscapeButton->labelsize(12);
			landscapeButton->labelfont(iMAP->boldFont);
			landscapeButton->labelcolor(FL_WHITE);
			landscapeButton->value(0);
			landscapeButton->callback((Fl_Callback*)squareLandscapeCallback);
//			landscapeButton->deactivate();
			landscapeButton->show();

			scaleLandscapeSlider = new Slider(10,70,w/2-15,20,"Scale");
			scaleLandscapeSlider->precision(2);
			scaleLandscapeSlider->labelfont(iMAP->boldFont);
			scaleLandscapeSlider->labelcolor(FL_WHITE);
			scaleLandscapeSlider->bounds(0,5.0);
			scaleLandscapeSlider->value(1.0);
			scaleLandscapeSlider->deactivate();
			scaleLandscapeSlider->show();

			landscapeAlphaSlider = new Slider(w/2+5,70,w/2-15,20,"Alpha");
			landscapeAlphaSlider->precision(2);
			landscapeAlphaSlider->bounds(0.01,1);
			landscapeAlphaSlider->labelfont(iMAP->boldFont);
			landscapeAlphaSlider->labelcolor(FL_WHITE);
			landscapeAlphaSlider->value(1.0);
			landscapeAlphaSlider->callback((Fl_Callback*)alphaSquareLandscapeCallback);
			landscapeAlphaSlider->deactivate();
			landscapeAlphaSlider->show();

			xLightPositionSlider = new Slider(10,105,w/2-15,20,"Light x");
			xLightPositionSlider->precision(2);
			xLightPositionSlider->bounds(-100,100);
			xLightPositionSlider->labelfont(iMAP->boldFont);
			xLightPositionSlider->labelcolor(FL_WHITE);
			xLightPositionSlider->value(0.0);
			xLightPositionSlider->deactivate();
			xLightPositionSlider->show();

			yLightPositionSlider = new Slider(10,140,w/2-15,20,"Light y");
			yLightPositionSlider->precision(2);
			yLightPositionSlider->bounds(-100,100);
			yLightPositionSlider->labelfont(iMAP->boldFont);
			yLightPositionSlider->labelcolor(FL_WHITE);
			yLightPositionSlider->value(0.0);
			yLightPositionSlider->deactivate();
			yLightPositionSlider->show();

			zLightPositionSlider = new Slider(10,175,w/2-15,20,"Light z");
			zLightPositionSlider->precision(2);
			zLightPositionSlider->labelfont(iMAP->boldFont);
			zLightPositionSlider->labelcolor(FL_WHITE);
			zLightPositionSlider->bounds(-100,100);
			zLightPositionSlider->value(100.0);
			zLightPositionSlider->deactivate();
			zLightPositionSlider->show();

			fogButton = new Fl_Check_Button(w/2+5,105,50,25,"Fog");
			fogButton->labelsize(12);
			fogButton->labelfont(iMAP->boldFont);
			fogButton->labelcolor(FL_WHITE);
			fogButton->value(0);
			fogButton->deactivate();
			fogButton->show();

			fogStartSlider = new Slider(w/2+5,140,w/2-15,20,"Fog Start");
			fogStartSlider->precision(2);
			fogStartSlider->bounds(0,100);
			fogStartSlider->labelfont(iMAP->boldFont);
			fogStartSlider->labelcolor(FL_WHITE);
			fogStartSlider->value(0.0);
			fogStartSlider->deactivate();
			fogStartSlider->show();

			fogEndSlider = new Slider(w/2+5,175,w/2-15,20,"Fog End");
			fogEndSlider->precision(2);
			fogEndSlider->labelfont(iMAP->boldFont);
			fogEndSlider->labelcolor(FL_WHITE);
			fogEndSlider->bounds(0,100);
			fogEndSlider->value(100.0);
			fogEndSlider->deactivate();
			fogEndSlider->show();

			landscapeAxisButton = new Fl_Check_Button(10,200,50,20,"Axes");
			landscapeAxisButton->labelsize(12);
			landscapeAxisButton->labelfont(iMAP->boldFont);
			landscapeAxisButton->labelcolor(FL_WHITE);
			landscapeAxisButton->value(1);
			landscapeAxisButton->callback((Fl_Callback*)axesCallback);
			landscapeAxisButton->deactivate();
			landscapeAxisButton->show();
		}
		landscapeGroup->deactivate();
		landscapeGroup->end();

		tabs->end();
		meshApplied = false;
	}

}

void SquareMesh::savePosterior() {
	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	switch(file->squareMeshGui->getPosteriorType()) {
	case 0: // diffusion
		native->title("Save Diffusion Posterior Data");
		native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
		native->filter("Diffusion Posterior File\t*.dpost\n");
	    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
	    native->preset_file(file->fileName);
		if (native->show() == 0) {
			// store filename
			sprintf(nativeFilename,"%s%s",native->filename(),".dpost");
			FILE * writeFile;

			writeFile = fopen (nativeFilename,"wb");

			// write header
			fprintf(writeFile,"File Name: %s\n\n"
							  "Number of Points: %i\n"
							  "Zone Centre: [%.3f,%.3f]\n"
							  "D_MAP: %f [um2/s]\n\n",
							  file->fileName,
							  getCount(getCurrentZoneX(),getCurrentZoneY()),
							  getCentroidX(getCurrentZoneX(),getCurrentZoneY()),
							  getCentroidY(getCurrentZoneX(),getCurrentZoneY()),
							  getDiffusion(getCurrentZoneX(),getCurrentZoneY()));

			for (int f = 0; f < file->squareMeshGui->posterioriSamples; f++) {
				fprintf(writeFile,"%f\t%f\n",file->squareMeshGui->minPosteriorSlider->value() + (double)f*file->squareMeshGui->posterioriIncrement,file->squareMeshGui->posterioriValues[f]);
			}
			fclose(writeFile);
		}
		break;
	case 1: // force
		native->title("Save Force Posterior Data");
		native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
		native->filter("Force Posterior File\t*.fpost\n");
	    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
	    native->preset_file(file->fileName);
		if (native->show() == 0) {
			// store filename
			sprintf(nativeFilename,"%s%s",native->filename(),".fpost");
			FILE * writeFile;

			writeFile = fopen (nativeFilename,"wb");

			// write header
			fprintf(writeFile,"File Name: %s\n\n"
							  "Number of Points: %i\n"
							  "Zone Centre: [%.3f,%.3f]\n"
							  "F_MAP: %f [pN]\n\n",
							  file->fileName,
							  getCount(getCurrentZoneX(),getCurrentZoneY()),
							  getCentroidX(getCurrentZoneX(),getCurrentZoneY()),
							  getCentroidY(getCurrentZoneX(),getCurrentZoneY()),
							  getForce(getCurrentZoneX(),getCurrentZoneY()));

			const double minF2 = getForce(getCurrentZoneX(),getCurrentZoneY()) - file->squareMeshGui->posterioriRange/2.0;

			for (int f = 0; f < file->squareMeshGui->posterioriSamples; f++) {
				fprintf(writeFile,"%f\t%f\n",minF2 + (double)f*file->squareMeshGui->posterioriIncrement,file->squareMeshGui->posterioriValues[f]);
			}
			fclose(writeFile);
		}
		break;
	case 2: // potential
		native->title("Save Potential Posterior Data");
		native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
		native->filter("Potential Posterior File\t*.vpost\n");
	    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
	    native->preset_file(file->fileName);
		if (native->show() == 0) {
			// store filename
			sprintf(nativeFilename,"%s%s",native->filename(),".vpost");
			FILE * writeFile;

			writeFile = fopen (nativeFilename,"wb");

			if (file->squareMeshGui->potentialReferenceButton->value()) {
				// write header
				const double deltaV = file->squareMesh->getPotential(file->squareMesh->getCurrentZone2X(),file->squareMesh->getCurrentZone2Y())
									  -file->squareMesh->getPotential(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
				fprintf(writeFile,"File Name: %s\n\n"
								  "Number of Points: %i\n"
								  "Zone 1 Centre: [%.3f,%.3f]\n"
						          "Zone 2 Centre: [%.3f,%.3f]\n"
								  "V_MAP (reference): %f [kT]\n\n",
								  file->fileName,
								  getCount(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY()),
								  getCount(file->squareMesh->getCurrentZone2X(),file->squareMesh->getCurrentZone2Y()),
								  getCentroidX(getCurrentZoneX(),getCurrentZoneY()),
								  getCentroidY(getCurrentZoneX(),getCurrentZoneY()),
								  deltaV);

				for (int f = 0; f < 2*file->squareMeshGui->posterioriSamples; f++) {
//					fprintf(writeFile,"%f\t%f\n",deltaV-file->squareMeshGui->maxPosteriorSlider->value() + (double)f*file->squareMeshGui->posterioriIncrement,file->squareMeshGui->posterioriValues[f]);
					fprintf(writeFile,"%f\t%f\n",deltaV-file->squareMeshGui->maxPosteriorSlider->value()/2.0 + (float)f/(2*file->squareMeshGui->posterioriSamples)*file->squareMeshGui->posterioriRange,file->squareMeshGui->posterioriValues[f]);
				}
			} else {
				// write header
				fprintf(writeFile,"File Name: %s\n\n"
								  "Number of Points: %i\n"
								  "Zone Centre: [%.3f,%.3f]\n"
								  "V_MAP (no reference): %f [kT]\n\n",
								  file->fileName,
								  getCount(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY()),
								  getCentroidX(getCurrentZoneX(),getCurrentZoneY()),
								  getCentroidY(getCurrentZoneX(),getCurrentZoneY()),
								  getPotential(getCurrentZoneX(),getCurrentZoneY()));

				for (int f = 0; f < file->squareMeshGui->posterioriSamples; f++) {
					fprintf(writeFile,"%f\t%f\n",file->squareMeshGui->minPosteriorSlider->value() + (double)f*file->squareMeshGui->posterioriIncrement,file->squareMeshGui->posterioriValues[f]);
				}
			}
			fclose(writeFile);
		}
		break;
	}
}

void samplePosteriorCallback(Fl_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {
	case 0:
		// find out which variable to sample
		file->squareMeshGui->samplePosterior(file->squareMeshGui->getPosteriorType());
		file->squareMeshGui->savePosteriorButton->activate();
		break;
	case 1:
		file->voronoiMeshGui->samplePosterior(file->voronoiMeshGui->getPosteriorType());
		file->voronoiMeshGui->savePosteriorButton->activate();
		break;
	case 2:
		file->treeMeshGui->samplePosterior(file->treeMeshGui->getPosteriorType());
		file->treeMeshGui->savePosteriorButton->activate();
		break;
	}
}

void maximumNeighbourDistanceCallback(Slider*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {
	case 0:
		file->squareMesh->setMaximumNeighbourDistance(w->value()/1000.0);
		file->squareMesh->updateNeighbours();
		break;
	case 1:
		file->voronoiMesh->setMaximumNeighbourDistance(w->value()/1000.0);
		file->voronoiMesh->updateNeighbours();
		break;
	case 2:
		file->treeMesh->setMaximumNeighbourDistance(w->value()/1000.0);
		file->treeMesh->updateNeighbours(file->treeMesh->quadTree);
		break;
	}
	Fl::redraw();
}

void roButtonCallback(Fl_Check_Button*w,void*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {
	case 0:
		if (w->value()) {
			file->squareMeshGui->roToleranceSlider->activate();
			file->squareMeshGui->roMaximumIterationsSlider->activate();
			file->squareMeshGui->roRadiusSlider->activate();
			if (file->squareMeshOverlay) { file->squareMesh->enableRo(); }
		}
		else {
			file->squareMeshGui->roToleranceSlider->deactivate();
			file->squareMeshGui->roMaximumIterationsSlider->deactivate();
			file->squareMeshGui->roRadiusSlider->deactivate();
			if (file->squareMeshOverlay) { file->squareMesh->disableRo(); }
		}
		break;
	case 1:
		if (w->value()) {
			file->voronoiMeshGui->roToleranceSlider->activate();
			file->voronoiMeshGui->roMaximumIterationsSlider->activate();
			file->voronoiMeshGui->roRadiusSlider->activate();
			if (file->voronoiMeshOverlay) { file->voronoiMesh->enableRo(); }
		}
		else {
			file->voronoiMeshGui->roToleranceSlider->deactivate();
			file->voronoiMeshGui->roMaximumIterationsSlider->deactivate();
			file->voronoiMeshGui->roRadiusSlider->deactivate();
			if (file->voronoiMeshOverlay) { file->voronoiMesh->disableRo(); }
		}
		break;
	case 2:
		if (w->value()) {
			file->treeMeshGui->roToleranceSlider->activate();
			file->treeMeshGui->roMaximumIterationsSlider->activate();
			file->treeMeshGui->roRadiusSlider->activate();
			if (file->treeMeshOverlay) { file->treeMesh->enableRo(); }
		}
		else {
			file->treeMeshGui->roToleranceSlider->deactivate();
			file->treeMeshGui->roMaximumIterationsSlider->deactivate();
			file->treeMeshGui->roRadiusSlider->deactivate();
			if (file->treeMeshOverlay) { file->treeMesh->disableRo(); }
		}
		break;
	}
}

void generateSquareLandscape() {
	if (file->squareMesh->landscapeTriangles != NULL) {
		delete [] file->squareMesh->landscapeTriangles;
		delete [] file->squareMesh->landscapeNormals;
		delete [] file->squareMesh->landscapeColors;
	}

	const int xZones = file->squareMesh->getXCells();
	const int yZones = file->squareMesh->getYCells();
	const int zones = (xZones)*(yZones);

	file->squareMesh->landscapeTriangles = new float[4*3*3*zones];
	file->squareMesh->landscapeNormals = new float[4*3*3*zones];
	float *normals3DViewTemp = new float[4*3*3*zones];
	file->squareMesh->landscapeColors = new float[4*3*4*zones];

	const float dx = (float)file->squareMesh->getDx();
	int i = 0;

	setLandscapeSquare();

//	if (iMAP->regularMeshGui->convolveLandscapeButton->value()) { convolveLandscapeRegular(3,7); }

	const float scale = file->squareMeshGui->scaleLandscapeSlider->value();
	const float alpha = file->squareMeshGui->landscapeAlphaSlider->value();

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

	// determine plot type
	float bl,br,tr,tl,mid;
	for (int a = 0; a < xZones; a++) {
		for (int b = 0; b < yZones; b++) {

			// set colours to zero (black) in beginning
			for (int p = 0; p < 48; p++) { file->squareMesh->landscapeColors[48*i+p] = 0.0; }
			for (int q = 0; q < 36; q++) { file->squareMesh->landscapeNormals[36*i+q] = 0.0; }

			if (a == 0 && b > 0 && b < yZones-1) {
				bl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape())/2.0;
				br = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b-1)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape())/4.0;
				tr = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b+1)->getLandscape()+
					  file->squareMesh->getCell(a,b+1)->getLandscape())/4.0;
				tl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b+1)->getLandscape())/2.0;
				mid = (bl+br+tr+tl)/4.0;
			}
			else if (a == xZones-1 && b > 0 && b < yZones-1) {
				bl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b-1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b)->getLandscape())/4.0;
				br = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape())/2.0;
				tr = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b+1)->getLandscape())/2.0;
				tl = (file->squareMesh->getCell(a,b)->getLandscape()+
				      file->squareMesh->getCell(a,b+1)->getLandscape()+
				      file->squareMesh->getCell(a-1,b+1)->getLandscape()+
				      file->squareMesh->getCell(a-1,b)->getLandscape())/4.0;
				mid = (bl+br+tr+tl)/4.0;
			}
			else if (b == 0 && a > 0 && a < xZones-1) {
				bl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a-1,b)->getLandscape())/2.0;
				br = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape())/2.0;
				tr = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b+1)->getLandscape()+
					  file->squareMesh->getCell(a,b+1)->getLandscape())/4.0;
				tl = (file->squareMesh->getCell(a,b)->getLandscape()+
				      file->squareMesh->getCell(a,b+1)->getLandscape()+
				      file->squareMesh->getCell(a-1,b+1)->getLandscape()+
				      file->squareMesh->getCell(a-1,b)->getLandscape())/4.0;
				mid = (bl+br+tr+tl)/4.0;
			}
			else if (b == yZones-1 && a > 0 && a < xZones-1) {
				bl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b-1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b)->getLandscape())/4.0;
				br = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b-1)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape())/4.0;
				tr = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape())/2.0;
				tl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a-1,b)->getLandscape())/2.0;
				mid = (bl+br+tr+tl)/4.0;
			}
			else if (a == 0 && b == 0) {
				bl = (file->squareMesh->getCell(a,b)->getLandscape());
				br = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape())/2.0;
				tr = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b+1)->getLandscape()+
					  file->squareMesh->getCell(a,b+1)->getLandscape())/4.0;
				tl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b+1)->getLandscape())/2.0;
				mid = (bl+br+tr+tl)/4.0;
			}
			else if (a == 0 && b == yZones-1) {
				bl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape())/2.0;
				br = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b-1)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape())/4.0;
				tr = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape())/2.0;
				tl = (file->squareMesh->getCell(a,b)->getLandscape());
				mid = (bl+br+tr+tl)/4.0;
			}
			else if (a == xZones-1 && b == 0) {
				bl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a-1,b)->getLandscape())/2.0;
				br = (file->squareMesh->getCell(a,b)->getLandscape());
				tr = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b+1)->getLandscape())/2.0;
				tl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b+1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b+1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b)->getLandscape())/4.0;
				mid = (bl+br+tr+tl)/4.0;
			}
			else if (a == xZones-1 && b == yZones-1) {
				bl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b-1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b)->getLandscape())/4.0;
				br = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape())/2.0;
				tr = (file->squareMesh->getCell(a,b)->getLandscape());
				tl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a-1,b)->getLandscape())/2.0;
				mid = (bl+br+tr+tl)/4.0;
			}
			else {
				bl = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b-1)->getLandscape()+
					  file->squareMesh->getCell(a-1,b)->getLandscape())/4.0;
				br = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b-1)->getLandscape()+
					  file->squareMesh->getCell(a,b-1)->getLandscape())/4.0;
				tr = (file->squareMesh->getCell(a,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b)->getLandscape()+
					  file->squareMesh->getCell(a+1,b+1)->getLandscape()+
					  file->squareMesh->getCell(a,b+1)->getLandscape())/4.0;
				tl = (file->squareMesh->getCell(a,b)->getLandscape()+
				      file->squareMesh->getCell(a,b+1)->getLandscape()+
				      file->squareMesh->getCell(a-1,b+1)->getLandscape()+
				      file->squareMesh->getCell(a-1,b)->getLandscape())/4.0;
				mid = (bl+br+tr+tl)/4.0;
			}

			/* BOTTOM TRIANGLE */
			// vertex 1
			file->squareMesh->landscapeTriangles[36*i+0] = (float)(a+0.5)*dx+xMin;
			file->squareMesh->landscapeTriangles[36*i+1] = (float)(b+0.5)*dx+yMin;
			file->squareMesh->landscapeTriangles[36*i+2] = bl;
			colormap(&file->squareMesh->landscapeColors[48*i],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+2],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+3] = FMIN(1.0,20.0*alpha*bl);
			// vertex 2
			file->squareMesh->landscapeTriangles[36*i+3] = (float)FMIN(((float)a+1.5)*dx+xMin,xMax);
			file->squareMesh->landscapeTriangles[36*i+4] = (float)(b+0.5)*dx+yMin;
			file->squareMesh->landscapeTriangles[36*i+5] = br;
			colormap(&file->squareMesh->landscapeColors[48*i+4],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+5],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+7] = FMIN(1.0,20.0*alpha*br);
			// vertex 3
			file->squareMesh->landscapeTriangles[36*i+6] = (float)FMIN(((float)a+1.0)*dx+xMin,xMax);
			file->squareMesh->landscapeTriangles[36*i+7] = (float)(b+1.0)*dx+yMin;
			file->squareMesh->landscapeTriangles[36*i+8] = mid;
			colormap(&file->squareMesh->landscapeColors[48*i+8],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+8],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+11] = FMIN(1.0,20.0*alpha*mid);

			/* RIGHT TRIANGLE */
			// vertex 1
			file->squareMesh->landscapeTriangles[36*i+9] = (float)FMIN(((float)a+1.5)*dx+xMin,xMax);
			file->squareMesh->landscapeTriangles[36*i+10] = (float)(b+0.5)*dx+yMin;
			file->squareMesh->landscapeTriangles[36*i+11] = br;
			colormap(&file->squareMesh->landscapeColors[48*i+12],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+11],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+15] = FMIN(1.0,20.0*alpha*br);
			// vertex 2
			file->squareMesh->landscapeTriangles[36*i+12] = (float)FMIN(((float)a+1.5)*dx+xMin,xMax);
			file->squareMesh->landscapeTriangles[36*i+13] = (float)FMIN(((float)b+1.5)*dx+yMin,yMax);
			file->squareMesh->landscapeTriangles[36*i+14] = tr;
			colormap(&file->squareMesh->landscapeColors[48*i+16],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+14],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+19] = FMIN(1.0,20.0*alpha*tr);
			// vertex 3
			file->squareMesh->landscapeTriangles[36*i+15] = (float)FMIN(((float)a+1.0)*dx+xMin,xMax);
			file->squareMesh->landscapeTriangles[36*i+16] = (float)(b+1.0)*dx+yMin;
			file->squareMesh->landscapeTriangles[36*i+17] = mid;
			colormap(&file->squareMesh->landscapeColors[48*i+20],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+17],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+23] = FMIN(1.0,20.0*alpha*mid);

			/* TOP TRIANGLE */
			// vertex 1
			file->squareMesh->landscapeTriangles[36*i+18] = (float)FMIN(((float)a+1.5)*dx+xMin,xMax);
			file->squareMesh->landscapeTriangles[36*i+19] = (float)FMIN(((float)b+1.5)*dx+yMin,yMax);
			file->squareMesh->landscapeTriangles[36*i+20] = tr;
			colormap(&file->squareMesh->landscapeColors[48*i+24],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+20],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+27] = FMIN(1.0,20.0*alpha*tr);
			// vertex 2
			file->squareMesh->landscapeTriangles[36*i+21] = (float)(a+0.5)*dx+xMin;
			file->squareMesh->landscapeTriangles[36*i+22] = (float)FMIN(((float)b+1.5)*dx+yMin,yMax);
			file->squareMesh->landscapeTriangles[36*i+23] = tl;
			colormap(&file->squareMesh->landscapeColors[48*i+28],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+23],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+31] = FMIN(1.0,20.0*alpha*tl);
			// vertex 3
			file->squareMesh->landscapeTriangles[36*i+24] = (float)FMIN(((float)a+1.0)*dx+xMin,xMax);
			file->squareMesh->landscapeTriangles[36*i+25] = (float)(b+1.0)*dx+yMin;
			file->squareMesh->landscapeTriangles[36*i+26] = mid;
			colormap(&file->squareMesh->landscapeColors[48*i+32],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+26],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+35] = FMIN(1.0,20.0*alpha*mid);

			/* LEFT TRIANGLE */
			// vertex 1
			file->squareMesh->landscapeTriangles[36*i+27] = (float)(a+0.5)*dx+xMin;
			file->squareMesh->landscapeTriangles[36*i+28] = (float)FMIN(((float)b+1.5)*dx+yMin,yMax);
			file->squareMesh->landscapeTriangles[36*i+29] = tl;
			colormap(&file->squareMesh->landscapeColors[48*i+36],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+29],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+39] = FMIN(1.0,20.0*alpha*tl);
			// vertex 2
			file->squareMesh->landscapeTriangles[36*i+30] = (float)(a+0.5)*dx+xMin;
			file->squareMesh->landscapeTriangles[36*i+31] = (float)(b+0.5)*dx+yMin;
			file->squareMesh->landscapeTriangles[36*i+32] = bl;
			colormap(&file->squareMesh->landscapeColors[48*i+40],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+32],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+43] = FMIN(1.0,20.0*alpha*bl);
			// vertex 3
			file->squareMesh->landscapeTriangles[36*i+33] = (float)FMIN(((float)a+1.0)*dx+xMin,xMax);
			file->squareMesh->landscapeTriangles[36*i+34] = (float)(b+1.0)*dx+yMin;
			file->squareMesh->landscapeTriangles[36*i+35] = mid;
			colormap(&file->squareMesh->landscapeColors[48*i+44],file->squareMeshGui->getColormap(),file->squareMesh->landscapeTriangles[36*i+35],file->cMin,file->cMax,file->flipColormap);
			file->squareMesh->landscapeColors[48*i+47] = FMIN(1.0,20.0*alpha*mid);

			// apply scaling on z coordinate
			file->squareMesh->landscapeTriangles[36*i+2] *= scale;
			file->squareMesh->landscapeTriangles[36*i+5] *= scale;
			file->squareMesh->landscapeTriangles[36*i+8] *= scale;
			file->squareMesh->landscapeTriangles[36*i+11] *= scale;
			file->squareMesh->landscapeTriangles[36*i+14] *= scale;
			file->squareMesh->landscapeTriangles[36*i+17] *= scale;
			file->squareMesh->landscapeTriangles[36*i+20] *= scale;
			file->squareMesh->landscapeTriangles[36*i+23] *= scale;
			file->squareMesh->landscapeTriangles[36*i+26] *= scale;
			file->squareMesh->landscapeTriangles[36*i+29] *= scale;
			file->squareMesh->landscapeTriangles[36*i+32] *= scale;
			file->squareMesh->landscapeTriangles[36*i+35] *= scale;

			/* BOTTOM TRIANGLE */
			crossProduct(&normals3DViewTemp[36*i+0],
					file->squareMesh->landscapeTriangles[36*i+6]-file->squareMesh->landscapeTriangles[36*i+0],
					file->squareMesh->landscapeTriangles[36*i+7]-file->squareMesh->landscapeTriangles[36*i+1],
					file->squareMesh->landscapeTriangles[36*i+8]-file->squareMesh->landscapeTriangles[36*i+2],
					file->squareMesh->landscapeTriangles[36*i+0]-file->squareMesh->landscapeTriangles[36*i+3],
					file->squareMesh->landscapeTriangles[36*i+1]-file->squareMesh->landscapeTriangles[36*i+4],
					file->squareMesh->landscapeTriangles[36*i+2]-file->squareMesh->landscapeTriangles[36*i+5]);

			crossProduct(&normals3DViewTemp[36*i+3],
					file->squareMesh->landscapeTriangles[36*i+0]-file->squareMesh->landscapeTriangles[36*i+3],
					file->squareMesh->landscapeTriangles[36*i+1]-file->squareMesh->landscapeTriangles[36*i+4],
					file->squareMesh->landscapeTriangles[36*i+2]-file->squareMesh->landscapeTriangles[36*i+5],
					file->squareMesh->landscapeTriangles[36*i+3]-file->squareMesh->landscapeTriangles[36*i+6],
					file->squareMesh->landscapeTriangles[36*i+4]-file->squareMesh->landscapeTriangles[36*i+7],
					file->squareMesh->landscapeTriangles[36*i+5]-file->squareMesh->landscapeTriangles[36*i+8]);

			crossProduct(&normals3DViewTemp[36*i+6],
					file->squareMesh->landscapeTriangles[36*i+3]-file->squareMesh->landscapeTriangles[36*i+6],
					file->squareMesh->landscapeTriangles[36*i+4]-file->squareMesh->landscapeTriangles[36*i+7],
					file->squareMesh->landscapeTriangles[36*i+5]-file->squareMesh->landscapeTriangles[36*i+8],
					file->squareMesh->landscapeTriangles[36*i+6]-file->squareMesh->landscapeTriangles[36*i+0],
					file->squareMesh->landscapeTriangles[36*i+7]-file->squareMesh->landscapeTriangles[36*i+1],
					file->squareMesh->landscapeTriangles[36*i+8]-file->squareMesh->landscapeTriangles[36*i+2]);

			/* RIGHT TRIANGLE */
			crossProduct(&normals3DViewTemp[36*i+9],
					file->squareMesh->landscapeTriangles[36*i+15]-file->squareMesh->landscapeTriangles[36*i+9],
					file->squareMesh->landscapeTriangles[36*i+16]-file->squareMesh->landscapeTriangles[36*i+10],
					file->squareMesh->landscapeTriangles[36*i+17]-file->squareMesh->landscapeTriangles[36*i+11],
					file->squareMesh->landscapeTriangles[36*i+9]-file->squareMesh->landscapeTriangles[36*i+12],
					file->squareMesh->landscapeTriangles[36*i+10]-file->squareMesh->landscapeTriangles[36*i+13],
					file->squareMesh->landscapeTriangles[36*i+11]-file->squareMesh->landscapeTriangles[36*i+14]);

			crossProduct(&normals3DViewTemp[36*i+12],
					file->squareMesh->landscapeTriangles[36*i+9]-file->squareMesh->landscapeTriangles[36*i+12],
					file->squareMesh->landscapeTriangles[36*i+10]-file->squareMesh->landscapeTriangles[36*i+13],
					file->squareMesh->landscapeTriangles[36*i+11]-file->squareMesh->landscapeTriangles[36*i+14],
					file->squareMesh->landscapeTriangles[36*i+12]-file->squareMesh->landscapeTriangles[36*i+15],
					file->squareMesh->landscapeTriangles[36*i+13]-file->squareMesh->landscapeTriangles[36*i+16],
					file->squareMesh->landscapeTriangles[36*i+14]-file->squareMesh->landscapeTriangles[36*i+17]);

			crossProduct(&normals3DViewTemp[36*i+15],
					file->squareMesh->landscapeTriangles[36*i+12]-file->squareMesh->landscapeTriangles[36*i+15],
					file->squareMesh->landscapeTriangles[36*i+13]-file->squareMesh->landscapeTriangles[36*i+16],
					file->squareMesh->landscapeTriangles[36*i+14]-file->squareMesh->landscapeTriangles[36*i+17],
					file->squareMesh->landscapeTriangles[36*i+15]-file->squareMesh->landscapeTriangles[36*i+9],
					file->squareMesh->landscapeTriangles[36*i+16]-file->squareMesh->landscapeTriangles[36*i+10],
					file->squareMesh->landscapeTriangles[36*i+17]-file->squareMesh->landscapeTriangles[36*i+11]);

			/* TOP TRIANGLE */
			crossProduct(&normals3DViewTemp[36*i+18],
					file->squareMesh->landscapeTriangles[36*i+24]-file->squareMesh->landscapeTriangles[36*i+18],
					file->squareMesh->landscapeTriangles[36*i+25]-file->squareMesh->landscapeTriangles[36*i+19],
					file->squareMesh->landscapeTriangles[36*i+26]-file->squareMesh->landscapeTriangles[36*i+20],
					file->squareMesh->landscapeTriangles[36*i+18]-file->squareMesh->landscapeTriangles[36*i+21],
					file->squareMesh->landscapeTriangles[36*i+19]-file->squareMesh->landscapeTriangles[36*i+22],
					file->squareMesh->landscapeTriangles[36*i+20]-file->squareMesh->landscapeTriangles[36*i+23]);

			crossProduct(&normals3DViewTemp[36*i+21],
					file->squareMesh->landscapeTriangles[36*i+18]-file->squareMesh->landscapeTriangles[36*i+21],
					file->squareMesh->landscapeTriangles[36*i+19]-file->squareMesh->landscapeTriangles[36*i+22],
					file->squareMesh->landscapeTriangles[36*i+20]-file->squareMesh->landscapeTriangles[36*i+23],
					file->squareMesh->landscapeTriangles[36*i+21]-file->squareMesh->landscapeTriangles[36*i+24],
					file->squareMesh->landscapeTriangles[36*i+22]-file->squareMesh->landscapeTriangles[36*i+25],
					file->squareMesh->landscapeTriangles[36*i+23]-file->squareMesh->landscapeTriangles[36*i+26]);

			crossProduct(&normals3DViewTemp[36*i+24],
					file->squareMesh->landscapeTriangles[36*i+21]-file->squareMesh->landscapeTriangles[36*i+24],
					file->squareMesh->landscapeTriangles[36*i+22]-file->squareMesh->landscapeTriangles[36*i+25],
					file->squareMesh->landscapeTriangles[36*i+23]-file->squareMesh->landscapeTriangles[36*i+26],
					file->squareMesh->landscapeTriangles[36*i+24]-file->squareMesh->landscapeTriangles[36*i+18],
					file->squareMesh->landscapeTriangles[36*i+25]-file->squareMesh->landscapeTriangles[36*i+19],
					file->squareMesh->landscapeTriangles[36*i+26]-file->squareMesh->landscapeTriangles[36*i+20]);

			/* LEFT TRIANGLE */
			crossProduct(&normals3DViewTemp[36*i+27],
					file->squareMesh->landscapeTriangles[36*i+33]-file->squareMesh->landscapeTriangles[36*i+27],
					file->squareMesh->landscapeTriangles[36*i+34]-file->squareMesh->landscapeTriangles[36*i+28],
					file->squareMesh->landscapeTriangles[36*i+35]-file->squareMesh->landscapeTriangles[36*i+29],
					file->squareMesh->landscapeTriangles[36*i+27]-file->squareMesh->landscapeTriangles[36*i+30],
					file->squareMesh->landscapeTriangles[36*i+28]-file->squareMesh->landscapeTriangles[36*i+31],
					file->squareMesh->landscapeTriangles[36*i+29]-file->squareMesh->landscapeTriangles[36*i+32]);

			crossProduct(&normals3DViewTemp[36*i+30],
					file->squareMesh->landscapeTriangles[36*i+27]-file->squareMesh->landscapeTriangles[36*i+30],
					file->squareMesh->landscapeTriangles[36*i+28]-file->squareMesh->landscapeTriangles[36*i+31],
					file->squareMesh->landscapeTriangles[36*i+29]-file->squareMesh->landscapeTriangles[36*i+32],
					file->squareMesh->landscapeTriangles[36*i+30]-file->squareMesh->landscapeTriangles[36*i+33],
					file->squareMesh->landscapeTriangles[36*i+31]-file->squareMesh->landscapeTriangles[36*i+34],
					file->squareMesh->landscapeTriangles[36*i+32]-file->squareMesh->landscapeTriangles[36*i+35]);

			crossProduct(&normals3DViewTemp[36*i+33],
					file->squareMesh->landscapeTriangles[36*i+30]-file->squareMesh->landscapeTriangles[36*i+33],
					file->squareMesh->landscapeTriangles[36*i+31]-file->squareMesh->landscapeTriangles[36*i+34],
					file->squareMesh->landscapeTriangles[36*i+32]-file->squareMesh->landscapeTriangles[36*i+35],
					file->squareMesh->landscapeTriangles[36*i+33]-file->squareMesh->landscapeTriangles[36*i+27],
					file->squareMesh->landscapeTriangles[36*i+34]-file->squareMesh->landscapeTriangles[36*i+28],
					file->squareMesh->landscapeTriangles[36*i+35]-file->squareMesh->landscapeTriangles[36*i+29]);
			i++;
		}
	}

	// smooth normals (only one normal value per vertex Ñ to help with lighting effects)
	int div = 0;
	for (int a = 0; a < 12*zones; a++) {
		div = 0;
		const float x = file->squareMesh->landscapeTriangles[3*a+0];
		const float y = file->squareMesh->landscapeTriangles[3*a+1];
		for (int aa = 0; aa < 12*zones; aa++) {
			if ( fabs(x-file->squareMesh->landscapeTriangles[3*aa+0]) < 0.0001 && fabs(y-file->squareMesh->landscapeTriangles[3*aa+1]) < 0.0001) {
				file->squareMesh->landscapeNormals[3*a+0] += normals3DViewTemp[3*aa+0];
				file->squareMesh->landscapeNormals[3*a+1] += normals3DViewTemp[3*aa+1];
				file->squareMesh->landscapeNormals[3*a+2] += normals3DViewTemp[3*aa+2];
				div++;
			}
		}
		if (div > 0) {
			file->squareMesh->landscapeNormals[3*a+0] /= div;
			file->squareMesh->landscapeNormals[3*a+1] /= div;
			file->squareMesh->landscapeNormals[3*a+2] /= div;
		} else {
			file->squareMesh->landscapeNormals[3*a+0] = 0.0;
			file->squareMesh->landscapeNormals[3*a+1] = 0.0;
			file->squareMesh->landscapeNormals[3*a+2] = 0.0;
		}
	}

	delete [] normals3DViewTemp;

}

void squareLandscapeCallback(Fl_Check_Button*w,void*) {
	if (w->value()) {
		file->squareMeshGui->scaleLandscapeSlider->activate();
		file->squareMeshGui->landscapeAlphaSlider->activate();
		file->squareMeshGui->xLightPositionSlider->activate();
		file->squareMeshGui->yLightPositionSlider->activate();
		file->squareMeshGui->zLightPositionSlider->activate();
		file->squareMeshGui->fogButton->activate();
		file->squareMeshGui->fogStartSlider->activate();
		file->squareMeshGui->fogEndSlider->activate();
		file->squareMeshGui->landscapeAxisButton->activate();

		generateSquareLandscape();
		file->landscapePlot = true;
		file->orthoLimit = 0.6*file->maxRange;

		file->xTranslate = file->yTranslate = file->zTranslate = 0.0;
		file->xRotate = file->yRotate = file->zRotate = 0;
		iMAP->fovy = 2*atan(file->maxRange/2/15)*180/PI;

	}
	else {
		file->squareMeshGui->scaleLandscapeSlider->deactivate();
		file->squareMeshGui->landscapeAlphaSlider->deactivate();
		file->squareMeshGui->xLightPositionSlider->deactivate();
		file->squareMeshGui->yLightPositionSlider->deactivate();
		file->squareMeshGui->zLightPositionSlider->deactivate();
		file->squareMeshGui->fogButton->deactivate();
		file->squareMeshGui->fogStartSlider->deactivate();
		file->squareMeshGui->fogEndSlider->deactivate();
		file->squareMeshGui->landscapeAxisButton->deactivate();

		file->landscapePlot = false;
	}
}

void updateSquareLandscapeCallback(Fl_Widget*w,void*) {
	generateSquareLandscape();
}

void crossProduct(float *v, float v1x, float v1y, float v1z, float v2x, float v2y, float v2z) {
	v[0] = v1y*v2z - v1z*v2y;
	v[1] = v1z*v2x - v1x*v2z;
	v[2] = v1x*v2y - v1y*v2x;
	const float length = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if (length != 0.0) {
		v[0] /= length;
		v[1] /= length;
		v[2] /= length;
	}
}

void crossProduct(float*n, float*v1,float*v2) {
	n[0] = v1[1]*v2[2]-v1[2]*v2[1];
	n[1] = v1[2]*v2[0]-v1[0]*v2[2];
	n[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

void setLandscapeSquare() {
	// assign data
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				if (file->squareMeshGui->overlayDiffusionButton->value()) {
					file->squareMesh->getCell(a,b)->setLandscape((float) ( file->squareMesh->getCell(a,b)->getDiffusion() - file->squareMesh->getDiffusionMin() ) /
																  (float) ( file->squareMesh->getDiffusionMax() - file->squareMesh->getDiffusionMin() ));
				}
				else if (file->squareMeshGui->overlayPotentialButton->value()) {
					file->squareMesh->getCell(a,b)->setLandscape((float) ( file->squareMesh->getCell(a,b)->getPotential() - file->squareMesh->getPotentialMin() ) /
																  (float) ( file->squareMesh->getPotentialMax() - file->squareMesh->getPotentialMin() ));
				}
				else if (file->squareMeshGui->overlayForceMagnitudeButton->value()) {
					file->squareMesh->getCell(a,b)->setLandscape((float) ( file->squareMesh->getCell(a,b)->getForceMagnitude() - file->squareMesh->getForceMin() ) /
																  (float) ( file->squareMesh->getForceMax() - file->squareMesh->getForceMin() ));
				}
				else if (file->squareMeshGui->overlayPointNumberButton->value()) {
					file->squareMesh->getCell(a,b)->setLandscape((float) ( file->squareMesh->getCell(a,b)->getCount() - file->squareMesh->getMinCell()->getCount() ) /
																  (float) ( file->squareMesh->getMaxCell()->getCount() - file->squareMesh->getMinCell()->getCount() ));
				}
				else { file->squareMesh->getCell(a,b)->setLandscape(0.0f); }
			}
			else { file->squareMesh->getCell(a,b)->setLandscape(0.0f); }
		}
	}
}

void setLandscapeQuadTree(QuadTree *tree) {

	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {

		if (tree->active()) {

			if (file->treeMeshGui->overlayDiffusionButton->value()) {
				tree->landscape = ((float) ( tree->getDiffusion() - file->treeMesh->getDiffusionMin() ) /
											 (float) ( file->treeMesh->getDiffusionMax() - file->treeMesh->getDiffusionMin() ));
			}
			else if (file->treeMeshGui->overlayPotentialButton->value()) {
				tree->landscape = ((float) ( tree->getPotential() - file->treeMesh->getPotentialMin() ) /
											 (float) ( file->treeMesh->getPotentialMax() - file->treeMesh->getPotentialMin() ));
			}
			else if (file->treeMeshGui->overlayForceMagnitudeButton->value()) {
				tree->landscape = ((float) ( tree->getForce() - file->treeMesh->getForceMin() ) /
											 (float) ( file->treeMesh->getForceMax() - file->treeMesh->getForceMin() ));
			}
			else if (file->treeMeshGui->overlayPointNumberButton->value()) {
				tree->landscape = ((float) ( tree->getCount() - file->treeMesh->minCount ) /
											 (float) ( file->treeMesh->maxCount - file->treeMesh->minCount ));
			}
			else { tree->landscape = 0.0; }

		}
		return;
	}

	setLandscapeQuadTree(tree->nw);
	setLandscapeQuadTree(tree->ne);
	setLandscapeQuadTree(tree->sw);
	setLandscapeQuadTree(tree->se);
	return;

}

void countQuadTreeVertices(QuadTree *tree, int *count) {

	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			if (tree->nBottomNeighboursLandscape > 0) {
				*count += 3*tree->nBottomNeighboursLandscape;
			}
			else { *count += 3; }
			if (tree->nTopNeighboursLandscape > 0) {
				*count += 3*tree->nTopNeighboursLandscape;
			}
			else { *count += 3; }
			if (tree->nRightNeighboursLandscape > 0) {
				*count += 3*tree->nRightNeighboursLandscape;
			}
			else { *count += 3; }
			if (tree->nLeftNeighboursLandscape > 0) {
				*count += 3*tree->nLeftNeighboursLandscape;
			}
			else { *count += 3; }
		}
		return;
	}

	countQuadTreeVertices(tree->nw,count);
	countQuadTreeVertices(tree->ne,count);
	countQuadTreeVertices(tree->sw,count);
	countQuadTreeVertices(tree->se,count);
	return;

}

void countQuadTreeBorderVertices(QuadTree *tree, int *count) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (!tree->active()) {
			int rightActiveNeighbours = 0,leftActiveNeighbours = 0,topActiveNeighbours = 0,bottomActiveNeighbours = 0;
			for (int t = 0; t < tree->nTopNeighboursLandscape; t++) { if (tree->getTopNeighbourLandscape(t)->active()) { topActiveNeighbours++; } }
			for (int r = 0; r < tree->nRightNeighboursLandscape; r++) { if (tree->getRightNeighbourLandscape(r)->active()) { rightActiveNeighbours++; } }
			for (int b = 0; b < tree->nBottomNeighboursLandscape; b++) { if (tree->getBottomNeighbourLandscape(b)->active()) { bottomActiveNeighbours++; } }
			for (int l = 0; l < tree->nLeftNeighboursLandscape; l++) { if (tree->getLeftNeighbourLandscape(l)->active()) { leftActiveNeighbours++; } }

			// verify that neighbours are active
			for (int i = 0; i < tree->nTopNeighboursLandscape; i++) {
				if (tree->getTopNeighbourLandscape(i)->active()) { *count += 3; }
				// only check one side to avoid duplicate triangles
//				else if (tree->getTopNeighbourLandscape(i)->getPower() == tree->getPower()) { *count += 3; }
				else { *count += 3; }
			}
			for (int i = 0; i < tree->nRightNeighboursLandscape; i++) {
				if (tree->getRightNeighbourLandscape(i)->active()) { *count += 3; }
				// only check one side to avoid duplicate triangles
//				else if (tree->getRightNeighbourLandscape(i)->getPower() == tree->getPower()) { *count += 3; }
				else { *count += 3; }
			}
			for (int i = 0; i < tree->nBottomNeighboursLandscape; i++) {
				if (tree->getBottomNeighbourLandscape(i)->active()) { *count += 3; }
			}
			for (int i = 0; i < tree->nLeftNeighboursLandscape; i++) {
				if (tree->getLeftNeighbourLandscape(i)->active()) { *count += 3; }
			}
		}
		return;
	}
	countQuadTreeBorderVertices(tree->nw,count);
	countQuadTreeBorderVertices(tree->ne,count);
	countQuadTreeBorderVertices(tree->sw,count);
	countQuadTreeBorderVertices(tree->se,count);
	return;
}

void setCornerLandscapesQuadTree(double x, double y, double *value, int *num, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			// bottom left corner
			if (x-0.0001 > tree->bounds.x && x-0.0001 < tree->bounds.x+tree->bounds.w &&
				y-0.0001 > tree->bounds.y && y-0.0001 < tree->bounds.y+tree->bounds.h) {
				*value += tree->landscape;
				*num += 1;
			}
			// bottom right corner
			else if (x+0.0001 > tree->bounds.x && x+0.0001 < tree->bounds.x+tree->bounds.w &&
					 y-0.0001 > tree->bounds.y && y-0.0001 < tree->bounds.y+tree->bounds.h) {
				*value += tree->landscape;
				*num += 1;
			}
			// top right corner
			else if (x+0.0001 > tree->bounds.x && x+0.0001 < tree->bounds.x+tree->bounds.w &&
					 y+0.0001 > tree->bounds.y && y+0.0001 < tree->bounds.y+tree->bounds.h) {
				*value += tree->landscape;
				*num += 1;
			}
			// top left corner
			else if (x-0.0001 > tree->bounds.x && x-0.0001 < tree->bounds.x+tree->bounds.w &&
					 y+0.0001 > tree->bounds.y && y+0.0001 < tree->bounds.y+tree->bounds.h) {
				*value += tree->landscape;
				*num += 1;
			}
		}
		return;
	}
	setCornerLandscapesQuadTree(x,y,value,num,tree->nw);
	setCornerLandscapesQuadTree(x,y,value,num,tree->ne);
	setCornerLandscapesQuadTree(x,y,value,num,tree->sw);
	setCornerLandscapesQuadTree(x,y,value,num,tree->se);
	return;
}

void findCommonVertex(float *vertices, QuadTree *tree1, QuadTree *tree2) {
//	vertices[0] = vertices[1] = vertices[2] = 0.0;
//	if (tree1->getPower() != tree2->getPower()) { return; }
//	if (fabs(tree1->getXMax()-tree2->getXMax()) < 0.0001) {
//		if (fabs(tree1->getYMax()-tree2->getYMin()) < 0.0001) {
//			vertices[0] =
//		}
//		else {
//
//		}
//	}
}

void averageQuadTreeLandscapeNormals(QuadTree *tree, float x, float y, float *norms, int *num) {
//	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
//		if (tree->active()) {
//			// bottom left corner
//			if (x-0.00001 > tree->bounds.x && x-0.00001 < tree->bounds.x+tree->bounds.w &&
//				y-0.00001 > tree->bounds.y && y-0.00001 < tree->bounds.y+tree->bounds.h) {
//				norms[0] += tree->landscape;
//				*num += 1;
//			}
//			// bottom right corner
//			else if (x+0.00001 > tree->bounds.x && x+0.00001 < tree->bounds.x+tree->bounds.w &&
//				y-0.00001 > tree->bounds.y && y-0.00001 < tree->bounds.y+tree->bounds.h) {
//				*value += tree->landscape;
//				*num += 1;
//			}
//			// top right corner
//			else if (x+0.00001 > tree->bounds.x && x+0.00001 < tree->bounds.x+tree->bounds.w &&
//				y+0.00001 > tree->bounds.y && y+0.00001 < tree->bounds.y+tree->bounds.h) {
//				*value += tree->landscape;
//				*num += 1;
//			}
//			// top left corner
//			else if (x-0.0001 > tree->bounds.x && x-0.0001 < tree->bounds.x+tree->bounds.w &&
//				y+0.00001 > tree->bounds.y && y+0.00001 < tree->bounds.y+tree->bounds.h) {
//				*value += tree->landscape;
//				*num += 1;
//			}
//		}
//		return;
//	}
//	averageQuadTreeLandscapeNormals(x,y,value,num,tree->nw);
//	averageQuadTreeLandscapeNormals(x,y,value,num,tree->ne);
//	averageQuadTreeLandscapeNormals(x,y,value,num,tree->sw);
//	averageQuadTreeLandscapeNormals(x,y,value,num,tree->se);
//	return;
}

void defineQuadTreeLandscape(QuadTree *tree,int *i,int *j) {

	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {

			/* DEFINE VERTICES */
			int w = 0;
			int v = 0;
			QuadTree *neighbour = NULL;
			double value = 0.0;
			int num = 0;
			const float alpha = file->treeMeshGui->landscapeAlphaSlider->value();
			const float scale = file->treeMeshGui->scaleLandscapeSlider->value();
			const int cm = file->treeMeshGui->getColormap();
			// bottom neighbours
			if (tree->nBottomNeighboursLandscape > 1) {
				for (int b = 0; b < tree->nBottomNeighboursLandscape; b++) {
					neighbour = tree->getBottomNeighbourLandscape(b);
					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = scale*tree->getLandscape(); w++;
//					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
//					} else {
//						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
//						v += 4;
//					}
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
					value = 0.0; num = 0;
					setCornerLandscapesQuadTree(neighbour->getXMin(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					switch (num) {
					case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
//					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
//					} else {
//						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
//						v += 4;
//					}
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
					value = 0.0; num = 0;
					setCornerLandscapesQuadTree(neighbour->getXMax(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					switch (num) {
					case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
//					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
//					} else {
//						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
//						v += 4;
//					}
				}
			} else {
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = scale*tree->getLandscape(); w++;
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXMin(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYMin(); w++;
				value = 0.0; num = 0;
				setCornerLandscapesQuadTree(tree->getXMin(),tree->getYMin(),&value,&num,file->treeMesh->quadTree);
				switch (num) {
				case 1:
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					break;
				default:
					file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
					break;
				}
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXMax(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYMin(); w++;
				value = 0.0; num = 0;
				setCornerLandscapesQuadTree(tree->getXMax(),tree->getYMin(),&value,&num,file->treeMesh->quadTree);
				switch (num) {
				case 1:
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					break;
				default:
					file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
					break;
				}
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
			}
			// right neighbours
			if (tree->nRightNeighboursLandscape > 1) {
				for (int r = 0; r < tree->nRightNeighboursLandscape; r++) {
					neighbour = tree->getRightNeighbourLandscape(r);
					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = scale*tree->getLandscape(); w++;
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
					value = 0.0; num = 0;
					setCornerLandscapesQuadTree(neighbour->getXMin(),neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					switch (num) {
					case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
					value = 0.0; num = 0;
					setCornerLandscapesQuadTree(neighbour->getXMin(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					switch (num) {
					case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
				}
			} else {
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = scale*tree->getLandscape(); w++;
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXMax(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYMin(); w++;
				value = 0.0; num = 0;
				setCornerLandscapesQuadTree(tree->getXMax(),tree->getYMin(),&value,&num,file->treeMesh->quadTree);
				switch (num) {
				case 1:
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					break;
				default:
					file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
					break;
				}
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXMax(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYMax(); w++;
				value = 0.0; num = 0;
				setCornerLandscapesQuadTree(tree->getXMax(),tree->getYMax(),&value,&num,file->treeMesh->quadTree);
				switch (num) {
				case 1:
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					break;
				default:
					file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
					break;
				}
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
			}
			// top neighbours
			if (tree->nTopNeighboursLandscape > 1) {
				for (int t = 0; t < tree->nTopNeighboursLandscape; t++) {
					neighbour = tree->getTopNeighbourLandscape(t);
					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = scale*tree->getLandscape(); w++;
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
					value = 0.0; num = 0;
					setCornerLandscapesQuadTree(neighbour->getXMax(),neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					switch (num) {
					case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
					value = 0.0; num = 0;
					setCornerLandscapesQuadTree(neighbour->getXMin(),neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					switch (num) {
					case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
				}
			} else {
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = scale*tree->getLandscape(); w++;
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXMax(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYMax(); w++;
				value = 0.0; num = 0;
				setCornerLandscapesQuadTree(tree->getXMax(),tree->getYMax(),&value,&num,file->treeMesh->quadTree);
				switch (num) {
				case 1:
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					break;
				default:
					file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
					break;
				}
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXMin(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYMax(); w++;
				value = 0.0; num = 0;
				setCornerLandscapesQuadTree(tree->getXMin(),tree->getYMax(),&value,&num,file->treeMesh->quadTree);
				switch (num) {
				case 1:
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					break;
				default:
					file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
					break;
				}
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
			}
			// left neighbours
			if (tree->nLeftNeighboursLandscape > 1) {
				for (int l = 0; l < tree->nLeftNeighboursLandscape; l++) {
					neighbour = tree->getLeftNeighbourLandscape(l);
					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = scale*tree->getLandscape(); w++;
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
					value = 0.0; num = 0;
					setCornerLandscapesQuadTree(neighbour->getXMax(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					switch (num) {
					case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
					value = 0.0; num = 0;
					setCornerLandscapesQuadTree(neighbour->getXMax(),neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					switch (num) {
					case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
				}
			} else {
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = scale*tree->getLandscape(); w++;
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXMin(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYMax(); w++;
				value = 0.0; num = 0;
				setCornerLandscapesQuadTree(tree->getXMin(),tree->getYMax(),&value,&num,file->treeMesh->quadTree);
				switch (num) {
				case 1:
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					break;
				default:
					file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
					break;
				}
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
				file->treeMesh->landscapeTriangles[*i+w] = tree->getXMin(); w++;
				file->treeMesh->landscapeTriangles[*i+w] = tree->getYMin(); w++;
				value = 0.0; num = 0;
				setCornerLandscapesQuadTree(tree->getXMin(),tree->getYMin(),&value,&num,file->treeMesh->quadTree);
				switch (num) {
				case 1:
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					break;
				default:
					file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
					break;
				}
				if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+v+3] = alpha;
					v += 4;
				} else {
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;
				}
			}

			// average the mid vector (to reduce pyramid effects)
			float mid = 0.0;
			int nmid = 0;
			for (int ww = 0; ww < w/9; ww++) {
				mid += file->treeMesh->landscapeTriangles[*i+9*ww+5] + file->treeMesh->landscapeTriangles[*i+9*ww+8];
				nmid += 2;
			}
			mid /= (float)nmid;
			for (int ww = 0; ww < w/9; ww++) {
				file->treeMesh->landscapeTriangles[*i+9*ww+2] = mid;
				if (mid > 0.0) {
					colormap(&file->treeMesh->landscapeColors[*j+12*ww],cm,mid/scale,file->cMin,file->cMax,file->flipColormap);
					file->treeMesh->landscapeColors[*j+12*ww+3] = alpha;
				} else {
					file->treeMesh->landscapeColors[*j+12*ww] = file->treeMesh->landscapeColors[*j+12*ww+1] = file->treeMesh->landscapeColors[*j+12*ww+2] = file->treeMesh->landscapeColors[*j+12*ww+3] = 0.0;
				}
			}

			// define normals
			float v01x,v01y,v01z;
			float v12x,v12y,v12z;
			float v20x,v20y,v20z;
			for (int ww = 0; ww < w/9; ww++ ) {
				v01x = file->treeMesh->landscapeTriangles[*i+9*ww+0]-file->treeMesh->landscapeTriangles[*i+9*ww+3];
				v01y = file->treeMesh->landscapeTriangles[*i+9*ww+1]-file->treeMesh->landscapeTriangles[*i+9*ww+4];
				v01z = file->treeMesh->landscapeTriangles[*i+9*ww+2]-file->treeMesh->landscapeTriangles[*i+9*ww+5];
				v12x = file->treeMesh->landscapeTriangles[*i+9*ww+3]-file->treeMesh->landscapeTriangles[*i+9*ww+6];
				v12y = file->treeMesh->landscapeTriangles[*i+9*ww+4]-file->treeMesh->landscapeTriangles[*i+9*ww+7];
				v12z = file->treeMesh->landscapeTriangles[*i+9*ww+5]-file->treeMesh->landscapeTriangles[*i+9*ww+8];
				v20x = file->treeMesh->landscapeTriangles[*i+9*ww+6]-file->treeMesh->landscapeTriangles[*i+9*ww+0];
				v20y = file->treeMesh->landscapeTriangles[*i+9*ww+7]-file->treeMesh->landscapeTriangles[*i+9*ww+1];
				v20z = file->treeMesh->landscapeTriangles[*i+9*ww+8]-file->treeMesh->landscapeTriangles[*i+9*ww+2];
				// 0 node
				crossProduct(&file->treeMesh->landscapeNormals[*i+9*ww+0],v20x,v20y,v20z,v01x,v01y,v01z);
				// 1 node
				crossProduct(&file->treeMesh->landscapeNormals[*i+9*ww+3],v01x,v01y,v01z,v12x,v12y,v12z);
				// 2 node
				crossProduct(&file->treeMesh->landscapeNormals[*i+9*ww+6],v12x,v12y,v12z,v20x,v20y,v20z);
			}

			*i += w;
			*j += v;

		}
		return;
	}

	defineQuadTreeLandscape(tree->nw,i,j);
	defineQuadTreeLandscape(tree->ne,i,j);
	defineQuadTreeLandscape(tree->sw,i,j);
	defineQuadTreeLandscape(tree->se,i,j);
	return;

}

void defineQuadTreeLandscapeBorders(QuadTree *tree,int *i,int *j) {

	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (!tree->active()) {
			int w = 0;
			int v = 0;
			double value = 0.0;
			int num = 0;
			QuadTree *neighbour = NULL;
			const float alpha = file->treeMeshGui->landscapeAlphaSlider->value();
			const float scale = file->treeMeshGui->scaleLandscapeSlider->value();
			const int cm = file->treeMeshGui->getColormap();
			bool right = false;
			bool top = false;
			int rightActiveNeighbours = 0,leftActiveNeighbours = 0,topActiveNeighbours = 0,bottomActiveNeighbours = 0;
			for (int t = 0; t < tree->nTopNeighboursLandscape; t++) { if (tree->getTopNeighbourLandscape(t)->active()) { topActiveNeighbours++; } }
			for (int r = 0; r < tree->nRightNeighboursLandscape; r++) { if (tree->getRightNeighbourLandscape(r)->active()) { rightActiveNeighbours++; } }
			for (int b = 0; b < tree->nBottomNeighboursLandscape; b++) { if (tree->getBottomNeighbourLandscape(b)->active()) { bottomActiveNeighbours++; } }
			for (int l = 0; l < tree->nLeftNeighboursLandscape; l++) { if (tree->getLeftNeighbourLandscape(l)->active()) { leftActiveNeighbours++; } }

			for (int t = 0; t < tree->nTopNeighboursLandscape; t++) {
				if (tree->getTopNeighbourLandscape(t)->active()) {
					neighbour = tree->getTopNeighbourLandscape(t);

					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					// colors
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;

					if (tree->getPower() > neighbour->getPower()) {
						file->treeMesh->landscapeTriangles[*i+w] = tree->getXMin(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(tree->getXMin(),neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					} else {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMin(),neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					}
					switch (num) {
					case 0: case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					// colors
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}

					if (tree->getPower() > neighbour->getPower()) {
						file->treeMesh->landscapeTriangles[*i+w] = tree->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree((float)tree->getXMax(),(float)neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					} else {
						value = 0.0; num = 0;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
						setCornerLandscapesQuadTree((float)neighbour->getXMax(),(float)neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					}
					switch (num) {
					case 0: case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					// colors
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}

				} else {
					neighbour = tree->getTopNeighbourLandscape(t);

					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					// colors
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;

					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					// colors
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;

					if (rightActiveNeighbours > 0 && right == false) {
						right = true;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(tree->getXMax(),tree->getYMax(),&value,&num,file->treeMesh->quadTree);
						if (num > 1) {
							file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
							// colors
							colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
							file->treeMesh->landscapeColors[*j+v+3] = alpha;
							v += 4;
						}
						else {
							file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
							// colors
							file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
							v += 4;
						}
					}
					else {
						file->treeMesh->landscapeTriangles[*i+w] = tree->getXMin(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(tree->getXMin(),tree->getYMax(),&value,&num,file->treeMesh->quadTree);
						if (num > 1) {
							file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
							// colors
							colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
							file->treeMesh->landscapeColors[*j+v+3] = alpha;
							v += 4;
						}
						else {
							file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
							// colors
							file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
							v += 4;
						}
					}
				}
			}
			for (int r = 0; r < tree->nRightNeighboursLandscape; r++) {
				if (tree->getRightNeighbourLandscape(r)->active()) {
					neighbour = tree->getRightNeighbourLandscape(r);

					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					// colors
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;

					if (tree->getPower() > neighbour->getPower()) {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMin(),tree->getYMax(),&value,&num,file->treeMesh->quadTree);
					} else {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMin(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					}
					switch (num) {
					case 0: case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					// colors
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}

					if (tree->getPower() > neighbour->getPower()) {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getYMin(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMin(),tree->getYMin(),&value,&num,file->treeMesh->quadTree);
					} else {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMin(),neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					}
					switch (num) {
					case 0: case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					// colors
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}

				} else {
					neighbour = tree->getRightNeighbourLandscape(r);

					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					// colors
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;

					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					// colors
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;

					if (topActiveNeighbours > 0 && top == false) {
						top = true;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(tree->getXMax(),tree->getYMax(),&value,&num,file->treeMesh->quadTree);
						if (num > 1) {
							file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
							// colors
							colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
							file->treeMesh->landscapeColors[*j+v+3] = alpha;
							v += 4;
						}
						else {
							file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
							// colors
							file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
							v += 4;
						}
					}
					else {
						file->treeMesh->landscapeTriangles[*i+w] = tree->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getYMin(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(tree->getXMax(),tree->getYMin(),&value,&num,file->treeMesh->quadTree);
						if (num > 1) {
							file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
							// colors
							colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
							file->treeMesh->landscapeColors[*j+v+3] = alpha;
							v += 4;
						}
						else {
							file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
							// colors
							file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
							v += 4;
						}
					}
				}
			}
			for (int b = 0; b < tree->nBottomNeighboursLandscape; b++) {
				if (tree->getBottomNeighbourLandscape(b)->active()) {
					neighbour = tree->getBottomNeighbourLandscape(b);

					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					// colors
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;

					if (tree->getPower() > neighbour->getPower()) {
						file->treeMesh->landscapeTriangles[*i+w] = tree->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(tree->getXMax(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					} else {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMax(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					}
					switch (num) {
					case 0: case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					// colors
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}

					if (tree->getPower() > neighbour->getPower()) {
						file->treeMesh->landscapeTriangles[*i+w] = tree->getXMin(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(tree->getXMin(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					} else {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMin(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMin(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					}
					switch (num) {
					case 0: case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					// colors
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
				}
			}
			for (int l = 0; l < tree->nLeftNeighboursLandscape; l++) {
				if (tree->getLeftNeighbourLandscape(l)->active()) {
					neighbour = tree->getLeftNeighbourLandscape(l);

					file->treeMesh->landscapeTriangles[*i+w] = tree->getXCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = tree->getYCentroid(); w++;
					file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
					// colors
					file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
					v += 4;

					if (tree->getPower() > neighbour->getPower()) {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMax(),tree->getYMax(),&value,&num,file->treeMesh->quadTree);
					} else {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMax(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMax(),neighbour->getYMax(),&value,&num,file->treeMesh->quadTree);
					}
					switch (num) {
					case 0: case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					// colors
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}

					if (tree->getPower() > neighbour->getPower()) {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = tree->getYMin(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMax(),tree->getYMin(),&value,&num,file->treeMesh->quadTree);
					} else {
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getXMax(); w++;
						file->treeMesh->landscapeTriangles[*i+w] = neighbour->getYMin(); w++;
						value = 0.0; num = 0;
						setCornerLandscapesQuadTree(neighbour->getXMax(),neighbour->getYMin(),&value,&num,file->treeMesh->quadTree);
					}
					switch (num) {
					case 0: case 1:
						file->treeMesh->landscapeTriangles[*i+w] = 0.0; w++;
						break;
					default:
						file->treeMesh->landscapeTriangles[*i+w] = scale*value/(float)num; w++;
						break;
					}
					// colors
					if (file->treeMesh->landscapeTriangles[*i+w-1] > 0.0) {
						colormap(&file->treeMesh->landscapeColors[*j+v],cm,file->treeMesh->landscapeTriangles[*i+w-1]/scale,file->cMin,file->cMax,file->flipColormap);
						file->treeMesh->landscapeColors[*j+v+3] = alpha;
						v += 4;
					} else {
						file->treeMesh->landscapeColors[*j+v] = file->treeMesh->landscapeColors[*j+v+1] = file->treeMesh->landscapeColors[*j+v+2] = file->treeMesh->landscapeColors[*j+v+3] = 0.0;
						v += 4;
					}
				}
			}

			// define normals
			float v01x,v01y,v01z;
			float v12x,v12y,v12z;
			float v20x,v20y,v20z;
			for (int ww = 0; ww < w/9; ww++ ) {
				v01x = file->treeMesh->landscapeTriangles[*i+9*ww+0]-file->treeMesh->landscapeTriangles[*i+9*ww+3];
				v01y = file->treeMesh->landscapeTriangles[*i+9*ww+1]-file->treeMesh->landscapeTriangles[*i+9*ww+4];
				v01z = file->treeMesh->landscapeTriangles[*i+9*ww+2]-file->treeMesh->landscapeTriangles[*i+9*ww+5];
				v12x = file->treeMesh->landscapeTriangles[*i+9*ww+3]-file->treeMesh->landscapeTriangles[*i+9*ww+6];
				v12y = file->treeMesh->landscapeTriangles[*i+9*ww+4]-file->treeMesh->landscapeTriangles[*i+9*ww+7];
				v12z = file->treeMesh->landscapeTriangles[*i+9*ww+5]-file->treeMesh->landscapeTriangles[*i+9*ww+8];
				v20x = file->treeMesh->landscapeTriangles[*i+9*ww+6]-file->treeMesh->landscapeTriangles[*i+9*ww+0];
				v20y = file->treeMesh->landscapeTriangles[*i+9*ww+7]-file->treeMesh->landscapeTriangles[*i+9*ww+1];
				v20z = file->treeMesh->landscapeTriangles[*i+9*ww+8]-file->treeMesh->landscapeTriangles[*i+9*ww+2];
				// 0 node
				crossProduct(&file->treeMesh->landscapeNormals[*i+9*ww+0],v20x,v20y,v20z,v01x,v01y,v01z);
				// 1 node
				crossProduct(&file->treeMesh->landscapeNormals[*i+9*ww+3],v01x,v01y,v01z,v12x,v12y,v12z);
				// 2 node
				crossProduct(&file->treeMesh->landscapeNormals[*i+9*ww+6],v12x,v12y,v12z,v20x,v20y,v20z);
			}

			*i += w;
			*j += v;
		}
		return;
	}

	defineQuadTreeLandscapeBorders(tree->nw,i,j);
	defineQuadTreeLandscapeBorders(tree->ne,i,j);
	defineQuadTreeLandscapeBorders(tree->sw,i,j);
	defineQuadTreeLandscapeBorders(tree->se,i,j);
	return;

}

void generateQuadTreeLandscape() {
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

	// count base vertices
	int nVertices = 0;
	countQuadTreeVertices(file->treeMesh->quadTree,&nVertices);
	file->treeMesh->landscapeVertices = nVertices;
	// count border vertices
	int nBorderVertices = 0;
	countQuadTreeBorderVertices(file->treeMesh->quadTree,&nBorderVertices);
	file->treeMesh->landscapeVertices += nBorderVertices;

	setLandscapeQuadTree(file->treeMesh->quadTree);

	if (file->treeMesh->landscapeTriangles != NULL) {
		delete [] file->treeMesh->landscapeTriangles;
		delete [] file->treeMesh->landscapeNormals;
		delete [] file->treeMesh->landscapeColors;
	}
	file->treeMesh->landscapeTriangles = new float[3*file->treeMesh->landscapeVertices];
	file->treeMesh->landscapeNormals = new float[3*file->treeMesh->landscapeVertices];
	file->treeMesh->landscapeColors = new float[4*file->treeMesh->landscapeVertices];

	int i = 0;
	int j = 0;

	defineQuadTreeLandscape(file->treeMesh->quadTree,&i,&j);
	defineQuadTreeLandscapeBorders(file->treeMesh->quadTree,&i,&j);

	float norms[3];
	int num = 0;
	float x,y;

	for (int g = 0; g < file->treeMesh->landscapeVertices; g++) {
		num = 0;
		norms[0] = norms[1] = norms[2] = 0.0;
		x = file->treeMesh->landscapeTriangles[3*g+0];
		y = file->treeMesh->landscapeTriangles[3*g+1];
		for (int gg = 0; gg < file->treeMesh->landscapeVertices; gg++) {
			if (fabs(x-file->treeMesh->landscapeTriangles[3*gg+0]) < 0.000001 && fabs(y-file->treeMesh->landscapeTriangles[3*gg+1]) < 0.000001) {
				norms[0] += file->treeMesh->landscapeNormals[3*gg+0];
				norms[1] += file->treeMesh->landscapeNormals[3*gg+1];
				norms[2] += file->treeMesh->landscapeNormals[3*gg+2];
				num++;
			}
		}
		if (num > 0) {
			file->treeMesh->landscapeNormals[3*g+0] = norms[0]/(float)num;
			file->treeMesh->landscapeNormals[3*g+1] = norms[1]/(float)num;
			file->treeMesh->landscapeNormals[3*g+2] = norms[2]/(float)num;
		}
	}

}

void quadTreeLandscapeCallback(Fl_Check_Button*w,void*) {
	if (w->value()) {
		file->treeMeshGui->scaleLandscapeSlider->activate();
		file->treeMeshGui->landscapeAlphaSlider->activate();
		file->treeMeshGui->xLightPositionSlider->activate();
		file->treeMeshGui->yLightPositionSlider->activate();
		file->treeMeshGui->zLightPositionSlider->activate();
		file->treeMeshGui->fogButton->activate();
		file->treeMeshGui->fogStartSlider->activate();
		file->treeMeshGui->fogEndSlider->activate();
		file->treeMeshGui->landscapeAxisButton->activate();

		generateQuadTreeLandscape();
		file->orthoLimit = 0.6*file->maxRange;
		file->landscapePlot = true;
		file->xTranslate = file->yTranslate = file->zTranslate = 0.0;
		file->xRotate = file->yRotate = file->zRotate = 0;
		iMAP->fovy = 2*atan(file->maxRange/2/15)*180/PI;
	}
	else {
		file->treeMeshGui->scaleLandscapeSlider->deactivate();
		file->treeMeshGui->landscapeAlphaSlider->deactivate();
		file->treeMeshGui->xLightPositionSlider->deactivate();
		file->treeMeshGui->yLightPositionSlider->deactivate();
		file->treeMeshGui->zLightPositionSlider->deactivate();
		file->treeMeshGui->fogButton->deactivate();
		file->treeMeshGui->fogStartSlider->deactivate();
		file->treeMeshGui->fogEndSlider->deactivate();
		file->treeMeshGui->landscapeAxisButton->deactivate();

		file->landscapePlot = false;
	}
}

void updateQuadTreeLandscapeCallback(Fl_Widget*w,void*) {
//	int i = 0;
//	int j = 0;
//	defineQuadTreeLandscape(file->treeMesh->quadTree,&i,&j);
	generateQuadTreeLandscape();
}

void alphaQuadTreeLandscapeCallback(Fl_Widget*w,void*v) {
	const float alpha = file->treeMeshGui->landscapeAlphaSlider->value();
	for (int a = 0; a < file->treeMesh->landscapeVertices; a++) {
		// change alpha on landscape
		if (file->treeMesh->landscapeColors[4*a+3] > 0.0) {
			file->treeMesh->landscapeColors[4*a+3] = alpha;
		}
	}
}

void alphaSquareLandscapeCallback(Fl_Widget*w,void*v) {
	int i = 0;
	const float alpha = file->squareMeshGui->landscapeAlphaSlider->value();
	const int xZones = file->squareMesh->getXCells();
	const int yZones = file->squareMesh->getYCells();

	for (int a = 0; a < xZones; a++) {
		for (int b = 0; b < yZones; b++) {
			// apply scaling on z coordinate
			if (file->squareMesh->landscapeTriangles[36*i+2] > 0.0) {
				file->squareMesh->landscapeColors[48*i+3] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+2]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+5] > 0.0) {
				file->squareMesh->landscapeColors[48*i+7] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+5]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+8] > 0.0) {
				file->squareMesh->landscapeColors[48*i+11] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+8]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+11] > 0.0) {
				file->squareMesh->landscapeColors[48*i+15] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+11]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+14] > 0.0) {
				file->squareMesh->landscapeColors[48*i+19] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+14]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+17] > 0.0) {
				file->squareMesh->landscapeColors[48*i+23] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+17]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+20] > 0.0) {
				file->squareMesh->landscapeColors[48*i+27] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+20]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+23] > 0.0) {
				file->squareMesh->landscapeColors[48*i+31] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+23]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+26] > 0.0) {
				file->squareMesh->landscapeColors[48*i+35] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+26]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+29] > 0.0) {
				file->squareMesh->landscapeColors[48*i+39] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+29]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+32] > 0.0) {
				file->squareMesh->landscapeColors[48*i+43] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+32]);
			}
			if (file->squareMesh->landscapeTriangles[36*i+35] > 0.0) {
				file->squareMesh->landscapeColors[48*i+47] = FMIN(1.0,20.0*alpha*file->squareMesh->landscapeTriangles[36*i+35]);
			}
			i++;
		}
	}
}

void saveTreeMeshCallback(Fl_Button*w,int*v) {
	file->treeMesh->exportMesh();
}

void smoothingPriorButtonCallback(Fl_Check_Button*w,int*v) {

	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {
		case 0:
			if (w->value()) {
				file->squareMeshGui->muSlider->activate();
				if (file->optimizationMode == 3) {
					file->squareMeshGui->lambdaSlider->activate();
				} else { file->squareMeshGui->lambdaSlider->deactivate(); }
				if (file->squareMeshApplied()) {
					file->squareMesh->enableSmoothingPrior();
				}
				file->squareMeshGui->roButton->activate();
				iMAP->smoothingPriorActive = true;
			} else {
				file->squareMeshGui->muSlider->deactivate();
				file->squareMeshGui->lambdaSlider->deactivate();
				if (file->squareMeshApplied()) {
					file->squareMesh->disableSmoothingPrior();
				}
				file->squareMeshGui->roButton->deactivate();
				file->squareMeshGui->roButton->value(0);
				iMAP->smoothingPriorActive = false;
			}
			break;
		case 1:
			if (w->value()) {
				file->voronoiMeshGui->muSlider->activate();
				if (file->optimizationMode == 3) {
					file->voronoiMeshGui->lambdaSlider->activate();
				} else { file->voronoiMeshGui->lambdaSlider->deactivate(); }
				if (file->voronoiMeshApplied()) {
					file->voronoiMesh->enableSmoothingPrior();
				}
				file->voronoiMeshGui->roButton->activate();
				iMAP->smoothingPriorActive = true;
			} else {
				file->voronoiMeshGui->muSlider->deactivate();
				file->voronoiMeshGui->lambdaSlider->deactivate();
				if (file->voronoiMeshApplied()) {
					file->voronoiMesh->disableSmoothingPrior();
				}
				file->voronoiMeshGui->roButton->deactivate();
				file->voronoiMeshGui->roButton->value(0);
				iMAP->smoothingPriorActive = false;
			}
			break;
		case 2:
			if (w->value()) {
				file->treeMeshGui->muSlider->activate();
				if (file->optimizationMode == 3) {
					file->treeMeshGui->lambdaSlider->activate();
				} else { file->treeMeshGui->lambdaSlider->deactivate(); }
				if (file->treeMeshApplied()) {
					file->treeMesh->enableSmoothingPrior();
				}
				file->treeMeshGui->roButton->activate();
				iMAP->smoothingPriorActive= true;
			} else {
				file->treeMeshGui->muSlider->deactivate();
				file->treeMeshGui->lambdaSlider->deactivate();
				if (file->treeMeshApplied()) {
					file->treeMesh->disableSmoothingPrior();
				}
				file->treeMeshGui->roButton->deactivate();
				file->treeMeshGui->roButton->value(0);
				iMAP->smoothingPriorActive = false;
			}
			break;
	}
	Fl::redraw();
}

void jeffreysPriorButtonCallback(Fl_Check_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {
		case 0:
			if (file->squareMeshApply) {
				if (w->value()) { file->squareMesh->enableJeffreysPrior(); }
				else { file->squareMesh->disableJeffreysPrior(); }
			}
			break;
		case 1:
			if (file->voronoiMeshOverlay) {
				if (w->value()) { file->voronoiMesh->enableJeffreysPrior();}
				else { file->voronoiMesh->disableJeffreysPrior(); }
			}
			break;
		case 2:
			if (file->treeMeshOverlay) {
				if (w->value()) { file->treeMesh->enableJeffreysPrior(); }
				else { file->treeMesh->disableJeffreysPrior(); }
			}
			break;
	}
	Fl::redraw();
}

void squareMeshWindowCallback(Fl_Widget *w,void*v) {
	file->squareMeshOverlay = false;
	file->squareMeshGui->resetButton->do_callback();
	file->squareMeshGui->hide();
	file->squareMeshGui->saveButton->deactivate();
	file->inferred = false;
	file->landscapePlot = false;
	iMAP->customSelectionInferenceGui->activate();
	if (file->squareMeshOverlay) {
		if ( file->squareMesh->randomizedOptimizationGui != NULL) {
			file->squareMesh->randomizedOptimizationGui->hide();
	//		delete file->squareMesh->randomizedOptimizationGui;
		}
	}
	delete file->squareMeshGui;
	file->squareMeshGui = NULL;
	iMAP->overlayAdjustment = true;
}

void minPointsSliderCallback(Slider*w,void*v) {
	if (file->squareMeshApplied()) {
		for (int a = 0; a < file->squareMesh->getXCells(); a++) {
			for (int b = 0; b < file->squareMesh->getYCells(); b++) {
				if (file->squareMesh->getCellCount(a,b) < (int)w->value()) {
					file->squareMesh->deactivateCell(a,b);
				}
				else {
					file->squareMesh->activateCell(a,b);
				}
				file->squareMeshGui->updateVariables(iMAP->activeZones);
				iMAP->overlayAdjustment = true;
			}
		}
	}
	if (file->voronoiMeshApplied()) {
		iMAP->activeZones = 0;
		for (int g = 0; g < file->voronoiMesh->getNumberOfClusters(); g++) {
			if (file->voronoiMesh->getCellCount(g) < (int)w->value()) {
				file->voronoiMesh->getCell(g)->deactivate();
			}
			else {
				file->voronoiMesh->getCell(g)->activate();
				iMAP->activeZones++;
			}
		}
		file->voronoiMeshGui->updateVariables(iMAP->activeZones);
	}
	if (file->treeMeshApplied()) {

	}
	Fl::redraw();
}

void optimizationChoiceCallback(Fl_Choice*w,int*v) {

	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	const int choice = w->value();
	file->optimizationMode = choice;

	switch(meshType) {
		case 0 : // square mesh
			switch(choice) {
				case 0: // D
					file->squareMeshGui->roButton->deactivate();
					file->squareMeshGui->roButton->value(0);
					file->squareMeshGui->roButton->do_callback();
					file->squareMeshGui->betaSlider->deactivate();
					file->squareMeshGui->polynomialOrderSlider->deactivate();
					file->squareMeshGui->priorGroup->activate();
					file->squareMeshGui->dPosterioriButton->activate();
					file->squareMeshGui->fPosteriorButton->deactivate();
					file->squareMeshGui->vPosteriorButton->deactivate();
					file->squareMeshGui->potentialReferenceButton->deactivate();
					file->squareMeshGui->overlayForceArrowsButton->copy_label("Arrows");
					file->squareMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
					file->squareMeshGui->smoothingPriorButton->activate();
					file->squareMeshGui->lambdaSlider->deactivate();
					file->squareMeshGui->muSlider->deactivate();
					file->squareMeshGui->smoothingPriorButton->value(0);
					file->squareMeshGui->smoothingPriorButton->do_callback();
					file->squareMeshGui->pauseButton->label("@||");
					break;
				case 1: // DF
					file->squareMeshGui->roButton->deactivate();
					file->squareMeshGui->roButton->value(0);
					file->squareMeshGui->roButton->do_callback();
					file->squareMeshGui->betaSlider->activate();
					file->squareMeshGui->polynomialOrderSlider->deactivate();
					file->squareMeshGui->priorGroup->activate();
					file->squareMeshGui->dPosterioriButton->activate();
					file->squareMeshGui->fPosteriorButton->activate();
					file->squareMeshGui->vPosteriorButton->deactivate();
					file->squareMeshGui->potentialReferenceButton->deactivate();
					file->squareMeshGui->overlayForceArrowsButton->copy_label("Arrows");
					file->squareMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
					file->squareMeshGui->smoothingPriorButton->activate();
					file->squareMeshGui->lambdaSlider->deactivate();
					file->squareMeshGui->muSlider->deactivate();
					file->squareMeshGui->smoothingPriorButton->value(0);
					file->squareMeshGui->smoothingPriorButton->do_callback();
					file->squareMeshGui->pauseButton->label("@||");
					break;
				case 2: // DDr
					file->squareMeshGui->roButton->deactivate();
					file->squareMeshGui->roButton->value(0);
					file->squareMeshGui->roButton->do_callback();
					file->squareMeshGui->betaSlider->activate();
					file->squareMeshGui->polynomialOrderSlider->deactivate();
					file->squareMeshGui->priorGroup->activate();
					file->squareMeshGui->dPosterioriButton->activate();
					file->squareMeshGui->fPosteriorButton->activate();
					file->squareMeshGui->vPosteriorButton->deactivate();
					file->squareMeshGui->potentialReferenceButton->deactivate();
					file->squareMeshGui->overlayForceArrowsButton->copy_label("Arrows");
					file->squareMeshGui->overlayForceMagnitudeButton->copy_label("Drift Magnitude");
					file->squareMeshGui->smoothingPriorButton->activate();
					file->squareMeshGui->lambdaSlider->deactivate();
					file->squareMeshGui->muSlider->deactivate();
					file->squareMeshGui->smoothingPriorButton->value(0);
					file->squareMeshGui->smoothingPriorButton->do_callback();
					file->squareMeshGui->pauseButton->label("@||");
					break;
				case 3: // DV
					file->squareMeshGui->smoothingPriorButton->activate();
					file->squareMeshGui->smoothingPriorButton->value(0);
					file->squareMeshGui->smoothingPriorButton->do_callback();
					file->squareMeshGui->roButton->value(0);
					file->squareMeshGui->roButton->activate();
					file->squareMeshGui->betaSlider->deactivate();
					file->squareMeshGui->polynomialOrderSlider->deactivate();
					file->squareMeshGui->priorGroup->activate();
					file->squareMeshGui->lambdaSlider->deactivate();
					file->squareMeshGui->muSlider->deactivate();
					file->squareMeshGui->dPosterioriButton->activate();
					file->squareMeshGui->fPosteriorButton->deactivate();
					file->squareMeshGui->vPosteriorButton->activate();
					file->squareMeshGui->potentialReferenceButton->activate();
					file->squareMeshGui->overlayForceArrowsButton->copy_label("Arrows");
					file->squareMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
					break;
				case 4: // Poly
					file->squareMeshGui->roButton->deactivate();
					file->squareMeshGui->roButton->value(0);
					file->squareMeshGui->roButton->do_callback();
					file->squareMeshGui->betaSlider->deactivate();
					file->squareMeshGui->polynomialOrderSlider->activate();
					file->squareMeshGui->priorGroup->activate();
					file->squareMeshGui->dPosterioriButton->activate();
					file->squareMeshGui->fPosteriorButton->activate();
					file->squareMeshGui->vPosteriorButton->deactivate();
					file->squareMeshGui->potentialReferenceButton->deactivate();
					file->squareMeshGui->overlayForceArrowsButton->copy_label("Arrows");
					file->squareMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
					file->squareMeshGui->smoothingPriorButton->deactivate();
					file->squareMeshGui->lambdaSlider->deactivate();
					file->squareMeshGui->muSlider->deactivate();
					file->squareMeshGui->smoothingPriorButton->value(0);
					file->squareMeshGui->smoothingPriorButton->do_callback();
					file->squareMeshGui->pauseButton->label("@||");
					break;
				default:
					break;
			}
			break;
		case 1: // voronoi mesh
			switch(choice) {
			case 0: // D
				file->voronoiMeshGui->roButton->deactivate();
				file->voronoiMeshGui->roButton->value(0);
				file->voronoiMeshGui->roButton->do_callback();
				file->voronoiMeshGui->betaSlider->deactivate();
				file->voronoiMeshGui->polynomialOrderSlider->deactivate();
				file->voronoiMeshGui->priorGroup->activate();
				file->voronoiMeshGui->dPosterioriButton->activate();
				file->voronoiMeshGui->fPosteriorButton->deactivate();
				file->voronoiMeshGui->vPosteriorButton->deactivate();
				file->voronoiMeshGui->potentialReferenceButton->deactivate();
				file->voronoiMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->voronoiMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
				file->voronoiMeshGui->smoothingPriorButton->activate();
				file->voronoiMeshGui->lambdaSlider->deactivate();
				file->voronoiMeshGui->muSlider->deactivate();
				file->voronoiMeshGui->smoothingPriorButton->value(0);
				file->voronoiMeshGui->smoothingPriorButton->do_callback();
				file->voronoiMeshGui->pauseButton->label("@||");
				break;
			case 1: // DF
				file->voronoiMeshGui->roButton->deactivate();
				file->voronoiMeshGui->roButton->value(0);
				file->voronoiMeshGui->roButton->do_callback();
				file->voronoiMeshGui->betaSlider->activate();
				file->voronoiMeshGui->polynomialOrderSlider->deactivate();
				file->voronoiMeshGui->priorGroup->activate();
				file->voronoiMeshGui->dPosterioriButton->activate();
				file->voronoiMeshGui->fPosteriorButton->activate();
				file->voronoiMeshGui->vPosteriorButton->deactivate();
				file->voronoiMeshGui->potentialReferenceButton->deactivate();
				file->voronoiMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->voronoiMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
				file->voronoiMeshGui->smoothingPriorButton->activate();
				file->voronoiMeshGui->lambdaSlider->deactivate();
				file->voronoiMeshGui->muSlider->deactivate();
				file->voronoiMeshGui->smoothingPriorButton->value(0);
				file->voronoiMeshGui->smoothingPriorButton->do_callback();
				file->voronoiMeshGui->pauseButton->label("@||");
				break;
			case 2: // DDr
				file->voronoiMeshGui->roButton->deactivate();
				file->voronoiMeshGui->roButton->value(0);
				file->voronoiMeshGui->roButton->do_callback();
				file->voronoiMeshGui->betaSlider->activate();
				file->voronoiMeshGui->polynomialOrderSlider->deactivate();
				file->voronoiMeshGui->priorGroup->activate();
				file->voronoiMeshGui->dPosterioriButton->activate();
				file->voronoiMeshGui->fPosteriorButton->activate();
				file->voronoiMeshGui->vPosteriorButton->deactivate();
				file->voronoiMeshGui->potentialReferenceButton->deactivate();
				file->voronoiMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->voronoiMeshGui->overlayForceMagnitudeButton->copy_label("Drift Magnitude");
				file->voronoiMeshGui->lambdaSlider->deactivate();
				file->voronoiMeshGui->muSlider->deactivate();
				file->voronoiMeshGui->smoothingPriorButton->activate();
				file->voronoiMeshGui->smoothingPriorButton->value(0);
				file->voronoiMeshGui->smoothingPriorButton->do_callback();
				file->voronoiMeshGui->pauseButton->label("@||");
				break;
			case 3: // DV
				file->voronoiMeshGui->smoothingPriorButton->activate();
				file->voronoiMeshGui->smoothingPriorButton->value(0);
				file->voronoiMeshGui->smoothingPriorButton->do_callback();
				file->voronoiMeshGui->roButton->value(0);
				file->voronoiMeshGui->roButton->activate();
				file->voronoiMeshGui->betaSlider->deactivate();
				file->voronoiMeshGui->polynomialOrderSlider->deactivate();
				file->voronoiMeshGui->priorGroup->activate();
				file->voronoiMeshGui->dPosterioriButton->activate();
				file->voronoiMeshGui->fPosteriorButton->deactivate();
				file->voronoiMeshGui->vPosteriorButton->activate();
				file->voronoiMeshGui->potentialReferenceButton->activate();
				file->voronoiMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->voronoiMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
				file->voronoiMeshGui->lambdaSlider->deactivate();
				file->voronoiMeshGui->muSlider->deactivate();
				break;
			case 4: // Poly
				file->voronoiMeshGui->roButton->deactivate();
				file->voronoiMeshGui->roButton->value(0);
				file->voronoiMeshGui->roButton->do_callback();
				file->voronoiMeshGui->betaSlider->deactivate();
				file->voronoiMeshGui->polynomialOrderSlider->activate();
				file->voronoiMeshGui->priorGroup->activate();
				file->voronoiMeshGui->dPosterioriButton->activate();
				file->voronoiMeshGui->fPosteriorButton->activate();
				file->voronoiMeshGui->vPosteriorButton->deactivate();
				file->voronoiMeshGui->potentialReferenceButton->deactivate();
				file->voronoiMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->voronoiMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
				file->voronoiMeshGui->smoothingPriorButton->deactivate();
				file->voronoiMeshGui->lambdaSlider->deactivate();
				file->voronoiMeshGui->muSlider->deactivate();
				file->voronoiMeshGui->smoothingPriorButton->value(0);
				file->voronoiMeshGui->smoothingPriorButton->do_callback();
				file->voronoiMeshGui->pauseButton->label("@||");
				break;
			default:
				break;
		}
			break;
		case 2: // quad-tree mesh
			switch(choice) {
			case 0: // D
				file->treeMeshGui->roButton->deactivate();
				file->treeMeshGui->roButton->value(0);
				file->treeMeshGui->roButton->do_callback();
				file->treeMeshGui->betaSlider->deactivate();
				file->treeMeshGui->polynomialOrderSlider->deactivate();
				file->treeMeshGui->priorGroup->activate();
				file->treeMeshGui->dPosteriorButton->activate();
				file->treeMeshGui->fPosteriorButton->deactivate();
				file->treeMeshGui->vPosteriorButton->deactivate();
				file->treeMeshGui->potentialReferenceButton->deactivate();
				file->treeMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->treeMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
				file->treeMeshGui->smoothingPriorButton->activate();
				file->treeMeshGui->lambdaSlider->deactivate();
				file->treeMeshGui->muSlider->deactivate();
				file->treeMeshGui->smoothingPriorButton->value(0);
				file->treeMeshGui->smoothingPriorButton->do_callback();
				file->treeMeshGui->pauseButton->label("@||");
				break;
			case 1: // DF
				file->treeMeshGui->roButton->deactivate();
				file->treeMeshGui->roButton->value(0);
				file->treeMeshGui->roButton->do_callback();
				file->treeMeshGui->betaSlider->activate();
				file->treeMeshGui->polynomialOrderSlider->deactivate();
				file->treeMeshGui->priorGroup->activate();
				file->treeMeshGui->dPosteriorButton->activate();
				file->treeMeshGui->fPosteriorButton->activate();
				file->treeMeshGui->vPosteriorButton->deactivate();
				file->treeMeshGui->potentialReferenceButton->deactivate();
				file->treeMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->treeMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
				file->treeMeshGui->smoothingPriorButton->activate();
				file->treeMeshGui->lambdaSlider->deactivate();
				file->treeMeshGui->muSlider->deactivate();
				file->treeMeshGui->smoothingPriorButton->value(0);
				file->treeMeshGui->smoothingPriorButton->do_callback();
				file->treeMeshGui->pauseButton->label("@||");
				break;
			case 2: // DDr
				file->treeMeshGui->roButton->deactivate();
				file->treeMeshGui->roButton->value(0);
				file->treeMeshGui->roButton->do_callback();
				file->treeMeshGui->betaSlider->activate();
				file->treeMeshGui->polynomialOrderSlider->deactivate();
				file->treeMeshGui->priorGroup->activate();
				file->treeMeshGui->dPosteriorButton->activate();
				file->treeMeshGui->fPosteriorButton->activate();
				file->treeMeshGui->vPosteriorButton->deactivate();
				file->treeMeshGui->potentialReferenceButton->deactivate();
				file->treeMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->treeMeshGui->overlayForceMagnitudeButton->copy_label("Drift Magnitude");
				file->treeMeshGui->lambdaSlider->deactivate();
				file->treeMeshGui->muSlider->deactivate();
				file->treeMeshGui->smoothingPriorButton->activate();
				file->treeMeshGui->smoothingPriorButton->value(0);
				file->treeMeshGui->smoothingPriorButton->do_callback();
				file->treeMeshGui->pauseButton->label("@||");
				break;
			case 3: // DV
				file->treeMeshGui->smoothingPriorButton->activate();
				file->treeMeshGui->smoothingPriorButton->value(0);
				file->treeMeshGui->smoothingPriorButton->do_callback();
				file->treeMeshGui->roButton->value(0);
				file->treeMeshGui->roButton->activate();
				file->treeMeshGui->betaSlider->deactivate();
				file->treeMeshGui->polynomialOrderSlider->deactivate();
				file->treeMeshGui->priorGroup->activate();
				file->treeMeshGui->dPosteriorButton->activate();
				file->treeMeshGui->fPosteriorButton->deactivate();
				file->treeMeshGui->vPosteriorButton->activate();
				file->treeMeshGui->potentialReferenceButton->activate();
				file->treeMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->treeMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
				file->treeMeshGui->lambdaSlider->deactivate();
				file->treeMeshGui->muSlider->deactivate();
				break;
			case 4: // Poly
				file->treeMeshGui->roButton->deactivate();
				file->treeMeshGui->roButton->value(0);
				file->treeMeshGui->roButton->do_callback();
				file->treeMeshGui->betaSlider->deactivate();
				file->treeMeshGui->polynomialOrderSlider->activate();
				file->treeMeshGui->priorGroup->activate();
				file->treeMeshGui->dPosteriorButton->activate();
				file->treeMeshGui->fPosteriorButton->activate();
				file->treeMeshGui->vPosteriorButton->deactivate();
				file->treeMeshGui->potentialReferenceButton->deactivate();
				file->treeMeshGui->overlayForceArrowsButton->copy_label("Arrows");
				file->treeMeshGui->overlayForceMagnitudeButton->copy_label("Force Magnitude");
				file->treeMeshGui->smoothingPriorButton->deactivate();
				file->treeMeshGui->lambdaSlider->deactivate();
				file->treeMeshGui->muSlider->deactivate();
				file->treeMeshGui->smoothingPriorButton->value(0);
				file->treeMeshGui->smoothingPriorButton->do_callback();
				file->treeMeshGui->pauseButton->label("@||");
				break;
			default:
				break;
		}
		break;
	}
	Fl::redraw();
}

void optimizationFunctionChoiceCallback(Fl_Choice*w,int*v) {
	file->optimizationFunction = w->value();
}

void saveSquareMeshCallback(Fl_Button*w,int*v) {
	file->squareMesh->exportMesh();
}

void SquareMeshGui::samplePosterior(int type) {

	switch(type) {
		case 0: // diffusion
		{
			switch(file->optimizationMode) {
				case 0: // D
					{
						posterioriSamples = (int) posterioriSampleNumberSlider->value();
						const double minD = minPosteriorSlider->value();
						const double maxD = maxPosteriorSlider->value();
						posterioriRange = maxD-minD;
						posterioriIncrement = (posterioriRange)/posterioriSampleNumberSlider->value();
						double logMax = -10000000.0;
						if (posterioriSamples > 5) {
							posterioriValues = new double[posterioriSamples];
							double sampler[1];
							// calculate diffusion values for each of sample positions
							for (int g = 0; g < posterioriSamples; g++) {
								sampler[0] = minD + (double)g*posterioriIncrement;
								posterioriValues[g] = -dPosteriorSquare(sampler);
								if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
							}
							// calculate log values
							for (int h = 0; h < posterioriSamples; h++) {
								posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
							}
							drawPosteriorPlot(type,0);
						}
					}
					break;
				case 1: // DF
					{
						posterioriSamples = (int) posterioriSampleNumberSlider->value();
						const double minD = minPosteriorSlider->value();
						const double maxD = maxPosteriorSlider->value();
						posterioriRange = maxD-minD;
						posterioriIncrement = (posterioriRange)/posterioriSampleNumberSlider->value();
						double logMax = -10000000.0;
						if (posterioriSamples > 5) {
							posterioriValues = new double[posterioriSamples];
							double sampler[3];
							sampler[0] = file->squareMesh->getForceX(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
							sampler[1] = file->squareMesh->getForceY(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
							// calculate diffusion values for each of sample positions
							for (int g = 0; g < posterioriSamples; g++) {
								sampler[2] = minD + (double)g*posterioriIncrement;
								posterioriValues[g] = -dfPosteriorSquare(sampler);
								if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
							}
							// calculate log values
							for (int h = 0; h < posterioriSamples; h++) {
								posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
							}
							drawPosteriorPlot(type,0);
						}
					}
					break;
				case 2: // DDr
					{
						posterioriSamples = (int) posterioriSampleNumberSlider->value();
						const double minD = minPosteriorSlider->value();
						const double maxD = maxPosteriorSlider->value();
						posterioriRange = maxD-minD;
						posterioriIncrement = (posterioriRange)/posterioriSampleNumberSlider->value();
						double logMax = -10000000.0;
						if (posterioriSamples > 5) {
							posterioriValues = new double[posterioriSamples];
							double sampler[3];
							sampler[0] = file->squareMesh->getForceX(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
							sampler[1] = file->squareMesh->getForceY(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
							// calculate diffusion values for each of sample positions
							for (int g = 0; g < posterioriSamples; g++) {
								sampler[2] = minD + (double)g*posterioriIncrement;
								posterioriValues[g] = -ddrPosteriorSquare(sampler);
								if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
							}
							// calculate log values
							for (int h = 0; h < posterioriSamples; h++) {
								posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
							}
							drawPosteriorPlot(type,0);
						}
					}
					break;
				case 3: // DV
					{
						posterioriSamples = (int) posterioriSampleNumberSlider->value();
						const double minD = minPosteriorSlider->value();
						const double maxD = maxPosteriorSlider->value();
						posterioriIncrement = (maxD-minD)/posterioriSampleNumberSlider->value();
						posterioriRange = maxD-minD;
						double logMax = -10000000.0;
						if (posterioriSamples > 5) {
							posterioriValues = new double[posterioriSamples];
							// copy DV values to array for posteriori evaluation
							int sampleIndex = file->squareMesh->getIdentifier(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
							double *sampler = new double[2*file->squareMesh->getTotalVariables()];
							// initialize potentials
							int i = 0;
							for (int a = 0; a < file->squareMesh->xCells; a++) {
								for (int b = 0; b < file->squareMesh->yCells; b++) {
									if ((file->squareMesh->active(a,b))) {
										sampler[2*i] = file->squareMesh->getDiffusion(a,b);
										sampler[2*i+1] = file->squareMesh->getPotential(a,b);
										i++;
									}
								}
							}
							// calculate diffusion values for each of sample positions
							for (int g = 0; g < posterioriSamples; g++) {
								sampler[2*sampleIndex] = minD + (double)g*posterioriIncrement;
								if (file->squareMesh->selectionMode()) { posterioriValues[g] = -dvPosteriorSquareSelection(sampler); }
								else { posterioriValues[g] = -dvPosteriorSquare(sampler); }
								if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
							}
							// calculate log values
							for (int h = 0; h < posterioriSamples; h++) {
								posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
							}
							drawPosteriorPlot(type,0);
							delete [] sampler;
						}
					}
					break;
			}
			break;
		}
		case 1: // force magnitude
		{
			if (file->optimizationMode == 1) { // only for DF mode
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				const double minF = minPosteriorSlider->value();
				const double maxF = maxPosteriorSlider->value();
				const double Fx = file->squareMesh->getForceX(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
				const double Fy = file->squareMesh->getForceY(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
				posterioriIncrement = (maxF-minF)/posterioriSampleNumberSlider->value();
				posterioriRange = maxF-minF;

				double logMax = -10000000.0;

				if (posterioriSamples > 5) {
					posterioriValues = new double[posterioriSamples];
					double sampler[3];
					sampler[2] = file->squareMesh->getDiffusion(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
					// calculate diffusion values for each of sample positions
					const double theta = atan2(Fy,Fx);
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[0] = (Fx-fabs(posterioriRange/2.0)*cos(theta)) + (double)g*posterioriIncrement*cos(theta);
						sampler[1] = (Fy-fabs(posterioriRange/2.0)*sin(theta)) + (double)g*posterioriIncrement*sin(theta);
						posterioriValues[g] = -dfPosteriorSquare(sampler);
						if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
					}

					// calculate log values
					for (int h = 0; h < posterioriSamples; h++) {
						posterioriValues[h] = exp(posterioriValues[h]-logMax);
					}
					drawPosteriorPlot(type,0);
				}
			}
			else if (file->optimizationMode == 2) { // only for DF mode
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				const double minF = minPosteriorSlider->value();
				const double maxF = maxPosteriorSlider->value();
				const double Fx = file->squareMesh->getForceX(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
				const double Fy = file->squareMesh->getForceY(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
				posterioriIncrement = (maxF-minF)/posterioriSampleNumberSlider->value();
				posterioriRange = maxF-minF;

				double logMax = -10000000.0;

				if (posterioriSamples > 5) {
					posterioriValues = new double[posterioriSamples];
					double sampler[3];
					sampler[2] = file->squareMesh->getDiffusion(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
					// calculate diffusion values for each of sample positions
					const double theta = atan2(Fy,Fx);
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[0] = (Fx-fabs(posterioriRange/2.0)*cos(theta)) + (double)g*posterioriIncrement*cos(theta);
						sampler[1] = (Fy-fabs(posterioriRange/2.0)*sin(theta)) + (double)g*posterioriIncrement*sin(theta);
						posterioriValues[g] = -ddrPosteriorSquare(sampler);
						if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
					}

					// calculate log values
					for (int h = 0; h < posterioriSamples; h++) {
						posterioriValues[h] = exp(posterioriValues[h]-logMax);
					}
					drawPosteriorPlot(type,0);
				}
			}
			break;
		}
		case 2: // potential
		{
			if (file->optimizationMode == 3 && this->potentialReferenceButton->value()) {
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				posterioriRange = fabs(maxPosteriorSlider->value());
				double logMax = -10000000.0;
				if (posterioriSamples > 5) {

					const int v0_id = file->squareMesh->getIdentifier(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
					const double v0 = file->squareMesh->getPotential(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
					const int v1_id = file->squareMesh->getIdentifier(file->squareMesh->getCurrentZone2X(),file->squareMesh->getCurrentZone2Y());
					const double v1 = file->squareMesh->getPotential(file->squareMesh->getCurrentZone2X(),file->squareMesh->getCurrentZone2Y());

					posterioriIncrement = posterioriRange/(double)posterioriSamples;

					posterioriValues = new double[2*posterioriSamples-1];

					// create delta matrix
					for (int s = 0; s < 2*posterioriSamples-1; s++) {
						posterioriValues[s] = 0.0;
					}

					double *sampler = new double[2*file->squareMesh->getTotalVariables()];

					// initialize potentials
					int i = 0;
					for (int a = 0; a < file->squareMesh->xCells; a++) {
						for (int b = 0; b < file->squareMesh->yCells; b++) {
							if ((file->squareMesh->active(a,b))) {
								sampler[2*i] = file->squareMesh->getDiffusion(a,b);
								sampler[2*i+1] = file->squareMesh->getPotential(a,b);
								i++;
							}
						}
					}

					double MAP_posterior;
					if (file->squareMesh->selectionMode()) { MAP_posterior = -dvPosteriorSquareSelection(sampler); }
					else { MAP_posterior = -dvPosteriorSquare(sampler); }

					// calculate potential values for each of sample positions
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[2*v0_id+1] = v0 - posterioriRange/2.0 + (double)g*posterioriIncrement;
						for (int h = 0; h < posterioriSamples; h++) {
							sampler[2*v1_id+1] = v1 - posterioriRange/2.0 + (double)h*posterioriIncrement;
							if (file->squareMesh->selectionMode()) {
								posterioriValues[ h-g+(posterioriSamples-1) ] += exp(-MAP_posterior-dvPosteriorSquareSelection(sampler));
							}
							else {
								posterioriValues[ h-g+(posterioriSamples-1) ] += exp(-MAP_posterior-dvPosteriorSquare(sampler));
							}
						}
					}

					// average the posteriori
					for (int e = 0; e < 2*posterioriSamples-1; e++) {
						posterioriValues[e] /= (posterioriSamples*posterioriSamples);
						if (logMax < posterioriValues[e]) { logMax = posterioriValues[e]; }
					}

					// calculate log values
					for (int h = 0; h < 2*posterioriSamples-1; h++) { posterioriValues[h] /= logMax; }

					drawPosteriorPlot(type,0);

					delete [] sampler;
				}
			} else if (file->optimizationMode == 3 && !this->potentialReferenceButton->value()) {
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				const double minV = minPosteriorSlider->value();
				const double maxV = maxPosteriorSlider->value();
				posterioriRange = maxV-minV;
				double logMax = -10000000.0;
				if (posterioriSamples > 5) {

					const int v0_id = file->squareMesh->getIdentifier(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());

					posterioriIncrement = posterioriRange/posterioriSampleNumberSlider->value();

					posterioriValues = new double[posterioriSamples];

					double *sampler = new double[2*file->squareMesh->getTotalVariables()];

					// initialize potentials
					int i = 0;
					for (int a = 0; a < file->squareMesh->xCells; a++) {
						for (int b = 0; b < file->squareMesh->yCells; b++) {
							if ((file->squareMesh->active(a,b))) {
								sampler[2*i] = file->squareMesh->getDiffusion(a,b);
								sampler[2*i+1] = file->squareMesh->getPotential(a,b);
								i++;
							}
						}
					}

					// calculate diffusion values for each of sample positions
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[2*v0_id+1] = minV + (double)g*posterioriIncrement;
						if (file->squareMesh->selectionMode()) { posterioriValues[g] = -dvPosteriorSquareSelection(sampler); }
						else { posterioriValues[g] = -dvPosteriorSquare(sampler); }
					}

					// average the posteriori
					for (int e = 0; e < posterioriSamples; e++) {
						if (logMax < posterioriValues[e]) { logMax = posterioriValues[e]; }
					}

					// calculate log values
					for (int h = 0; h < posterioriSamples; h++) {
						posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
					}
					drawPosteriorPlot(type,1);
					delete [] sampler;
				}

			}
			break;
		}
	}
}

void SquareMeshGui::drawPosteriorPlot(int type,int sup) {

	switch(type) {
		case 0: // diffusion
		{
			const double minD = minPosteriorSlider->value();
			const double maxD = maxPosteriorSlider->value();

			posterioriRange = maxD - minD;
			posterioriIncrement = (maxD-minD)/posterioriSampleNumberSlider->value();

			if (posterioriRange > 0.0 && file->inferred) {
				window->make_current();

				const int w = 320;
				const int h = 195;

	        	fl_line_style(FL_SOLID, 1, 0);

				fl_color(FL_WHITE); fl_rectf(170,35,w,h);
				fl_color(FL_BLACK); fl_rect(170,35,w,h);

				fl_line(170,190,170+w-1,190);
//				fl_font(iMAP->normalFont,10); fl_draw("[um\262/s]",170+w-43,225);

		        int pCurrent,pOld;
		        const double spacing = ((double)(w-1)/(double)posterioriSamples);
		        char xLabel [50];

		        for (int q = 1; q < posterioriSamples; q++) {
		        	// draw data points and connectors
		        	pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
		        	pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

		        	fl_color(FL_RED);
		        	fl_line_style(FL_SOLID, 2, 0);
					fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
		        }
		        fl_color(FL_BLACK);
		        fl_font(iMAP->normalFont,10);
		        // draw x-axis labels
		        for (int t = 1; t < 9; t++) {
					sprintf(xLabel,"%.3f",minD+t*(posterioriSamples/9.0)*posterioriIncrement);
					fl_draw(90,xLabel,177+t*w/9,227);
		        	fl_line_style(FL_DASH, 1, 0);
					fl_line(174+t*w/9,55,174+t*w/9,190);
		        }
				// draw MAP label
				sprintf(xLabel,"MAP: %.3f [um\262/s]",file->squareMesh->getDiffusion(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY()));
		        fl_font(iMAP->normalFont,10);
				fl_draw(xLabel,170+2,45);
			}
		}
		break;
		case 1: // force
		{
			const double minF = minPosteriorSlider->value();
			const double maxF = maxPosteriorSlider->value();
			const double Fx = file->squareMesh->getForceX(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
			const double Fy = file->squareMesh->getForceY(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());
			posterioriIncrement = (maxF-minF)/posterioriSampleNumberSlider->value();
			posterioriRange = maxF-minF;

			if (posterioriRange > 0.0 && file->inferred) {
				window->make_current();

				const int w = 320;
				const int h = 195;

	        	fl_line_style(FL_SOLID, 1, 0);

				fl_color(FL_WHITE); fl_rectf(170,35,w,h);
				fl_color(FL_BLACK); fl_rect(170,35,w,h);

				fl_line(170,190,170+w-1,190);
//				fl_font(iMAP->normalFont,10); fl_draw("[pN]",170+w-25,225);

		        int pCurrent,pOld;
		        const double spacing = ((double)(w-1)/(double)posterioriSamples);
		        char xLabel [50];

		        for (int q = 1; q < posterioriSamples; q++) {
		        	// draw data points and connectors
		        	pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
		        	pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

		        	fl_color(FL_RED);
		        	fl_line_style(FL_SOLID, 2, 0);
					fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
		        }
		        fl_color(FL_BLACK);
		        fl_font(iMAP->normalFont,10);
		        // draw x-axis labels

		        const double minF2 = sqrt(Fx*Fx+Fy*Fy) - posterioriRange/2.0;

		        for (int t = 1; t < 9; t++) {
					sprintf(xLabel,"%.1f",minF2+t*(posterioriSamples/9.0)*posterioriIncrement);
					fl_draw(90,xLabel,177+t*w/9,222);
		        	fl_line_style(FL_DASH, 1, 0);
					fl_line(174+t*w/9,55,174+t*w/9,190);
		        }
				// draw MAP label
				sprintf(xLabel,"MAP: %.1f [pN]",file->squareMesh->getForce(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY()));
		        fl_font(iMAP->normalFont,10);
				fl_draw(xLabel,170+2,45);
			}

		}
			break;
		case 2: // potential
		{
			switch(sup) {
				case 0:
				{
					const double maxSample = fabs(maxPosteriorSlider->value())/2.0;

					if (posterioriRange > 0.0 && file->inferred) {
						window->make_current();

						const int w = 320;
						const int h = 195;

						fl_line_style(FL_SOLID, 1, 0);

						fl_color(FL_WHITE); fl_rectf(170,35,w,h);
						fl_color(FL_BLACK); fl_rect(170,35,w,h);

						fl_line(170,190,170+w-1,190);
//						fl_font(iMAP->normalFont,10); fl_draw("[kT]",170+w-25,225);

						// draw MAP label
						const double deltaV = file->squareMesh->getPotential(file->squareMesh->getCurrentZone2X(),file->squareMesh->getCurrentZone2Y())
											  -file->squareMesh->getPotential(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY());

						int pCurrent,pOld;
						const double spacing = ((double)(w-1)/(double)(2*posterioriSamples));
						char xLabel [50];

						for (int q = 1; q < 2*posterioriSamples; q++) {
							// draw data points and connectors
							pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
							pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

							fl_color(FL_RED);
							fl_line_style(FL_SOLID, 2, 0);
							fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
						}
						fl_color(FL_BLACK);
						fl_font(iMAP->normalFont,10);
						// draw x-axis labels
						for (int t = 1; t < 9; t++) {
							sprintf(xLabel,"%.2f",deltaV-maxSample+(float)t/9.0*posterioriRange);
							fl_draw(90,xLabel,177+t*w/9,227);
							fl_line_style(FL_DASH, 1, 0);
							fl_line(174+t*w/9,55,174+t*w/9,190);
						}

						sprintf(xLabel,"MAP (reference): %.2f [kT]",deltaV);
						fl_font(iMAP->normalFont,10);
						fl_draw(xLabel,170+2,45);
					}

				}
					break;
				case 1:
				{
					const double minSample = minPosteriorSlider->value();
					const double maxSample = maxPosteriorSlider->value();

					posterioriRange = maxSample - minSample;

					if (posterioriRange > 0.0 && file->inferred) {
						window->make_current();

						const int w = 320;
						const int h = 195;

			        	fl_line_style(FL_SOLID, 1, 0);

						fl_color(FL_WHITE); fl_rectf(170,35,w,h);
						fl_color(FL_BLACK); fl_rect(170,35,w,h);

						fl_line(170,190,170+w-1,190);
//						fl_font(iMAP->normalFont,10); fl_draw("[kT]",170+w-25,225);

				        int pCurrent,pOld;
				        const double spacing = ((double)(w-1)/(double)posterioriSamples);
				        char xLabel [80];

				        for (int q = 1; q < posterioriSamples; q++) {
				        	// draw data points and connectors
				        	pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
				        	pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

				        	fl_color(FL_RED);
				        	fl_line_style(FL_SOLID, 2, 0);
							fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
				        }
				        fl_color(FL_BLACK);
				        fl_font(iMAP->normalFont,10);
				        // draw x-axis labels
				        for (int t = 1; t < 9; t++) {
							sprintf(xLabel,"%.2f",minSample+t*(posterioriSamples/9.0)*posterioriIncrement);
							fl_draw(90,xLabel,177+t*w/9,227);
				        	fl_line_style(FL_DASH, 1, 0);
							fl_line(174+t*w/9,55,174+t*w/9,190);
				        }
						// draw MAP label
						sprintf(xLabel,"MAP: %.2f [kT] (no reference)",file->squareMesh->getPotential(file->squareMesh->getCurrentZoneX(),file->squareMesh->getCurrentZoneY()));
				        fl_font(iMAP->normalFont,10);
						fl_draw(xLabel,170+2,45);
					}

				}
				break;
			}
			break;
		}
	}

}

void VoronoiMeshGui::show() {
	window->show();
	int cellNumber;
	if (iMAP->selectionButtonPressed) {	cellNumber = (int)((float)file->selection->cell.count/40.0f); }
	else { cellNumber = (int)((float)file->localizationCountVoronoi/40.0f); }
	clusteringChoice->value(0);
	maxIterationsSlider->value(50);
	polynomialOrderSlider->value(2);
	betaSlider->value(2.0f);
	inferenceModeChoice->value(0);
}

VoronoiMeshGui::VoronoiMeshGui() {

	if (iMAP->fileNumber > 0) {

		clickInformationButton = NULL;
		colormap = 9;

		const int w = 500;
		const int h = 300;

		window = new Fl_Window(w,h,"Voronoi Tessellation");
		window->callback(voronoiMeshWindowCallback);
		window->color(iMAP->bgColor);
		window->begin();
		window->set_non_modal();

		// Permanent interface buttons
		{
			overlayPointNumberButton = new Fl_Check_Button(5,245,65,20,"Points");
			overlayPointNumberButton->labelsize(12);
			overlayPointNumberButton->value(1);
			overlayPointNumberButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayPointNumberButton->labelfont(iMAP->normalFont);
			overlayPointNumberButton->labelcolor(FL_WHITE);
			overlayPointNumberButton->callback((Fl_Callback*)pointNumberOverlayCallback,(int*)1);
			overlayPointNumberButton->show();

			overlayDiffusionButton = new Fl_Check_Button(70,245,165,20,"Diffusion Coefficient");
			overlayDiffusionButton->labelsize(12);
			overlayDiffusionButton->callback((Fl_Callback*)diffusionOverlayCallback,(int*)1);
			overlayDiffusionButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayDiffusionButton->labelfont(iMAP->normalFont);
			overlayDiffusionButton->labelcolor(iMAP->bgColor);
			overlayDiffusionButton->deactivate();
			overlayDiffusionButton->show();

			overlayDiffusionLogButton = new Fl_Check_Button(70,260,165,20,"Log");
			overlayDiffusionLogButton->labelsize(12);
			overlayDiffusionLogButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayDiffusionLogButton->callback((Fl_Callback*)diffusionLogOverlayCallback,(int*)1);
			overlayDiffusionLogButton->labelfont(iMAP->normalFont);
			overlayDiffusionLogButton->labelcolor(iMAP->bgColor);
			overlayDiffusionLogButton->deactivate();
			overlayDiffusionLogButton->show();

			overlayForceMagnitudeButton = new Fl_Check_Button(240,245,135,20,"Force Magnitude");
			overlayForceMagnitudeButton->labelsize(12);
			overlayForceMagnitudeButton->callback((Fl_Callback*)forceMagnitudeOverlayCallback,(int*)1);
			overlayForceMagnitudeButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayForceMagnitudeButton->labelfont(iMAP->normalFont);
			overlayForceMagnitudeButton->labelcolor(iMAP->bgColor);
			overlayForceMagnitudeButton->deactivate();
			overlayForceMagnitudeButton->show();

			overlayForceArrowsButton = new Fl_Check_Button(240,260,135,20,"Arrows");
			overlayForceArrowsButton->labelsize(12);
			overlayForceArrowsButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayForceArrowsButton->labelfont(iMAP->normalFont);
			overlayForceArrowsButton->labelcolor(iMAP->bgColor);
			overlayForceArrowsButton->deactivate();
			overlayForceArrowsButton->show();

			overlayPotentialButton = new Fl_Check_Button(365,245,140,20,"Potential Energy");
			overlayPotentialButton->labelsize(12);
			overlayPotentialButton->callback((Fl_Callback*)potentialOverlayCallback,(int*)1);
			overlayPotentialButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayPotentialButton->labelfont(iMAP->normalFont);
			overlayPotentialButton->labelcolor(iMAP->bgColor);
			overlayPotentialButton->deactivate();
			overlayPotentialButton->show();

			variablesBox = new Fl_Box(0,275,w,25,"0 Zones");
			variablesBox->labelsize(14);
			variablesBox->align(FL_ALIGN_CENTER);
			variablesBox->labelfont(iMAP->boldFont);
			variablesBox->labelcolor(FL_WHITE);
			variablesBox->show();
		}

		tabs = new Fl_Tabs(0,0,w,h-60);
		tabs->labelsize(12);
		tabs->box(FL_BORDER_FRAME);
		tabs->labelcolor(FL_WHITE);
		tabs->labelfont(iMAP->normalFont);
		tabs->begin();

		inferenceGroup = new Fl_Group(0,25,w,h-85,"Inference");
		inferenceGroup->labelsize(12);
		inferenceGroup->labelcolor(FL_DARK1);
		inferenceGroup->labelfont(iMAP->normalFont);
		inferenceGroup->box(FL_BORDER_FRAME);
		inferenceGroup->begin();
		{
			cellsSlider = new Slider(10,45,w-20,20,"Number of Zones");
			cellsSlider->precision(0);
			cellsSlider->labelsize(12);
			cellsSlider->labelcolor(FL_WHITE);
			cellsSlider->labelfont(iMAP->boldFont);
			int cellNumber;
			if (iMAP->selectionButtonPressed) {	cellNumber = (int)((float)file->selection->cell.count/80.0f); }
			else { cellNumber = (int)((float)file->localizationCountVoronoi/80.0f); }
			cellsSlider->value(cellNumber);
			cellsSlider->callback((Fl_Callback*)nullCallback);
			cellsSlider->bounds(cellNumber/2,1.5*cellNumber);
			cellsSlider->show();

			minPointsSlider = new Slider(10,85,w-20,20,"Minimum Points / Zone");
			minPointsSlider->precision(0);
			minPointsSlider->labelsize(12);
			minPointsSlider->value(20);
			minPointsSlider->labelcolor(FL_WHITE);
			minPointsSlider->labelfont(iMAP->boldFont);
			minPointsSlider->callback((Fl_Callback*)minPointsSliderCallback);
			minPointsSlider->show();
			minPointsSlider->deactivate();

			clusteringChoice = new Fl_Choice(90,125,90,20,"Clustering");
			clusteringChoice->menu(clusteringMenu);
			clusteringChoice->textsize(12);
			clusteringChoice->textfont(iMAP->boldFont);
			clusteringChoice->labelsize(12);
			clusteringChoice->labelfont(iMAP->boldFont);
			clusteringChoice->labelcolor(FL_WHITE);
			clusteringChoice->callback((Fl_Callback*)clusteringChoiceCallback);
			clusteringChoice->show();

			distanceButton = new DistanceButton(205,125,145,20);
			distanceButton->show();

			maxIterationsSlider = new Slider(360,125,130,20,"Maximum Iterations");
			maxIterationsSlider->precision(0);
			maxIterationsSlider->labelsize(12);
			maxIterationsSlider->labelfont(iMAP->boldFont);
			maxIterationsSlider->labelcolor(FL_WHITE);
			maxIterationsSlider->value(50);
			maxIterationsSlider->bounds(10,100);
			maxIterationsSlider->show();

			inferenceModeChoice = new Fl_Choice(10,170,(w-24)/2,20,"Mode");
			inferenceModeChoice->menu(optmizationMenu);
			inferenceModeChoice->align(FL_ALIGN_TOP);
			inferenceModeChoice->textsize(11);
			inferenceModeChoice->textfont(iMAP->boldFont);
			inferenceModeChoice->labelcolor(FL_WHITE);
			inferenceModeChoice->labelfont(iMAP->boldFont);
//			inferenceModeChoice->box(FL_BORDER_FRAME);
			inferenceModeChoice->labelsize(11);
			inferenceModeChoice->color(FL_GRAY,FL_GRAY);
			inferenceModeChoice->callback((Fl_Callback*)optimizationChoiceCallback,(int*)1);
			inferenceModeChoice->deactivate();
			inferenceModeChoice->show();

			localizationPrecisionSlider = new Slider(15+(w-24)/2,170,(w-24)/2,20,"Localization Precision [nm]");
			localizationPrecisionSlider->labelsize(12);
			localizationPrecisionSlider->labelcolor(FL_WHITE);
			localizationPrecisionSlider->labelfont(iMAP->boldFont);
			localizationPrecisionSlider->precision(0);
			localizationPrecisionSlider->value(30);
			localizationPrecisionSlider->bounds(0,100);
			localizationPrecisionSlider->show();


			applyButton = new Fl_Button(10,200,(w-40)/6,30,"Apply");
			applyButton->labelsize(12);
			applyButton->labelfont(iMAP->boldFont);
			applyButton->callback((Fl_Callback*)applyMeshCallback,(int*)1);
			applyButton->show();

			const int buttonWidth = 76;

			resetButton = new Fl_Button(10+buttonWidth+5,200,buttonWidth,30,"Reset");
			resetButton->labelsize(12);
			resetButton->callback((Fl_Callback*)resetMeshCallback,(int*)1);
			resetButton->labelfont(iMAP->boldFont);
			resetButton->show();
			resetButton->deactivate();

			inferButton = new Fl_Button(10+2*buttonWidth+10,200,buttonWidth,30,"Infer");
			inferButton->labelsize(12);
			inferButton->labelfont(iMAP->boldFont);
			inferButton->callback((Fl_Callback*)inferMeshCallback,(int*)1);
			inferButton->show();
			inferButton->deactivate();

			pauseButton = new Fl_Button(10+3*buttonWidth+15,200,buttonWidth,30,"Stop");
			pauseButton->labelsize(14);
			pauseButton->label("@||");
//			pauseButton->labelcolor(iMAP->bgColor);
			pauseButton->callback((Fl_Callback*)pauseCalculationButtonCallback,(int*)1);
			pauseButton->deactivate();
			pauseButton->show();

			stopButton = new Fl_Button(10+4*buttonWidth+20,200,buttonWidth,30,"Stop");
			stopButton->labelsize(12);
			stopButton->label("@square");
			stopButton->labelcolor(FL_RED);
			stopButton->callback((Fl_Callback*)stopCalculationButtonCallback,(int*)1);
			stopButton->deactivate();
			stopButton->show();

			saveButton = new Fl_Button(10+5*buttonWidth+25,200,buttonWidth,30,"Save");
			saveButton->labelsize(12);
			saveButton->labelfont(iMAP->boldFont);
			saveButton->callback((Fl_Callback*)saveVoronoiMeshCallback);
			saveButton->show();
			saveButton->deactivate();

		}
		inferenceGroup->end();

		annotationsGroup = new Fl_Group(0,25,w,h-85,"Overlay");
		annotationsGroup->labelsize(12);
		annotationsGroup->labelcolor(FL_DARK1);
		annotationsGroup->box(FL_BORDER_FRAME);
		annotationsGroup->labelfont(iMAP->normalFont);
		annotationsGroup->begin();
		{
			// font options
			localizationNumberLabelButton = new Fl_Check_Button(10,45,150,20,"Localization Number");
			localizationNumberLabelButton->labelsize(12);
			localizationNumberLabelButton->labelcolor(FL_WHITE);
			localizationNumberLabelButton->labelfont(iMAP->boldFont);
			localizationNumberLabelButton->value(1);
			localizationNumberLabelButton->show();

			fontScaleSlider = new Slider(175,45,w-185,20,"Font Scale");
			fontScaleSlider->labelsize(12);
			fontScaleSlider->precision(2);
			fontScaleSlider->bounds(0.5f,1.5f);
			fontScaleSlider->value(1.0f);
			fontScaleSlider->labelfont(iMAP->boldFont);
			fontScaleSlider->labelcolor(FL_WHITE);
			fontScaleSlider->show();

			// grid options
			displayGridButton = new Fl_Check_Button(10,85,105,20,"Display Grid");
			displayGridButton->labelsize(12);
			displayGridButton->labelcolor(FL_WHITE);
			displayGridButton->labelfont(iMAP->boldFont);
			displayGridButton->value(1);
			displayGridButton->show();

			gridColorButton = new Fl_Button(120,85,85,20,"Grid Color");
			gridColorButton->labelsize(12);
			gridColorButton->labelfont(iMAP->boldFont);
			gridColorButton->callback((Fl_Callback*)gridColorButtonCallback,(int*)1);
			gridColorButton->show();

			gridRGB[0] = gridRGB[1] = gridRGB[2] = 1.0;

			gridAlphaSlider = new Slider(215,85,w-225,20,"Grid Alpha");
			gridAlphaSlider->labelsize(12);
			gridAlphaSlider->precision(2);
			gridAlphaSlider->labelcolor(FL_WHITE);
			gridAlphaSlider->labelfont(iMAP->boldFont);
			gridAlphaSlider->value(1.0f);
			gridAlphaSlider->bounds(0.0,1.0);
			gridAlphaSlider->show();

			// color overlay options
			overlayButton = new Fl_Check_Button(10,125,115,20,"Color Overlay");
			overlayButton->labelsize(12);
			overlayButton->labelcolor(FL_WHITE);
			overlayButton->labelfont(iMAP->boldFont);
			overlayButton->value(1);
			overlayButton->show();

			overlayColorChoice = new Fl_Choice(155,125,70,20,"Map");
			overlayColorChoice->labelsize(12);
			overlayColorChoice->menu(colormapItemsVoronoi);
			overlayColorChoice->textsize(12);
			overlayColorChoice->textfont(iMAP->boldFont);
			overlayColorChoice->labelcolor(FL_WHITE);
			overlayColorChoice->labelfont(iMAP->boldFont);
			overlayColorChoice->callback((Fl_Callback*)colorOverlayCallback,(int*)1);
			overlayColorChoice->value(colormap); // jet by default
			overlayColorChoice->show();

			flipColormapButton = new Fl_Check_Button(225,125,40,20,"Flip");
			flipColormapButton->callback((Fl_Callback*)flipColormapCallback,(int*)1);
			flipColormapButton->labelsize(12);
			flipColormapButton->labelcolor(FL_WHITE);
			flipColormapButton->labelfont(iMAP->boldFont);
			flipColormapButton->value(0);
			flipColormapButton->show();

			overlayAlphaSlider = new Slider(280,125,w-290,20,"Overlay Alpha");
			overlayAlphaSlider->labelsize(12);
			overlayAlphaSlider->precision(2);
			overlayAlphaSlider->bounds(0.0f,1.0f);
			overlayAlphaSlider->value(0.5f);
			overlayAlphaSlider->labelfont(iMAP->boldFont);
			overlayAlphaSlider->labelcolor(FL_WHITE);
			overlayAlphaSlider->callback((Fl_Callback*)colorOverlayCallback,(int*)1);
			overlayAlphaSlider->show();

			// spot visualization options
			spotVisualizationButton = new Fl_Check_Button(10,165,135,20,"Spot Visualization");
			spotVisualizationButton->labelsize(12);
			spotVisualizationButton->labelcolor(FL_WHITE);
			spotVisualizationButton->labelfont(iMAP->boldFont);
			spotVisualizationButton->value(0);
			spotVisualizationButton->show();

			spotVisualizationScaleSlider = new Slider(160,165,w-170,20,"Scale");
			spotVisualizationScaleSlider->labelsize(12);
			spotVisualizationScaleSlider->precision(2);
			spotVisualizationScaleSlider->bounds(0.0f,5.0f);
			spotVisualizationScaleSlider->labelcolor(FL_WHITE);
			spotVisualizationScaleSlider->labelfont(iMAP->boldFont);
			spotVisualizationScaleSlider->value(1.0f);
			spotVisualizationScaleSlider->show();

			// hover information
			infoOverlayButton = new Fl_Check_Button(10,190,140,20,"Hover Information");
			infoOverlayButton->labelsize(12);
			infoOverlayButton->value(0);
			infoOverlayButton->labelfont(iMAP->boldFont);
			infoOverlayButton->labelcolor(FL_WHITE);
			infoOverlayButton->deactivate();
			infoOverlayButton->show();

			neighbourDistanceViewButton = new Fl_Check_Button(10,210,195,20,"View Neighbor Connections");
			neighbourDistanceViewButton->labelcolor(FL_WHITE);
			neighbourDistanceViewButton->labelfont(iMAP->boldFont);
			neighbourDistanceViewButton->labelsize(12);
			neighbourDistanceViewButton->value(1);
			neighbourDistanceViewButton->show();

		}
		annotationsGroup->deactivate();
		annotationsGroup->end();

		advancedGroup = new Fl_Group(0,25,w,h-85,"Advanced");
		advancedGroup->labelsize(12);
		advancedGroup->labelcolor(FL_DARK1);
		advancedGroup->box(FL_BORDER_FRAME);
		advancedGroup->labelfont(iMAP->normalFont);
		advancedGroup->begin();
		{

			roButton = new Fl_Check_Button(10,30,180,20,"Randomized Optimization");
			roButton->labelcolor(FL_WHITE);
			roButton->labelfont(iMAP->boldFont);
			roButton->labelsize(12);
			roButton->value(0);
			roButton->callback((Fl_Callback*)roButtonCallback,(int*)1);
			roButton->deactivate();
			roButton->show();

			roRadiusSlider = new Slider(10,65,(w-20)/3-5,20,"Selection Radius [nm]");
			roRadiusSlider->labelsize(12);
			roRadiusSlider->labelcolor(FL_WHITE);
			roRadiusSlider->labelfont(iMAP->boldFont);
			roRadiusSlider->precision(0);
			roRadiusSlider->value(500);
			roRadiusSlider->bounds(100,2000);
			roRadiusSlider->deactivate();
			roRadiusSlider->show();

			roToleranceSlider = new Slider(12+(w-20)/3,65,w/3-10,20,"Cost Tolerance [%]");
			roToleranceSlider->labelsize(12);
			roToleranceSlider->labelcolor(FL_WHITE);
			roToleranceSlider->labelfont(iMAP->boldFont);
			roToleranceSlider->precision(4);
			roToleranceSlider->value(0.001);
			roToleranceSlider->bounds(0.0,1.0);
			roToleranceSlider->callback((Fl_Callback*)randomizedOptimizationToleranceCallback,(int*)1);
			roToleranceSlider->deactivate();
			roToleranceSlider->show();

			roMaximumIterationsSlider = new Slider(14+2*(w-20)/3,65,w/3-10,20,"Maximum Iterations");
			roMaximumIterationsSlider->labelsize(12);
			roMaximumIterationsSlider->labelcolor(FL_WHITE);
			roMaximumIterationsSlider->labelfont(iMAP->boldFont);
			roMaximumIterationsSlider->precision(0);
			roMaximumIterationsSlider->value(5);
			roMaximumIterationsSlider->bounds(3,10);
			roMaximumIterationsSlider->deactivate();
			roMaximumIterationsSlider->show();

			neighbourDistanceSlider = new Slider(10,105,w-20,20,"Maximum Neighbor Distance [nm]");
			neighbourDistanceSlider->labelsize(12);
			neighbourDistanceSlider->labelcolor(FL_WHITE);
			neighbourDistanceSlider->labelfont(iMAP->boldFont);
			neighbourDistanceSlider->precision(1);
			neighbourDistanceSlider->value(1000.0);
			neighbourDistanceSlider->bounds(50,2000);
			neighbourDistanceSlider->callback((Fl_Callback*)maximumNeighbourDistanceCallback,(int*)1);
			neighbourDistanceSlider->show();

			Fl_Box *dfLabel = new Fl_Box(15,145,100,20,"(D,F) Inference:");
			dfLabel->labelsize(12);
			dfLabel->labelcolor(FL_WHITE);
			dfLabel->labelfont(iMAP->boldFont);
			dfLabel->show();

			betaSlider = new Slider(130,145,w-140,20,"Potential Energy Penalization (Beta)");
			betaSlider->precision(1);
			betaSlider->labelcolor(FL_WHITE);
			betaSlider->labelfont(iMAP->boldFont);
			betaSlider->bounds(0.0f,5.0f);
			betaSlider->value(2.0f);
			betaSlider->deactivate();
			betaSlider->show();

			Fl_Box *polynomialLabel = new Fl_Box(10,185,220,20,"Polynomial Potential Inference:");
			polynomialLabel->labelsize(12);
			polynomialLabel->labelcolor(FL_WHITE);
			polynomialLabel->labelfont(iMAP->boldFont);
			polynomialLabel->show();

			polynomialOrderSlider = new Slider(235,185,w-245,20,"Polynomial Order");
			polynomialOrderSlider->precision(0);
			polynomialOrderSlider->labelsize(12);
			polynomialOrderSlider->labelcolor(FL_WHITE);
			polynomialOrderSlider->value(2);
			polynomialOrderSlider->labelfont(iMAP->boldFont);
			polynomialOrderSlider->bounds(2,8);
			polynomialOrderSlider->show();
			polynomialOrderSlider->deactivate();

		}
		advancedGroup->deactivate();
		advancedGroup->end();

		priorGroup = new Fl_Group(0,25,w,h-85,"Prior");
		priorGroup->labelsize(12);
//		priorGroup->color(FL_WHITE,FL_WHITE);
		priorGroup->labelcolor(FL_DARK1);
		priorGroup->box(FL_BORDER_FRAME);
		priorGroup->labelfont(iMAP->normalFont);
		priorGroup->begin();
		{

			jeffreysPriorButton = new Fl_Check_Button(10,30,180,25,"Enable Jeffreys Prior");
			jeffreysPriorButton->labelsize(12);
			jeffreysPriorButton->labelcolor(FL_WHITE);
			jeffreysPriorButton->labelfont(iMAP->boldFont);
			jeffreysPriorButton->value(1);
			jeffreysPriorButton->callback((Fl_Callback*)jeffreysPriorButtonCallback,(int*)1);
			jeffreysPriorButton->show();

			smoothingPriorButton = new Fl_Check_Button(10,50,180,25,"Enable Smoothing Prior");
			smoothingPriorButton->labelsize(12);
			smoothingPriorButton->labelcolor(FL_WHITE);
			smoothingPriorButton->labelfont(iMAP->boldFont);
			smoothingPriorButton->value(0);
			smoothingPriorButton->callback((Fl_Callback*)smoothingPriorButtonCallback,(int*)1);
			smoothingPriorButton->show();
//			smoothingPriorButton->deactivate();

			lambdaSlider = new Slider(10,85,(w-30)/2,20,"V Penalisation (lambda)");
			lambdaSlider->precision(3);
			lambdaSlider->labelcolor(FL_WHITE);
			lambdaSlider->labelfont(iMAP->boldFont);
			lambdaSlider->labelsize(12);
			lambdaSlider->value(0.1);
			lambdaSlider->bounds(0.0,10);
			lambdaSlider->deactivate();
			lambdaSlider->show();

			muSlider = new Slider((w-30)/2+15,85,(w-30)/2,20,"D Penalization (mu)");
			muSlider->precision(3);
			muSlider->labelsize(12);
			muSlider->labelcolor(FL_WHITE);
			muSlider->labelfont(iMAP->boldFont);
			muSlider->value(0.1);
			muSlider->bounds(0.0,10.0);
			muSlider->deactivate();
			muSlider->show();

		}
//		priorGroup->deactivate();
		priorGroup->end();

		posterioriGroup = new Fl_Group(0,25,w,h-85,"Posterior");
		posterioriGroup->labelsize(12);
		posterioriGroup->labelcolor(FL_DARK1);
		posterioriGroup->box(FL_BORDER_FRAME);
		posterioriGroup->labelfont(iMAP->normalFont);
		posterioriGroup->begin();
		{
			dPosterioriButton = new Fl_Round_Button(10,30,155,20,"Diffusion");
			dPosterioriButton->labelfont(1);
			dPosterioriButton->labelsize(12);
			dPosterioriButton->labelfont(iMAP->boldFont);
			dPosterioriButton->type(FL_TOGGLE_BUTTON);
			dPosterioriButton->labelcolor(FL_WHITE);
			dPosterioriButton->callback((Fl_Callback*)dPosteriorCallback);
			dPosterioriButton->value(0);
			dPosterioriButton->deactivate();
			dPosterioriButton->show();

			fPosteriorButton = new Fl_Round_Button(10,50,155,20,"Force");
			fPosteriorButton->labelfont(iMAP->boldFont);
			fPosteriorButton->labelsize(12);
			fPosteriorButton->type(FL_TOGGLE_BUTTON);
			fPosteriorButton->labelcolor(FL_WHITE);
			fPosteriorButton->callback((Fl_Callback*)fPosterioriCallback);
			fPosteriorButton->value(0);
			fPosteriorButton->deactivate();
			fPosteriorButton->show();

			vPosteriorButton = new Fl_Round_Button(10,70,155,20,"Potential");
			vPosteriorButton->labelfont(iMAP->boldFont);
			vPosteriorButton->labelsize(12);
			vPosteriorButton->type(FL_TOGGLE_BUTTON);
			vPosteriorButton->labelcolor(FL_WHITE);
			vPosteriorButton->callback((Fl_Callback*)vPosterioriCallback);
			vPosteriorButton->value(0);
			vPosteriorButton->deactivate();
			vPosteriorButton->show();

			potentialReferenceButton = new Fl_Button(95,70,70,20,"Reference");
			potentialReferenceButton->labelfont(iMAP->boldFont);
			potentialReferenceButton->labelsize(12);
			potentialReferenceButton->type(FL_TOGGLE_BUTTON);
			potentialReferenceButton->callback((Fl_Callback*)potentialReferenceCallback);
			potentialReferenceButton->value(0);
			potentialReferenceButton->deactivate();
			potentialReferenceButton->show();

			minPosteriorSlider = new Slider(10,105,155,20,"Minimum");
			minPosteriorSlider->labelsize(12);
			minPosteriorSlider->labelfont(iMAP->boldFont);
			minPosteriorSlider->precision(3);
			minPosteriorSlider->labelcolor(FL_WHITE);
			minPosteriorSlider->value(0.001f);
			minPosteriorSlider->bounds(0.001,1.0);
			minPosteriorSlider->deactivate();
			minPosteriorSlider->show();

			maxPosteriorSlider = new Slider(10,140,155,20,"Maximum");
			maxPosteriorSlider->labelsize(12);
			maxPosteriorSlider->labelfont(iMAP->boldFont);
			maxPosteriorSlider->precision(3);
			maxPosteriorSlider->labelcolor(FL_WHITE);
			maxPosteriorSlider->value(1.0f);
			maxPosteriorSlider->bounds(0.001,1.0);
			maxPosteriorSlider->deactivate();
			maxPosteriorSlider->show();

			posterioriSampleNumberSlider = new Slider(10,175,155,20,"Samples");
			posterioriSampleNumberSlider->labelsize(12);
			posterioriSampleNumberSlider->labelfont(iMAP->boldFont);
			posterioriSampleNumberSlider->precision(0);
			posterioriSampleNumberSlider->labelcolor(FL_WHITE);
			posterioriSampleNumberSlider->bounds(100,1000);
			posterioriSampleNumberSlider->value(100);
			posterioriSampleNumberSlider->deactivate();
			posterioriSampleNumberSlider->show();

			samplePosteriorButton = new Fl_Button(10,205,75,25,"Sample");
			samplePosteriorButton->labelfont(iMAP->boldFont);
			samplePosteriorButton->labelsize(12);
			samplePosteriorButton->callback((Fl_Callback*)samplePosteriorCallback,(int*)1);
//			samplePosterioriButton->deactivate();
			samplePosteriorButton->show();

			savePosteriorButton = new Fl_Button(90,205,75,25,"Save");
			savePosteriorButton->labelfont(iMAP->boldFont);
			savePosteriorButton->labelsize(12);
			savePosteriorButton->callback((Fl_Callback*)savePosteriorCallback,(int*)1);
			savePosteriorButton->deactivate();
			savePosteriorButton->show();

			// draw box
			Fl_Box *posterioriBox = new Fl_Box(170,35,320,195);
			posterioriBox->color(FL_WHITE);

		}
		posterioriGroup->deactivate();
		posterioriGroup->end();

		simulationGroup = new Fl_Group(0,25,w,h-85,"Simulation");
		simulationGroup->labelsize(12);
		simulationGroup->labelcolor(FL_DARK1);
		simulationGroup->box(FL_BORDER_FRAME);
		simulationGroup->labelfont(iMAP->normalFont);
		simulationGroup->begin();
		{
			numberOfTrajectoriesSlider = new Slider(10,45,w-20,20,"Number of Trajectories");
			numberOfTrajectoriesSlider->precision(0);
			numberOfTrajectoriesSlider->labelfont(iMAP->boldFont);
			numberOfTrajectoriesSlider->labelcolor(FL_WHITE);
			numberOfTrajectoriesSlider->bounds(100,1000);
			numberOfTrajectoriesSlider->value(500);
//			numberOfTrajectoriesSlider->callback((Fl_Callback*)simulationWidgetsCallback,(int*)1);
			numberOfTrajectoriesSlider->activate();

			deltaTSlider = new Slider(10,85,(w-30)/2,20,"Delta [ms]");
			deltaTSlider->labelsize(12);
			deltaTSlider->precision(0);
			deltaTSlider->labelfont(iMAP->boldFont);
			deltaTSlider->labelcolor(FL_WHITE);
			deltaTSlider->value(25);
			deltaTSlider->bounds(0,100);
//			deltaTSlider->callback((Fl_Callback*)simulationWidgetsCallback,(int*)1);
			deltaTSlider->activate();

			timeStepsSlider = new Slider((w-30)/2+20,85,(w-30)/2,20,"Maximum Time Steps");
			timeStepsSlider->labelsize(12);
			timeStepsSlider->labelfont(iMAP->boldFont);
			timeStepsSlider->labelcolor(FL_WHITE);
			timeStepsSlider->precision(0);
			timeStepsSlider->value(25);
			timeStepsSlider->bounds(0,100);
//			timeStepsSlider->callback((Fl_Callback*)simulationWidgetsCallback,(int*)0);
			timeStepsSlider->activate();

			generateTrajectoriesButton = new Fl_Button(w/2-50,115,100,40,"Save\nTrajectories");
			generateTrajectoriesButton->labelsize(12);
			generateTrajectoriesButton->labelfont(iMAP->boldFont);
			generateTrajectoriesButton->callback((Fl_Callback*)generateTrajectoriesCallback,(int*)1);
			generateTrajectoriesButton->show();
			generateTrajectoriesButton->activate();
//
//			saveTrajectoriesButton = new Fl_Button(100,115,80,40,"Save\nTrajectories");
//			saveTrajectoriesButton->labelsize(12);
//			saveTrajectoriesButton->labelfont(iMAP->boldFont);
//			saveTrajectoriesButton->callback((Fl_Callback*)saveTrajectoriesCallback,(int*)0);
//			saveTrajectoriesButton->show();
//			saveTrajectoriesButton->deactivate();
		}
		simulationGroup->deactivate();
		simulationGroup->end();

		landscapeGroup = new Fl_Group(0,25,w,h-85,"Landscape");
		landscapeGroup->labelsize(12);
		landscapeGroup->labelcolor(FL_DARK1);
		landscapeGroup->box(FL_BORDER_FRAME);
		landscapeGroup->labelfont(iMAP->normalFont);
		landscapeGroup->begin();
		{
			landscapeButton = new Fl_Check_Button(10,30,120,25,"View Landscape");
			landscapeButton->labelsize(12);
			landscapeButton->labelfont(iMAP->boldFont);
			landscapeButton->labelcolor(FL_WHITE);
			landscapeButton->value(0);
			landscapeButton->callback((Fl_Callback*)voronoiLandscapeCallback);
//			landscapeButton->deactivate();
			landscapeButton->show();

			scaleLandscapeSlider = new Slider(10,70,w/2-15,20,"Scale");
			scaleLandscapeSlider->precision(2);
			scaleLandscapeSlider->labelfont(iMAP->boldFont);
			scaleLandscapeSlider->labelcolor(FL_WHITE);
			scaleLandscapeSlider->bounds(0,5.0);
			scaleLandscapeSlider->value(1.0);
			scaleLandscapeSlider->deactivate();
			scaleLandscapeSlider->show();

			landscapeAlphaSlider = new Slider(w/2+5,70,w/2-15,20,"Alpha");
			landscapeAlphaSlider->precision(2);
			landscapeAlphaSlider->bounds(0.01,1);
			landscapeAlphaSlider->labelfont(iMAP->boldFont);
			landscapeAlphaSlider->labelcolor(FL_WHITE);
			landscapeAlphaSlider->value(1.0);
			landscapeAlphaSlider->callback((Fl_Callback*)alphaVoronoiLandscapeCallback);
			landscapeAlphaSlider->deactivate();
			landscapeAlphaSlider->show();

			xLightPositionSlider = new Slider(10,105,w/2-15,20,"Light x");
			xLightPositionSlider->precision(2);
			xLightPositionSlider->bounds(-100,100);
			xLightPositionSlider->labelfont(iMAP->boldFont);
			xLightPositionSlider->labelcolor(FL_WHITE);
			xLightPositionSlider->value(0.0);
			xLightPositionSlider->deactivate();
			xLightPositionSlider->show();

			yLightPositionSlider = new Slider(10,140,w/2-15,20,"Light y");
			yLightPositionSlider->precision(2);
			yLightPositionSlider->bounds(-100,100);
			yLightPositionSlider->labelfont(iMAP->boldFont);
			yLightPositionSlider->labelcolor(FL_WHITE);
			yLightPositionSlider->value(0.0);
			yLightPositionSlider->deactivate();
			yLightPositionSlider->show();

			zLightPositionSlider = new Slider(10,175,w/2-15,20,"Light z");
			zLightPositionSlider->precision(2);
			zLightPositionSlider->labelfont(iMAP->boldFont);
			zLightPositionSlider->labelcolor(FL_WHITE);
			zLightPositionSlider->bounds(-100,100);
			zLightPositionSlider->value(100.0);
			zLightPositionSlider->deactivate();
			zLightPositionSlider->show();

			fogButton = new Fl_Check_Button(w/2+5,105,50,25,"Fog");
			fogButton->labelsize(12);
			fogButton->labelfont(iMAP->boldFont);
			fogButton->labelcolor(FL_WHITE);
			fogButton->value(0);
			fogButton->deactivate();
			fogButton->show();

			fogStartSlider = new Slider(w/2+5,140,w/2-15,20,"Fog Start");
			fogStartSlider->precision(2);
			fogStartSlider->bounds(0,100);
			fogStartSlider->labelfont(iMAP->boldFont);
			fogStartSlider->labelcolor(FL_WHITE);
			fogStartSlider->value(0.0);
			fogStartSlider->deactivate();
			fogStartSlider->show();

			fogEndSlider = new Slider(w/2+5,175,w/2-15,20,"Fog End");
			fogEndSlider->precision(2);
			fogEndSlider->labelfont(iMAP->boldFont);
			fogEndSlider->labelcolor(FL_WHITE);
			fogEndSlider->bounds(0,100);
			fogEndSlider->value(100.0);
			fogEndSlider->deactivate();
			fogEndSlider->show();

			landscapeAxisButton = new Fl_Check_Button(10,200,50,20,"Axes");
			landscapeAxisButton->labelsize(12);
			landscapeAxisButton->labelfont(iMAP->boldFont);
			landscapeAxisButton->labelcolor(FL_WHITE);
			landscapeAxisButton->value(1);
			landscapeAxisButton->callback((Fl_Callback*)axesCallback);
			landscapeAxisButton->deactivate();
			landscapeAxisButton->show();
		}
		landscapeGroup->deactivate();
		landscapeGroup->end();

		tabs->end();

		meshApplied = false;
	}

}

void nullCallback(Fl_Widget*w,void*v) { Fl::redraw(); }

void voronoiMeshWindowCallback(Fl_Widget *w,void*v) {
	file->voronoiMeshOverlay = false;
	file->voronoiMeshGui->resetButton->do_callback();
	file->voronoiMeshGui->hide();
	file->voronoiMeshGui->saveButton->deactivate();
	file->inferred = false;
	file->landscapePlot = false;
	iMAP->customSelectionInferenceGui->activate();
	if (file->voronoiMeshOverlay) {
		if ( file->voronoiMesh->randomizedOptimizationGui != NULL) {
			file->voronoiMesh->randomizedOptimizationGui->hide();
	//		delete file->voronoiMesh->randomizedOptimizationGui;
		}
	}
	delete file->voronoiMeshGui;
	file->voronoiMeshGui = NULL;
	iMAP->overlayAdjustment = true;
}

void clusteringChoiceCallback(Fl_Choice*w,int*v) { Fl::redraw(); }

void colormapChoiceVoronoiCallback(Fl_Choice*w,int*v) {
	char colorString [5];
	sprintf(colorString,"%i",v);
	const int color = atoi(colorString);
	float *rgb = new float[3];

	if (file->voronoiMesh->selectionMode()) {
		for (int e = 0; e < file->voronoiMesh->selection.count; e++) {
			rgb = colormap(color,file->voronoiMesh->getCell(file->voronoiMesh->getClusterIndex(e))->getRandomValue(),file->cMin,file->cMax,false);
			file->voronoiMesh->colorArray[4*e] = rgb[0];
			file->voronoiMesh->colorArray[4*e+1] = rgb[1];
			file->voronoiMesh->colorArray[4*e+2] = rgb[2];
			file->voronoiMesh->colorArray[4*e+3] = 0.5f;
		}
	} else {
		for (int e = 0; e < file->localizationCount; e++) {
			rgb = colormap(color,file->voronoiMesh->getCell(file->voronoiMesh->getClusterIndex(e))->getRandomValue(),file->cMin,file->cMax,false);
			file->voronoiMesh->colorArray[4*e] = rgb[0];
			file->voronoiMesh->colorArray[4*e+1] = rgb[1];
			file->voronoiMesh->colorArray[4*e+2] = rgb[2];
			file->voronoiMesh->colorArray[4*e+3] = 0.5f;
		}
	}
	file->voronoiMeshGui->setColormap(color);
	iMAP->overlayAdjustment = true;

	if (file->landscapePlot) { generateVoronoiLandscape(); }

	Fl::redraw();
}

void colormapChoiceSquareCallback(Fl_Choice*w,int*v) {

	char colorString [5];
	sprintf(colorString,"%i",v);
	const int color = atoi(colorString);

	file->squareMeshGui->setColormap(color);
	iMAP->overlayAdjustment = true;

	if (file->landscapePlot) { generateSquareLandscape(); }

	Fl::redraw();
}

void colormapChoiceQuadTreeCallback(Fl_Choice*w,int*v) {
	char colorString [5];
	sprintf(colorString,"%i",v);
	const int color = atoi(colorString);

	file->treeMeshGui->setColormap(color);
	iMAP->overlayAdjustment = true;

	if (file->landscapePlot) { generateQuadTreeLandscape(); }

	Fl::redraw();
}

void saveVoronoiMeshCallback(Fl_Button*w,int*v) {
	file->voronoiMesh->exportMesh();
}

void VoronoiMeshGui::samplePosterior(int type) {

	switch(type) {
		case 0: // diffusion
		{
			switch(file->optimizationMode) {
			case 0: // D
				{
					posterioriSamples = (int) posterioriSampleNumberSlider->value();
					const double minD = minPosteriorSlider->value();
					const double maxD = maxPosteriorSlider->value();
					posterioriRange = fabs(maxD-minD);
					posterioriIncrement = posterioriRange/(float)posterioriSamples;
					double logMax = -10000000.0;
					if (posterioriSamples > 5) {
						posterioriValues = new double[posterioriSamples];
						double sampler[1];
						// calculate diffusion values for each of sample positions
						for (int g = 0; g < posterioriSamples; g++) {
							sampler[0] = minD + (double)g*posterioriIncrement;
							posterioriValues[g] = -dPosteriorVoronoi(sampler);

							if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
						}
						// calculate log values
						for (int h = 0; h < posterioriSamples; h++) {
							posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
//							fprintf(stderr,"SAMPLE\t%f\t%f\n",minD + (double)h*posterioriIncrement,posterioriValues[h]);
						}
//						fprintf(stderr,"SAMPLE\trange %f\tinc %f\tsamples %i\n",posterioriRange,posterioriIncrement,posterioriSamples);
						drawPosteriorPlot(type,0);
					}
				}
				break;
			case 1: // DF
				{
					posterioriSamples = (int) posterioriSampleNumberSlider->value();
					const double minD = minPosteriorSlider->value();
					const double maxD = maxPosteriorSlider->value();
					posterioriRange = fabs(maxD-minD);
					posterioriIncrement = posterioriRange/(float)posterioriSamples;
					double logMax = -10000000.0;
					if (posterioriSamples > 5) {
						posterioriValues = new double[posterioriSamples];
						double sampler[3];
//						fprintf(stderr,"SAMPLE %i\n",file->voronoiMesh->getCurrentZone());
						sampler[0] = file->voronoiMesh->getForceX(file->voronoiMesh->getCurrentZone());
						sampler[1] = file->voronoiMesh->getForceY(file->voronoiMesh->getCurrentZone());
						// calculate diffusion values for each of sample positions
//						fprintf(stderr,"current zone : %i\n",file->voronoiMesh->getCurrentZone());
						for (int g = 0; g < posterioriSamples; g++) {
							sampler[2] = minD + (double)g*posterioriIncrement;
							posterioriValues[g] = -dfPosteriorVoronoi(sampler);

							if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
						}
						// calculate log values
						for (int h = 0; h < posterioriSamples; h++) {
							posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
//							fprintf(stderr,"SAMPLE\t%f\t%f\n",minD + (double)h*posterioriIncrement,posterioriValues[h]);
						}
//						fprintf(stderr,"SAMPLE\trange %f\tinc %f\tsamples %i\n",posterioriRange,posterioriIncrement,posterioriSamples);
						drawPosteriorPlot(type,0);
					}
				}
				break;
			case 2: // DDr
				{
					posterioriSamples = (int) posterioriSampleNumberSlider->value();
					const double minD = minPosteriorSlider->value();
					const double maxD = maxPosteriorSlider->value();
					posterioriRange = fabs(maxD-minD);
					posterioriIncrement = posterioriRange/(float)posterioriSamples;
					double logMax = -10000000.0;
					if (posterioriSamples > 5) {
						posterioriValues = new double[posterioriSamples];
						double sampler[3];
//						fprintf(stderr,"SAMPLE %i\n",file->voronoiMesh->getCurrentZone());
						sampler[0] = file->voronoiMesh->getForceX(file->voronoiMesh->getCurrentZone());
						sampler[1] = file->voronoiMesh->getForceY(file->voronoiMesh->getCurrentZone());
						// calculate diffusion values for each of sample positions
//						fprintf(stderr,"current zone : %i\n",file->voronoiMesh->getCurrentZone());
						for (int g = 0; g < posterioriSamples; g++) {
							sampler[2] = minD + (double)g*posterioriIncrement;
							posterioriValues[g] = -ddrPosteriorVoronoi(sampler);

							if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
						}
						// calculate log values
						for (int h = 0; h < posterioriSamples; h++) {
							posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
//							fprintf(stderr,"SAMPLE\t%f\t%f\n",minD + (double)h*posterioriIncrement,posterioriValues[h]);
						}
//						fprintf(stderr,"SAMPLE\trange %f\tinc %f\tsamples %i\n",posterioriRange,posterioriIncrement,posterioriSamples);
						drawPosteriorPlot(type,0);
					}
				}
				break;
			case 3: // DV
				{
					posterioriSamples = (int) posterioriSampleNumberSlider->value();
					const double minD = minPosteriorSlider->value();
					const double maxD = maxPosteriorSlider->value();
					posterioriIncrement = (maxD-minD)/posterioriSampleNumberSlider->value();
					posterioriRange = maxD-minD;
					double logMax = -10000000.0;
					if (posterioriSamples > 5) {
						posterioriValues = new double[posterioriSamples];
						// copy DV values to array for posteriori evaluation
						int sampleIndex = file->voronoiMesh->getIdentifier(file->voronoiMesh->getCurrentZone());
						double *sampler = new double[2*file->voronoiMesh->getTotalVariables()];
						// initialize potentials
						int i = 0;
						for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
							if ((file->voronoiMesh->active(a))) {
								sampler[2*i] = file->voronoiMesh->getDiffusion(a);
								sampler[2*i+1] = file->voronoiMesh->getPotential(a);
								i++;
							}
						}
						// calculate diffusion values for each of sample positions
						for (int g = 0; g < posterioriSamples; g++) {
							sampler[2*sampleIndex] = minD + (double)g*posterioriIncrement;

							if (file->voronoiMesh->selectionMode()) { posterioriValues[g] = -dvPosteriorVoronoiSelection(sampler); }
							else { posterioriValues[g] = -dvPosteriorVoronoi(sampler); }

							if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
						}
						// calculate log values
						for (int h = 0; h < posterioriSamples; h++) {
							posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
						}
						drawPosteriorPlot(type,0);
						delete [] sampler;
					}
				}
				break;
			}
			break;
		}
		case 1: // force
		{
			if (file->optimizationMode == 1) { // DF
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				const double minF = minPosteriorSlider->value();
				const double maxF = maxPosteriorSlider->value();
				const double Fx = file->voronoiMesh->getForceX(file->voronoiMesh->getCurrentZone());
				const double Fy = file->voronoiMesh->getForceY(file->voronoiMesh->getCurrentZone());
				posterioriIncrement = (maxF-minF)/posterioriSampleNumberSlider->value();
				posterioriRange = maxF-minF;

				double logMax = -10000000.0;
				if (posterioriSamples > 5) {
					posterioriValues = new double[posterioriSamples];
					double sampler[3];
					sampler[2] = file->voronoiMesh->getDiffusion(file->voronoiMesh->getCurrentZone());
					// calculate diffusion values for each of sample positions
					const double theta = atan2(Fy,Fx);
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[0] = (Fx-fabs(posterioriRange/2.0)*cos(theta)) + (double)g*posterioriIncrement*cos(theta);
						sampler[1] = (Fy-fabs(posterioriRange/2.0)*sin(theta)) + (double)g*posterioriIncrement*sin(theta);
						posterioriValues[g] = -dfPosteriorVoronoi(sampler);
						if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
					}

					// calculate log values
					for (int h = 0; h < posterioriSamples; h++) {
						posterioriValues[h] = exp(posterioriValues[h]-logMax);
					}
					drawPosteriorPlot(type,0);
	//				saveDPosterioriButton->activate();
				}
			}
			else if (file->optimizationMode == 2) { // DDr
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				const double minF = minPosteriorSlider->value();
				const double maxF = maxPosteriorSlider->value();
				const double Fx = file->voronoiMesh->getForceX(file->voronoiMesh->getCurrentZone());
				const double Fy = file->voronoiMesh->getForceY(file->voronoiMesh->getCurrentZone());
				posterioriIncrement = (maxF-minF)/posterioriSampleNumberSlider->value();
				posterioriRange = maxF-minF;

				double logMax = -10000000.0;
				if (posterioriSamples > 5) {
					posterioriValues = new double[posterioriSamples];
					double sampler[3];
					sampler[2] = file->voronoiMesh->getDiffusion(file->voronoiMesh->getCurrentZone());
					// calculate diffusion values for each of sample positions
					const double theta = atan2(Fy,Fx);
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[0] = (Fx-fabs(posterioriRange/2.0)*cos(theta)) + (double)g*posterioriIncrement*cos(theta);
						sampler[1] = (Fy-fabs(posterioriRange/2.0)*sin(theta)) + (double)g*posterioriIncrement*sin(theta);
						posterioriValues[g] = -ddrPosteriorVoronoi(sampler);
						if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
					}

					// calculate log values
					for (int h = 0; h < posterioriSamples; h++) {
						posterioriValues[h] = exp(posterioriValues[h]-logMax);
					}
					drawPosteriorPlot(type,0);
	//				saveDPosterioriButton->activate();
				}
			}
			break;
		}
		case 2: // potential
		{
			if (file->optimizationMode == 3 && this->potentialReferenceButton->value()) {
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				posterioriRange = fabs(maxPosteriorSlider->value());
				double logMax = -10000000.0;
				if (posterioriSamples > 5) {

					const int v0_id = file->voronoiMesh->getIdentifier(file->voronoiMesh->getCurrentZone());
					const double v0 = file->voronoiMesh->getPotential(file->voronoiMesh->getCurrentZone());
					const int v1_id = file->voronoiMesh->getIdentifier(file->voronoiMesh->getCurrentZone2());
					const double v1 = file->voronoiMesh->getPotential(file->voronoiMesh->getCurrentZone2());

					posterioriIncrement = posterioriRange/posterioriSampleNumberSlider->value();

					posterioriValues = new double[2*posterioriSamples];

					// create delta matrix
					for (int s = 0; s < 2*posterioriSamples; s++) {
						posterioriValues[s] = 0.0;
					}

					double *sampler = new double[2*file->voronoiMesh->getTotalVariables()];

					// initialize potentials
					int i = 0;
					for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
						if ((file->voronoiMesh->active(a))) {
							sampler[2*i] = file->voronoiMesh->getDiffusion(a);
							sampler[2*i+1] = file->voronoiMesh->getPotential(a);
							i++;
						}
					}

					double MAP_posterior;
					if (file->voronoiMesh->selectionMode()) { MAP_posterior = -dvPosteriorVoronoiSelection(sampler); }
					else { MAP_posterior = -dvPosteriorVoronoi(sampler); }

					// calculate potential values for each of sample positions
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[2*v0_id+1] = v0 - posterioriRange/2.0 + (double)g*posterioriIncrement;
						for (int h = 0; h < posterioriSamples; h++) {
							sampler[2*v1_id+1] = v1 - posterioriRange/2.0 + (double)h*posterioriIncrement;
							if (file->voronoiMesh->selectionMode()) {
								posterioriValues[ h-g+(posterioriSamples-1) ] += exp(-MAP_posterior-dvPosteriorVoronoiSelection(sampler));
							}
							else {
								posterioriValues[ h-g+(posterioriSamples-1) ] += exp(-MAP_posterior-dvPosteriorVoronoi(sampler));
							}
						}
					}

					// average the posteriori
					for (int e = 0; e < 2*posterioriSamples-1; e++) {
						posterioriValues[e] /= (posterioriSamples*posterioriSamples);
						if (logMax < posterioriValues[e]) { logMax = posterioriValues[e]; }
					}

					// calculate log values
					for (int h = 0; h < 2*posterioriSamples-1; h++) { posterioriValues[h] /= logMax; }

					drawPosteriorPlot(type,0);
					delete [] sampler;
				}
			} else if (file->optimizationMode == 3 && !this->potentialReferenceButton->value()) {
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				const double minV = minPosteriorSlider->value();
				const double maxV = maxPosteriorSlider->value();
				posterioriRange = maxV-minV;
				double logMax = -10000000.0;
				if (posterioriSamples > 5) {

					const int v0_id = file->voronoiMesh->getIdentifier(file->voronoiMesh->getCurrentZone());

					posterioriIncrement = posterioriRange/posterioriSampleNumberSlider->value();

					posterioriValues = new double[posterioriSamples];

					double *sampler = new double[2*file->voronoiMesh->getTotalVariables()];

					// initialize potentials
					int i = 0;
					for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
						if ((file->voronoiMesh->active(a))) {
							sampler[2*i] = file->voronoiMesh->getDiffusion(a);
							sampler[2*i+1] = file->voronoiMesh->getPotential(a);
							i++;
						}
					}

					// calculate diffusion values for each of sample positions
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[2*v0_id+1] = minV + (double)g*posterioriIncrement;
						if (file->voronoiMesh->selectionMode()) { posterioriValues[g] = -dvPosteriorVoronoiSelection(sampler); }
						else { posterioriValues[g] = -dvPosteriorVoronoi(sampler); }
					}

					// average the posteriori
					for (int e = 0; e < posterioriSamples; e++) {
						if (logMax < posterioriValues[e]) { logMax = posterioriValues[e]; }
					}

					// calculate log values
					for (int h = 0; h < posterioriSamples; h++) {
						posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
					}
					drawPosteriorPlot(type,1);
					delete [] sampler;
				}

			}
			break;
		}
	}
}

void VoronoiMeshGui::drawPosteriorPlot(int type,int sup) {

	switch(type) {
		case 0: // diffusion
		{
			const double minD = minPosteriorSlider->value();
			const double maxD = maxPosteriorSlider->value();

			posterioriRange = fabs(maxD - minD);
			posterioriIncrement = posterioriRange/(float)posterioriSamples;

			if (posterioriRange > 0.0 && file->inferred) {
				window->make_current();

				const int w = 320;
				const int h = 195;

	        	fl_line_style(FL_SOLID, 1, 0);

				fl_color(FL_WHITE); fl_rectf(170,35,w,h);
				fl_color(FL_BLACK); fl_rect(170,35,w,h);

				fl_line(170,190,170+w-1,190);
//				fl_font(iMAP->normalFont,10); fl_draw("[um\262/s]",170+w-43,225);

		        int pCurrent,pOld;
		        const double spacing = ((double)(w-1)/(double)posterioriSamples);
		        char xLabel [50];

//		        fprintf(stderr,"\n\n");
		        for (int q = 1; q < posterioriSamples; q++) {
		        	// draw data points and connectors
		        	pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
		        	pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

		        	fl_color(FL_RED);
		        	fl_line_style(FL_SOLID, 2, 0);
					fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
//					fprintf(stderr,"DRAW\t%f\t%f\n",minD + (double)q*posterioriIncrement,posterioriValues[q]);
		        }


		        fl_color(FL_BLACK);
		        fl_font(iMAP->normalFont,10);
		        // draw x-axis labels
		        for (int t = 1; t < 9; t++) {
					sprintf(xLabel,"%.3f",minD+(float)t*((float)posterioriSamples/9.0)*posterioriIncrement);
					fl_draw(90,xLabel,177+t*w/9,227);
		        	fl_line_style(FL_DASH, 1, 0);
					fl_line(174+t*w/9,55,174+t*w/9,190);
		        }
				// draw MAP label
				sprintf(xLabel,"MAP: %.3f [um\262/s]",file->voronoiMesh->getDiffusion(file->voronoiMesh->getCurrentZone()));
		        fl_font(iMAP->normalFont,10);
				fl_draw(xLabel,170+2,45);
			}
			break;
		}
		case 1: // force
		{
			const double minF = minPosteriorSlider->value();
			const double maxF = maxPosteriorSlider->value();

			const double Fx = file->voronoiMesh->getForceX(file->voronoiMesh->getCurrentZone());
			const double Fy = file->voronoiMesh->getForceY(file->voronoiMesh->getCurrentZone());
			posterioriIncrement = (maxF-minF)/posterioriSampleNumberSlider->value();
			posterioriRange = maxF-minF;

			if (posterioriRange > 0.0 && file->inferred) {
				window->make_current();

				const int w = 320;
				const int h = 195;

	        	fl_line_style(FL_SOLID, 1, 0);

				fl_color(FL_WHITE); fl_rectf(170,35,w,h);
				fl_color(FL_BLACK); fl_rect(170,35,w,h);

				fl_line(170,190,170+w-1,190);
//				fl_font(iMAP->normalFont,10); fl_draw("[pN]",170+w-25,225);

		        int pCurrent,pOld;
		        const double spacing = ((double)(w-1)/(double)posterioriSamples);
		        char xLabel [50];

		        for (int q = 1; q < posterioriSamples; q++) {
		        	// draw data points and connectors
		        	pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
		        	pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

		        	fl_color(FL_RED);
		        	fl_line_style(FL_SOLID, 2, 0);
					fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
		        }
		        fl_color(FL_BLACK);
		        fl_font(iMAP->normalFont,10);
		        // draw x-axis labels
		        const double minF2 = sqrt(Fx*Fx+Fy*Fy) - posterioriRange/2.0;

		        for (int t = 1; t < 9; t++) {
					sprintf(xLabel,"%.1f",minF2+t*(posterioriSamples/9.0)*posterioriIncrement);
					fl_draw(90,xLabel,177+t*w/9,222);
		        	fl_line_style(FL_DASH, 1, 0);
					fl_line(174+t*w/9,55,174+t*w/9,190);
		        }
				// draw MAP label
				sprintf(xLabel,"MAP: %.1f [pN]",file->voronoiMesh->getForce(file->voronoiMesh->getCurrentZone()));
		        fl_font(iMAP->normalFont,10);
				fl_draw(xLabel,170+2,45);
			}

			break;
		}
		case 2: // potential
		{
			switch(sup) {
				case 0:
				{
					const double maxSample = fabs(maxPosteriorSlider->value())/2.0;

					if (posterioriRange > 0.0 && file->inferred) {
						window->make_current();

						const int w = 320;
						const int h = 195;

						fl_line_style(FL_SOLID, 1, 0);

						fl_color(FL_WHITE); fl_rectf(170,35,w,h);
						fl_color(FL_BLACK); fl_rect(170,35,w,h);

						fl_line(170,190,170+w-1,190);
//						fl_font(iMAP->normalFont,10); fl_draw("[kT]",170+w-25,225);

						const double deltaV = file->voronoiMesh->getPotential(file->voronoiMesh->getCurrentZone())
											  -file->voronoiMesh->getPotential(file->voronoiMesh->getCurrentZone2());

						int pCurrent,pOld;
						const double spacing = ((double)(w-1)/(double)(2*posterioriSamples));
						char xLabel [50];

						for (int q = 1; q < 2*posterioriSamples; q++) {
							// draw data points and connectors
							pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
							pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

							fl_color(FL_RED);
							fl_line_style(FL_SOLID, 2, 0);
							fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
						}
						fl_color(FL_BLACK);
						fl_font(iMAP->normalFont,10);
						// draw x-axis labels
						for (int t = 1; t < 9; t++) {
							sprintf(xLabel,"%.2f",deltaV-maxSample+(float)t/9.0*posterioriRange);
							fl_draw(90,xLabel,177+t*w/9,227);
							fl_line_style(FL_DASH, 1, 0);
							fl_line(174+t*w/9,55,174+t*w/9,190);
						}
						// draw MAP label
						sprintf(xLabel,"MAP: %.2f [kT]",deltaV);
						fl_font(iMAP->normalFont,10);
						fl_draw(xLabel,170+2,45);
					}

					break;
				}
				case 1:
				{
					const double minSample = minPosteriorSlider->value();
					const double maxSample = maxPosteriorSlider->value();

					posterioriRange = maxSample - minSample;

					if (posterioriRange > 0.0 && file->inferred) {
						window->make_current();

						const int w = 320;
						const int h = 195;

			        	fl_line_style(FL_SOLID, 1, 0);

						fl_color(FL_WHITE); fl_rectf(170,35,w,h);
						fl_color(FL_BLACK); fl_rect(170,35,w,h);

						fl_line(170,190,170+w-1,190);
//						fl_font(iMAP->normalFont,10); fl_draw("[kT]",170+w-25,225);

				        int pCurrent,pOld;
				        const double spacing = ((double)(w-1)/(double)posterioriSamples);
				        char xLabel [80];

				        for (int q = 1; q < posterioriSamples; q++) {
				        	// draw data points and connectors
				        	pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
				        	pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

				        	fl_color(FL_RED);
				        	fl_line_style(FL_SOLID, 2, 0);
							fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
				        }
				        fl_color(FL_BLACK);
				        fl_font(iMAP->normalFont,10);
				        // draw x-axis labels
				        for (int t = 1; t < 9; t++) {
							sprintf(xLabel,"%.2f",minSample+t*(posterioriSamples/9.0)*posterioriIncrement);
							fl_draw(90,xLabel,177+t*w/9,227);
				        	fl_line_style(FL_DASH, 1, 0);
							fl_line(174+t*w/9,55,174+t*w/9,190);
				        }
						// draw MAP label
						sprintf(xLabel,"MAP: %.2f [kT] (no reference)",file->voronoiMesh->getPotential(file->voronoiMesh->getCurrentZone()));
				        fl_font(iMAP->normalFont,10);
						fl_draw(xLabel,170+2,45);
					}

					break;
				}
			}
			break;
		}
	}

}

void VoronoiMesh::savePosterior() {
	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	switch(file->voronoiMeshGui->getPosteriorType()) {
	case 0: // diffusion
		native->title("Save Diffusion Posterior Data");
		native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
		native->filter("Diffusion Posterior File\t*.dpost\n");
	    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
	    native->preset_file(file->fileName);
		if (native->show() == 0) {
			// store filename
			sprintf(nativeFilename,"%s%s",native->filename(),".dpost");
			FILE * writeFile;

			writeFile = fopen (nativeFilename,"wb");

			// write header
			fprintf(writeFile,"File Name: %s\n\n"
							  "Number of Points: %i\n"
							  "Zone Centre: [%.3f,%.3f]\n"
							  "D_MAP: %f [um2/s]\n\n",
							  file->fileName,
							  getCount(getCurrentZone()),
							  getCentroidX(getCurrentZone()),
							  getCentroidY(getCurrentZone()),
							  getDiffusion(getCurrentZone()));

			for (int f = 0; f < file->voronoiMeshGui->posterioriSamples; f++) {
				fprintf(writeFile,"%f\t%f\n",file->voronoiMeshGui->minPosteriorSlider->value() + (double)f*file->voronoiMeshGui->posterioriIncrement,file->voronoiMeshGui->posterioriValues[f]);
			}
			fclose(writeFile);
		}
		break;
	case 1: // force
		native->title("Save Force Posterior Data");
		native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
		native->filter("Force Posterior File\t*.fpost\n");
	    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
	    native->preset_file(file->fileName);
		if (native->show() == 0) {
			// store filename
			sprintf(nativeFilename,"%s%s",native->filename(),".fpost");
			FILE * writeFile;

			writeFile = fopen (nativeFilename,"wb");

			// write header
			fprintf(writeFile,"File Name: %s\n\n"
							  "Number of Points: %i\n"
							  "Zone Centre: [%.3f,%.3f]\n"
							  "F_MAP: %f [pN]\n\n",
							  file->fileName,
							  getCount(getCurrentZone()),
							  getCentroidX(getCurrentZone()),
							  getCentroidY(getCurrentZone()),
							  getForce(getCurrentZone()));

			const double minF2 = getForce(getCurrentZone()) - file->voronoiMeshGui->posterioriRange/2.0;

			for (int f = 0; f < file->voronoiMeshGui->posterioriSamples; f++) {
				fprintf(writeFile,"%f\t%f\n",minF2 + (double)f*file->voronoiMeshGui->posterioriIncrement,file->voronoiMeshGui->posterioriValues[f]);
			}
			fclose(writeFile);
		}
		break;
	case 2: // potential
		native->title("Save Potential Posterior Data");
		native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
		native->filter("Potential Posterior File\t*.vpost\n");
	    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
	    native->preset_file(file->fileName);
		if (native->show() == 0) {
			// store filename
			sprintf(nativeFilename,"%s%s",native->filename(),".vpost");
			FILE * writeFile;

			writeFile = fopen (nativeFilename,"wb");

			if (file->voronoiMeshGui->potentialReferenceButton->value()) {
				// write header
				const double deltaV = file->voronoiMesh->getPotential(file->voronoiMesh->getCurrentZone2())
									  -file->voronoiMesh->getPotential(file->voronoiMesh->getCurrentZone());
				fprintf(writeFile,"File Name: %s\n\n"
								  "Number of Points: %i\n"
						  	  	  "Zone 1 Centre: [%.3f,%.3f]\n"
						  	  	  "Zone 2 Centre: [%.3f,%.3f]\n"
								  "V_MAP (reference): %f [kT]\n\n",
								  file->fileName,
								  getCount(file->voronoiMesh->getCurrentZone()),
								  getCentroidX(file->voronoiMesh->getCurrentZone()),getCentroidY(file->voronoiMesh->getCurrentZone()),
								  getCentroidX(file->voronoiMesh->getCurrentZone2()),getCentroidY(file->voronoiMesh->getCurrentZone2()),
								  deltaV);

				for (int f = 0; f < 2*file->voronoiMeshGui->posterioriSamples; f++) {
					fprintf(writeFile,"%f\t%f\n",deltaV-file->voronoiMeshGui->maxPosteriorSlider->value()/2.0 + (float)f/(2*file->voronoiMeshGui->posterioriSamples)*file->voronoiMeshGui->posterioriRange,file->voronoiMeshGui->posterioriValues[f]);
				}
			} else {
				// write header
				fprintf(writeFile,"File Name: %s\n\n"
								  "Number of Points: %i\n"
								  "Zone Centre: [%.3f,%.3f]\n"
								  "V_MAP: %f [kT]\n\n",
								  file->fileName,
								  getCount(getCurrentZone()),
								  getCentroidX(getCurrentZone()),
								  getCentroidY(getCurrentZone()),
								  getPotential(getCurrentZone()));

				for (int f = 0; f < file->voronoiMeshGui->posterioriSamples; f++) {
					fprintf(writeFile,"%f\t%f\n",file->voronoiMeshGui->minPosteriorSlider->value() + (double)f*file->voronoiMeshGui->posterioriIncrement,file->voronoiMeshGui->posterioriValues[f]);
				}
			}
			fclose(writeFile);
		}
		break;
	}
}

void voronoiLandscapeCallback(Fl_Check_Button*w,void*) {
	if (w->value()) {
		file->voronoiMeshGui->scaleLandscapeSlider->activate();
		file->voronoiMeshGui->landscapeAlphaSlider->activate();
		file->voronoiMeshGui->xLightPositionSlider->activate();
		file->voronoiMeshGui->yLightPositionSlider->activate();
		file->voronoiMeshGui->zLightPositionSlider->activate();
		file->voronoiMeshGui->fogButton->activate();
		file->voronoiMeshGui->fogStartSlider->activate();
		file->voronoiMeshGui->fogEndSlider->activate();
		file->voronoiMeshGui->landscapeAxisButton->activate();

		generateVoronoiLandscape();
		file->landscapePlot = true;
		file->orthoLimit = 0.6*file->maxRange;
		file->xTranslate = file->yTranslate = file->zTranslate = 0.0;
		file->xRotate = file->yRotate = file->zRotate = 0;
		iMAP->fovy = 2*atan(file->maxRange/2/15)*180/PI;
	}
	else {
		file->voronoiMeshGui->scaleLandscapeSlider->deactivate();
		file->voronoiMeshGui->landscapeAlphaSlider->deactivate();
		file->voronoiMeshGui->xLightPositionSlider->deactivate();
		file->voronoiMeshGui->yLightPositionSlider->deactivate();
		file->voronoiMeshGui->zLightPositionSlider->deactivate();
		file->voronoiMeshGui->fogButton->deactivate();
		file->voronoiMeshGui->fogStartSlider->deactivate();
		file->voronoiMeshGui->fogEndSlider->deactivate();
		file->voronoiMeshGui->landscapeAxisButton->deactivate();

		file->landscapePlot = false;
	}
}

void updateVoronoiLandscapeCallback(Fl_Widget*w,void*) {
	generateVoronoiLandscape();
}

void alphaVoronoiLandscapeCallback(Fl_Widget*w,void*v) {
	const float alpha = file->voronoiMeshGui->landscapeAlphaSlider->value();
	for (int i = 0; i < file->voronoiMesh->landscapeVertices/3; i++) {
		if (file->voronoiMesh->landscapeColors[12*i+3] > 0.0) {
			file->voronoiMesh->landscapeColors[12*i+3] = alpha;
		}
		if (file->voronoiMesh->landscapeColors[12*i+7] > 0.0) {
			file->voronoiMesh->landscapeColors[12*i+7] = alpha;
		}
		if (file->voronoiMesh->landscapeColors[12*i+11] > 0.0) {
			file->voronoiMesh->landscapeColors[12*i+11] = alpha;
		}
	}
}

void setLandscapeVoronoi() {

	// assign data
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			if (file->voronoiMeshGui->overlayDiffusionButton->value()) {
				file->voronoiMesh->getCell(a)->setLandscape((float) ( file->voronoiMesh->getCell(a)->getDiffusion() - file->voronoiMesh->getDiffusionMin() ) /
															  (float) ( file->voronoiMesh->getDiffusionMax() - file->voronoiMesh->getDiffusionMin() ));
			}
			else if (file->voronoiMeshGui->overlayPotentialButton->value()) {
				file->voronoiMesh->getCell(a)->setLandscape((float) ( file->voronoiMesh->getCell(a)->getPotential() - file->voronoiMesh->getPotentialMin() ) /
															  (float) ( file->voronoiMesh->getPotentialMax() - file->voronoiMesh->getPotentialMin() ));
			}
			else if (file->voronoiMeshGui->overlayForceMagnitudeButton->value()) {
				file->voronoiMesh->getCell(a)->setLandscape((float) ( file->voronoiMesh->getCell(a)->getForceMagnitude() - file->voronoiMesh->getForceMin() ) /
															  (float) ( file->voronoiMesh->getForceMax() - file->voronoiMesh->getForceMin() ));
			}
			else if (file->voronoiMeshGui->overlayPointNumberButton->value()) {
				file->voronoiMesh->getCell(a)->setLandscape((float) ( file->voronoiMesh->getCell(a)->getCount() - file->voronoiMesh->getMinCell()->getCount() ) /
															  (float) ( file->voronoiMesh->getMaxCell()->getCount() - file->voronoiMesh->getMinCell()->getCount() ));
			}
			else { file->voronoiMesh->getCell(a)->setLandscape(0.0f); }
		}
		else { file->voronoiMesh->getCell(a)->setLandscape(0.0f); }
	}

}

void generateVoronoiLandscape() {

	float xMin,xMax,yMin,yMax;
	if (file->voronoiMesh->selectionMode()) {
		xMin = (float)file->voronoiMesh->selection.xMin;
		xMax = (float)file->voronoiMesh->selection.xMax;
		yMin = (float)file->voronoiMesh->selection.yMin;
		yMax = (float)file->voronoiMesh->selection.yMax;
	} else {
		xMin = (float)file->xMin;
		xMax = (float)file->xMax;
		yMin = (float)file->yMin;
		yMax = (float)file->yMax;
	}

	setLandscapeVoronoi();

	int nVertices = 0;
	// count number of vertices needed (3*nVertices for each cell)
	for (int e = 0; e < file->voronoiMesh->getNumberOfClusters(); e++) {
		nVertices += 3*file->voronoiMesh->getCell(e)->nVertices;
	}

	file->voronoiMesh->landscapeVertices = nVertices;

	// allocate memory for rendering variables
	if (file->voronoiMesh->landscapeTriangles != NULL) {
		delete [] file->voronoiMesh->landscapeTriangles;
		delete [] file->voronoiMesh->landscapeNormals;
		delete [] file->voronoiMesh->landscapeColors;
	}
	file->voronoiMesh->landscapeTriangles = new float[3*nVertices];
	file->voronoiMesh->landscapeNormals = new float[3*nVertices];
	file->voronoiMesh->landscapeColors = new float[4*nVertices];

	const float scale = file->voronoiMeshGui->scaleLandscapeSlider->value();
	const float alpha = file->voronoiMeshGui->landscapeAlphaSlider->value();

	// populate vertices, normals, colours
	int i = 0;
	int jj = 0;
	int nVerts;
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		nVerts = file->voronoiMesh->getCell(a)->nVertices;
		for (int k = 0; k < nVerts; k++) {
			// corner 1
			file->voronoiMesh->landscapeTriangles[9*i+0] = file->voronoiMesh->getCell(a)->vertices[2*k+0];
			file->voronoiMesh->landscapeTriangles[9*i+1] = file->voronoiMesh->getCell(a)->vertices[2*k+1];
			file->voronoiMesh->landscapeTriangles[9*i+2] = scale*file->voronoiMesh->getCell(a)->getLandscape();
//			if (file->voronoiMesh->getCell(a)->distanceToCentroid[k] < file->voronoiMesh->getPerimeterMax()/4.0) {
				file->voronoiMesh->landscapeColors[12*jj+3]	= alpha;
//			} else {
//				file->voronoiMesh->landscapeColors[12*jj+3] = 0.0;
//			}
			// corner 2
			if (k+1 == nVerts) {
				file->voronoiMesh->landscapeTriangles[9*i+3] = file->voronoiMesh->getCell(a)->vertices[0];
				file->voronoiMesh->landscapeTriangles[9*i+4] = file->voronoiMesh->getCell(a)->vertices[1];
				file->voronoiMesh->landscapeTriangles[9*i+5] = scale*file->voronoiMesh->getCell(a)->getLandscape();
//				if (file->voronoiMesh->getCell(a)->distanceToCentroid[0] < file->voronoiMesh->getPerimeterMax()/4.0) {
					file->voronoiMesh->landscapeColors[12*jj+7]	= alpha;
//				} else {
//					file->voronoiMesh->landscapeColors[12*jj+7] = 0.0;
//				}
			} else {
				file->voronoiMesh->landscapeTriangles[9*i+3] = file->voronoiMesh->getCell(a)->vertices[2*(k+1)+0];
				file->voronoiMesh->landscapeTriangles[9*i+4] = file->voronoiMesh->getCell(a)->vertices[2*(k+1)+1];
				file->voronoiMesh->landscapeTriangles[9*i+5] = scale*file->voronoiMesh->getCell(a)->getLandscape();
//				if (file->voronoiMesh->getCell(a)->distanceToCentroid[k+1] < file->voronoiMesh->getPerimeterMax()/4.0) {
					file->voronoiMesh->landscapeColors[12*jj+7]	= alpha;
//				} else {
//					file->voronoiMesh->landscapeColors[12*jj+7] = 0.0;
//				}
			}
			// centroid
			file->voronoiMesh->landscapeTriangles[9*i+6] = file->voronoiMesh->getCell(a)->getXCentre();
			file->voronoiMesh->landscapeTriangles[9*i+7] = file->voronoiMesh->getCell(a)->getYCentre();
			file->voronoiMesh->landscapeTriangles[9*i+8] = scale*file->voronoiMesh->getCell(a)->getLandscape();
			if (file->voronoiMesh->active(a)) { file->voronoiMesh->landscapeColors[12*jj+11] = alpha; }
			else { file->voronoiMesh->landscapeColors[12*jj+11] = 0.0; }
			i++;
//			file->voronoiMesh->landscapeColors[12*jj+3] = (1.0-normArea);
//			file->voronoiMesh->landscapeColors[12*jj+7] = (1.0-normArea);
//			file->voronoiMesh->landscapeColors[12*jj+11] = (1.0-normArea);
			jj++;
		}
	}

	// pass to average z coordinates at cell boundaries
	float zAvg = 0.0;
	int num = 0;
	float *zCoordsTemp = new float[nVertices];
	i = 0;
	for (int k = 0; k < nVertices; k++) {
		zAvg = 0.0;
		num = 0;
		for (int j = 0; j < nVertices; j++) {
			const float x = file->voronoiMesh->landscapeTriangles[3*j+0];
			const float y = file->voronoiMesh->landscapeTriangles[3*j+1];
			if (fabs(file->voronoiMesh->landscapeTriangles[3*k+0]-x) < 0.001 && fabs(file->voronoiMesh->landscapeTriangles[3*k+1]-y) < 0.001) {
				zAvg += file->voronoiMesh->landscapeTriangles[3*j+2];
				num++;
			}
		}
		// store coordinates in temporary vector
		zCoordsTemp[k] = zAvg/(float)num;
	}
	// copy stored coordinates to actual rendering vector
	for (int k = 0; k < nVertices; k++) {
		file->voronoiMesh->landscapeTriangles[3*k+2] = zCoordsTemp[k];
	}
	delete [] zCoordsTemp;

	// smooth centroid coordinate (to suppress pyramid effects) by averaging z values of vertices
	i = 0;
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		zAvg = 0.0;
		nVerts = file->voronoiMesh->getCell(a)->nVertices;
		for (int k = 0; k < nVerts; k++) {
			// corner 1
			zAvg += file->voronoiMesh->landscapeTriangles[9*i+2];
			// corner 2
			zAvg += file->voronoiMesh->landscapeTriangles[9*i+5];
			// centroid
			i++;
		}
		i -= nVerts;
		for (int k = 0; k < nVerts; k++) {
			file->voronoiMesh->landscapeTriangles[9*i+8] = zAvg/(float)(2*nVerts);
			i++;
		}
	}

	// calculate normals
	float v02[3];
	float v10[3];
	float v21[3];
	i = 0;
	for (int j = 0; j < nVertices/3; j++) {

		// calculate normals
		v02[0] = file->voronoiMesh->landscapeTriangles[9*j+6] - file->voronoiMesh->landscapeTriangles[9*j+0];
		v02[1] = file->voronoiMesh->landscapeTriangles[9*j+7] - file->voronoiMesh->landscapeTriangles[9*j+1];
		v02[2] = file->voronoiMesh->landscapeTriangles[9*j+8] - file->voronoiMesh->landscapeTriangles[9*j+2];

		v10[0] = file->voronoiMesh->landscapeTriangles[9*j+0] - file->voronoiMesh->landscapeTriangles[9*j+3];
		v10[1] = file->voronoiMesh->landscapeTriangles[9*j+1] - file->voronoiMesh->landscapeTriangles[9*j+4];
		v10[2] = file->voronoiMesh->landscapeTriangles[9*j+2] - file->voronoiMesh->landscapeTriangles[9*j+5];

		v21[0] = file->voronoiMesh->landscapeTriangles[9*j+3] - file->voronoiMesh->landscapeTriangles[9*j+6];
		v21[1] = file->voronoiMesh->landscapeTriangles[9*j+4] - file->voronoiMesh->landscapeTriangles[9*j+7];
		v21[2] = file->voronoiMesh->landscapeTriangles[9*j+5] - file->voronoiMesh->landscapeTriangles[9*j+8];

		// vertex 0
		crossProduct(&file->voronoiMesh->landscapeNormals[9*j+0],v02[0],v02[1],v02[2],v10[0],v10[1],v10[2]);
		// vertex 1
		crossProduct(&file->voronoiMesh->landscapeNormals[9*j+3],v10[0],v10[1],v10[2],v21[0],v21[1],v21[2]);
		// vertex 2
		crossProduct(&file->voronoiMesh->landscapeNormals[9*j+6],v10[0],v10[1],v10[2],v21[0],v21[1],v21[2]);

		// assign colours
		// vertex 1
		if (file->voronoiMesh->landscapeTriangles[9*j+2] > 0.0) {
			colormap(&file->voronoiMesh->landscapeColors[12*j+0],file->voronoiMeshGui->getColormap(),file->voronoiMesh->landscapeTriangles[9*j+2]/scale,file->cMin,file->cMax,file->flipColormap);
			file->voronoiMesh->landscapeColors[12*j+3] *= alpha;
//			file->voronoiMesh->landscapeColors[12*j+3] *= 3.0*file->voronoiMesh->landscapeTriangles[9*j+2]/scale;
		} else {
			file->voronoiMesh->landscapeColors[12*j+0]=file->voronoiMesh->landscapeColors[12*j+1]=file->voronoiMesh->landscapeColors[12*j+2]=file->voronoiMesh->landscapeColors[12*j+3]=0.0;
		}
		// vertex 2
		if (file->voronoiMesh->landscapeTriangles[9*j+5] > 0.0) {
			colormap(&file->voronoiMesh->landscapeColors[12*j+4],file->voronoiMeshGui->getColormap(),file->voronoiMesh->landscapeTriangles[9*j+5]/scale,file->cMin,file->cMax,file->flipColormap);
			file->voronoiMesh->landscapeColors[12*j+7] *= alpha;
//			file->voronoiMesh->landscapeColors[12*j+7] *= 3.0*file->voronoiMesh->landscapeTriangles[9*j+5]/scale;
		} else {
			file->voronoiMesh->landscapeColors[12*j+4]=file->voronoiMesh->landscapeColors[12*j+5]=file->voronoiMesh->landscapeColors[12*j+6]=file->voronoiMesh->landscapeColors[12*j+7]=0.0;
		}
		// vertex 3
		if (file->voronoiMesh->landscapeTriangles[9*j+8] > 0.0) {
			colormap(&file->voronoiMesh->landscapeColors[12*j+8],file->voronoiMeshGui->getColormap(),file->voronoiMesh->landscapeTriangles[9*j+8]/scale,file->cMin,file->cMax,file->flipColormap);
			file->voronoiMesh->landscapeColors[12*j+11] *= alpha;
//			file->voronoiMesh->landscapeColors[12*j+11] *= 3.0*file->voronoiMesh->landscapeTriangles[9*j+8]/scale;
		} else {
			file->voronoiMesh->landscapeColors[12*j+8]=file->voronoiMesh->landscapeColors[12*j+9]=file->voronoiMesh->landscapeColors[12*j+10]=file->voronoiMesh->landscapeColors[12*j+11]=0.0;
		}
	}

	// average normals
	// pass to average z coordinates at cell boundaries
	float norm[3];
	num = 0;
	i = 0;
	for (int k = 0; k < nVertices; k++) {
		norm[0] = norm[1] = norm[2] = 0.0;
		num = 0;
		for (int j = 0; j < nVertices; j++) {
			const float x = file->voronoiMesh->landscapeTriangles[3*j+0];
			const float y = file->voronoiMesh->landscapeTriangles[3*j+1];
			if (fabs(file->voronoiMesh->landscapeTriangles[3*k+0]-x) < 0.001 && fabs(file->voronoiMesh->landscapeTriangles[3*k+1]-y) < 0.001) {
				norm[0] += file->voronoiMesh->landscapeNormals[3*j+0];
				norm[1] += file->voronoiMesh->landscapeNormals[3*j+1];
				norm[2] += file->voronoiMesh->landscapeNormals[3*j+2];
				num++;
			}
		}
		// update rendering normal
		file->voronoiMesh->landscapeNormals[3*k+0] = norm[0]/(float)num;
		file->voronoiMesh->landscapeNormals[3*k+1] = norm[1]/(float)num;
		file->voronoiMesh->landscapeNormals[3*k+2] = norm[2]/(float)num;
	}

}

void TreeMeshGui::samplePosterior(int type) {

	switch(type) {
		case 0: // diffusion
		{
			switch(file->optimizationMode) {
			case 0: // D
				{
					posterioriSamples = (int) posterioriSampleNumberSlider->value();
					const double minD = minPosteriorSlider->value();
					const double maxD = maxPosteriorSlider->value();
					posterioriIncrement = (maxD-minD)/posterioriSampleNumberSlider->value();
					posterioriRange = maxD-minD;
					double logMax = -10000000.0;
					file->treeMesh->setCurrentZone(file->treeMesh->selectedQuadLeaf);
					if (posterioriSamples > 5) {
						posterioriValues = new double[posterioriSamples];
						double sampler[1];
						// calculate diffusion values for each of sample positions
						for (int g = 0; g < posterioriSamples; g++) {
							sampler[0] = minD + (double)g*posterioriIncrement;
							posterioriValues[g] = -dPosteriorQuadTree(sampler);
							if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
						}
						// calculate log values
						for (int h = 0; h < posterioriSamples; h++) {
							posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
//								fprintf(stderr,"SAMPLE\t%f\t%f\n",minD + (double)h*posterioriIncrement,posterioriValues[h]);
						}
						drawPosteriorPlot(type,0);
		//				saveDPosterioriButton->activate();
					}
				}
				break;
			case 1: // DF
				{
					posterioriSamples = (int) posterioriSampleNumberSlider->value();
					const double minD = minPosteriorSlider->value();
					const double maxD = maxPosteriorSlider->value();
					posterioriIncrement = (maxD-minD)/posterioriSampleNumberSlider->value();
					posterioriRange = maxD-minD;
					double logMax = -10000000.0;
					file->treeMesh->setCurrentZone(file->treeMesh->selectedQuadLeaf);
					if (posterioriSamples > 5) {
						posterioriValues = new double[posterioriSamples];
						double sampler[3];
						sampler[0] = file->treeMesh->getCurrentQuadZone()->xForce;
						sampler[1] = file->treeMesh->getCurrentQuadZone()->yForce;
						// calculate diffusion values for each of sample positions
						for (int g = 0; g < posterioriSamples; g++) {
							sampler[2] = minD + (double)g*posterioriIncrement;
							posterioriValues[g] = -dfPosteriorQuadTree(sampler);
							if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
						}
						// calculate log values
						for (int h = 0; h < posterioriSamples; h++) {
							posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
//								fprintf(stderr,"SAMPLE\t%f\t%f\n",minD + (double)h*posterioriIncrement,posterioriValues[h]);
						}
						drawPosteriorPlot(type,0);
		//				saveDPosterioriButton->activate();
					}
				}
				break;
			case 2: // DDr
				{
					posterioriSamples = (int) posterioriSampleNumberSlider->value();
					const double minD = minPosteriorSlider->value();
					const double maxD = maxPosteriorSlider->value();
					posterioriIncrement = (maxD-minD)/posterioriSampleNumberSlider->value();
					posterioriRange = maxD-minD;
					double logMax = -10000000.0;
					file->treeMesh->setCurrentZone(file->treeMesh->selectedQuadLeaf);
					if (posterioriSamples > 5) {
						posterioriValues = new double[posterioriSamples];
						double sampler[3];
						sampler[0] = file->treeMesh->getCurrentQuadZone()->xForce;
						sampler[1] = file->treeMesh->getCurrentQuadZone()->yForce;
						// calculate diffusion values for each of sample positions
						for (int g = 0; g < posterioriSamples; g++) {
							sampler[2] = minD + (double)g*posterioriIncrement;
							posterioriValues[g] = -ddrPosteriorQuadTree(sampler);
							if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
						}
						// calculate log values
						for (int h = 0; h < posterioriSamples; h++) {
							posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
//								fprintf(stderr,"SAMPLE\t%f\t%f\n",minD + (double)h*posterioriIncrement,posterioriValues[h]);
						}
						drawPosteriorPlot(type,0);
		//				saveDPosterioriButton->activate();
					}
				}
				break;
				case 3: // DV
					{
						posterioriSamples = (int) posterioriSampleNumberSlider->value();
						const double minD = minPosteriorSlider->value();
						const double maxD = maxPosteriorSlider->value();
						posterioriIncrement = (maxD-minD)/posterioriSampleNumberSlider->value();
						posterioriRange = maxD-minD;
						double logMax = -10000000.0;
						if (posterioriSamples > 5) {
							posterioriValues = new double[posterioriSamples];
							// copy DV values to array for posteriori evaluation
							int sampleIndex = file->treeMesh->selectedQuadLeaf->getIdentifier();
							double *sampler = new double[2*file->treeMesh->getTotalVariables()];
							// initialize potentials
							for (int i = 0; i < file->treeMesh->getTotalVariables(); i++) {
								file->treeMesh->getTree(file->treeMesh->quadTree,i);
								sampler[2*i] = file->treeMesh->identifiedQuadLeaf->getDiffusion();
								sampler[2*i+1] = file->treeMesh->identifiedQuadLeaf->getPotential();
							}
							// calculate diffusion values for each of sample positions
							for (int g = 0; g < posterioriSamples; g++) {
								sampler[2*sampleIndex] = minD + (double)g*posterioriIncrement;

								if (file->treeMesh->selectionMode()) { posterioriValues[g] = -dvPosteriorQuadTreeSelection(sampler); }
								else { posterioriValues[g] = -dvPosteriorQuadTree(sampler); }

								if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
							}
							// calculate log values
							for (int h = 0; h < posterioriSamples; h++) {
								posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
							}
							drawPosteriorPlot(type,0);
							delete [] sampler;
						}
					}
					break;
				}
			}
			break;
		case 1: // force
		{
			if (file->optimizationMode == 1) { // DF
				file->treeMesh->setCurrentZone(file->treeMesh->selectedQuadLeaf);

				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				const double minF = minPosteriorSlider->value();
				const double maxF = maxPosteriorSlider->value();
				const double Fx = file->treeMesh->getCurrentQuadZone()->xForce;
				const double Fy = file->treeMesh->getCurrentQuadZone()->yForce;

				posterioriIncrement = (maxF-minF)/posterioriSampleNumberSlider->value();
				posterioriRange = maxF-minF;

				double logMax = -10000000.0;
				if (posterioriSamples > 5) {
					posterioriValues = new double[posterioriSamples];
					double sampler[3];
					sampler[2] = file->treeMesh->selectedQuadLeaf->getDiffusion();
					// calculate diffusion values for each of sample positions
					const double theta = atan2(Fy,Fx);
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[0] = (Fx-fabs(posterioriRange/2.0)*cos(theta)) + (double)g*posterioriIncrement*cos(theta);
						sampler[1] = (Fy-fabs(posterioriRange/2.0)*sin(theta)) + (double)g*posterioriIncrement*sin(theta);
						posterioriValues[g] = -dfPosteriorQuadTree(sampler);
						if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
					}

					// calculate log values
					for (int h = 0; h < posterioriSamples; h++) {
						posterioriValues[h] = exp(posterioriValues[h]-logMax);
					}
					drawPosteriorPlot(type,0);
	//				saveDPosterioriButton->activate();
				}
			}
			else if (file->optimizationMode == 2) { // DF
				file->treeMesh->setCurrentZone(file->treeMesh->selectedQuadLeaf);

				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				const double minF = minPosteriorSlider->value();
				const double maxF = maxPosteriorSlider->value();
				const double Fx = file->treeMesh->getCurrentQuadZone()->xForce;
				const double Fy = file->treeMesh->getCurrentQuadZone()->yForce;

				posterioriIncrement = (maxF-minF)/posterioriSampleNumberSlider->value();
				posterioriRange = maxF-minF;

				double logMax = -10000000.0;
				if (posterioriSamples > 5) {
					posterioriValues = new double[posterioriSamples];
					double sampler[3];
					sampler[2] = file->treeMesh->selectedQuadLeaf->getDiffusion();
					// calculate diffusion values for each of sample positions
					const double theta = atan2(Fy,Fx);
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[0] = (Fx-fabs(posterioriRange/2.0)*cos(theta)) + (double)g*posterioriIncrement*cos(theta);
						sampler[1] = (Fy-fabs(posterioriRange/2.0)*sin(theta)) + (double)g*posterioriIncrement*sin(theta);
						posterioriValues[g] = -dfPosteriorQuadTree(sampler);
						if (logMax < posterioriValues[g]) { logMax = posterioriValues[g]; }
					}

					// calculate log values
					for (int h = 0; h < posterioriSamples; h++) {
						posterioriValues[h] = exp(posterioriValues[h]-logMax);
					}
					drawPosteriorPlot(type,0);
	//				saveDPosterioriButton->activate();
				}
			}
		}
			break;
		case 2: // potential
		{
			if (file->optimizationMode == 3 && this->potentialReferenceButton->value()) {
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				posterioriRange = fabs(maxPosteriorSlider->value());
				double logMax = -10000000.0;
				if (posterioriSamples > 5) {

					const int v0_id = file->treeMesh->selectedQuadLeaf->getIdentifier();
					const double v0 = file->treeMesh->selectedQuadLeaf->getPotential();
					const int v1_id = file->treeMesh->selectedQuadLeaf2->getIdentifier();
					const double v1 = file->treeMesh->selectedQuadLeaf2->getPotential();

					posterioriIncrement = posterioriRange/posterioriSampleNumberSlider->value();

					posterioriValues = new double[2*posterioriSamples];

					// create delta matrix
					for (int s = 0; s < 2*posterioriSamples; s++) {
						posterioriValues[s] = 0.0;
					}

					double *sampler = new double[2*file->treeMesh->getTotalVariables()];

					// initialize potentials
					for (int i = 0; i < file->treeMesh->getTotalVariables(); i++) {
						file->treeMesh->getTree(file->treeMesh->quadTree,i);
						sampler[2*i] = file->treeMesh->identifiedQuadLeaf->getDiffusion();
						sampler[2*i+1] = file->treeMesh->identifiedQuadLeaf->getPotential();
					}

					double MAP_posterior;
					if (file->treeMesh->selectionMode()) { MAP_posterior = -dvPosteriorQuadTreeSelection(sampler); }
					else { MAP_posterior = -dvPosteriorQuadTree(sampler); }

					// calculate potential values for each of sample positions
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[2*v0_id+1] = v0 - posterioriRange/2.0 + (double)g*posterioriIncrement;
						for (int h = 0; h < posterioriSamples; h++) {
							sampler[2*v1_id+1] = v1 - posterioriRange/2.0 + (double)h*posterioriIncrement;
							if (file->treeMesh->selectionMode()) {
								posterioriValues[ h-g+(posterioriSamples-1) ] += exp(-MAP_posterior-dvPosteriorQuadTreeSelection(sampler));
							}
							else {
								posterioriValues[ h-g+(posterioriSamples-1) ] += exp(-MAP_posterior-dvPosteriorQuadTree(sampler));
							}
						}
					}

					// average the posteriori
					for (int e = 0; e < 2*posterioriSamples; e++) {
						posterioriValues[e] /= (posterioriSamples*posterioriSamples);
						if (logMax < posterioriValues[e]) { logMax = posterioriValues[e]; }
					}

					// calculate log values
					for (int h = 0; h < 2*posterioriSamples-1; h++) { posterioriValues[h] /= logMax; }

					drawPosteriorPlot(type,0);
					delete [] sampler;
				}
			} else if (file->optimizationMode == 3 && !this->potentialReferenceButton->value()) {
				posterioriSamples = (int) posterioriSampleNumberSlider->value();
				const double minV = minPosteriorSlider->value();
				const double maxV = maxPosteriorSlider->value();
				posterioriRange = maxV-minV;
				double logMax = -10000000.0;
				if (posterioriSamples > 5) {

					const int v0_id = file->treeMesh->selectedQuadLeaf->getIdentifier();

					posterioriIncrement = posterioriRange/posterioriSampleNumberSlider->value();

					posterioriValues = new double[posterioriSamples];

					double *sampler = new double[2*file->treeMesh->getTotalVariables()];

					// initialize potentials
					for (int i = 0; i < file->treeMesh->getTotalVariables(); i++) {
						file->treeMesh->getTree(file->treeMesh->quadTree,i);
						sampler[2*i] = file->treeMesh->identifiedQuadLeaf->getDiffusion();
						sampler[2*i+1] = file->treeMesh->identifiedQuadLeaf->getPotential();
					}

					// calculate diffusion values for each of sample positions
					for (int g = 0; g < posterioriSamples; g++) {
						sampler[2*v0_id+1] = minV + (double)g*posterioriIncrement;
						if (file->treeMesh->selectionMode()) { posterioriValues[g] = -dvPosteriorQuadTreeSelection(sampler); }
						else { posterioriValues[g] = -dvPosteriorQuadTree(sampler); }
					}

					// average the posteriori
					for (int e = 0; e < posterioriSamples; e++) {
						if (logMax < posterioriValues[e]) { logMax = posterioriValues[e]; }
					}

					// calculate log values
					for (int h = 0; h < posterioriSamples; h++) {
						posterioriValues[h] = exp(posterioriValues[h]-logMax) >= 0.0 ? exp(posterioriValues[h]-logMax) : 0.0;
					}
					drawPosteriorPlot(type,1);
					delete [] sampler;
				}

			}
		}
		break;
	}
}

void TreeMeshGui::drawPosteriorPlot(int type,int sup) {

	switch(type) {
		case 0: // diffusion
		{
			const double minD = minPosteriorSlider->value();
			const double maxD = maxPosteriorSlider->value();

			posterioriRange = maxD - minD;
			posterioriIncrement = (maxD-minD)/posterioriSampleNumberSlider->value();

			if (posterioriRange > 0.0 && file->inferred) {
				window->make_current();

				const int w = 320;
				const int h = 195;

	        	fl_line_style(FL_SOLID, 1, 0);

				fl_color(FL_WHITE); fl_rectf(170,35,w,h);
				fl_color(FL_BLACK); fl_rect(170,35,w,h);

				fl_line(170,190,170+w-1,190);
//				fl_font(iMAP->normalFont,10); fl_draw("[um\262/s]",170+w-43,225);

		        int pCurrent,pOld;
		        const double spacing = ((double)(w-1)/(double)posterioriSamples);
		        char xLabel [50];

		        for (int q = 1; q < posterioriSamples; q++) {
		        	// draw data points and connectors
		        	pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
		        	pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

		        	fl_color(FL_RED);
		        	fl_line_style(FL_SOLID, 2, 0);
					fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
//					fprintf(stderr,"DRAW\t%f\t%f\n",minD + (double)q*posterioriIncrement,posterioriValues[q]);
		        }
		        fl_color(FL_BLACK);
		        fl_font(iMAP->normalFont,10);
		        // draw x-axis labels
		        for (int t = 1; t < 9; t++) {
					sprintf(xLabel,"%.3f",minD+(float)t*((float)posterioriSamples/9.0)*posterioriIncrement);
					fl_draw(90,xLabel,177+t*w/9,227);
		        	fl_line_style(FL_DASH, 1, 0);
					fl_line(174+t*w/9,55,174+t*w/9,190);
		        }

				// draw MAP label
				sprintf(xLabel,"MAP: %.3f [um\262/s]",file->treeMesh->selectedQuadLeaf->getDiffusion());
		        fl_font(iMAP->normalFont,10);
				fl_draw(xLabel,170+2,45);
			}
		}
		break;
		case 1: // force
		{
			const double minF = minPosteriorSlider->value();
			const double maxF = maxPosteriorSlider->value();

			const double Fx = file->treeMesh->selectedQuadLeaf->getForceX();
			const double Fy = file->treeMesh->selectedQuadLeaf->getForceY();
			posterioriIncrement = (maxF-minF)/posterioriSampleNumberSlider->value();
			posterioriRange = maxF-minF;

			if (posterioriRange > 0.0 && file->inferred) {
				window->make_current();

				const int w = 320;
				const int h = 195;

	        	fl_line_style(FL_SOLID, 1, 0);

				fl_color(FL_WHITE); fl_rectf(170,35,w,h);
				fl_color(FL_BLACK); fl_rect(170,35,w,h);

				fl_line(170,190,170+w-1,190);
//				fl_font(iMAP->normalFont,10); fl_draw("[pN]",170+w-25,225);

		        int pCurrent,pOld;
		        const double spacing = ((double)(w-1)/(double)posterioriSamples);
		        char xLabel [50];

		        for (int q = 1; q < posterioriSamples; q++) {
		        	// draw data points and connectors
		        	pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
		        	pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

		        	fl_color(FL_RED);
		        	fl_line_style(FL_SOLID, 2, 0);
					fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
		        }
		        fl_color(FL_BLACK);
		        fl_font(iMAP->normalFont,10);
		        // draw x-axis labels
		        const double minF2 = sqrt(Fx*Fx+Fy*Fy) - posterioriRange/2.0;
		        for (int t = 1; t < 9; t++) {
					sprintf(xLabel,"%.1f",minF2+t*(posterioriSamples/9.0)*posterioriIncrement);
					fl_draw(90,xLabel,177+t*w/9,222);
		        	fl_line_style(FL_DASH, 1, 0);
					fl_line(174+t*w/9,55,174+t*w/9,190);
		        }
				// draw MAP label
				sprintf(xLabel,"MAP: %.1f [pN]",file->treeMesh->selectedQuadLeaf->getForce());
		        fl_font(iMAP->normalFont,10);
				fl_draw(xLabel,170+2,45);
			}
		}
			break;
		case 2:
		{
			switch(sup) {
				case 0:
				{
					const double maxSample = fabs(maxPosteriorSlider->value())/2.0;

					if (posterioriRange > 0.0 && file->inferred) {
						window->make_current();

						const int w = 320;
						const int h = 195;

						fl_line_style(FL_SOLID, 1, 0);

						fl_color(FL_WHITE); fl_rectf(170,35,w,h);
						fl_color(FL_BLACK); fl_rect(170,35,w,h);

						fl_line(170,190,170+w-1,190);
//						fl_font(iMAP->normalFont,10); fl_draw("[kT]",170+w-25,225);

						const double deltaV = file->treeMesh->selectedQuadLeaf->getPotential()
											  -file->treeMesh->selectedQuadLeaf2->getPotential();

						int pCurrent,pOld;
						const double spacing = ((double)(w-1)/(double)(2*posterioriSamples));
						char xLabel [50];

						for (int q = 1; q < 2*posterioriSamples; q++) {
							// draw data points and connectors
							pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
							pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

							fl_color(FL_RED);
							fl_line_style(FL_SOLID, 2, 0);
							fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
						}
						fl_color(FL_BLACK);
						fl_font(iMAP->normalFont,10);
						// draw x-axis labels
						for (int t = 1; t < 9; t++) {
							sprintf(xLabel,"%.2f",deltaV-maxSample+(float)t/9.0*posterioriRange);
							fl_draw(90,xLabel,177+t*w/9,227);
							fl_line_style(FL_DASH, 1, 0);
							fl_line(174+t*w/9,55,174+t*w/9,190);
						}
						// draw MAP label
						sprintf(xLabel,"MAP (reference): %.2f [kT]",deltaV);
						fl_font(iMAP->normalFont,10);
						fl_draw(xLabel,170+2,45);
					}
				}
					break;
				case 1:
				{
					const double minSample = minPosteriorSlider->value();
					const double maxSample = maxPosteriorSlider->value();

					posterioriRange = maxSample - minSample;

					if (posterioriRange > 0.0 && file->inferred) {
						window->make_current();

						const int w = 320;
						const int h = 195;

			        	fl_line_style(FL_SOLID, 1, 0);

						fl_color(FL_WHITE); fl_rectf(170,35,w,h);
						fl_color(FL_BLACK); fl_rect(170,35,w,h);

						fl_line(170,190,170+w-1,190);
//						fl_font(iMAP->normalFont,10); fl_draw("[kT]",170+w-25,225);

				        int pCurrent,pOld;
				        const double spacing = ((double)(w-1)/(double)posterioriSamples);
				        char xLabel [80];

				        for (int q = 1; q < posterioriSamples; q++) {
				        	// draw data points and connectors
				        	pOld = (int)( (posterioriValues[q-1])*((double)h-60.0));
				        	pCurrent = (int)( (posterioriValues[q])*((double)h-60.0));

				        	fl_color(FL_RED);
				        	fl_line_style(FL_SOLID, 2, 0);
							fl_line(171+(q-1)*spacing,h-6-pOld,171+q*spacing,h-6-pCurrent);
				        }
				        fl_color(FL_BLACK);
				        fl_font(iMAP->normalFont,10);
				        // draw x-axis labels
				        for (int t = 1; t < 9; t++) {
							sprintf(xLabel,"%.2f",minSample+t*(posterioriSamples/9.0)*posterioriIncrement);
							fl_draw(90,xLabel,177+t*w/9,227);
				        	fl_line_style(FL_DASH, 1, 0);
							fl_line(174+t*w/9,55,174+t*w/9,190);
				        }
						// draw MAP label
						sprintf(xLabel,"MAP: %.2f [kT] (no reference)",file->treeMesh->selectedQuadLeaf->getPotential());
				        fl_font(iMAP->normalFont,10);
						fl_draw(xLabel,170+2,45);
					}

				}
				break;
			}
		}
		break;
	}

}

void TreeMesh::savePosterior() {
	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	switch(file->treeMeshGui->getPosteriorType()) {
	case 0: // diffusion
		native->title("Save Diffusion Posterior Data");
		native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
		native->filter("Diffusion Posterior File\t*.dpost\n");
	    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
	    native->preset_file(file->fileName);
		if (native->show() == 0) {
			// store filename
			sprintf(nativeFilename,"%s%s",native->filename(),".dpost");
			FILE * writeFile;

			writeFile = fopen (nativeFilename,"wb");

			// write header
			fprintf(writeFile,"File Name: %s\n\n"
							  "Number of Points: %i\n"
							  "Zone Centre: [%.3f,%.3f]\n"
							  "D_MAP: %f [um2/s]\n\n",
							  file->fileName,
							  selectedQuadLeaf->getCount(),
							  selectedQuadLeaf->getXCentroid(),
							  selectedQuadLeaf->getYCentroid(),
							  selectedQuadLeaf->getDiffusion());

			for (int f = 0; f < file->treeMeshGui->posterioriSamples; f++) {
				fprintf(writeFile,"%f\t%f\n",file->treeMeshGui->minPosteriorSlider->value() + (double)f*file->treeMeshGui->posterioriIncrement,file->treeMeshGui->posterioriValues[f]);
			}
			fclose(writeFile);
		}
		break;
	case 1: // force
		native->title("Save Force Posterior Data");
		native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
		native->filter("Force Posterior File\t*.fpost\n");
	    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
	    native->preset_file(file->fileName);
		if (native->show() == 0) {
			// store filename
			sprintf(nativeFilename,"%s%s",native->filename(),".fpost");
			FILE * writeFile;

			writeFile = fopen (nativeFilename,"wb");

			// write header
			fprintf(writeFile,"File Name: %s\n\n"
							  "Number of Points: %i\n"
							  "Zone Centre: [%.3f,%.3f]\n"
							  "F_MAP: %f [pN]\n\n",
							  file->fileName,
							  selectedQuadLeaf->getCount(),
							  selectedQuadLeaf->getXCentroid(),
							  selectedQuadLeaf->getYCentroid(),
							  selectedQuadLeaf->getForce());

			const double minF2 = selectedQuadLeaf->getForce() - file->treeMeshGui->posterioriRange/2.0;

			for (int f = 0; f < file->treeMeshGui->posterioriSamples; f++) {
				fprintf(writeFile,"%f\t%f\n",minF2 + (double)f*file->treeMeshGui->posterioriIncrement,file->treeMeshGui->posterioriValues[f]);
			}
			fclose(writeFile);
		}
		break;
	case 2: // potential
		native->title("Save Potential Posterior Data");
		native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
		native->filter("Potential Posterior File\t*.vpost\n");
	    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
	    native->preset_file(file->fileName);
		if (native->show() == 0) {
			// store filename
			sprintf(nativeFilename,"%s%s",native->filename(),".vpost");
			FILE * writeFile;

			writeFile = fopen (nativeFilename,"wb");

			if (file->treeMeshGui->potentialReferenceButton->value()) {
				// write header
				const double deltaV = selectedQuadLeaf2->getPotential() - selectedQuadLeaf->getPotential();
				fprintf(writeFile,"File Name: %s\n\n"
								  "Number of Points: %i\n"
						  	  	  "Zone 1 Centre: [%.3f,%.3f]\n"
						  	  	  "Zone 2 Centre: [%.3f,%.3f]\n"
								  "V_MAP (reference): %f [kT]\n\n",
								  file->fileName,
								  selectedQuadLeaf->getCount(),
								  selectedQuadLeaf->getXCentroid(),selectedQuadLeaf->getYCentroid(),
								  selectedQuadLeaf2->getXCentroid(),selectedQuadLeaf2->getYCentroid(),
								  deltaV);

				for (int f = 0; f < 2*file->treeMeshGui->posterioriSamples; f++) {
//					fprintf(writeFile,"%f\t%f\n",deltaV-file->squareMeshGui->maxPosteriorSlider->value() + (double)f*file->squareMeshGui->posterioriIncrement,file->squareMeshGui->posterioriValues[f]);
					fprintf(writeFile,"%f\t%f\n",deltaV-file->treeMeshGui->maxPosteriorSlider->value()/2.0 + (float)f/(2*file->treeMeshGui->posterioriSamples)*file->treeMeshGui->posterioriRange,file->treeMeshGui->posterioriValues[f]);
				}
			} else {
				// write header
				fprintf(writeFile,"File Name: %s\n\n"
								  "Number of Points: %i\n"
								  "Zone Centre: [%.3f,%.3f]\n"
								  "V_MAP: %f [kT]\n\n",
								  file->fileName,
								  selectedQuadLeaf->getCount(),
								  selectedQuadLeaf->getXCentroid(),
								  selectedQuadLeaf->getYCentroid(),
								  selectedQuadLeaf->getPotential());
				for (int f = 0; f < file->treeMeshGui->posterioriSamples; f++) {
					fprintf(writeFile,"%f\t%f\n",file->treeMeshGui->minPosteriorSlider->value() + (double)f*file->treeMeshGui->posterioriIncrement,file->treeMeshGui->posterioriValues[f]);
				}
			}
			fclose(writeFile);
		}
		break;
	}
}

void diffusionOverlayCallback(Fl_Check_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	w->labelcolor(FL_WHITE);

	switch(meshType) {
	case 0:
		if (w->value()) {
			file->squareMeshGui->overlayPotentialButton->value(0);
			file->squareMeshGui->overlayForceMagnitudeButton->value(0);
			file->squareMeshGui->overlayPointNumberButton->value(0);
			file->squareMeshGui->overlayDiffusionLogButton->activate();
			if (file->squareMeshGui->landscapeButton->value()) {
				file->squareMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	case 1:
		if (w->value()) {
			file->voronoiMeshGui->overlayPotentialButton->value(0);
			file->voronoiMeshGui->overlayForceMagnitudeButton->value(0);
			file->voronoiMeshGui->overlayPointNumberButton->value(0);
			file->voronoiMeshGui->overlayDiffusionLogButton->activate();
			if (file->voronoiMeshGui->landscapeButton->value()) {
				file->voronoiMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	case 2:
		if (w->value()) {
			file->treeMeshGui->overlayPotentialButton->value(0);
			file->treeMeshGui->overlayForceMagnitudeButton->value(0);
			file->treeMeshGui->overlayPointNumberButton->value(0);
			file->treeMeshGui->overlayDiffusionLogButton->activate();
			if (file->treeMeshGui->landscapeButton->value()) {
				file->treeMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	}
	iMAP->overlayAdjustment = true;
	Fl::redraw();
}

void diffusionLogOverlayCallback(Fl_Check_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	w->labelcolor(FL_WHITE);

	switch(meshType) {
	case 0:
		if (w->value()) {
			file->squareMeshGui->overlayPotentialButton->value(0);
			file->squareMeshGui->overlayForceMagnitudeButton->value(0);
			file->squareMeshGui->overlayPointNumberButton->value(0);
			if (file->squareMeshGui->landscapeButton->value()) {
				file->squareMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	case 1:
		if (w->value()) {
			file->voronoiMeshGui->overlayPotentialButton->value(0);
			file->voronoiMeshGui->overlayForceMagnitudeButton->value(0);
			file->voronoiMeshGui->overlayPointNumberButton->value(0);
			if (file->voronoiMeshGui->landscapeButton->value()) {
				file->voronoiMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	case 2:
		if (w->value()) {
			file->treeMeshGui->overlayPotentialButton->value(0);
			file->treeMeshGui->overlayForceMagnitudeButton->value(0);
			file->treeMeshGui->overlayPointNumberButton->value(0);
			if (file->treeMeshGui->landscapeButton->value()) {
				file->treeMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	}

	iMAP->overlayAdjustment = true;
	Fl::redraw();
}

void potentialOverlayCallback(Fl_Check_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	w->labelcolor(FL_WHITE);

	switch(meshType) {
	case 0:
		if (w->value()) {
			file->squareMeshGui->overlayDiffusionButton->value(0);
			file->squareMeshGui->overlayForceMagnitudeButton->value(0);
			file->squareMeshGui->overlayPointNumberButton->value(0);
			file->squareMeshGui->overlayDiffusionLogButton->deactivate();
			if (file->squareMeshGui->landscapeButton->value()) {
				file->squareMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	case 1:
		if (w->value()) {
			file->voronoiMeshGui->overlayDiffusionButton->value(0);
			file->voronoiMeshGui->overlayForceMagnitudeButton->value(0);
			file->voronoiMeshGui->overlayPointNumberButton->value(0);
			file->voronoiMeshGui->overlayDiffusionLogButton->deactivate();
			if (file->voronoiMeshGui->landscapeButton->value()) {
				file->voronoiMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	case 2:
		if (w->value()) {
			file->treeMeshGui->overlayDiffusionButton->value(0);
			file->treeMeshGui->overlayForceMagnitudeButton->value(0);
			file->treeMeshGui->overlayPointNumberButton->value(0);
			file->treeMeshGui->overlayDiffusionLogButton->deactivate();
			if (file->treeMeshGui->landscapeButton->value()) {
				file->treeMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	}
	iMAP->overlayAdjustment = true;
	Fl::redraw();
}

void forceMagnitudeOverlayCallback(Fl_Check_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	w->labelcolor(FL_WHITE);

	switch(meshType) {
	case 0:
		if (w->value()) {
			file->squareMeshGui->overlayDiffusionButton->value(0);
			file->squareMeshGui->overlayPotentialButton->value(0);
			file->squareMeshGui->overlayPointNumberButton->value(0);
			file->squareMeshGui->overlayDiffusionLogButton->deactivate();
			if (file->squareMeshGui->landscapeButton->value()) {
				file->squareMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	case 1:
		if (w->value()) {
			file->voronoiMeshGui->overlayDiffusionButton->value(0);
			file->voronoiMeshGui->overlayPotentialButton->value(0);
			file->voronoiMeshGui->overlayPointNumberButton->value(0);
			file->voronoiMeshGui->overlayDiffusionLogButton->deactivate();
			if (file->voronoiMeshGui->landscapeButton->value()) {
				file->voronoiMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	case 2:
		if (w->value()) {
			file->treeMeshGui->overlayDiffusionButton->value(0);
			file->treeMeshGui->overlayPotentialButton->value(0);
			file->treeMeshGui->overlayPointNumberButton->value(0);
			file->treeMeshGui->overlayDiffusionLogButton->deactivate();
			if (file->treeMeshGui->landscapeButton->value()) {
				file->treeMeshGui->landscapeButton->do_callback();
			}
		}
		break;
	}
	iMAP->overlayAdjustment = true;
	Fl::redraw();
}

void pointNumberOverlayCallback(Fl_Check_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	w->labelcolor(FL_WHITE);

	switch(meshType) {
		case 0:
			if (w->value()) {
				file->squareMeshGui->overlayDiffusionButton->value(0);
				file->squareMeshGui->overlayPotentialButton->value(0);
				file->squareMeshGui->overlayForceMagnitudeButton->value(0);
				file->squareMeshGui->overlayDiffusionLogButton->deactivate();
				if (file->squareMeshGui->landscapeButton->value()) {
					file->squareMeshGui->landscapeButton->do_callback();
				}
			}
			break;
		case 1:
			if (w->value()) {
				file->voronoiMeshGui->overlayDiffusionButton->value(0);
				file->voronoiMeshGui->overlayPotentialButton->value(0);
				file->voronoiMeshGui->overlayForceMagnitudeButton->value(0);
				file->voronoiMeshGui->overlayDiffusionLogButton->deactivate();
				if (file->voronoiMeshGui->landscapeButton->value()) {
					file->voronoiMeshGui->landscapeButton->do_callback();
				}
			}
			break;
		case 2:
			if (w->value()) {
				file->treeMeshGui->overlayDiffusionButton->value(0);
				file->treeMeshGui->overlayPotentialButton->value(0);
				file->treeMeshGui->overlayForceMagnitudeButton->value(0);
				file->treeMeshGui->overlayDiffusionLogButton->deactivate();
				if (file->treeMeshGui->landscapeButton->value()) {
					file->treeMeshGui->landscapeButton->do_callback();
				}
			}
			break;
	}
	iMAP->overlayAdjustment = true;
	Fl::redraw();
}

void ColorBar::draw() {
	if (iMAP->fileNumber > 0) {
		char label [70];
		int colorBoxes = 90;
		unsigned char r,g,b;

		if (file->squareMeshOverlay || file->voronoiMeshOverlay || file->treeMeshOverlay) {
			if (file->squareMeshApplied() || file->voronoiMeshApplied() || file->treeMeshApplied()) {
				int cm;
				switch (file->meshType) {
				case 0: // square mesh
					cm = file->squareMeshGui->getColormap();
					break;
				case 1: // voronoi mesh
					cm = file->voronoiMeshGui->getColormap();
					break;
				case 2: // quad-tree mesh
					cm = file->treeMeshGui->getColormap();
					break;
				}

				for (int u = 0; u <= colorBoxes; u++) {
					const float* rgb = colormap(cm,(float)u/(float)colorBoxes,file->cMin,file->cMax,file->flipColormap);
					// draw colorbar
					r = 255*(rgb[0]);
					g = 255*(rgb[1]);
					b = 255*(rgb[2]);
					fl_color(r,g,b); fl_rectf(x(),y()-u*3.38,w(),h());
					// determine color overlay type
					if (u%15==0) {
						int overlayType = 0;
						switch (file->meshType) {
						case 0: // square mesh
							if (file->squareMeshGui->overlayForceMagnitudeButton->value() && file->optimizationMode == 2) { overlayType = 0; }
							else if (file->squareMeshGui->overlayForceMagnitudeButton->value()) { overlayType = 1; }
							else if (file->squareMeshGui->overlayDiffusionButton->value() && file->squareMeshGui->overlayDiffusionLogButton->value()) { overlayType = 3; }
							else if (file->squareMeshGui->overlayDiffusionButton->value()) { overlayType = 2; }
							else if (file->squareMeshGui->overlayPotentialButton->value()) { overlayType = 4; }
							else if (file->squareMeshGui->overlayPointNumberButton->value()) { overlayType = 5; }
							break;
						case 1: // voronoi mesh
							if (file->voronoiMeshGui->overlayForceMagnitudeButton->value() && file->optimizationMode == 2) { overlayType = 0; }
							else if (file->voronoiMeshGui->overlayForceMagnitudeButton->value()) { overlayType = 1; }
							else if (file->voronoiMeshGui->overlayDiffusionButton->value() && file->voronoiMeshGui->overlayDiffusionLogButton->value()) { overlayType = 3; }
							else if (file->voronoiMeshGui->overlayDiffusionButton->value()) { overlayType = 2; }
							else if (file->voronoiMeshGui->overlayPotentialButton->value()) { overlayType = 4; }
							else if (file->voronoiMeshGui->overlayPointNumberButton->value()) { overlayType = 5; }
							break;
						case 2: // quadtree mesh
							if (file->treeMeshGui->overlayForceMagnitudeButton->value() && file->optimizationMode == 2) { overlayType = 0; }
							else if (file->treeMeshGui->overlayForceMagnitudeButton->value()) { overlayType = 1; }
							else if (file->treeMeshGui->overlayDiffusionButton->value() && file->treeMeshGui->overlayDiffusionLogButton->value()) { overlayType = 3; }
							else if (file->treeMeshGui->overlayDiffusionButton->value()) { overlayType = 2; }
							else if (file->treeMeshGui->overlayPotentialButton->value()) { overlayType = 4; }
							else if (file->treeMeshGui->overlayPointNumberButton->value()) { overlayType = 5; }
							break;
						}

						switch(overlayType) {
						case 0: // drift magnitude
							if (file->meshType == 0) {
								sprintf(label,"%.3f um/s",(file->squareMesh->getForceMin()+(file->squareMesh->getForceMax()-file->squareMesh->getForceMin())*u/colorBoxes));
							} else if (file->meshType == 1) {
								sprintf(label,"%.3f um/s",(file->voronoiMesh->getForceMin()+(file->voronoiMesh->getForceMax()-file->voronoiMesh->getForceMin())*u/colorBoxes));
							} else if (file->meshType == 2) {
								sprintf(label,"%.3f um/s",(file->treeMesh->getForceMin()+(file->treeMesh->getForceMax()-file->treeMesh->getForceMin())*u/colorBoxes));
							}
							fl_draw(label,x()+55,y()-u*3.38+15);
							break;
						case 1: // force magnitude
							if (file->meshType == 0) {
								sprintf(label,"%.3f pN",(file->squareMesh->getForceMin()+(file->squareMesh->getForceMax()-file->squareMesh->getForceMin())*u/colorBoxes));
							} else if (file->meshType == 1) {
								sprintf(label,"%.3f pN",(file->voronoiMesh->getForceMin()+(file->voronoiMesh->getForceMax()-file->voronoiMesh->getForceMin())*u/colorBoxes));
							} else if (file->meshType == 2) {
								sprintf(label,"%.3f pN",(file->treeMesh->getForceMin()+(file->treeMesh->getForceMax()-file->treeMesh->getForceMin())*u/colorBoxes));
							}
							fl_draw(label,x()+55,y()-u*3.38+15);
							break;
						case 2: // diffusion coefficient
							if (file->meshType == 0) {
								sprintf(label,"%.3f um\262/s",(file->squareMesh->getDiffusionMin()+(file->squareMesh->getDiffusionMax()-file->squareMesh->getDiffusionMin())*u/colorBoxes));
							} else if (file->meshType == 1) {
								sprintf(label,"%.3f um\262/s",(file->voronoiMesh->getDiffusionMin()+(file->voronoiMesh->getDiffusionMax()-file->voronoiMesh->getDiffusionMin())*u/colorBoxes));
							} else if (file->meshType == 2) {
								sprintf(label,"%.3f um\262/s",(file->treeMesh->getDiffusionMin()+(file->treeMesh->getDiffusionMax()-file->treeMesh->getDiffusionMin())*u/colorBoxes));
							}
							fl_draw(label,x()+55,y()-u*3.38+15);
							break;
						case 3: // log diffusion coefficient
							if (file->meshType == 0) {
								sprintf(label,"%.3f log\num\262/s",(file->squareMesh->getDiffusionLogMin()+(file->squareMesh->getDiffusionLogMax()-file->squareMesh->getDiffusionLogMin())*u/colorBoxes));
							} else if (file->meshType == 1) {
								sprintf(label,"%.3f log\num\262/s",(file->voronoiMesh->getDiffusionLogMin()+(file->voronoiMesh->getDiffusionLogMax()-file->voronoiMesh->getDiffusionLogMin())*u/colorBoxes));
							} else if (file->meshType == 2) {
								sprintf(label,"%.3f log\num\262/s",(file->treeMesh->getDiffusionLogMin()+(file->treeMesh->getDiffusionLogMax()-file->treeMesh->getDiffusionLogMin())*u/colorBoxes));
							}
							fl_draw(label,x()+55,y()-u*3.38+15);
							break;
						case 4: // potential energy
							if (file->meshType == 0) {
								sprintf(label,"%.3f kT",(file->squareMesh->getPotentialMin()+(file->squareMesh->getPotentialMax()-file->squareMesh->getPotentialMin())*u/colorBoxes));
							} else if (file->meshType == 1) {
								sprintf(label,"%.3f kT",(file->voronoiMesh->getPotentialMin()+(file->voronoiMesh->getPotentialMax()-file->voronoiMesh->getPotentialMin())*u/colorBoxes));
							} else if (file->meshType == 2) {
								sprintf(label,"%.3f kT",(file->treeMesh->getPotentialMin()+(file->treeMesh->getPotentialMax()-file->treeMesh->getPotentialMin())*u/colorBoxes));
							}
							fl_draw(label,x()+55,y()-u*3.38+15);
							break;
						case 5: // point number
							if (file->meshType == 0) {
								sprintf(label,"%i",(file->squareMesh->getMinCell()->getCount()+(file->squareMesh->getMaxCell()->getCount()-file->squareMesh->getMinCell()->getCount())*u/colorBoxes));
							} else if (file->meshType == 1) {
								sprintf(label,"%i",(file->voronoiMesh->getMinCell()->getCount()+(file->voronoiMesh->getMaxCell()->getCount()-file->voronoiMesh->getMinCell()->getCount())*u/colorBoxes));
							} else if (file->meshType == 2) {
								sprintf(label,"%i",(file->treeMesh->getMinQuadTree()->getCount()+(file->treeMesh->getMaxQuadTree()->getCount()-file->treeMesh->getMinQuadTree()->getCount())*u/colorBoxes));
							}
							fl_draw(label,x()+55,y()-u*3.38+15);
							break;
						default:
							break;
						}
					}
				}
			}
		}
		fl_color(FL_WHITE); fl_rect(x(),y()-305,w(),325);
	}
	Fl_Box::draw();
}

void flipColormapCallback(Fl_Check_Button*w,int*v) {
	if (w->value()) {
		file->flipColormap = true;
	} else {
		file->flipColormap = false;
	}

	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {
	case 0: // square mesh
		if (file->landscapePlot) { generateSquareLandscape(); }
		break;
	case 1: // voronoi mesh
		if (file->landscapePlot) { generateVoronoiLandscape(); }
		break;
	case 2: // quad-tree mesh
		if (file->landscapePlot) { generateQuadTreeLandscape(); }
		break;
	}
}

SingleTrajectoryInferenceGui::SingleTrajectoryInferenceGui(char* filename) {

	const int h = 210;
	const int w = 400;

	minSampledDiffusion = 1.0e-3;
	maxSampledDiffusion = 1.0;
	diffusionRange = maxSampledDiffusion-minSampledDiffusion;
	diffusionSamples = 100;
	diffusionRange = maxSampledDiffusion-minSampledDiffusion;
	valuesSampled = false;

	startIndex = 0;
	endIndex = file->localizationCount-1;

	window = new Fl_Window(w,h,"Single Trajectory Inference");
	window->begin();
	window->color(iMAP->bgColor);
	window->set_non_modal();
	window->callback((Fl_Callback*)singleTrajectoryInferenceWindowCallback);

	tabs = new Fl_Tabs(0,0,w,h);
	tabs->labelsize(12);
	tabs->box(FL_BORDER_FRAME);
	tabs->labelcolor(FL_WHITE);
	tabs->labelfont(iMAP->normalFont);
	tabs->begin();

	inferenceGroup = new Fl_Group(0,25,w,h,"Inference");
	inferenceGroup->labelsize(12);
	inferenceGroup->labelcolor(FL_DARK1);
	inferenceGroup->labelfont(iMAP->boldFont);
	inferenceGroup->box(FL_BORDER_FRAME);
	inferenceGroup->begin();
	{
		filenameBox = new Fl_Box(10,30,w-20,20,filename);
		filenameBox->labelsize(12);
		filenameBox->align(FL_ALIGN_CENTER);
		filenameBox->labelcolor(FL_WHITE);
		filenameBox->labelfont(iMAP->boldFont);
		filenameBox->show();

		char label[100];

		sprintf(label,"Number of Points = %i",file->localizationCount);
		pointsBox = new Fl_Box(10,60,100,20,label);
		pointsBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
		pointsBox->labelsize(12);
		pointsBox->labelcolor(FL_WHITE);
		pointsBox->labelfont(iMAP->boldFont);
		pointsBox->show();

		noiseSigmaSlider = new Slider(w/2+5,65,(w-30)/2,20,"Noise Sigma [nm]");
		noiseSigmaSlider->precision(0);
		noiseSigmaSlider->labelsize(12);
		noiseSigmaSlider->labelcolor(FL_WHITE);
		noiseSigmaSlider->labelfont(iMAP->boldFont);
		noiseSigmaSlider->bounds(0.0,100.0);
		noiseSigmaSlider->value(30);
		noiseSigmaSlider->show();

		trajectoryStartSlider = new Slider(10,105,(w-30)/2,20,"Trajectory Start [ms]");
		trajectoryStartSlider->precision(3);
		trajectoryStartSlider->labelsize(12);
		trajectoryStartSlider->bounds(file->tMin,file->tMax);
		trajectoryStartSlider->value(file->tMin);
		trajectoryStartSlider->callback((Fl_Callback*)singleTrajectoryStartCallback);
		trajectoryStartSlider->labelcolor(FL_WHITE);
		trajectoryStartSlider->labelfont(iMAP->boldFont);
		trajectoryStartSlider->show();

		trajectoryEndSlider = new Slider(w/2+5,105,(w-30)/2,20,"Trajectory End [ms]");
		trajectoryEndSlider->precision(3);
		trajectoryEndSlider->labelsize(12);
		trajectoryEndSlider->bounds(file->tMin,file->tMax);
		trajectoryEndSlider->value(file->tMax);
		trajectoryEndSlider->labelcolor(FL_WHITE);
		trajectoryEndSlider->labelfont(iMAP->boldFont);
		trajectoryEndSlider->callback((Fl_Callback*)singleTrajectoryEndCallback);
		trajectoryEndSlider->show();

		inferenceButton = new Fl_Button(10,130,(w-30)/2,25,"Infer");
		inferenceButton->labelsize(12);
		inferenceButton->labelcolor(FL_BLACK);
		inferenceButton->labelfont(iMAP->boldFont);
		inferenceButton->callback((Fl_Callback*)singleTrajectoryInferenceCallback);
		inferenceButton->show();

		saveButton = new Fl_Button(w/2+5,130,(w-30)/2,25,"Save");
		saveButton->labelsize(12);
		saveButton->labelcolor(FL_BLACK);
		saveButton->labelfont(iMAP->boldFont);
		saveButton->callback((Fl_Callback*)singleTrajectorySaveCallback);
		saveButton->show();

		diffusionBox = new Fl_Box(10,160,(w-30)/2,20,"D = - [um\262/s]");
		diffusionBox->align(FL_ALIGN_CENTER);
		diffusionBox->labelsize(12);
		diffusionBox->labelcolor(FL_WHITE);
		diffusionBox->labelfont(iMAP->boldFont);
		diffusionBox->show();

		xForceBox = new Fl_Box(w/2,160,(w-30)/2,20,"Fx = - [pN]");
		xForceBox->align(FL_ALIGN_CENTER);
		xForceBox->labelsize(12);
		xForceBox->labelcolor(FL_WHITE);
		xForceBox->labelfont(iMAP->boldFont);
		xForceBox->show();

		yForceBox = new Fl_Box(w/2,180,(w-30)/2,20,"Fy = - [pN]");
		yForceBox->align(FL_ALIGN_CENTER);
		yForceBox->labelsize(12);
		yForceBox->labelcolor(FL_WHITE);
		yForceBox->labelfont(iMAP->boldFont);
		yForceBox->show();

		forceMagnitudeBox = new Fl_Box(10,180,(w-30)/2,20,"||F|| = - [pN]");
		forceMagnitudeBox->align(FL_ALIGN_CENTER);
		forceMagnitudeBox->labelsize(12);
		forceMagnitudeBox->labelcolor(FL_WHITE);
		forceMagnitudeBox->labelfont(iMAP->boldFont);
		forceMagnitudeBox->show();
	}
	inferenceGroup->end();

	diffusionGroup = new Fl_Group(0,25,w,h,"Diffusion Posterior");
	diffusionGroup->labelsize(12);
	diffusionGroup->labelcolor(FL_DARK1);
	diffusionGroup->labelfont(iMAP->boldFont);
	diffusionGroup->box(FL_BORDER_FRAME);
	diffusionGroup->begin();
	{
		minDiffusionSlider = new Slider(290,45,100,20,"Minimum [um\262/s]");
		minDiffusionSlider->precision(3);
		minDiffusionSlider->value(1.0e-3);
		minDiffusionSlider->bounds(1.0e-3,1.0);
		minDiffusionSlider->labelcolor(FL_WHITE);
		minDiffusionSlider->labelfont(iMAP->boldFont);
		minDiffusionSlider->show();

		maxDiffusionSlider = new Slider(290,80,100,20,"Maximum [um\262/s]");
		maxDiffusionSlider->precision(3);
		maxDiffusionSlider->value(1.0);
		maxDiffusionSlider->bounds(1.0e-3,1.0);
		maxDiffusionSlider->labelcolor(FL_WHITE);
		maxDiffusionSlider->labelfont(iMAP->boldFont);
		maxDiffusionSlider->show();

		diffusionSampleNumberSlider = new Slider(290,115,100,20,"Samples");
		diffusionSampleNumberSlider->precision(0);
		diffusionSampleNumberSlider->value(10000);
		diffusionSampleNumberSlider->bounds(100,10000);
		diffusionSampleNumberSlider->labelcolor(FL_WHITE);
		diffusionSampleNumberSlider->labelfont(iMAP->boldFont);
		diffusionSampleNumberSlider->show();

		sampleDiffusionButton = new Fl_Button(290,140,100,35,"Sample\nPosteriori");
		sampleDiffusionButton->labelcolor(FL_BLACK);
		sampleDiffusionButton->labelfont(iMAP->boldFont);
		sampleDiffusionButton->labelsize(12);
		sampleDiffusionButton->callback((Fl_Callback*)singleTrajectorySampleDiffusionCallback);
		sampleDiffusionButton->deactivate();
		sampleDiffusionButton->show();

		saveDPosterioriButton = new Fl_Button(290,180,100,20,"Save");
		saveDPosterioriButton->labelcolor(FL_BLACK);
		saveDPosterioriButton->labelfont(iMAP->boldFont);
		saveDPosterioriButton->labelsize(12);
		saveDPosterioriButton->callback((Fl_Callback*)singleTrajectorySampleDiffusionSaveCallback);
		saveDPosterioriButton->deactivate();
		saveDPosterioriButton->show();
	}
	diffusionGroup->end();

	forceGroup = new Fl_Group(0,25,w,h,"Force Posterior");
	forceGroup->labelsize(12);
	forceGroup->labelcolor(FL_DARK1);
	forceGroup->labelfont(iMAP->boldFont);
	forceGroup->box(FL_BORDER_FRAME);
	forceGroup->begin();
	{
		minForceSlider = new Slider(290,45,100,20,"Minimum [pN]");
		minForceSlider->precision(1);
		minForceSlider->value(-10);
		minForceSlider->bounds(-10,10);
		minForceSlider->labelcolor(FL_WHITE);
		minForceSlider->labelfont(iMAP->boldFont);
		minForceSlider->show();

		maxForceSlider = new Slider(290,80,100,20,"Maximum [pN]");
		maxForceSlider->precision(1);
		maxForceSlider->value(10);
		maxForceSlider->bounds(-10,10);
		maxForceSlider->labelcolor(FL_WHITE);
		maxForceSlider->labelfont(iMAP->boldFont);
		maxForceSlider->show();

		forceSampleNumberSlider = new Slider(290,115,100,20,"Samples");
		forceSampleNumberSlider->precision(0);
		forceSampleNumberSlider->value(100);
		forceSampleNumberSlider->bounds(100,10000);
		forceSampleNumberSlider->labelcolor(FL_WHITE);
		forceSampleNumberSlider->labelfont(iMAP->boldFont);
		forceSampleNumberSlider->show();

		sampleForceButton = new Fl_Button(290,140,100,35,"Sample\nPosteriori");
		sampleForceButton->labelcolor(FL_BLACK);
		sampleForceButton->labelfont(iMAP->boldFont);
		sampleForceButton->labelsize(12);
		sampleForceButton->callback((Fl_Callback*)singleTrajectorySampleForceCallback);
		sampleForceButton->deactivate();
		sampleForceButton->show();

		saveForcePosterioriButton = new Fl_Button(290,180,100,20,"Save");
		saveForcePosterioriButton->labelcolor(FL_BLACK);
		saveForcePosterioriButton->labelfont(iMAP->boldFont);
		saveForcePosterioriButton->labelsize(12);
		saveForcePosterioriButton->callback((Fl_Callback*)singleTrajectorySampleForceSaveCallback);
		saveForcePosterioriButton->deactivate();
		saveForcePosterioriButton->show();
	}
	forceGroup->end();

	tabs->end();

	window->end();
	window->show();

	iMAP->singleTrajectoryInferenceMode = true;
}

//void batchTrajectoryInferenceWindowCallback(Fl_Window*w,void*) {
//	iMAP->batchTrajectoryInferenceGui->hide();
//	iMAP->batchTrajectoryInferenceMode = false;
//}

void singleTrajectoryInferenceWindowCallback(Fl_Window*w,void*) {
	iMAP->singleTrajectoryInferenceGui->hide();
	iMAP->singleTrajectoryInferenceMode = false;
	iMAP->drawTrajectoriesButton->value(1);
	iMAP->drawTrajectoriesButton->do_callback();
	iMAP->overlayAdjustment = true;
}

void singleTrajectoryStartCallback(Slider*w,void*) {
	iMAP->singleTrajectoryInferenceGui->setStart(w->value());
}

void singleTrajectoryEndCallback(Slider*w,void*) {
	iMAP->singleTrajectoryInferenceGui->setEnd(w->value());
}

void SingleTrajectoryInferenceGui::setStart(double start) {
	for (int c = 0; c < file->localizationCount; c++) {
		if (file->tPointer[c]-start < 0.0001 && c <= endIndex) {
			// mark detections to be blacked out
			startIndex = c;
		}
	}
	char label[100];
	sprintf(label,"Number of Points = %i",endIndex-startIndex);
	pointsBox->copy_label(label);
}

void SingleTrajectoryInferenceGui::setEnd(double end) {
	for (int c = file->localizationCount-1; c >= 0; c--) {
		if (file->tPointer[c]-end > 0.0001 && c >= startIndex) {
			// mark detections to be blacked out
			endIndex = c;
		}
	}
	char label[100];
	sprintf(label,"Number of Points = %i",endIndex-startIndex);
	pointsBox->copy_label(label);
}

void SingleTrajectoryInferenceGui::draw() {
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glLineWidth(2.0);
	glPushMatrix();
		updateView();
		glBegin(GL_LINES);
		float *rgb = new float[3];

		for (int p = startIndex; p < endIndex; p++) {
			rgb = colormap(9,(double)(p-startIndex+1)/(double)(endIndex-startIndex),file->cMin,file->cMax,false);
			glColor4f(rgb[0],rgb[1],rgb[2],1.0);
			glVertex2f(file->xPointer[p],file->yPointer[p]);
			glVertex2f(file->xPointer[p+1],file->yPointer[p+1]);
		}
		glEnd();
	glPopMatrix();

	glDisable(GL_BLEND);
}

void singleTrajectoryInferenceCallback(Fl_Button*,void*) {
	iMAP->singleTrajectoryInferenceGui->infer();
}

void SingleTrajectoryInferenceGui::infer() {
	const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
	const double D_eff_y = file->averageDy*file->averageDy/file->averageDt; // valeur initial de diffusion

	optimizationArray[0] = 0.5*(D_eff_x + D_eff_y);
	optimizationArray[1] = optimizationArray[2] = 0.0;

	double fret[1];
	int iterations[1];

	// perform inference
	dfpmin(optimizationArray, 3, GTOL, iterations, fret, singleTrajectoryPosterior, (dfunc));

	char label[100];

	sprintf(label,"D = %f [um\262/s]",optimizationArray[0]);
	diffusionBox->labelcolor(FL_WHITE);
	diffusionBox->labelfont(iMAP->boldFont);
	diffusionBox->copy_label(label);

	sprintf(label,"Fx = %f [pN]",optimizationArray[1]);
	xForceBox->labelcolor(FL_WHITE);
	xForceBox->labelfont(iMAP->boldFont);
	xForceBox->copy_label(label);

	sprintf(label,"Fy = %f [pN]",optimizationArray[2]);
	yForceBox->labelcolor(FL_WHITE);
	yForceBox->labelfont(iMAP->boldFont);
	yForceBox->copy_label(label);

	sprintf(label,"||F|| = %f [pN]",sqrt(optimizationArray[1]*optimizationArray[1]+optimizationArray[2]*optimizationArray[2]));
	forceMagnitudeBox->labelfont(iMAP->boldFont);
	forceMagnitudeBox->labelcolor(FL_WHITE);
	forceMagnitudeBox->copy_label(label);

	// adjust sliders in posterior window so that MAP is centred
//	this->minDiffusionSlider->value(optimizationArray[0]*0.2);
//	this->minDiffusionSlider->bounds(optimizationArray[0]*0.2,optimizationArray[0]*5);
//	this->maxDiffusionSlider->value(optimizationArray[0]*5);
//	this->maxDiffusionSlider->bounds(optimizationArray[0]*0.2,optimizationArray[0]*5);

	this->minDiffusionSlider->value(0.00001);
	this->minDiffusionSlider->bounds(0.00001,optimizationArray[0]*10);
	this->maxDiffusionSlider->value(optimizationArray[0]*10);
	this->maxDiffusionSlider->bounds(0.00001,optimizationArray[0]*10);

	const double min = 0.0;
	const double max = 10.0;

	this->minForceSlider->value(min);
	this->minForceSlider->bounds(min,max);
	this->minForceSlider->deactivate();
	this->maxForceSlider->value(max);
	this->maxForceSlider->bounds(min,max);

	this->sampleDiffusionButton->activate();
	this->sampleForceButton->activate();
	valuesSampled = true;
}

void singleTrajectorySaveCallback(Fl_Button*,void*) {
	iMAP->singleTrajectoryInferenceGui->save();
}

void SingleTrajectoryInferenceGui::save() {
	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Single Trajectory Inference Data");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Single Trajectory File\t*.straj\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".straj");
		FILE * writeFile;

		writeFile = fopen (nativeFilename,"wb");

		// write header
		fprintf(writeFile,"File Name: %s\n\n"
						  "Number of Points: %i\n"
						  "Duration: [%.3f,%.3f]\n"
						  "Noise Sigma: %.0f [nm]\n"
						  "Diffusion Coefficient: %f [um^2/s]\n"
						  "x-Force: %f [pN]\n"
						  "y-Force: %f [pN]\n"
						  "Force Magnitude: %f [pN]\n",
						  file->fileName,
						  endIndex-startIndex,
						  file->tPointer[startIndex],file->tPointer[endIndex],
						  noiseSigmaSlider->value(),
						  optimizationArray[0],
						  optimizationArray[1],
						  optimizationArray[2],
						  sqrt(optimizationArray[1]*optimizationArray[1] + optimizationArray[2]*optimizationArray[2]));

		fclose(writeFile);
	}
}

void singleTrajectorySampleDiffusionCallback(Fl_Button*,void*) {
	iMAP->singleTrajectoryInferenceGui->sampleDiffusionPosteriori();
}

void SingleTrajectoryInferenceGui::sampleDiffusionPosteriori() {
	diffusionSamples = (int) diffusionSampleNumberSlider->value();
	const double minD = minDiffusionSlider->value();
	const double maxD = maxDiffusionSlider->value();
	diffusionIncrement = (maxD-minD)/diffusionSampleNumberSlider->value();
	diffusionRange = maxD-minD;
	double logMax = -10000000.0;
	if (diffusionSamples > 5) {
		diffusionValues = new double[diffusionSamples];
		double sampler[3];
		sampler[1] = optimizationArray[1];
		sampler[2] = optimizationArray[2];
		// calculate diffusion values for each of sample positions
		for (int g = 0; g < diffusionSamples; g++) {
			sampler[0] = minD + (double)g*diffusionIncrement;
			diffusionValues[g] = -singleTrajectoryPosterior(sampler);
			if (logMax < diffusionValues[g]) { logMax = diffusionValues[g]; }
		}
		// calculate log values
		for (int h = 0; h < diffusionSamples; h++) {
			diffusionValues[h] = exp(diffusionValues[h]-logMax);
		}
		drawDiffusionPosterioriPlot();
		saveDPosterioriButton->activate();
	}
}

void SingleTrajectoryInferenceGui::drawDiffusionPosterioriPlot() {
	minSampledDiffusion = minDiffusionSlider->value();
	maxSampledDiffusion = maxDiffusionSlider->value();

	diffusionRange = maxSampledDiffusion - minSampledDiffusion;

	if (diffusionRange > 0.0 && valuesSampled) {
		window->make_current();

		const int w = 270;
		const int h = 165;

		char xLabel [50];

		fl_color(FL_WHITE); fl_rectf(10,35,w,h);
		fl_color(FL_BLACK); fl_rect(10,35,w,h);

		fl_line(10,163,279,163);
		sprintf(xLabel,"MAP: %.3f [um\262/s]",this->optimizationArray[0]);
		fl_font(iMAP->normalFont,10); fl_draw(xLabel,14,47);

        int pCurrent,pOld;
        const double spacing = ((double)(w-1)/(double)diffusionSamples);

        for (int q = 1; q < diffusionSamples; q++) {
        	// draw data points and connectors
        	pOld = (int)( (diffusionValues[q-1])*((double)h-60.0) );
        	pCurrent = (int)( (diffusionValues[q])*((double)h-60.0) );

        	fl_color(FL_RED);
        	fl_line_style(FL_SOLID, 2, 0);
			fl_line(11+(q-1)*spacing,h-3-pOld,11+q*spacing,h-3-pCurrent);

        }
        fl_color(FL_BLACK);
        fl_font(iMAP->normalFont,10);
        for (int t = 1; t < 8; t++) {
			sprintf(xLabel,"%.3f",minSampledDiffusion+(t)*(diffusionSamples/8.0)*diffusionIncrement);
			fl_draw(90,xLabel,20+t*w/8,196);

        	fl_line_style(FL_DASH, 1, 0);
			fl_line(17+t*w/8,55,17+t*w/8,162);
        }
	}
}

void singleTrajectorySampleDiffusionSaveCallback(Fl_Button*,void*) {
	iMAP->singleTrajectoryInferenceGui->saveDiffusionPosteriori();
}

void SingleTrajectoryInferenceGui::saveDiffusionPosteriori() {


	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Diffusion Posteriori Data");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Diffusion Posteriori\t*.dpost\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".dpost");
		FILE * writeFile;

		writeFile = fopen (nativeFilename,"wb");

		// write data
        for (int q = 0; q < diffusionSamples; q++) {
        	fprintf(writeFile,"%f\t%f\n",
        			minSampledDiffusion+q*diffusionIncrement,
        			diffusionValues[q]);
        }
		fclose(writeFile);
	}

}

void singleTrajectorySampleForceCallback(Fl_Button*,void*) {
	iMAP->singleTrajectoryInferenceGui->sampleForcePosteriori();
}

void SingleTrajectoryInferenceGui::sampleForcePosteriori() {
	forceSamples = (int) forceSampleNumberSlider->value();
	const double minFx = minForceSlider->value();
	const double maxFx = maxForceSlider->value();

	const double Fx = this->optimizationArray[1];
	const double Fy = this->optimizationArray[2];

	forceRange = fabs(maxFx-minFx);
	forceIncrement = forceRange/forceSampleNumberSlider->value();
	double logMax = -10000000.0;

	if (forceSamples > 5) {
		forceValues = new double[forceSamples];
		double sampler[3];
		sampler[0] = optimizationArray[0]; // diffusion
		// calculate diffusion values for each of sample positions
		const double theta = atan2(Fy,Fx);
		for (int g = 0; g < forceSamples; g++) {
			sampler[1] = (Fx-fabs(forceRange/2.0)*cos(theta)) + (double)g*forceIncrement*cos(theta);
			sampler[2] = (Fy-fabs(forceRange/2.0)*sin(theta)) + (double)g*forceIncrement*sin(theta);

			forceValues[g] = -singleTrajectoryPosterior(sampler);
			if (logMax < forceValues[g]) { logMax = forceValues[g]; }
		}
		// calculate log values
		for (int h = 0; h < forceSamples; h++) {
			forceValues[h] = exp(forceValues[h]-logMax) >= 0.0 ? exp(forceValues[h]-logMax) : 0.0;
		}
		drawForcePosterioriPlot();
		saveForcePosterioriButton->activate();
	}
}

void SingleTrajectoryInferenceGui::drawForcePosterioriPlot() {
	minSampledForce= minForceSlider->value();
	maxSampledForce= maxForceSlider->value();

	forceRange = maxSampledForce - minSampledForce;

	if (forceRange > 0.0 && valuesSampled) {
		window->make_current();

		const int w = 270;
		const int h = 165;

        char xLabel [50];

		fl_color(FL_WHITE); fl_rectf(10,35,w,h);
		fl_color(FL_BLACK); fl_rect(10,35,w,h);

		const double forceMag = sqrt(optimizationArray[1]*optimizationArray[1]+optimizationArray[2]*optimizationArray[2]);

		fl_line(10,163,279,163);
		sprintf(xLabel,"MAP: %.2f [pN]",forceMag);
		fl_font(iMAP->normalFont,10); fl_draw(xLabel,14,47);

        int pCurrent,pOld;
        const double spacing = ((double)(w-1)/(double)forceSamples);

        for (int q = 1; q < forceSamples; q++) {
        	// draw data points and connectors
        	pOld = (int)( (forceValues[q-1])*((double)h-60.0));
        	pCurrent = (int)( (forceValues[q])*((double)h-60.0));

        	fl_color(FL_RED);
        	fl_line_style(FL_SOLID, 2, 0);
			fl_line(11+(q-1)*spacing,h-3-pOld,11+q*spacing,h-3-pCurrent);

        }
        fl_color(FL_BLACK);
        fl_font(iMAP->normalFont,10);

		const double minF2 = forceMag - forceRange/2.0;

        for (int t = 1; t < 8; t++) {
			sprintf(xLabel,"%.2f",minF2+(t)*(forceSamples/8.0)*forceIncrement);
			fl_draw(90,xLabel,20+t*w/8,196);
        	fl_line_style(FL_DASH, 1, 0);
			fl_line(17+t*w/8,55,17+t*w/8,162);
        }
	}
}

void singleTrajectorySampleForceSaveCallback(Fl_Button*,void*) {
	iMAP->singleTrajectoryInferenceGui->saveForcePosteriori();
}

void SingleTrajectoryInferenceGui::saveForcePosteriori() {

	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Fx Posteriori Data");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Fx Posteriori\t*.fpost\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".fxpost");
		FILE * writeFile;

		writeFile = fopen (nativeFilename,"wb");

		const double forceMag = sqrt(optimizationArray[1]*optimizationArray[1]+optimizationArray[2]*optimizationArray[2]);
		const double minF2 = forceMag - forceRange/2.0;

		// write data
        for (int q = 0; q < forceSamples; q++) {
        	fprintf(writeFile,"%f\t%f\n",
        			minF2+q*forceIncrement,
        			forceValues[q]);
        }
		fclose(writeFile);
	}

}

//BatchTrajectoryInferenceGui::BatchTrajectoryInferenceGui() {
//	nPointsBins = 40;
//	nDiffusionBins = 40;
//	currentTrack = 0;
//	activeFiles = file->numberOfFiles;
//	diffusionRange = 0.0;
//	pointsRange = file->getMaxPoints() - file->getMinPoints();
//	minPoints = file->getMinPoints();
//	maxPoints = file->getMaxPoints();
//	inferred = false;
//
//	this->pointsBinTotals = NULL;
//	this->diffusionBinTotals = NULL;
//
//	const int h = 500;
//	const int w = 400;
//	char label[100];
//
//	window = new Fl_Window(w,h,"Batch Trajectory Inference");
//	window->begin();
//	window->set_non_modal();
//	window->callback((Fl_Callback*)batchTrajectoryInferenceWindowCallback);
//
//	tabs = new Fl_Tabs(0,0,w,h);
//	tabs->labelsize(12);
//	tabs->begin();
//
//	inferenceGroup = new Fl_Group(10,25,w,h,"Inference");
//	inferenceGroup->labelsize(12);
//	inferenceGroup->labelfont(1);
//	inferenceGroup->begin();
//	{
//		sprintf(label,"%i of %i Files Selected",file->numberOfFiles,file->numberOfFiles);
//		filenameBox = new Fl_Box(10,28,w-20,20,label);
//		filenameBox->labelsize(12);
//		filenameBox->labelfont(1);
//		filenameBox->align(FL_ALIGN_CENTER);
//		filenameBox->show();
//
//		minPointsSlider = new Slider(10,230,(w-30)/2,20,"Minimum Points");
//		minPointsSlider->precision(0);
//		minPointsSlider->labelsize(12);
//		minPointsSlider->bounds(file->getMinPoints(),file->getMaxPoints());
//		minPointsSlider->value(file->getMinPoints());
//		minPointsSlider->callback((Fl_Callback*)batchMinPointsSliderCallback);
//		minPointsSlider->labelfont(1);
//		minPointsSlider->show();
//
//		maxPointsSlider = new Slider(w/2+5,230,(w-30)/2,20,"Maximum Points");
//		maxPointsSlider->precision(0);
//		maxPointsSlider->labelsize(12);
//		maxPointsSlider->bounds(file->getMinPoints(),file->getMaxPoints());
//		maxPointsSlider->value(file->getMaxPoints());
//		maxPointsSlider->callback((Fl_Callback*)batchMaxPointsSliderCallback);
//		maxPointsSlider->labelfont(1);
//		maxPointsSlider->show();
//
//		inferButton = new Fl_Button(10,260,(w-30)/2,25,"Infer");
//		inferButton->labelsize(12);
//		inferButton->labelfont(1);
//		inferButton->callback((Fl_Callback*)batchTrajectoryInferenceCallback);
//		inferButton->show();
//
//		saveButton = new Fl_Button(w/2+5,260,(w-30)/2,25,"Save");
//		saveButton->labelsize(12);
//		saveButton->labelfont(1);
//		saveButton->callback((Fl_Callback*)batchTrajectorySaveCallback);
//		saveButton->deactivate();
//		saveButton->show();
//
//		minDiffusionSlider = new Slider(10,470,(w-30)/2,20,"Minimum Diffusion [um\262/s]");
//		minDiffusionSlider->precision(3);
//		minDiffusionSlider->labelsize(12);
//		minDiffusionSlider->labelfont(1);
//		minDiffusionSlider->show();
//		minDiffusionSlider->deactivate();
//
//		maxDiffusionSlider = new Slider(w/2+5,470,(w-30)/2,20,"Maximum Diffusion [um\262/s]");
//		maxDiffusionSlider->precision(3);
//		maxDiffusionSlider->labelsize(12);
//		maxDiffusionSlider->labelfont(1);
//		maxDiffusionSlider->show();
//		maxDiffusionSlider->deactivate();
//	}
//	inferenceGroup->end();
//
//	advancedGroup = new Fl_Group(10,25,w,h,"Advanced");
//	advancedGroup->labelsize(12);
//	advancedGroup->labelfont(1);
//	advancedGroup->begin();
//	{
//		noiseSigmaSlider = new Slider(10,45,(w-30)/2,20,"Noise Sigma [nm]");
//		noiseSigmaSlider->precision(0);
//		noiseSigmaSlider->labelsize(12);
//		noiseSigmaSlider->bounds(0.0,100.0);
//		noiseSigmaSlider->value(30);
//		noiseSigmaSlider->show();
//	}
//	advancedGroup->end();
//
//	tabs->end();
//	window->end();
//	window->show();
//}
//
//int BatchTrajectoryInferenceGui::getActiveFiles() {
//	// count number of active files
//	activeFiles = 0;
//	for (int a = 0; a < file->numberOfFiles; a++) {
//		if (file->tracks[a].n >= (int)minPointsSlider->value() && file->tracks[a].n <= (int)maxPointsSlider->value()) {
//			activeFiles++;
//		}
//	}
//	return activeFiles;
//}
//
//void batchTrajectoryPointsSliderCallback(Slider*w,void*) {
//	iMAP->batchTrajectoryInferenceGui->getActiveFiles();
//}
//
//void batchTrajectoryInferenceCallback(Fl_Button*w,void*) {
//	iMAP->batchTrajectoryInferenceGui->infer();
//}
//
//void batchMinPointsSliderCallback(Slider*w,void*) {
//	iMAP->batchTrajectoryInferenceGui->minPoints = (int)w->value();
//
//	char label[100];
//	sprintf(label,"%i of %i Files Selected",iMAP->batchTrajectoryInferenceGui->getActiveFiles(),file->numberOfFiles);
//	iMAP->batchTrajectoryInferenceGui->filenameBox->copy_label(label);
//
//	iMAP->batchTrajectoryInferenceGui->drawPointsHistogram();
//
//	iMAP->batchTrajectoryInferenceGui->saveButton->deactivate();
//}
//
//void batchMaxPointsSliderCallback(Slider*w,void*) {
//	iMAP->batchTrajectoryInferenceGui->maxPoints = (int)w->value();
//
//	char label[100];
//	sprintf(label,"%i of %i Files Selected",iMAP->batchTrajectoryInferenceGui->getActiveFiles(),file->numberOfFiles);
//	iMAP->batchTrajectoryInferenceGui->filenameBox->copy_label(label);
//
//	iMAP->batchTrajectoryInferenceGui->drawPointsHistogram();
//
//	iMAP->batchTrajectoryInferenceGui->saveButton->deactivate();
//}
//
//void batchMinDiffusionSliderCallback(Slider*w,void*) {
//	if (iMAP->batchTrajectoryInferenceGui->inferred) {
//		iMAP->batchTrajectoryInferenceGui->minDiffusion = w->value();
//		iMAP->batchTrajectoryInferenceGui->drawDiffusionHistogram();
//	}
//}
//
//void batchMaxDiffusionSliderCallback(Slider*w,void*) {
//	if (iMAP->batchTrajectoryInferenceGui->inferred) {
//		iMAP->batchTrajectoryInferenceGui->maxDiffusion = w->value();
//		iMAP->batchTrajectoryInferenceGui->drawDiffusionHistogram();
//	}
//}
//
//void BatchTrajectoryInferenceGui::drawPointsHistogram() {
//	pointsRange = maxPoints - minPoints;
//	if (pointsRange > 0) {
//		// create bins
//		if (pointsBinTotals != NULL) { delete [] pointsBinTotals; }
//		pointsBinTotals = new int[nPointsBins];
//
//		for (int d = 0; d < nPointsBins; d++) {
//			pointsBinTotals[d] = 0;
//		}
//
//		// accumulate trajectory point totals
//		int total = 0;
//		for (int f = 0; f < file->numberOfFiles; f++) {
//			if (file->tracks[f].n <= maxPoints && file->tracks[f].n >= minPoints) {
//				const int idx = (int)( (float)(file->tracks[f].n-minPoints)/(float)(pointsRange)*(float)(nPointsBins) )-1;
//				pointsBinTotals[idx]++;
//				total++;
//			}
//		}
//		if (total > 1) {
//            int maxBin = -1000000;
//            int maxBinIndex;
//            for (int j = 0; j < nPointsBins; j++) {
//                    if (maxBin < pointsBinTotals[j]) {
//                            maxBin = pointsBinTotals[j];
//                            maxBinIndex = j;
//                    }
//            }
//
//            window->make_current();
//            fl_color(FL_WHITE); fl_rectf(10,50,400-20,160);
//            fl_color(FL_BLACK); fl_rect(10,50,400-20,160);
//
//            fl_line(10,180,400-11,180);
//            fl_font(1,12); fl_draw("Number of Points",14,63);
//
//            const int binWidth = (int)((float)(400-20)/(float)nPointsBins);
//            int binLength;
//            const int h = 160;
//            for (int n = 0; n < nPointsBins; n++) {
//				binLength = (h-70)*pointsBinTotals[n]/maxBin;
//
//				fl_color((unsigned char)(255*0.3),
//								 (unsigned char)(255*1.0),
//								 (unsigned char)(255*0.3));
//
//				fl_rectf(22+n*binWidth,51+h-binLength-30,binWidth,binLength);
//				fl_color(FL_BLACK); fl_rect(22+n*binWidth,51+h-binLength-30,binWidth+1,binLength);
//
//				if (n%5 == 0 || n == nPointsBins-1) {
//						char points [30];
//						sprintf(points,"%i",(int)((float)minPoints+(float)((float)n/(float)nPointsBins)*(float)pointsRange));
//						fl_font(0,10); fl_draw(90,points,27+n*binWidth,203);
//				}
//
//				if (n == maxBinIndex) {
//						char maxBinString [30];
//						sprintf(maxBinString,"%i",pointsBinTotals[maxBinIndex]);
//						fl_font(0,10); fl_draw(45,maxBinString,24+n*binWidth,88);
//				}
//            }
//		}
//	}
//}
//
//void BatchTrajectoryInferenceGui::drawDiffusionHistogram() {
//
//	minDiffusion = minDiffusionSlider->value();
//	maxDiffusion = maxDiffusionSlider->value();
//
//	diffusionRange = maxDiffusion - minDiffusion;
//
//	if (diffusionRange > 0.01) {
//		// create bins
//		if (diffusionBinTotals != NULL) { delete [] diffusionBinTotals; }
//		diffusionBinTotals = new int[nDiffusionBins];
//
//		for (int d = 0; d < nDiffusionBins; d++) {
//			diffusionBinTotals[d] = 0;
//		}
//
//		// accumulate trajectory point totals
//		int total = 0;
//		for (int f = 0; f < activeFiles; f++) {
//			if (diffusionValues[f] <= maxDiffusion && diffusionValues[f] >= minDiffusion) {
//				const int idx = (int)( (diffusionValues[f]-minDiffusion)/diffusionRange*(float)(nDiffusionBins) )-1;
//				diffusionBinTotals[idx]++;
//				total++;
//			}
//		}
//		if (total > 1) {
//            int maxBin = -1000000;
//            int maxBinIndex;
//            for (int j = 0; j < nDiffusionBins; j++) {
//                    if (maxBin < diffusionBinTotals[j]) {
//                            maxBin = diffusionBinTotals[j];
//                            maxBinIndex = j;
//                    }
//            }
//
//            window->make_current();
//            fl_color(FL_WHITE); fl_rectf(10,290,400-20,160);
//            fl_color(FL_BLACK); fl_rect(10,290,400-20,160);
//
//            fl_line(10,420,400-11,420);
//            fl_font(1,12); fl_draw("Diffusion Coefficient [um\262/s]",14,303);
//
//            const int binWidth = (int)((float)(400-20)/(float)nDiffusionBins);
//            int binLength;
//            const int h = 160;
//            for (int n = 0; n < nDiffusionBins; n++) {
//				binLength = (h-70)*diffusionBinTotals[n]/maxBin;
//
//				fl_color((unsigned char)(255*0.3),
//								 (unsigned char)(255*0.3),
//								 (unsigned char)(255*1.0));
//
//				fl_rectf(22+n*binWidth,291+h-binLength-30,binWidth,binLength);
//				fl_color(FL_BLACK); fl_rect(22+n*binWidth,291+h-binLength-30,binWidth+1,binLength);
//
//				if (n%5 == 0 || n == nDiffusionBins-1) {
//						char points [30];
//						sprintf(points,"%.3f",minDiffusion+((float)n/(float)nDiffusionBins)*diffusionRange);
//						fl_font(0,10); fl_draw(90,points,27+n*binWidth,446);
//				}
//
//				if (n == maxBinIndex) {
//						char maxBinString [30];
//						sprintf(maxBinString,"%i",diffusionBinTotals[maxBinIndex]);
//						fl_font(0,10); fl_draw(45,maxBinString,24+n*binWidth,328);
//				}
//            }
//
//		}
//	}
//
//}
//
//void BatchTrajectoryInferenceGui::infer() {
//
//	double fret[1];
//	int iterations[1];
//
//	diffusionValues = new double[activeFiles];
//	xForceValues = new double[activeFiles];
//	yForceValues = new double[activeFiles];
//
//	int a = 0;
////	double wtime = omp_get_wtime ();
////	#pragma omp for
//	for (int b = 0; b < file->numberOfFiles; b++) {
//		if (file->tracks[b].n >= (int)minPointsSlider->value() && file->tracks[b].n <= (int)maxPointsSlider->value()) {
//
//			currentTrack = b;
//
//			const double D_eff_x = file->tracks[b].averageDx*file->tracks[b].averageDx/file->tracks[b].averageDt;
//			const double D_eff_y = file->tracks[b].averageDy*file->tracks[b].averageDy/file->tracks[b].averageDt;
//
//			optimizationArray[0] = 0.5*(D_eff_x + D_eff_y);
//			optimizationArray[1] = optimizationArray[2] = 0.0;
//
//			// perform inference
//			dfpmin(optimizationArray, 3, GTOL, iterations, fret, batchTrajectoryPosterior, (dfunc));
//
//			diffusionValues[a] = optimizationArray[0];
//			xForceValues[a] = optimizationArray[1];
//			yForceValues[a] = optimizationArray[2];
//
//			a++;
//		}
//	}
////	wtime = omp_get_wtime () - wtime;
////	fprintf(stderr,"Calculation Time = %f\n",wtime);
//
//	minDiffusion = 10000000.0;
//	maxDiffusion = -10000000.0;
//
//	for (int h = 0; h < activeFiles; h++) {
//		if (minDiffusion > diffusionValues[h]) { minDiffusion = diffusionValues[h]; }
//		if (maxDiffusion < diffusionValues[h]) { maxDiffusion = diffusionValues[h]; }
//	}
//
//	// update diffusion gui & histogram
//	inferred = true;
//	minDiffusionSlider->activate();
//	minDiffusionSlider->bounds(minDiffusion,maxDiffusion);
//	minDiffusionSlider->value(minDiffusion);
//	minDiffusionSlider->callback((Fl_Callback*)batchMinDiffusionSliderCallback);
//	maxDiffusionSlider->activate();
//	maxDiffusionSlider->bounds(minDiffusion,maxDiffusion);
//	maxDiffusionSlider->value(maxDiffusion);
//	maxDiffusionSlider->callback((Fl_Callback*)batchMaxDiffusionSliderCallback);
//	saveButton->activate();
//
//	drawDiffusionHistogram();
//}
//
//void batchTrajectorySaveCallback(Fl_Button*,void*) {
//	iMAP->batchTrajectoryInferenceGui->save();
//}
//
//void BatchTrajectoryInferenceGui::save() {
//
//	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
//	char nativeFilename[FILENAME_MAX];
//
//	native->title("Save Batch Trajectory Inference Data");
//	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
//	native->filter("Batch Trajectory File\t*.btraj\n");
//    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
//    native->preset_file(file->fileName);
//
//	if (native->show() == 0) {
//		// store filename
//		sprintf(nativeFilename,"%s%s",native->filename(),".btraj");
//		FILE * writeFile;
//
//		writeFile = fopen (nativeFilename,"wb");
//
//		// write header
//		fprintf(writeFile,"%i Files\n\n"
//						  "Points Range [%i,%i]\n"
//						  "Diffusion Range [%f,%f] [um^2/s]\n\n",
//						  getActiveFiles(),
//						  minPoints,maxPoints,
//						  minDiffusion,maxDiffusion);
//
//		// write mesh data
////		fprintf(writeFile,"Points\tD\t\tFx\t\tFy\t\t||F||\t\tFile Name\n");
//		fprintf(writeFile,"Points\tD\t\tFx\t\tFy\t\t||F||\n");
//
//		int a = 0;
//		for (int b = 0; b < file->numberOfFiles; b++) {
//			if (file->tracks[b].n >= (int)minPointsSlider->value() && file->tracks[b].n <= (int)maxPointsSlider->value()) {
////				fprintf(writeFile,"%i\t%f\t%f\t%f\t%f\t%s\n",
////								  file->tracks[b].n,
////								  this->diffusionValues[a],
////								  this->xForceValues[a],
////								  this->yForceValues[a],
////								  sqrt(xForceValues[a]*xForceValues[a]+yForceValues[a]*yForceValues[a]),
////								  file->fileNameList[b]
////								  );
////				a++;
//				fprintf(writeFile,"%i\t%f\t%f\t%f\t%f\n",
//								  file->tracks[b].n,
//								  this->diffusionValues[a],
//								  this->xForceValues[a],
//								  this->yForceValues[a],
//								  sqrt(xForceValues[a]*xForceValues[a]+yForceValues[a]*yForceValues[a])
//								  );
//				a++;
//			}
//		}
//
//		fclose (writeFile);
//	}
//
//}

TreeMeshGui::TreeMeshGui() {

	if (iMAP->fileNumber > 0) {

		colormap = 9;

		const int w = 500;
		const int h = 300;

		window = new Fl_Window(w,h,"Quad-Tree Meshing");
		window->callback(treeMeshWindowCallback);
		window->color(iMAP->bgColor);
		window->begin();
		window->set_non_modal();

		// Permanent interface buttons
		{
			overlayPointNumberButton = new Fl_Check_Button(5,245,65,20,"Points");
			overlayPointNumberButton->labelsize(12);
			overlayPointNumberButton->value(1);
			overlayPointNumberButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayPointNumberButton->labelfont(iMAP->normalFont);
			overlayPointNumberButton->labelcolor(FL_WHITE);
			overlayPointNumberButton->callback((Fl_Callback*)pointNumberOverlayCallback,(int*)2);
			overlayPointNumberButton->show();

			overlayDiffusionButton = new Fl_Check_Button(70,245,165,20,"Diffusion Coefficient");
			overlayDiffusionButton->labelsize(12);
			overlayDiffusionButton->callback((Fl_Callback*)diffusionOverlayCallback,(int*)2);
			overlayDiffusionButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayDiffusionButton->labelfont(iMAP->normalFont);
			overlayDiffusionButton->labelcolor(iMAP->bgColor);
			overlayDiffusionButton->deactivate();
			overlayDiffusionButton->show();

			overlayDiffusionLogButton = new Fl_Check_Button(70,260,165,20,"Log");
			overlayDiffusionLogButton->labelsize(12);
			overlayDiffusionLogButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayDiffusionLogButton->labelfont(iMAP->normalFont);
			overlayDiffusionLogButton->labelcolor(iMAP->bgColor);
			overlayDiffusionLogButton->callback((Fl_Callback*)diffusionLogOverlayCallback,(int*)2);
			overlayDiffusionLogButton->deactivate();
			overlayDiffusionLogButton->show();

			overlayForceMagnitudeButton = new Fl_Check_Button(240,245,135,20,"Force Magnitude");
			overlayForceMagnitudeButton->labelsize(12);
			overlayForceMagnitudeButton->callback((Fl_Callback*)forceMagnitudeOverlayCallback,(int*)2);
			overlayForceMagnitudeButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayForceMagnitudeButton->labelfont(iMAP->normalFont);
			overlayForceMagnitudeButton->labelcolor(iMAP->bgColor);
			overlayForceMagnitudeButton->deactivate();
			overlayForceMagnitudeButton->show();

			overlayForceArrowsButton = new Fl_Check_Button(240,260,135,20,"Arrows");
			overlayForceArrowsButton->labelsize(12);
			overlayForceArrowsButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayForceArrowsButton->labelfont(iMAP->normalFont);
			overlayForceArrowsButton->labelcolor(iMAP->bgColor);
			overlayForceArrowsButton->deactivate();
			overlayForceArrowsButton->show();

			overlayPotentialButton = new Fl_Check_Button(365,245,140,20,"Potential Energy");
			overlayPotentialButton->labelsize(12);
			overlayPotentialButton->callback((Fl_Callback*)potentialOverlayCallback,(int*)2);
			overlayPotentialButton->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
			overlayPotentialButton->labelfont(iMAP->normalFont);
			overlayPotentialButton->labelcolor(iMAP->bgColor);
			overlayPotentialButton->deactivate();
			overlayPotentialButton->show();

			variablesBox = new Fl_Box(0,275,w,25,"0 Zones");
			variablesBox->labelsize(14);
			variablesBox->align(FL_ALIGN_CENTER);
			variablesBox->labelfont(iMAP->boldFont);
			variablesBox->labelcolor(FL_WHITE);
			variablesBox->show();
		}

		tabs = new Fl_Tabs(0,0,w,h-60);
		tabs->labelsize(12);
//		tabs->color(FL_DARK1,FL_DARK1);
		tabs->box(FL_BORDER_FRAME);
		tabs->labelcolor(FL_WHITE);
		tabs->labelfont(iMAP->normalFont);
		tabs->begin();

		inferenceGroup = new Fl_Group(0,25,w,h-85,"Inference");
		inferenceGroup->labelsize(12);
//		inferenceGroup->color(FL_WHITE,iMAP->bgColor);
		inferenceGroup->labelcolor(FL_DARK1);
		inferenceGroup->labelfont(iMAP->boldFont);
		inferenceGroup->box(FL_BORDER_FRAME);
		inferenceGroup->begin();
		{

			minCapacitySlider = new Slider(10,45,(w-20)/2-5,20,"Minimum Capacity");
			minCapacitySlider->precision(0);
			minCapacitySlider->labelsize(12);
			minCapacitySlider->value(100);
			minCapacitySlider->labelfont(iMAP->boldFont);
			minCapacitySlider->labelcolor(FL_WHITE);
			minCapacitySlider->bounds(10,150);
			minCapacitySlider->show();

			minLeafPowerSlider = new Slider(10+(w-20)/2+5,45,(w-20)/2-5,20,"Minimum Leaf Power");
			minLeafPowerSlider->precision(0);
			minLeafPowerSlider->bounds(0,1);
			minLeafPowerSlider->value(0);
			minLeafPowerSlider->callback((Fl_Callback*)updateTreeThresholdsCallback);
			minLeafPowerSlider->labelfont(iMAP->boldFont);
			minLeafPowerSlider->labelcolor(FL_WHITE);
			minLeafPowerSlider->deactivate();
			minLeafPowerSlider->show();

			const double avgStep = 1000.*sqrt(file->averageDx*file->averageDx+file->averageDy*file->averageDy);
			minSideSizeSlider = new Slider(10,85,(w-20)/2-5,20,"Minimum Side Size [nm]");
			minSideSizeSlider->precision(0);
			minSideSizeSlider->labelsize(12);
			minSideSizeSlider->value(avgStep);
			minSideSizeSlider->labelfont(iMAP->boldFont);
			minSideSizeSlider->labelcolor(FL_WHITE);
			minSideSizeSlider->bounds(avgStep,2.0*avgStep);
			minSideSizeSlider->show();

			minPointsSlider = new Slider(10,125,w-20,20,"Minimum Points / Zone");
			minPointsSlider->precision(0);
			minPointsSlider->labelsize(12);
			minPointsSlider->value(20);
			minPointsSlider->bounds(2,100);
			minPointsSlider->labelfont(iMAP->boldFont);
			minPointsSlider->callback((Fl_Callback*)updateTreeThresholdsCallback);
			minPointsSlider->labelcolor(FL_WHITE);
			minPointsSlider->deactivate();
			minPointsSlider->show();

			inferenceModeChoice = new Fl_Choice(10,170,(w-24)/2,20,"Mode");
			inferenceModeChoice->menu(optmizationMenu);
			inferenceModeChoice->align(FL_ALIGN_TOP);
			inferenceModeChoice->textsize(11);
			inferenceModeChoice->textfont(iMAP->boldFont);
			inferenceModeChoice->labelcolor(FL_WHITE);
			inferenceModeChoice->labelfont(iMAP->boldFont);
//			inferenceModeChoice->box(FL_BORDER_FRAME);
			inferenceModeChoice->labelsize(11);
			inferenceModeChoice->color(FL_GRAY,FL_GRAY);
			inferenceModeChoice->callback((Fl_Callback*)optimizationChoiceCallback,(int*)2);
			inferenceModeChoice->deactivate();
			inferenceModeChoice->show();

			localizationPrecisionSlider = new Slider(15+(w-24)/2,170,(w-24)/2,20,"Localization Precision [nm]");
			localizationPrecisionSlider->labelsize(12);
			localizationPrecisionSlider->labelcolor(FL_WHITE);
			localizationPrecisionSlider->labelfont(iMAP->boldFont);
			localizationPrecisionSlider->precision(0);
			localizationPrecisionSlider->value(30);
			localizationPrecisionSlider->bounds(0,100);
			localizationPrecisionSlider->show();

			applyButton = new Fl_Button(10,200,(w-40)/6,30,"Apply");
			applyButton->labelsize(11);
			applyButton->labelfont(iMAP->boldFont);
//			applyButton->color(iMAP->bgColor);
//			applyButton->labelcolor(FL_WHITE);
			applyButton->callback((Fl_Callback*)applyMeshCallback,(int*)2);
			applyButton->show();

			const int buttonWidth = 76;
			resetButton = new Fl_Button(10+buttonWidth+5,200,buttonWidth,30,"Reset");
			resetButton->labelsize(11);
			resetButton->labelfont(iMAP->boldFont);
//			resetButton->color(iMAP->bgColor);
//			resetButton->labelcolor(FL_WHITE);
			resetButton->callback((Fl_Callback*)resetMeshCallback,(int*)2);
			resetButton->show();
			resetButton->deactivate();

			inferButton = new Fl_Button(10+2*buttonWidth+10,200,buttonWidth,30,"Infer");
			inferButton->labelsize(11);
			inferButton->labelfont(iMAP->boldFont);
//			inferButton->color(iMAP->bgColor);
//			inferButton->labelcolor(FL_WHITE);
			inferButton->callback((Fl_Callback*)inferMeshCallback,(int*)2);
			inferButton->show();
			inferButton->deactivate();

			pauseButton = new Fl_Button(10+3*buttonWidth+15,200,buttonWidth,30,"Pause");
			pauseButton->labelsize(14);
//			pauseButton->color(iMAP->bgColor);
			pauseButton->label("@||");
			pauseButton->callback((Fl_Callback*)pauseCalculationButtonCallback,(int*)2);
			pauseButton->deactivate();
			pauseButton->show();

			stopButton = new Fl_Button(10+4*buttonWidth+20,200,buttonWidth,30,"Stop");
			stopButton->labelsize(12);
			stopButton->label("@square");
//			stopButton->color(iMAP->bgColor);
			stopButton->labelcolor(FL_RED);
			stopButton->callback((Fl_Callback*)stopCalculationButtonCallback,(int*)2);
			stopButton->deactivate();
			stopButton->show();

			saveButton = new Fl_Button(10+5*buttonWidth+25,200,buttonWidth,30,"Save");
			saveButton->labelsize(12);
//			stopButton->color(iMAP->bgColor);
			saveButton->labelfont(iMAP->boldFont);
//			saveButton->labelcolor(FL_WHITE);
			saveButton->callback((Fl_Callback*)saveTreeMeshCallback);
			saveButton->deactivate();
			saveButton->show();
		}
		inferenceGroup->end();

		annotationsGroup = new Fl_Group(0,25,w,h-85,"Overlay");
		annotationsGroup->labelsize(12);
//		annotationsGroup->color(FL_WHITE,FL_WHITE);
		annotationsGroup->labelcolor(FL_DARK1);
		annotationsGroup->box(FL_BORDER_FRAME);
		annotationsGroup->labelfont(iMAP->normalFont);
		annotationsGroup->begin();
		{
			// font options
			localizationNumberLabelButton = new Fl_Check_Button(10,45,150,20,"Localization Number");
			localizationNumberLabelButton->labelsize(12);
			localizationNumberLabelButton->labelcolor(FL_WHITE);
			localizationNumberLabelButton->labelfont(iMAP->boldFont);
			localizationNumberLabelButton->value(1);
			localizationNumberLabelButton->show();

			fontScaleSlider = new Slider(175,45,w-185,20,"Font Scale");
			fontScaleSlider->labelsize(12);
			fontScaleSlider->precision(2);
			fontScaleSlider->bounds(0.5f,1.5f);
			fontScaleSlider->value(1.0f);
			fontScaleSlider->labelfont(iMAP->boldFont);
			fontScaleSlider->labelcolor(FL_WHITE);
			fontScaleSlider->show();

			// grid options
			displayGridButton = new Fl_Check_Button(10,85,105,20,"Display Grid");
			displayGridButton->labelsize(12);
			displayGridButton->labelcolor(FL_WHITE);
			displayGridButton->labelfont(iMAP->boldFont);
			displayGridButton->value(1);
			displayGridButton->show();

			gridColorButton = new Fl_Button(120,85,85,20,"Grid Color");
			gridColorButton->labelsize(12);
			gridColorButton->labelfont(iMAP->boldFont);
			gridColorButton->callback((Fl_Callback*)gridColorButtonCallback,(int*)2);
			gridColorButton->show();

			gridRGB[0] = gridRGB[1] = gridRGB[2] = 1.0;

			gridAlphaSlider = new Slider(215,85,w-225,20,"Grid Alpha");
			gridAlphaSlider->labelsize(12);
			gridAlphaSlider->precision(2);
			gridAlphaSlider->labelcolor(FL_WHITE);
			gridAlphaSlider->labelfont(iMAP->boldFont);
			gridAlphaSlider->value(1.0f);
			gridAlphaSlider->bounds(0.0,1.0);
			gridAlphaSlider->show();

			// color overlay options
			overlayButton = new Fl_Check_Button(10,125,115,20,"Color Overlay");
			overlayButton->labelsize(12);
			overlayButton->labelcolor(FL_WHITE);
			overlayButton->labelfont(iMAP->boldFont);
			overlayButton->value(1);
			overlayButton->show();

			overlayColorChoice = new Fl_Choice(155,125,70,20,"Map");
			overlayColorChoice->labelsize(12);
			overlayColorChoice->menu(colormapItemsQuadTree);
			overlayColorChoice->textsize(12);
			overlayColorChoice->labelcolor(FL_WHITE);
			overlayColorChoice->textfont(iMAP->boldFont);
			overlayColorChoice->labelfont(iMAP->boldFont);
			overlayColorChoice->callback((Fl_Callback*)colorOverlayCallback,(int*)2);
			overlayColorChoice->value(colormap); // jet by default
			overlayColorChoice->show();

			flipColormapButton = new Fl_Check_Button(225,125,40,20,"Flip");
			flipColormapButton->callback((Fl_Callback*)flipColormapCallback,(int*)2);
			flipColormapButton->labelsize(12);
			flipColormapButton->labelcolor(FL_WHITE);
			flipColormapButton->labelfont(iMAP->boldFont);
			flipColormapButton->value(0);
			flipColormapButton->show();

			overlayAlphaSlider = new Slider(280,125,w-290,20,"Overlay Alpha");
			overlayAlphaSlider->labelsize(12);
			overlayAlphaSlider->precision(2);
			overlayAlphaSlider->bounds(0.0f,1.0f);
			overlayAlphaSlider->value(0.5f);
			overlayAlphaSlider->labelfont(iMAP->boldFont);
			overlayAlphaSlider->labelcolor(FL_WHITE);
			overlayAlphaSlider->callback((Fl_Callback*)colorOverlayCallback,(int*)2);
			overlayAlphaSlider->show();

			// spot visualization options
			spotVisualizationButton = new Fl_Check_Button(10,165,135,20,"Spot Visualization");
			spotVisualizationButton->labelsize(12);
			spotVisualizationButton->labelcolor(FL_WHITE);
			spotVisualizationButton->labelfont(iMAP->boldFont);
			spotVisualizationButton->value(0);
			spotVisualizationButton->show();

			spotVisualizationScaleSlider = new Slider(160,165,w-170,20,"Scale");
			spotVisualizationScaleSlider->labelsize(12);
			spotVisualizationScaleSlider->precision(2);
			spotVisualizationScaleSlider->bounds(0.0f,5.0f);
			spotVisualizationScaleSlider->labelcolor(FL_WHITE);
			spotVisualizationScaleSlider->labelfont(iMAP->boldFont);
			spotVisualizationScaleSlider->value(1.0f);
			spotVisualizationScaleSlider->show();

			// hover information
			infoOverlayButton = new Fl_Check_Button(10,190,140,20,"Hover Information");
			infoOverlayButton->labelsize(12);
			infoOverlayButton->value(0);
			infoOverlayButton->labelfont(iMAP->boldFont);
			infoOverlayButton->labelcolor(FL_WHITE);
			infoOverlayButton->deactivate();
			infoOverlayButton->show();

			neighbourDistanceViewButton = new Fl_Check_Button(10,210,195,20,"View Neighbor Connections");
			neighbourDistanceViewButton->labelcolor(FL_WHITE);
			neighbourDistanceViewButton->labelfont(iMAP->boldFont);
			neighbourDistanceViewButton->labelsize(12);
			neighbourDistanceViewButton->value(1);
			neighbourDistanceViewButton->show();
		}
		annotationsGroup->deactivate();
		annotationsGroup->end();

		advancedGroup = new Fl_Group(0,25,w,h-85,"Advanced");
		advancedGroup->labelsize(12);
//		advancedGroup->color(FL_WHITE,FL_WHITE);
		advancedGroup->labelcolor(FL_DARK1);
		advancedGroup->box(FL_BORDER_FRAME);
		advancedGroup->labelfont(iMAP->normalFont);
		advancedGroup->begin();
		{

			roButton = new Fl_Check_Button(10,30,180,20,"Randomized Optimization");
			roButton->labelcolor(FL_WHITE);
			roButton->labelfont(iMAP->boldFont);
			roButton->labelsize(12);
			roButton->value(0);
			roButton->callback((Fl_Callback*)roButtonCallback,(int*)2);
			roButton->deactivate();
			roButton->show();

			roRadiusSlider = new Slider(10,65,(w-20)/3-5,20,"Selection Radius [nm]");
			roRadiusSlider->labelsize(12);
			roRadiusSlider->labelcolor(FL_WHITE);
			roRadiusSlider->labelfont(iMAP->boldFont);
			roRadiusSlider->precision(0);
			roRadiusSlider->value(500);
			roRadiusSlider->bounds(100,2000);
			roRadiusSlider->deactivate();
			roRadiusSlider->show();

			roToleranceSlider = new Slider(12+(w-20)/3,65,w/3-10,20,"Cost Tolerance [%]");
			roToleranceSlider->labelsize(12);
			roToleranceSlider->labelcolor(FL_WHITE);
			roToleranceSlider->labelfont(iMAP->boldFont);
			roToleranceSlider->precision(4);
			roToleranceSlider->value(0.001);
			roToleranceSlider->bounds(0.0,1.0);
			roToleranceSlider->callback((Fl_Callback*)randomizedOptimizationToleranceCallback,(int*)2);
			roToleranceSlider->deactivate();
			roToleranceSlider->show();

			roMaximumIterationsSlider = new Slider(14+2*(w-20)/3,65,w/3-10,20,"Maximum Iterations");
			roMaximumIterationsSlider->labelsize(12);
			roMaximumIterationsSlider->labelcolor(FL_WHITE);
			roMaximumIterationsSlider->labelfont(iMAP->boldFont);
			roMaximumIterationsSlider->precision(0);
			roMaximumIterationsSlider->value(5);
			roMaximumIterationsSlider->bounds(3,10);
			roMaximumIterationsSlider->deactivate();
			roMaximumIterationsSlider->show();

			neighbourDistanceSlider = new Slider(10,105,w-20,20,"Maximum Neighbor Distance [nm]");
			neighbourDistanceSlider->labelsize(12);
			neighbourDistanceSlider->labelcolor(FL_WHITE);
			neighbourDistanceSlider->labelfont(iMAP->boldFont);
			neighbourDistanceSlider->precision(1);
			neighbourDistanceSlider->value(1000.0);
			neighbourDistanceSlider->bounds(50,2000);
			neighbourDistanceSlider->callback((Fl_Callback*)maximumNeighbourDistanceCallback,(int*)2);
			neighbourDistanceSlider->show();

			Fl_Box *dfLabel = new Fl_Box(15,145,100,20,"(D,F) Inference:");
			dfLabel->labelsize(12);
			dfLabel->labelcolor(FL_WHITE);
			dfLabel->labelfont(iMAP->boldFont);
			dfLabel->show();

			betaSlider = new Slider(130,145,w-140,20,"Potential Energy Penalization (Beta)");
			betaSlider->precision(1);
			betaSlider->labelcolor(FL_WHITE);
			betaSlider->labelfont(iMAP->boldFont);
			betaSlider->bounds(0.0f,5.0f);
			betaSlider->value(2.0f);
			betaSlider->deactivate();
			betaSlider->show();

			Fl_Box *polynomialLabel = new Fl_Box(10,185,220,20,"Polynomial Potential Inference:");
			polynomialLabel->labelsize(12);
			polynomialLabel->labelcolor(FL_WHITE);
			polynomialLabel->labelfont(iMAP->boldFont);
			polynomialLabel->show();

			polynomialOrderSlider = new Slider(235,185,w-245,20,"Polynomial Order");
			polynomialOrderSlider->precision(0);
			polynomialOrderSlider->labelsize(12);
			polynomialOrderSlider->labelcolor(FL_WHITE);
			polynomialOrderSlider->value(2);
			polynomialOrderSlider->labelfont(iMAP->boldFont);
			polynomialOrderSlider->bounds(2,8);
			polynomialOrderSlider->show();
			polynomialOrderSlider->deactivate();

		}
		advancedGroup->deactivate();
		advancedGroup->end();

		priorGroup = new Fl_Group(0,25,w,h-85,"Prior");
		priorGroup->labelsize(12);
//		priorGroup->color(FL_WHITE,FL_WHITE);
		priorGroup->labelcolor(FL_DARK1);
		priorGroup->box(FL_BORDER_FRAME);
		priorGroup->labelfont(iMAP->normalFont);
		priorGroup->begin();
		{

			jeffreysPriorButton = new Fl_Check_Button(10,30,180,25,"Enable Jeffreys Prior");
			jeffreysPriorButton->labelsize(12);
			jeffreysPriorButton->labelcolor(FL_WHITE);
			jeffreysPriorButton->labelfont(iMAP->boldFont);
			jeffreysPriorButton->value(1);
			jeffreysPriorButton->callback((Fl_Callback*)jeffreysPriorButtonCallback,(int*)2);
			jeffreysPriorButton->show();

			smoothingPriorButton = new Fl_Check_Button(10,50,180,25,"Enable Smoothing Prior");
			smoothingPriorButton->labelsize(12);
			smoothingPriorButton->labelcolor(FL_WHITE);
			smoothingPriorButton->labelfont(iMAP->boldFont);
			smoothingPriorButton->value(0);
			smoothingPriorButton->callback((Fl_Callback*)smoothingPriorButtonCallback,(int*)2);
			smoothingPriorButton->show();
//			smoothingPriorButton->deactivate();

			lambdaSlider = new Slider(10,85,(w-30)/2,20,"V Penalisation (lambda)");
			lambdaSlider->precision(3);
			lambdaSlider->labelcolor(FL_WHITE);
			lambdaSlider->labelfont(iMAP->boldFont);
			lambdaSlider->labelsize(12);
			lambdaSlider->value(0.1);
			lambdaSlider->bounds(0.0,10);
			lambdaSlider->deactivate();
			lambdaSlider->show();

			muSlider = new Slider((w-30)/2+15,85,(w-30)/2,20,"D Penalization (mu)");
			muSlider->precision(3);
			muSlider->labelsize(12);
			muSlider->labelcolor(FL_WHITE);
			muSlider->labelfont(iMAP->boldFont);
			muSlider->value(0.1);
			muSlider->bounds(0.0,10.0);
			muSlider->deactivate();
			muSlider->show();

		}
//		priorGroup->deactivate();
		priorGroup->end();

		posterioriGroup = new Fl_Group(0,25,w,h-85,"Posterior");
		posterioriGroup->labelsize(12);
//		posterioriGroup->color(FL_WHITE,FL_WHITE);
		posterioriGroup->labelcolor(FL_DARK1);
		posterioriGroup->box(FL_BORDER_FRAME);
		posterioriGroup->labelfont(iMAP->normalFont);
		posterioriGroup->begin();
		{
			dPosteriorButton = new Fl_Round_Button(10,30,155,20,"Diffusion");
			dPosteriorButton->labelfont(1);
			dPosteriorButton->labelsize(12);
			dPosteriorButton->labelfont(iMAP->boldFont);
			dPosteriorButton->type(FL_TOGGLE_BUTTON);
			dPosteriorButton->labelcolor(FL_WHITE);
			dPosteriorButton->callback((Fl_Callback*)dPosteriorCallback);
			dPosteriorButton->value(0);
			dPosteriorButton->deactivate();
			dPosteriorButton->show();

			fPosteriorButton = new Fl_Round_Button(10,50,155,20,"Force");
			fPosteriorButton->labelfont(iMAP->boldFont);
			fPosteriorButton->labelsize(12);
			fPosteriorButton->type(FL_TOGGLE_BUTTON);
			fPosteriorButton->labelcolor(FL_WHITE);
			fPosteriorButton->callback((Fl_Callback*)fPosterioriCallback);
			fPosteriorButton->value(0);
			fPosteriorButton->deactivate();
			fPosteriorButton->show();

			vPosteriorButton = new Fl_Round_Button(10,70,155,20,"Potential");
			vPosteriorButton->labelfont(iMAP->boldFont);
			vPosteriorButton->labelsize(12);
			vPosteriorButton->type(FL_TOGGLE_BUTTON);
			vPosteriorButton->labelcolor(FL_WHITE);
			vPosteriorButton->callback((Fl_Callback*)vPosterioriCallback);
			vPosteriorButton->value(0);
			vPosteriorButton->deactivate();
			vPosteriorButton->show();

			potentialReferenceButton = new Fl_Button(95,70,70,20,"Reference");
			potentialReferenceButton->labelfont(iMAP->boldFont);
			potentialReferenceButton->labelsize(12);
			potentialReferenceButton->type(FL_TOGGLE_BUTTON);
			potentialReferenceButton->callback((Fl_Callback*)potentialReferenceCallback);
			potentialReferenceButton->value(0);
			potentialReferenceButton->deactivate();
			potentialReferenceButton->show();

			minPosteriorSlider = new Slider(10,105,155,20,"Minimum");
			minPosteriorSlider->labelsize(12);
			minPosteriorSlider->labelfont(iMAP->boldFont);
			minPosteriorSlider->precision(3);
			minPosteriorSlider->labelcolor(FL_WHITE);
			minPosteriorSlider->value(0.001f);
			minPosteriorSlider->bounds(0.001,1.0);
			minPosteriorSlider->deactivate();
			minPosteriorSlider->show();

			maxPosteriorSlider = new Slider(10,140,155,20,"Maximum");
			maxPosteriorSlider->labelsize(12);
			maxPosteriorSlider->labelfont(iMAP->boldFont);
			maxPosteriorSlider->precision(3);
			maxPosteriorSlider->labelcolor(FL_WHITE);
			maxPosteriorSlider->value(1.0f);
			maxPosteriorSlider->bounds(0.001,1.0);
			maxPosteriorSlider->deactivate();
			maxPosteriorSlider->show();

			posterioriSampleNumberSlider = new Slider(10,175,155,20,"Samples");
			posterioriSampleNumberSlider->labelsize(12);
			posterioriSampleNumberSlider->labelfont(iMAP->boldFont);
			posterioriSampleNumberSlider->precision(0);
			posterioriSampleNumberSlider->labelcolor(FL_WHITE);
			posterioriSampleNumberSlider->bounds(100,1000);
			posterioriSampleNumberSlider->value(100);
			posterioriSampleNumberSlider->deactivate();
			posterioriSampleNumberSlider->show();

			samplePosteriorButton = new Fl_Button(10,205,75,25,"Sample");
			samplePosteriorButton->labelfont(iMAP->boldFont);
			samplePosteriorButton->labelsize(12);
			samplePosteriorButton->callback((Fl_Callback*)samplePosteriorCallback,(int*)2);
//			samplePosterioriButton->deactivate();
			samplePosteriorButton->show();

			savePosteriorButton = new Fl_Button(90,205,75,25,"Save");
			savePosteriorButton->labelfont(iMAP->boldFont);
			savePosteriorButton->labelsize(12);
			savePosteriorButton->callback((Fl_Callback*)savePosteriorCallback,(int*)2);
			savePosteriorButton->deactivate();
			savePosteriorButton->show();

			// draw box
			Fl_Box *posterioriBox = new Fl_Box(170,35,320,195);
			posterioriBox->color(FL_WHITE);

		}
		posterioriGroup->deactivate();
		posterioriGroup->end();

		simulationGroup = new Fl_Group(0,25,w,h-85,"Simulation");
		simulationGroup->labelsize(12);
//		simulationGroup->color(FL_WHITE,FL_WHITE);
		simulationGroup->labelcolor(FL_DARK1);
		simulationGroup->box(FL_BORDER_FRAME);
		simulationGroup->labelfont(iMAP->normalFont);
		simulationGroup->begin();
		{
			numberOfTrajectoriesSlider = new Slider(10,45,w-20,20,"Number of Trajectories");
			numberOfTrajectoriesSlider->precision(0);
			numberOfTrajectoriesSlider->labelfont(iMAP->boldFont);
			numberOfTrajectoriesSlider->labelcolor(FL_WHITE);
			numberOfTrajectoriesSlider->bounds(100,1000);
			numberOfTrajectoriesSlider->value(500);
//			numberOfTrajectoriesSlider->callback((Fl_Callback*)simulationWidgetsCallback,(int*)2);
			numberOfTrajectoriesSlider->activate();

			deltaTSlider = new Slider(10,85,(w-30)/2,20,"Delta [ms]");
			deltaTSlider->labelsize(12);
			deltaTSlider->precision(0);
			deltaTSlider->labelfont(iMAP->boldFont);
			deltaTSlider->labelcolor(FL_WHITE);
			deltaTSlider->value(25);
			deltaTSlider->bounds(0,100);
//			deltaTSlider->callback((Fl_Callback*)simulationWidgetsCallback,(int*)2);
			deltaTSlider->activate();

			timeStepsSlider = new Slider((w-30)/2+20,85,(w-30)/2,20,"Maximum Time Steps");
			timeStepsSlider->labelsize(12);
			timeStepsSlider->labelfont(iMAP->boldFont);
			timeStepsSlider->labelcolor(FL_WHITE);
			timeStepsSlider->precision(0);
			timeStepsSlider->value(25);
			timeStepsSlider->bounds(0,100);
//			timeStepsSlider->callback((Fl_Callback*)simulationWidgetsCallback,(int*)2);
			timeStepsSlider->activate();

			generateTrajectoriesButton = new Fl_Button(w/2-50,115,100,40,"Save\nTrajectories");
			generateTrajectoriesButton->labelsize(12);
			generateTrajectoriesButton->labelfont(iMAP->boldFont);
			generateTrajectoriesButton->callback((Fl_Callback*)generateTrajectoriesCallback,(int*)2);
			generateTrajectoriesButton->show();
			generateTrajectoriesButton->activate();
//
//			saveTrajectoriesButton = new Fl_Button(100,115,80,40,"Save\nTrajectories");
//			saveTrajectoriesButton->labelsize(12);
//			saveTrajectoriesButton->labelfont(iMAP->boldFont);
//			saveTrajectoriesButton->callback((Fl_Callback*)saveTrajectoriesCallback,(int*)0);
//			saveTrajectoriesButton->show();
//			saveTrajectoriesButton->deactivate();
		}
		simulationGroup->deactivate();
		simulationGroup->end();

		landscapeGroup = new Fl_Group(0,25,w,h-85,"Landscape");
		landscapeGroup->labelsize(12);
//		landscapeGroup->color(FL_WHITE,FL_WHITE);
		landscapeGroup->labelcolor(FL_DARK1);
		landscapeGroup->box(FL_BORDER_FRAME);
		landscapeGroup->labelfont(iMAP->normalFont);
		landscapeGroup->begin();
		{
			landscapeButton = new Fl_Check_Button(10,30,120,25,"View Landscape");
			landscapeButton->labelsize(12);
			landscapeButton->labelfont(iMAP->boldFont);
			landscapeButton->labelcolor(FL_WHITE);
			landscapeButton->value(0);
			landscapeButton->callback((Fl_Callback*)quadTreeLandscapeCallback);
//			landscapeButton->deactivate();
			landscapeButton->show();

			scaleLandscapeSlider = new Slider(10,70,w/2-15,20,"Scale");
			scaleLandscapeSlider->precision(2);
			scaleLandscapeSlider->labelfont(iMAP->boldFont);
			scaleLandscapeSlider->labelcolor(FL_WHITE);
			scaleLandscapeSlider->bounds(0,5.0);
			scaleLandscapeSlider->value(1.0);
			scaleLandscapeSlider->deactivate();
			scaleLandscapeSlider->show();

			landscapeAlphaSlider = new Slider(w/2+5,70,w/2-15,20,"Alpha");
			landscapeAlphaSlider->precision(2);
			landscapeAlphaSlider->bounds(0.01,1);
			landscapeAlphaSlider->labelfont(iMAP->boldFont);
			landscapeAlphaSlider->labelcolor(FL_WHITE);
			landscapeAlphaSlider->value(1.0);
			landscapeAlphaSlider->callback((Fl_Callback*)alphaQuadTreeLandscapeCallback);
			landscapeAlphaSlider->deactivate();
			landscapeAlphaSlider->show();

			xLightPositionSlider = new Slider(10,105,w/2-15,20,"Light x");
			xLightPositionSlider->precision(2);
			xLightPositionSlider->bounds(-100,100);
			xLightPositionSlider->labelfont(iMAP->boldFont);
			xLightPositionSlider->labelcolor(FL_WHITE);
			xLightPositionSlider->value(0.0);
			xLightPositionSlider->deactivate();
			xLightPositionSlider->show();

			yLightPositionSlider = new Slider(10,140,w/2-15,20,"Light y");
			yLightPositionSlider->precision(2);
			yLightPositionSlider->bounds(-100,100);
			yLightPositionSlider->labelfont(iMAP->boldFont);
			yLightPositionSlider->labelcolor(FL_WHITE);
			yLightPositionSlider->value(0.0);
			yLightPositionSlider->deactivate();
			yLightPositionSlider->show();

			zLightPositionSlider = new Slider(10,175,w/2-15,20,"Light z");
			zLightPositionSlider->precision(2);
			zLightPositionSlider->labelfont(iMAP->boldFont);
			zLightPositionSlider->labelcolor(FL_WHITE);
			zLightPositionSlider->bounds(-100,100);
			zLightPositionSlider->value(100.0);
			zLightPositionSlider->deactivate();
			zLightPositionSlider->show();

			fogButton = new Fl_Check_Button(w/2+5,105,50,25,"Fog");
			fogButton->labelsize(12);
			fogButton->labelfont(iMAP->boldFont);
			fogButton->labelcolor(FL_WHITE);
			fogButton->value(0);
			fogButton->deactivate();
			fogButton->show();

			fogStartSlider = new Slider(w/2+5,140,w/2-15,20,"Fog Start");
			fogStartSlider->precision(2);
			fogStartSlider->bounds(0,100);
			fogStartSlider->labelfont(iMAP->boldFont);
			fogStartSlider->labelcolor(FL_WHITE);
			fogStartSlider->value(0.0);
			fogStartSlider->deactivate();
			fogStartSlider->show();

			fogEndSlider = new Slider(w/2+5,175,w/2-15,20,"Fog End");
			fogEndSlider->precision(2);
			fogEndSlider->labelfont(iMAP->boldFont);
			fogEndSlider->labelcolor(FL_WHITE);
			fogEndSlider->bounds(0,100);
			fogEndSlider->value(100.0);
			fogEndSlider->deactivate();
			fogEndSlider->show();

			landscapeAxisButton = new Fl_Check_Button(10,200,50,20,"Axes");
			landscapeAxisButton->labelsize(12);
			landscapeAxisButton->labelfont(iMAP->boldFont);
			landscapeAxisButton->labelcolor(FL_WHITE);
			landscapeAxisButton->value(1);
			landscapeAxisButton->callback((Fl_Callback*)axesCallback);
			landscapeAxisButton->deactivate();
			landscapeAxisButton->show();
		}
		landscapeGroup->deactivate();
		landscapeGroup->end();

		tabs->end();
		window->end();

		meshApplied = false;
	}
}

void updateTreeThresholdsCallback(Slider*w,void*) {
	file->treeMesh->totalVariables = 0;
	updateTreeThresholds(file->treeMesh->quadTree);
	file->treeMesh->setMinNumberOfPointsPerCell((int)file->treeMeshGui->minPointsSlider->value());
	file->treeMesh->setMinLeafPower((int)file->treeMeshGui->minLeafPowerSlider->value());
	file->treeMesh->setMinCapacity((int)file->treeMeshGui->minCapacitySlider->value());
	iMAP->overlayAdjustment = true;
	file->treeMeshGui->updateVariables(file->treeMesh->totalVariables);
}

void updateTreeThresholds(QuadTree *tree) {

	// only count points in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {

		if (tree->getCount() > 1 && tree->getCount() >= file->treeMeshGui->getMinPoints() && tree->power >= file->treeMeshGui->getMinLeafPower()) {
			tree->activate();
			if (file->treeMesh->maxCount < tree->getCount()) {
				file->treeMesh->maxCount = tree->getCount();
				file->treeMesh->setMaxTree(tree);
			}
			if (file->treeMesh->minCount > tree->getCount()) {
				file->treeMesh->minCount = tree->getCount();
				file->treeMesh->setMinTree(tree);
			}
			file->treeMesh->totalVariables++;
		} else {
			tree->deactivate();
		}

		return;
	}
	updateTreeThresholds(tree->nw);
	updateTreeThresholds(tree->ne);
	updateTreeThresholds(tree->sw);
	updateTreeThresholds(tree->se);
	return;
}

void treeMeshWindowCallback(Fl_Widget *w,void*v) {
	file->treeMeshOverlay = false;
	file->treeMeshGui->resetButton->do_callback();
	file->treeMeshGui->hide();
	file->treeMeshGui->saveButton->deactivate();
	file->inferred = false;
	file->landscapePlot = false;
	iMAP->customSelectionInferenceGui->activate();
	if (file->treeMeshOverlay) {
		if ( file->treeMesh->randomizedOptimizationGui != NULL) {
			file->treeMesh->randomizedOptimizationGui->hide();
	//		delete file->treeMesh->randomizedOptimizationGui;
		}
	}
	delete file->treeMeshGui;
	file->treeMeshGui = NULL;
	iMAP->overlayAdjustment = true;
}

void generateTrajectoriesCallback(Fl_Widget*,int*v) {

	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {
	case 0:
		if (file->simulation != NULL) { delete file->simulation; }
		file->simulation = new Simulation(file->squareMesh);
		file->simulation->createTrajectories(file->squareMesh);
		file->simulation->saveTrajectories();
		break;
	case 1:
		if (file->simulation != NULL) { delete file->simulation; }
		file->simulation = new Simulation(file->voronoiMesh);
		file->simulation->createTrajectories(file->voronoiMesh);
		file->simulation->saveTrajectories();
		break;
	case 2:
		if (file->simulation != NULL) { delete file->simulation; }
		file->simulation = new Simulation(file->treeMesh);
		file->simulation->createTrajectories(file->treeMesh);
		file->simulation->saveTrajectories();
		break;
	}

}

void saveTrajectoriesCallback(Fl_Widget*,int*v) {

	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {
	case 0:
		file->simulation->saveTrajectories();
		break;
	case 1:
		file->simulation->saveTrajectories();
		break;
	case 2:
		file->simulation->saveTrajectories();
		break;
	}

}

void selectionButtonCallback(Fl_Button*w, void*) {
	if (w->value()) {
		iMAP->selectionButtonPressed = true;
		iMAP->customSelectionInferenceGui->inferDFButton->activate();
		iMAP->customSelectionInferenceGui->inferDDrButton->activate();
	} else {
		iMAP->selectionButtonPressed = false;
		iMAP->customSelectionInferenceGui->inferDFButton->deactivate();
		iMAP->customSelectionInferenceGui->inferDDrButton->deactivate();
	}
}

void customSelectionInferDFButtonCallback(Fl_Button*w,void*) {
	file->meshType = -1;
	file->updateIntervalData();
	file->selection->inferDF();
	file->reloadOriginalData();
	iMAP->customSelectionInferenceGui->DFinferred = true;
	iMAP->customSelectionInferenceGui->DDrinferred = false;
}

void customSelectionInferDDrButtonCallback(Fl_Button*w,void*) {
	file->meshType = -1;
	file->updateIntervalData();
	file->selection->inferDDr();
	file->reloadOriginalData();
	iMAP->customSelectionInferenceGui->DFinferred = false;
	iMAP->customSelectionInferenceGui->DDrinferred = true;
}

void customSelectionInferenceWindowCallback(Fl_Window*w,void*) {
	iMAP->selectionMade = false;
	iMAP->customSelectionInferenceGui->selectionButton->value(0);
	iMAP->selectionButtonPressed = false;
}

void customSelectionSaveCallback(Fl_Button*,void*) {
	iMAP->customSelectionInferenceGui->save();
}

void CustomSelectionInferenceGui::disableSelection() {
	this->selectionButton->value(0);
	iMAP->selectionButtonPressed = false;
}

void CustomSelectionInferenceGui::save() {

	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Custom Selection Inference Data");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Custom Selection File\t*.cstraj\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".straj");
		FILE * writeFile;

		writeFile = fopen (nativeFilename,"wb");

		// write header
		fprintf(writeFile,"File Name: %s\n\n"
						  "Number of Points: %i\n"
						  "Area = %f [um^2]\n"
						  "Perimeter %f [um]\n"
						  "Noise Sigma: %.0f [nm]\n"
						  "Diffusion Coefficient: %f [um^2/s]\n"
						  "x-Force: %f [pN]\n"
						  "y-Force: %f [pN]\n"
						  "Force Magnitude: %f [pN]\n",
						  file->fileName,
						  file->selection->cell.count,
						  file->selection->cell.area,
						  file->selection->cell.perimeter,
						  noiseSigmaSlider->value(),
						  optimizationArray[0],
						  optimizationArray[1],
						  optimizationArray[2],
						  sqrt(optimizationArray[1]*optimizationArray[1] + optimizationArray[2]*optimizationArray[2]));

//		fprintf(writeFile,"Vertices\n");
//		for (int v = 0; v < file->selection->n; v++) {
//			fprintf(writeFile,"%f\t%f\n",file->selection->vertices[2*v],file->selection->vertices[2*v+1]);
//		}

		fclose(writeFile);
	}
}

void customSelectionSampleDiffusionCallback(Fl_Button*w,void*) {
	iMAP->customSelectionInferenceGui->sampleDiffusionPosteriori();
}

void CustomSelectionInferenceGui::sampleDiffusionPosteriori() {
	diffusionSamples = (int) diffusionSampleNumberSlider->value();
	const double minD = minDiffusionSlider->value();
	const double maxD = maxDiffusionSlider->value();
	diffusionIncrement = (maxD-minD)/diffusionSampleNumberSlider->value();
	diffusionRange = maxD-minD;
	double logMax = -10000000.0;
	if (diffusionSamples > 5) {
		diffusionValues = new double[diffusionSamples];
		double sampler[3];
		sampler[1] = optimizationArray[1];
		sampler[2] = optimizationArray[2];
		// calculate diffusion values for each of sample positions
		for (int g = 0; g < diffusionSamples; g++) {
			sampler[0] = minD + (double)g*diffusionIncrement;
			diffusionValues[g] = -customSelectionDFPosterior(sampler);
			if (logMax < diffusionValues[g]) { logMax = diffusionValues[g]; }
		}
		// calculate log values
		for (int h = 0; h < diffusionSamples; h++) {
			diffusionValues[h] = exp(diffusionValues[h]-logMax);
		}
		drawDiffusionPosterioriPlot();
		saveDPosterioriButton->activate();
	}
}

void customSelectionSampleDiffusionSaveCallback(Fl_Button*w,void*) {
	iMAP->customSelectionInferenceGui->saveDiffusionPosteriori();
}

void CustomSelectionInferenceGui::drawDiffusionPosterioriPlot() {
	minSampledDiffusion = minDiffusionSlider->value();
	maxSampledDiffusion = maxDiffusionSlider->value();

	diffusionRange = maxSampledDiffusion - minSampledDiffusion;

	if (diffusionRange > 0.0 && valuesSampled) {
//		window->make_current();

		const int w = 270;
		const int h = 165;

		fl_color(FL_WHITE); fl_rectf(10,35,w,h);
		fl_color(FL_BLACK); fl_rect(10,35,w,h);

		fl_line(10,163,279,163);
		fl_font(1,12); fl_draw("[um\262/s]",14,50);

        int pCurrent,pOld;
        const double spacing = ((double)(w-1)/(double)diffusionSamples);
        char xLabel [10];

        for (int q = 1; q < diffusionSamples; q++) {
        	// draw data points and connectors
        	pOld = (int)( (diffusionValues[q-1])*((double)h-60.0) );
        	pCurrent = (int)( (diffusionValues[q])*((double)h-60.0) );

        	fl_color(FL_RED);
        	fl_line_style(FL_SOLID, 2, 0);
			fl_line(11+(q-1)*spacing,h-3-pOld,11+q*spacing,h-3-pCurrent);

        }
        fl_color(FL_BLACK);
        fl_font(0,12);
        for (int t = 0; t < 5; t++) {
			sprintf(xLabel,"%.3f",minSampledDiffusion+t*(diffusionSamples/5.0)*diffusionIncrement);
			fl_draw(90,xLabel,20+t*w/5,196);
        }
		sprintf(xLabel,"%.3f",minSampledDiffusion+(diffusionSamples)*diffusionIncrement);
		fl_draw(90,xLabel,w+5,196);
	}
}

void CustomSelectionInferenceGui::saveDiffusionPosteriori() {

	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Diffusion Posteriori Data");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Diffusion Posteriori\t*.dpost\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".dpost");
		FILE * writeFile;

		writeFile = fopen (nativeFilename,"wb");

		// write header
//		fprintf(writeFile,"Diffusion\tPosteriori\n");

		// write data
        for (int q = 0; q < diffusionSamples; q++) {
        	fprintf(writeFile,"%f\t%f\n",
        			minSampledDiffusion+q*diffusionIncrement,
        			diffusionValues[q]);
        }
		fclose(writeFile);
	}

}

void customSelectionSampleFxCallback(Fl_Button*w,void*) {
	iMAP->customSelectionInferenceGui->sampleFxPosteriori();
}

void CustomSelectionInferenceGui::sampleFxPosteriori() {
	fxSamples = (int) fxSampleNumberSlider->value();
	const double minFx = minFxSlider->value();
	const double maxFx = maxFxSlider->value();
	fxIncrement = (maxFx-minFx)/fxSampleNumberSlider->value();
	fxRange = maxFx-minFx;
	double logMax = -10000000.0;
	if (fxSamples > 5) {
		fxValues = new double[fxSamples];
		double sampler[3];
		sampler[0] = optimizationArray[0];
		sampler[2] = optimizationArray[2];
		// calculate diffusion values for each of sample positions
		for (int g = 0; g < fxSamples; g++) {
			sampler[1] = minFx + (double)g*fxIncrement;
			fxValues[g] = -customSelectionDFPosterior(sampler);
			if (logMax < fxValues[g]) { logMax = fxValues[g]; }
		}
		// calculate log values
		for (int h = 0; h < fxSamples; h++) {
			fxValues[h] = exp(fxValues[h]-logMax);
		}
		drawFxPosterioriPlot();
		saveFxPosterioriButton->activate();
	}
}

void CustomSelectionInferenceGui::drawFxPosterioriPlot() {
	minSampledFx= minFxSlider->value();
	maxSampledFx= maxFxSlider->value();

	fxRange = maxSampledFx - minSampledFx;

	if (fxRange > 0.0 && valuesSampled) {
//		window->make_current();

		const int w = 270;
		const int h = 165;

		fl_color(FL_WHITE); fl_rectf(10,35,w,h);
		fl_color(FL_BLACK); fl_rect(10,35,w,h);

		fl_line(10,163,279,163);
		fl_font(1,12); fl_draw("[pN]",14,50);

        int pCurrent,pOld;
        const double spacing = ((double)(w-1)/(double)fxSamples);
        char xLabel [10];

        for (int q = 1; q < fxSamples; q++) {
        	// draw data points and connectors
        	pOld = (int)( (fxValues[q-1])*((double)h-60.0));
        	pCurrent = (int)( (fxValues[q])*((double)h-60.0));

        	fl_color(FL_RED);
        	fl_line_style(FL_SOLID, 2, 0);
			fl_line(11+(q-1)*spacing,h-3-pOld,11+q*spacing,h-3-pCurrent);

        }
        fl_color(FL_BLACK);
        fl_font(0,12);
        for (int t = 0; t < 5; t++) {
			sprintf(xLabel,"%.2f",minSampledFx+t*(fxSamples/5.0)*fxIncrement);
			fl_draw(90,xLabel,20+t*w/5,196);
        }
		sprintf(xLabel,"%.2f",minSampledFx+(fxSamples)*fxIncrement);
		fl_draw(90,xLabel,w+5,196);
	}
}

void customSelectionSampleFxSaveCallback(Fl_Button*w,void*) {
	iMAP->customSelectionInferenceGui->saveFxPosteriori();
}

void CustomSelectionInferenceGui::saveFxPosteriori() {

	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Fx Posteriori Data");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Fx Posteriori\t*.fxpost\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".fxpost");
		FILE * writeFile;

		writeFile = fopen (nativeFilename,"wb");

		// write header
//		fprintf(writeFile,"Fx\t\tPosteriori\n");

		// write data
        for (int q = 0; q < fxSamples; q++) {
        	fprintf(writeFile,"%f\t%f\n",
        			minSampledFx+q*fxIncrement,
        			fxValues[q]);
        }
		fclose(writeFile);
	}

}

void customSelectionSampleFyCallback(Fl_Button*w,void*) {
	iMAP->customSelectionInferenceGui->sampleFyPosteriori();
}

void CustomSelectionInferenceGui::sampleFyPosteriori() {
	fySamples = (int) fySampleNumberSlider->value();
	const double minFy = minFySlider->value();
	const double maxFy = maxFySlider->value();
	fyIncrement = (maxFy-minFy)/fySampleNumberSlider->value();
	fyRange = maxFy-minFy;
	double logMax = -10000000.0;
	if (fySamples > 5) {
		fyValues = new double[fySamples];
		double sampler[3];
		sampler[0] = optimizationArray[0];
		sampler[1] = optimizationArray[1];
		// calculate diffusion values for each of sample positions
		for (int g = 0; g < fySamples; g++) {
			sampler[2] = minFy + (double)g*fyIncrement;
			fyValues[g] = -customSelectionDFPosterior(sampler);
			if (logMax < fyValues[g]) { logMax = fyValues[g]; }
		}
		// calculate log values
		for (int h = 0; h < fySamples; h++) {
			fyValues[h] = exp(fyValues[h]-logMax);
		}
		drawFyPosterioriPlot();
		saveFyPosterioriButton->activate();
	}
}

void CustomSelectionInferenceGui::drawFyPosterioriPlot() {
	minSampledFy= minFySlider->value();
	maxSampledFy= maxFySlider->value();

	fyRange = maxSampledFy - minSampledFy;

	if (fyRange > 0.0 && valuesSampled) {
//		window->make_current();

		const int w = 270;
		const int h = 165;

		fl_color(FL_WHITE); fl_rectf(10,35,w,h);
		fl_color(FL_BLACK); fl_rect(10,35,w,h);

		fl_line(10,163,279,163);
		fl_font(1,12); fl_draw("[pN]",14,50);

        int pCurrent,pOld;
        const double spacing = ((double)(w-1)/(double)fySamples);
        char xLabel [10];

        for (int q = 1; q < fySamples; q++) {
        	// draw data points and connectors
        	pOld = (int)( (fyValues[q-1])*((double)h-60.0));
        	pCurrent = (int)( (fyValues[q])*((double)h-60.0));

        	fl_color(FL_RED);
        	fl_line_style(FL_SOLID, 2, 0);
			fl_line(11+(q-1)*spacing,h-3-pOld,11+q*spacing,h-3-pCurrent);

        }
        fl_color(FL_BLACK);
        fl_font(0,12);
        for (int t = 0; t < 5; t++) {
			sprintf(xLabel,"%.2f",minSampledFy+t*(fySamples/5.0)*fyIncrement);
			fl_draw(90,xLabel,20+t*w/5,196);
        }
		sprintf(xLabel,"%.2f",minSampledFy+(fySamples)*fyIncrement);
		fl_draw(90,xLabel,w+5,196);
	}
}

void customSelectionSampleFySaveCallback(Fl_Button*w,void*) {
	iMAP->customSelectionInferenceGui->saveFyPosteriori();
}

void CustomSelectionInferenceGui::saveFyPosteriori() {

	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Fy Posteriori Data");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Fy Posteriori\t*.fypost\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".fypost");
		FILE * writeFile;

		writeFile = fopen (nativeFilename,"wb");

		// write data
        for (int q = 0; q < fySamples; q++) {
        	fprintf(writeFile,"%f\t%f\n",
        			minSampledFx+q*fyIncrement,
        			fyValues[q]);
        }
		fclose(writeFile);
	}

}

void gridColorButtonCallback(Fl_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch (meshType) {
		case 0:
	        fl_color_chooser("Color:",
	                         file->squareMeshGui->gridRGB[0],
	                         file->squareMeshGui->gridRGB[1],
	                         file->squareMeshGui->gridRGB[2],
	                         0);
			break;
		case 1:
	        fl_color_chooser("Color:",
	                         file->voronoiMeshGui->gridRGB[0],
	                         file->voronoiMeshGui->gridRGB[1],
	                         file->voronoiMeshGui->gridRGB[2],
	                         0);
			break;
		case 2:
	        fl_color_chooser("Color:",
	                         file->treeMeshGui->gridRGB[0],
	                         file->treeMeshGui->gridRGB[1],
	                         file->treeMeshGui->gridRGB[2],
	                         0);
			break;
	}
	Fl::redraw();
}

void simulationWidgetsCallback(Fl_Widget*w,int*v) {

	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch (meshType) {
		case 0:
			file->squareMeshGui->saveTrajectoriesButton->deactivate();
			break;
		case 1:
			file->voronoiMeshGui->saveTrajectoriesButton->deactivate();
			break;
		case 2:
			file->treeMeshGui->saveTrajectoriesButton->deactivate();
			break;
	}
}

void randomizedOptimizationToleranceCallback(Slider*w,int*v) {
	file->randomizedOptimizationTolerance = (float)w->value()/100.0;
	Fl::redraw();
}

void stopCalculationButtonCallback(Fl_Button*w,int*v) {
	iMAP->stopCalculation = true;
	iMAP->pauseCalculation = false;

	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	w->deactivate();

	switch (meshType) {
		case 0:
			file->squareMeshGui->inferButton->activate();
			file->squareMeshGui->pauseButton->label("@||");
			file->squareMeshGui->pauseButton->deactivate();
			file->squareMeshGui->stopButton->activate();
			file->squareMeshGui->resetButton->activate();
			file->squareMesh->offsetPotentials();
			break;
		case 1:
			file->voronoiMeshGui->inferButton->activate();
			file->voronoiMeshGui->pauseButton->label("@||");
			file->voronoiMeshGui->pauseButton->deactivate();
			file->voronoiMeshGui->stopButton->activate();
			file->voronoiMeshGui->resetButton->activate();
			file->voronoiMesh->offsetPotentials();
			break;
		case 2:
			file->treeMeshGui->inferButton->activate();
			file->treeMeshGui->pauseButton->label("@||");
			file->treeMeshGui->pauseButton->deactivate();
			file->treeMeshGui->stopButton->activate();
			file->treeMeshGui->resetButton->activate();
			file->treeMesh->offsetPotentials();
			break;
	}

	w->deactivate();

	Fl::redraw();
}

void axesCallback(Fl_Check_Button*w,void*) {
	if (w->value()) { iMAP->landscapeAxesEnable = true; }
	else { iMAP->landscapeAxesEnable = false; }
}

void pauseCalculationButtonCallback(Fl_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	if (iMAP->pauseCalculation == true) {
		w->label("@||");
		switch(meshType) {
		case 0:
			file->squareMeshGui->inferButton->do_callback();
			file->squareMeshGui->saveButton->activate();
			file->squareMeshGui->stopButton->activate();
			break;
		case 1:
			file->voronoiMeshGui->inferButton->do_callback();
			file->voronoiMeshGui->saveButton->activate();
			file->voronoiMeshGui->stopButton->activate();
			break;
		case 2:
			file->treeMeshGui->inferButton->do_callback();
			file->treeMeshGui->saveButton->activate();
			file->treeMeshGui->stopButton->activate();
			break;
		}
		iMAP->pauseCalculation = false;
		return;
	}

	switch(meshType) {
	case 0:
		file->squareMeshGui->inferButton->deactivate();
		file->squareMeshGui->saveButton->activate();
		file->squareMeshGui->stopButton->activate();
		break;
	case 1:
		file->voronoiMeshGui->inferButton->deactivate();
		file->voronoiMeshGui->saveButton->activate();
		file->voronoiMeshGui->stopButton->activate();
		break;
	case 2:
		file->treeMeshGui->inferButton->deactivate();
		file->treeMeshGui->saveButton->activate();
		file->treeMeshGui->stopButton->activate();
		break;
	}

	if (iMAP->pauseCalculation == false) {
		w->label("@>");
		iMAP->pauseCalculation = true;
	}
	Fl::redraw();
}

void customSelectionNoiseSigmaSliderCallback(Fl_Widget*w,void*) {
	Fl::redraw();
}

void inferMeshCallback(Fl_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	// clear text buffer
	iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
	iMAP->updateDisplayButton->value(1);
	iMAP->updateDisplayButton->do_callback();

	w->deactivate();

	switch(meshType) {
	case 0: // square mesh
		file->squareMeshGui->resetButton->deactivate();
		file->squareMeshGui->applyButton->deactivate();
		file->squareMeshGui->pauseButton->activate();
		file->squareMeshGui->stopButton->deactivate();
		file->squareMeshGui->saveButton->deactivate();
		//file->squareMeshGui->minPointsSlider->do_callback();
		file->squareMeshGui->neighbourDistanceViewButton->value(0);
		w->deactivate();
		file->squareMeshGui->landscapeButton->value(0);
		file->squareMeshGui->landscapeButton->do_callback();
		file->squareMeshGui->posterioriGroup->activate();
		file->squareMesh->infer();
		file->squareMeshGui->infoOverlayButton->activate();
		file->squareMeshGui->simulationGroup->activate();
		file->squareMeshGui->resetButton->activate();
		file->squareMeshGui->applyButton->deactivate();
		if (iMAP->pauseCalculation == true) {
			file->squareMeshGui->resetButton->deactivate();
			file->squareMeshGui->inferButton->deactivate();
			file->squareMeshGui->stopButton->activate();
		} else {
			w->activate();
		}
		break;
	case 1: // voronoi mesh
		file->voronoiMeshGui->resetButton->deactivate();
		file->voronoiMeshGui->applyButton->deactivate();
		file->voronoiMeshGui->pauseButton->activate();
		file->voronoiMeshGui->stopButton->deactivate();
		file->voronoiMeshGui->saveButton->deactivate();
		file->voronoiMeshGui->posterioriGroup->activate();
		file->voronoiMeshGui->neighbourDistanceViewButton->value(0);
		//file->voronoiMeshGui->minPointsSlider->do_callback();
		w->deactivate();
		file->voronoiMeshGui->landscapeButton->value(0);
		file->voronoiMeshGui->landscapeButton->do_callback();
		file->voronoiMesh->infer();
		file->voronoiMeshGui->infoOverlayButton->activate();
		file->voronoiMeshGui->simulationGroup->activate();
		file->voronoiMeshGui->resetButton->activate();
		file->voronoiMeshGui->applyButton->deactivate();
		if (iMAP->pauseCalculation == true) {
			file->voronoiMeshGui->resetButton->deactivate();
			file->voronoiMeshGui->inferButton->deactivate();
			file->voronoiMeshGui->stopButton->activate();
		} else {
			w->activate();
		}
		break;
	case 2: // tree mesh
		file->treeMeshGui->resetButton->deactivate();
		file->treeMeshGui->applyButton->deactivate();
		file->treeMeshGui->pauseButton->activate();
		file->treeMeshGui->stopButton->deactivate();
		file->treeMeshGui->saveButton->deactivate();
		file->treeMeshGui->neighbourDistanceViewButton->value(0);
		//file->treeMeshGui->minPointsSlider->do_callback();
		w->deactivate();
		file->treeMeshGui->landscapeButton->value(0);
		file->treeMeshGui->landscapeButton->do_callback();
		file->treeMesh->infer();
		file->treeMeshGui->posterioriGroup->activate();
		file->treeMeshGui->infoOverlayButton->activate();
		file->treeMeshGui->simulationGroup->activate();
		file->treeMeshGui->resetButton->activate();
		file->treeMeshGui->applyButton->deactivate();
		if (iMAP->pauseCalculation == true) {
			file->treeMeshGui->resetButton->deactivate();
			file->treeMeshGui->inferButton->deactivate();
			file->treeMeshGui->stopButton->activate();
		} else {
			w->activate();
		}
		break;
	}
	Fl::redraw();
}

void resetMeshCallback(Fl_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	w->deactivate();

	switch(meshType) {
	case 0: // square mesh
		file->squareMeshGui->gridRGB[0] = 1.0;
		file->squareMeshGui->gridRGB[1] = 1.0;
		file->squareMeshGui->gridRGB[2] = 1.0;
		file->squareMeshGui->localizationNumberLabelButton->value(1);
		file->unapplySquareMesh();
		file->squareMeshGui->resetSliders();
		file->squareMeshGui->applyButton->activate();
		file->squareMeshGui->inferButton->deactivate();
		file->squareMeshGui->pauseButton->deactivate();
		file->squareMeshGui->saveButton->deactivate();
		file->squareMeshGui->stopButton->deactivate();
		file->squareMeshGui->inferenceModeChoice->deactivate();
		file->squareMeshGui->polynomialOrderSlider->activate();
		file->squareMeshGui->overlayForceArrowsButton->deactivate();
		file->squareMeshGui->overlayForceArrowsButton->labelcolor(iMAP->bgColor);
		file->squareMeshGui->overlayDiffusionButton->deactivate();
		file->squareMeshGui->overlayDiffusionButton->labelcolor(iMAP->bgColor);
		file->squareMeshGui->overlayDiffusionLogButton->deactivate();
		file->squareMeshGui->overlayDiffusionLogButton->labelcolor(iMAP->bgColor);
		file->squareMeshGui->overlayPotentialButton->deactivate();
		file->squareMeshGui->overlayPotentialButton->labelcolor(iMAP->bgColor);
		file->squareMeshGui->overlayForceMagnitudeButton->deactivate();
		file->squareMeshGui->overlayForceMagnitudeButton->labelcolor(iMAP->bgColor);
		file->squareMeshGui->landscapeButton->deactivate();
		file->squareMeshGui->pauseButton->label("@||");
		file->squareMeshGui->infoOverlayButton->value(0);
		file->squareMeshGui->infoOverlayButton->deactivate();
		file->squareMeshGui->roButton->value(0);
		file->squareMeshGui->roButton->do_callback();
		file->squareMeshGui->neighbourDistanceViewButton->value(1);
		file->squareMeshGui->inferenceModeChoice->deactivate();
		file->squareMeshGui->jeffreysPriorButton->value(1);
		file->squareMeshGui->jeffreysPriorButton->do_callback();
		file->squareMeshGui->smoothingPriorButton->value(0);
		file->squareMeshGui->smoothingPriorButton->do_callback();
		file->squareMeshGui->landscapeButton->value(0);
		file->squareMeshGui->landscapeButton->do_callback();
		iMAP->overlayAdjustment = true;
		break;
	case 1: // voronoi mesh
		file->voronoiMeshOverlay = false;
		file->voronoiMeshGui->localizationNumberLabelButton->value(1);
		file->voronoiMeshGui->gridRGB[0] = 1.0;
		file->voronoiMeshGui->gridRGB[1] = 1.0;
		file->voronoiMeshGui->gridRGB[2] = 1.0;
		file->unapplyVoronoiMesh();
		file->voronoiMeshGui->resetSliders();
		file->voronoiMeshGui->applyButton->activate();
		file->voronoiMeshGui->inferButton->deactivate();
		file->voronoiMeshGui->saveButton->deactivate();
		file->voronoiMeshGui->stopButton->deactivate();
		file->voronoiMeshGui->inferenceModeChoice->deactivate();
		file->voronoiMeshGui->polynomialOrderSlider->activate();
		file->voronoiMeshGui->distanceButton->activate();
		file->voronoiMeshGui->clusteringChoice->activate();
		file->voronoiMeshGui->cellsSlider->activate();
		file->voronoiMeshGui->landscapeButton->deactivate();
		file->voronoiMeshGui->overlayForceArrowsButton->deactivate();
		file->voronoiMeshGui->overlayForceArrowsButton->labelcolor(iMAP->bgColor);
		file->voronoiMeshGui->overlayDiffusionButton->deactivate();
		file->voronoiMeshGui->overlayDiffusionButton->labelcolor(iMAP->bgColor);
		file->voronoiMeshGui->overlayDiffusionLogButton->deactivate();
		file->voronoiMeshGui->overlayDiffusionLogButton->labelcolor(iMAP->bgColor);
		file->voronoiMeshGui->overlayPotentialButton->deactivate();
		file->voronoiMeshGui->overlayPotentialButton->labelcolor(iMAP->bgColor);
		file->voronoiMeshGui->overlayForceMagnitudeButton->deactivate();
		file->voronoiMeshGui->overlayForceMagnitudeButton->labelcolor(iMAP->bgColor);
		file->voronoiMeshGui->pauseButton->label("@||");
		file->voronoiMeshGui->infoOverlayButton->value(0);
		file->voronoiMeshGui->infoOverlayButton->deactivate();
		file->voronoiMeshGui->roButton->value(0);
		file->voronoiMeshGui->roButton->do_callback();
		file->voronoiMeshGui->neighbourDistanceViewButton->value(1);
		file->voronoiMeshGui->inferenceModeChoice->deactivate();
		file->voronoiMeshGui->jeffreysPriorButton->value(1);
		file->voronoiMeshGui->jeffreysPriorButton->do_callback();
		file->voronoiMeshGui->smoothingPriorButton->value(0);
		file->voronoiMeshGui->smoothingPriorButton->do_callback();

		file->voronoiMeshGui->landscapeButton->value(0);
		file->voronoiMeshGui->landscapeButton->do_callback();

		file->voronoiMeshGui->maxIterationsSlider->activate();
		break;
	case 2: // tree mesh
		file->treeMeshOverlay = false;
		file->unapplyTreeMesh();
		file->treeMeshGui->resetSliders();
		file->treeMeshGui->localizationNumberLabelButton->value(1);
		file->treeMeshGui->gridRGB[0] = 1.0;
		file->treeMeshGui->gridRGB[1] = 1.0;
		file->treeMeshGui->gridRGB[2] = 1.0;
		file->treeMeshGui->applyButton->activate();
		file->treeMeshGui->inferButton->deactivate();
		file->treeMeshGui->saveButton->deactivate();
		file->treeMeshGui->stopButton->deactivate();
		file->treeMeshGui->inferenceModeChoice->deactivate();
		file->treeMeshGui->polynomialOrderSlider->activate();
		file->treeMeshGui->overlayForceArrowsButton->deactivate();
		file->treeMeshGui->overlayForceArrowsButton->labelcolor(iMAP->bgColor);
		file->treeMeshGui->overlayDiffusionButton->deactivate();
		file->treeMeshGui->overlayDiffusionButton->labelcolor(iMAP->bgColor);
		file->treeMeshGui->overlayDiffusionLogButton->deactivate();
		file->treeMeshGui->overlayDiffusionLogButton->labelcolor(iMAP->bgColor);
		file->treeMeshGui->overlayPotentialButton->deactivate();
		file->treeMeshGui->overlayPotentialButton->labelcolor(iMAP->bgColor);
		file->treeMeshGui->overlayForceMagnitudeButton->deactivate();
		file->treeMeshGui->overlayForceMagnitudeButton->labelcolor(iMAP->bgColor);
		file->treeMeshGui->minPointsSlider->activate();
		file->treeMeshGui->landscapeButton->value(0);
		file->treeMeshGui->landscapeButton->deactivate();
		file->treeMeshOverlay = false;
		file->treeMeshGui->minLeafPowerSlider->deactivate();
		file->treeMeshGui->pauseButton->label("@||");
		file->treeMeshGui->infoOverlayButton->value(0);
		file->treeMeshGui->infoOverlayButton->deactivate();
		file->treeMeshGui->minSideSizeSlider->activate();
		file->treeMeshGui->roButton->value(0);
		file->treeMeshGui->roButton->do_callback();
		file->treeMeshGui->neighbourDistanceViewButton->value(1);
		file->treeMeshGui->inferenceModeChoice->deactivate();
		file->treeMeshGui->jeffreysPriorButton->value(1);
		file->treeMeshGui->jeffreysPriorButton->do_callback();
		file->treeMeshGui->smoothingPriorButton->value(0);
		file->treeMeshGui->smoothingPriorButton->do_callback();

		file->treeMeshGui->landscapeButton->value(0);
		file->treeMeshGui->landscapeButton->do_callback();

		break;
	}

	iMAP->startIntervalSlider->activate();
	iMAP->endIntervalSlider->activate();

	iMAP->customSelectionInferenceGui->activate();
	file->inferred = false;

	file->reloadOriginalData();

	w->activate();

	Fl::redraw();
}

void applyMeshCallback(Fl_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	if (file->updateIntervalData()) {

		w->deactivate();

		if (iMAP->selectionButtonPressed && file->selection == NULL) {
			iMAP->customSelectionInferenceGui->selectionButton->value(0);
			iMAP->customSelectionInferenceGui->selectionButton->do_callback();
			iMAP->selectionButtonPressed = false;
		}

		switch(meshType) {
		case 0: // square mesh
			file->squareMeshGui->inferButton->deactivate();
			file->squareMeshGui->resetButton->deactivate();
			file->squareMeshGui->inferenceModeChoice->activate();
			if (file->squareMesh != NULL) { delete file->squareMesh; }
			if (!iMAP->selectionButtonPressed) {
				file->squareMesh = new SquareMesh((double)file->squareMeshGui->getDx());
				file->squareMeshGui->gridRGB[0] = 1.0;
				file->squareMeshGui->gridRGB[1] = 1.0;
				file->squareMeshGui->gridRGB[2] = 1.0;
				file->squareMeshGui->localizationNumberLabelButton->value(0);
			} else {
				file->squareMesh = new SquareMesh((double)file->squareMeshGui->getDx(),file->selection->cell);
				file->squareMeshGui->gridRGB[0] = 1.0;
				file->squareMeshGui->gridRGB[1] = 1.0;
				file->squareMeshGui->gridRGB[2] = 1.0;
				file->squareMeshGui->localizationNumberLabelButton->value(0);
			}
			file->squareMeshGui->applySliders();
			file->squareMeshGui->minPointsSlider->bounds(2,file->squareMesh->getMaxCell()->getCount());
			file->squareMeshGui->minPointsSlider->value(20);
			file->applySquareMesh();
			file->squareMeshGui->landscapeButton->activate();
			iMAP->overlayAdjustment = true;
			file->squareMeshOverlay = true;
			file->meshType = 0;
			file->squareMeshGui->inferButton->activate();
			file->squareMeshGui->resetButton->activate();
			file->squareMeshGui->neighbourDistanceSlider->do_callback();
			break;
		case 1: // voronoi mesh
			file->voronoiMeshGui->inferButton->deactivate();
			file->voronoiMeshGui->resetButton->deactivate();
			file->voronoiMeshGui->inferenceModeChoice->activate();
			if (file->voronoiMesh != NULL) { delete file->voronoiMesh; }
			if (!iMAP->selectionButtonPressed) {
				file->voronoiMesh = new VoronoiMesh();
				file->voronoiMeshGui->gridRGB[0] = 1.0;
				file->voronoiMeshGui->gridRGB[1] = 1.0;
				file->voronoiMeshGui->gridRGB[2] = 1.0;
				file->voronoiMeshGui->localizationNumberLabelButton->value(0);
				file->applyVoronoiMesh();
				file->voronoiMesh->createClusters();
			} else {
				file->voronoiMesh = new VoronoiMesh(file->selection->cell);
				file->voronoiMeshGui->gridRGB[0] = 1.0;
				file->voronoiMeshGui->gridRGB[1] = 1.0;
				file->voronoiMeshGui->gridRGB[2] = 1.0;
				file->voronoiMeshGui->localizationNumberLabelButton->value(0);
				file->applyVoronoiMesh();
				file->voronoiMesh->createClustersSelection();
			}
			file->voronoiMeshGui->applySliders();
			file->voronoiMeshGui->landscapeButton->activate();
			file->voronoiMeshOverlay = true;
			iMAP->overlayAdjustment = true;
			file->voronoiMeshGui->minPointsSlider->bounds(2,file->voronoiMesh->getMaxCell()->getCount());
			file->voronoiMeshGui->minPointsSlider->value(20);
			file->meshType = 1;
			file->voronoiMeshGui->updateVariables(file->voronoiMesh->getTotalVariables());
			file->voronoiMeshGui->inferButton->activate();
			file->voronoiMeshGui->resetButton->activate();
			file->voronoiMeshGui->maxIterationsSlider->deactivate();
			file->voronoiMeshGui->neighbourDistanceSlider->do_callback();
			break;
		case 2: // tree mesh
			file->treeMeshGui->inferButton->deactivate();
			file->treeMeshGui->resetButton->deactivate();
			file->treeMeshGui->inferenceModeChoice->activate();
			if (file->treeMesh != NULL) { delete file->treeMesh; }

			if (!iMAP->selectionButtonPressed) {
				file->treeMesh = new TreeMesh();
			} else {
				file->treeMesh = new TreeMesh(file->selection->cell);
			}

			if (!iMAP->selectionButtonPressed) {
				file->treeMesh->quadTree = new QuadTree(Rect(file->xMin, file->yMin, file->maxRange, file->maxRange),
															file->treeMeshGui->getCount(),
															file->treeMeshGui->minSideSizeSlider->value()/1000.0);
				generateTree(file->treeMesh->quadTree);
			} else {
				file->treeMesh->quadTree = new QuadTree(Rect(file->treeMesh->selection.xMin, file->treeMesh->selection.yMin, FMAX(file->treeMesh->selection.xRange,file->treeMesh->selection.yRange),FMAX(file->treeMesh->selection.xRange,file->treeMesh->selection.yRange)),
															file->treeMeshGui->getCount(),
															file->treeMeshGui->minSideSizeSlider->value()/1000.0);
				generateSelectionTree(file->treeMesh->quadTree);
			}

			file->treeMeshOverlay = true;

			// determine whether leaf should be active or not
			file->treeMesh->minCount = 1000000;
			file->treeMesh->maxCount = -10000000;
			file->treeMesh->totalVariables = 0;
			file->treeMesh->applyTreeMesh(file->treeMesh->quadTree);

			// find neighbours
			file->treeMesh->resetNeighbourCount(file->treeMesh->quadTree);
			file->treeMesh->findNeighbours(file->treeMesh->quadTree);
			file->treeMesh->assignNeighbours(file->treeMesh->quadTree);
			file->treeMeshGui->minLeafPowerSlider->bounds(0,file->treeMesh->maxQuadTreePower);
			file->treeMeshGui->gridRGB[0] = 1.0;
			file->treeMeshGui->gridRGB[1] = 1.0;
			file->treeMeshGui->gridRGB[2] = 1.0;
			file->treeMeshGui->localizationNumberLabelButton->value(0);

			iMAP->overlayAdjustment = true;
			file->treeMeshGui->applySliders();
			file->applyTreeMesh();
			file->treeMeshGui->landscapeButton->activate();
			file->meshType = 2;
			file->treeMeshGui->updateVariables(file->treeMesh->totalVariables);

			file->treeMeshGui->inferButton->activate();
			file->treeMeshGui->resetButton->activate();
			file->treeMeshGui->neighbourDistanceSlider->do_callback();
			break;
		}
		iMAP->customSelectionInferenceGui->disableSelection();
		iMAP->customSelectionInferenceGui->deactivate();
		iMAP->startIntervalSlider->deactivate();
		iMAP->endIntervalSlider->deactivate();

	}
	else { fl_alert("Invalid interval selection."); }
	Fl::redraw();
}

void savePosteriorCallback(Fl_Button*w,int*v) {
	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {

	case 0:
		file->squareMesh->savePosterior();
		break;
	case 1:
		file->voronoiMesh->savePosterior();
		break;
	case 2:
		file->treeMesh->savePosterior();
		break;
	}
	Fl::redraw();
}

CustomSelectionInferenceGui::CustomSelectionInferenceGui(int x, int y, int w, int h) : Fl_Group(x,y,w,h) {

	DFinferred = false;
	DDrinferred = false;

	this->begin();

	surroundingBox = new Fl_Box(x,y,w,h);
	surroundingBox->box(FL_BORDER_FRAME);
	surroundingBox->show();

	selectionButton = new Fl_Button(x+5,y+5,w/2-6,25,"Make Selection");
	selectionButton->type(FL_TOGGLE_BUTTON);
	selectionButton->labelsize(12);
	selectionButton->labelfont(iMAP->boldFont);
	selectionButton->align(FL_ALIGN_CENTER);
	selectionButton->callback((Fl_Callback*)selectionButtonCallback);
	selectionButton->show();

	inferDFButton = new Fl_Button(x+5,y+35,w/4-5,25,"Infer (D,F)");
	inferDFButton->labelsize(10);
	inferDFButton->labelfont(iMAP->boldFont);
	inferDFButton->callback((Fl_Callback*)customSelectionInferDFButtonCallback);
	inferDFButton->deactivate();
	inferDFButton->show();

	inferDDrButton = new Fl_Button(x+5+w/4,y+35,w/4-5,25,"Infer (D,Drift)");
	inferDDrButton->labelsize(10);
	inferDDrButton->labelfont(iMAP->boldFont);
	inferDDrButton->callback((Fl_Callback*)customSelectionInferDDrButtonCallback);
	inferDDrButton->deactivate();
	inferDDrButton->show();

	noiseSigmaSlider = new Slider(x+5,y+80,w/2-10,20,"Localization Precision [nm]");
	noiseSigmaSlider->precision(0);
	noiseSigmaSlider->labelsize(12);
	noiseSigmaSlider->labelfont(iMAP->boldFont);
	noiseSigmaSlider->labelcolor(FL_WHITE);
	noiseSigmaSlider->bounds(0.0,100.0);
	noiseSigmaSlider->value(30);
	noiseSigmaSlider->callback((Fl_Callback*)customSelectionNoiseSigmaSliderCallback);
	noiseSigmaSlider->show();

	pointsLabelBox = new Fl_Box(x+w/2+5,y+5,w/4-5,20,"Points");
	pointsLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	pointsLabelBox->labelcolor(FL_WHITE);
	pointsLabelBox->color(iMAP->bgColor);
	pointsLabelBox->labelfont(iMAP->boldFont);
	pointsLabelBox->labelsize(12);
	pointsLabelBox->show();

	pointsBox = new Fl_Box(x+3*w/4+5,y+5,w/4-10,20,"0");
	pointsBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	pointsBox->labelcolor(FL_WHITE);
	pointsBox->color(iMAP->bgColor);
	pointsBox->labelfont(iMAP->normalFont);
	pointsBox->labelsize(12);
	pointsBox->show();

	diffusionLabelBox = new Fl_Box(x+w/2+5,y+24,w/4-5,20,"Diffusion [um\262/s]");
	diffusionLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	diffusionLabelBox->labelcolor(FL_WHITE);
	diffusionLabelBox->color(iMAP->bgColor);
	diffusionLabelBox->labelfont(iMAP->boldFont);
	diffusionLabelBox->labelsize(12);
	diffusionLabelBox->show();

	diffusionBox = new Fl_Box(x+3*w/4+25,y+24,w/4-30,20,"0");
	diffusionBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	diffusionBox->labelcolor(FL_WHITE);
	diffusionBox->color(iMAP->bgColor);
	diffusionBox->labelfont(iMAP->normalFont);
	diffusionBox->labelsize(12);
	diffusionBox->show();

	forceLabelBox = new Fl_Box(x+w/2+5,y+43,w/4-5,20,"Force [pN]");
	forceLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	forceLabelBox->labelcolor(FL_WHITE);
	forceLabelBox->color(iMAP->bgColor);
	forceLabelBox->labelfont(iMAP->boldFont);
	forceLabelBox->labelsize(12);
	forceLabelBox->show();

	forceBox = new Fl_Box(x+3*w/4+15,y+43,w/4-20,20,"{0,0}");
	forceBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	forceBox->labelcolor(FL_WHITE);
	forceBox->color(iMAP->bgColor);
	forceBox->labelfont(iMAP->normalFont);
	forceBox->labelsize(12);
	forceBox->show();

	forceMagnitudeLabelBox = new Fl_Box(x+w/2+5,y+62,w/4-5,20,"Force Magnitude [pN]");
	forceMagnitudeLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	forceMagnitudeLabelBox->labelcolor(FL_WHITE);
	forceMagnitudeLabelBox->color(iMAP->bgColor);
	forceMagnitudeLabelBox->labelfont(iMAP->boldFont);
	forceMagnitudeLabelBox->labelsize(12);
	forceMagnitudeLabelBox->show();

	forceMagnitudeBox = new Fl_Box(x+3*w/4+15,y+62,w/4-20,20,"0");
	forceMagnitudeBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	forceMagnitudeBox->labelcolor(FL_WHITE);
	forceMagnitudeBox->color(iMAP->bgColor);
	forceMagnitudeBox->labelfont(iMAP->normalFont);
	forceMagnitudeBox->labelsize(12);
	forceMagnitudeBox->show();

	areaLabelBox = new Fl_Box(x+w/2+5,y+81,w/4+20,20,"Area [um\262]");
	areaLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	areaLabelBox->labelcolor(FL_WHITE);
	areaLabelBox->color(iMAP->bgColor);
	areaLabelBox->labelfont(iMAP->boldFont);
	areaLabelBox->labelsize(12);
	areaLabelBox->show();

	areaBox = new Fl_Box(x+3*w/4+30,y+81,w/4-35,20,"0");
	areaBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	areaBox->labelcolor(FL_WHITE);
	areaBox->color(iMAP->bgColor);
	areaBox->labelfont(iMAP->normalFont);
	areaBox->labelsize(12);
	areaBox->show();

	this->end();
}

DistanceButton::DistanceButton(int x, int y, int w, int h, const char* l) : Fl_Group(x,y,w,h,l) {
	begin();
	Fl_Group* db_group = new Fl_Group(x,y,w,h,l);
	labelBox = new Fl_Box(x,y,60,20,"Distance Type:");
	labelBox->labelsize(12);
	labelBox->labelfont(iMAP->boldFont);
	labelBox->labelcolor(FL_WHITE);
	labelBox->show();
	l1Button = new Fl_Round_Button(x+80,y,35,20,"L1");
	l1Button->labelsize(12);
	l1Button->labelfont(iMAP->boldFont);
	l1Button->labelcolor(FL_WHITE);
	l1Button->type(102);
	l2Button = new Fl_Round_Button(x+115,y,35,20,"L2");
	l2Button->labelsize(12);
	l2Button->labelfont(iMAP->boldFont);
	l2Button->labelcolor(FL_WHITE);
	l2Button->type(102);
	l2Button->value(1);
	db_group->end();
	end();
	show();
}

Slider::Slider(int x, int y, int w, int h, const char *l) : Fl_Group(x,y,w,h,l) {
	this->labelsize(12);

	int in_w = 50;
	input = new Fl_Float_Input(x, y, in_w, h);
	input->callback(Input_CB, (void*)this);
	input->when(FL_WHEN_CHANGED);
	input->labelsize(12);
	input->textfont(iMAP->normalFont);
	input->textcolor(FL_WHITE);
	input->box(FL_BORDER_FRAME);
	input->labelcolor(FL_WHITE);
	input->color(FL_WHITE);
	input->textsize(10);

	slider = new Fl_Slider(x+in_w, y, w-in_w, h);
	slider->type(1);
	slider->box(FL_BORDER_FRAME);
	slider->callback(Slider_CB, (void*)this);
	slider->align(FL_ALIGN_TOP);
	slider->labelfont(iMAP->boldFont);
	slider->labelcolor(FL_WHITE);
	slider->color(FL_GREEN);
	slider->labelsize(12);

	bounds(0, 1);     // some usable default
	value(0);          // some usable default
	decimalPoints = 0;
	end();             // close the group
}

FileInfo::FileInfo(int x,int y,int w,int h) : Fl_Group(x,y,w,h) {

	this->begin();

	surroundingBox = new Fl_Box(x,y,w,h);
	surroundingBox->box(FL_BORDER_FRAME);
	surroundingBox->color(FL_WHITE);
	surroundingBox->show();

	filesLabelBox = new Fl_Box(x+5,y+5,40,18);
	filesLabelBox->labelsize(12);
	filesLabelBox->label("Trajectories ");
	filesLabelBox->labelfont(iMAP->boldFont);
	filesLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	filesLabelBox->labelcolor(FL_WHITE);
	filesLabelBox->color(FL_BLACK);
	filesLabelBox->show();

	filesBox = new Fl_Box(x+100,y,100,25);
	filesBox->labelsize(12);
	filesBox->label("0");
	filesBox->labelfont(iMAP->normalFont);
	filesBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	filesBox->labelcolor(FL_WHITE);
	filesBox->color(FL_BLACK);
	filesBox->show();

	durationLabelBox = new Fl_Box(x+5,y+20,100,25);
	durationLabelBox->labelsize(12);
	durationLabelBox->label("Duration [s] ");
	durationLabelBox->labelfont(iMAP->boldFont);
	durationLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	durationLabelBox->labelcolor(FL_WHITE);
	durationLabelBox->color(FL_BLACK);
	durationLabelBox->show();

	durationBox = new Fl_Box(x+100,y+20,100,25);
	durationBox->labelsize(12);
	durationBox->label("0.0");
	durationBox->labelfont(iMAP->normalFont);
	durationBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	durationBox->labelcolor(FL_WHITE);
	durationBox->color(FL_BLACK);
	durationBox->show();

	acquisitionTimeLabelBox = new Fl_Box(x+5,y+40,100,25);
	acquisitionTimeLabelBox->labelsize(12);
	acquisitionTimeLabelBox->label("Acquisition Time [ms] ");
	acquisitionTimeLabelBox->labelfont(iMAP->boldFont);
	acquisitionTimeLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	acquisitionTimeLabelBox->labelcolor(FL_WHITE);
	acquisitionTimeLabelBox->color(FL_BLACK);
	acquisitionTimeLabelBox->show();

	acquisitionTimeBox = new Fl_Box(x+100,y+40,100,25);
	acquisitionTimeBox->labelsize(12);
	acquisitionTimeBox->label("0");
	acquisitionTimeBox->labelfont(iMAP->normalFont);
	acquisitionTimeBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	acquisitionTimeBox->labelcolor(FL_WHITE);
	acquisitionTimeBox->color(FL_BLACK);
	acquisitionTimeBox->show();

	totalDetectionsLabelBox = new Fl_Box(x+210,y,100,25);
	totalDetectionsLabelBox->labelsize(12);
	totalDetectionsLabelBox->label("Points");
	totalDetectionsLabelBox->labelfont(iMAP->boldFont);
	totalDetectionsLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	totalDetectionsLabelBox->labelcolor(FL_WHITE);
	totalDetectionsLabelBox->color(FL_BLACK);
	totalDetectionsLabelBox->show();

	totalDetectionsBox = new Fl_Box(x+305,y,100,25);
	totalDetectionsBox->labelsize(12);
	totalDetectionsBox->label("0");
	totalDetectionsBox->labelfont(iMAP->normalFont);
	totalDetectionsBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	totalDetectionsBox->labelcolor(FL_WHITE);
	totalDetectionsBox->color(FL_BLACK);
	totalDetectionsBox->show();

	dimensionsLabelBox = new Fl_Box(x+210,y+20,100,25);
	dimensionsLabelBox->labelsize(12);
	dimensionsLabelBox->label("Dimensions [um] ");
	dimensionsLabelBox->labelfont(iMAP->boldFont);
	dimensionsLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	dimensionsLabelBox->labelcolor(FL_WHITE);
	dimensionsLabelBox->color(FL_BLACK);
	dimensionsLabelBox->show();

	dimensionsBox = new Fl_Box(x+305,y+20,100,25);
	dimensionsBox->labelsize(12);
	dimensionsBox->label("0.0 x 0.0");
	dimensionsBox->labelfont(iMAP->normalFont);
	dimensionsBox->labelcolor(FL_WHITE);
	dimensionsBox->color(FL_BLACK);
	dimensionsBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	dimensionsBox->show();

	averageStepLabelBox = new Fl_Box(x+210,y+40,100,25);
	averageStepLabelBox->labelsize(12);
	averageStepLabelBox->label("Average Step [nm] ");
	averageStepLabelBox->labelfont(iMAP->boldFont);
	averageStepLabelBox->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
	averageStepLabelBox->labelcolor(FL_WHITE);
	averageStepLabelBox->color(FL_BLACK);
	averageStepLabelBox->show();

	averageStepBox = new Fl_Box(x+305,y+40,100,25);
	averageStepBox->labelsize(12);
	averageStepBox->label("0");
	averageStepBox->labelfont(iMAP->normalFont);
	averageStepBox->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
	averageStepBox->labelcolor(FL_WHITE);
	averageStepBox->color(FL_BLACK);
	averageStepBox->show();

	this->end();
}

void colorOverlayCallback(Fl_Widget*w,int*v) {

	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch(meshType) {
		case 0: // square mesh
			file->squareMeshGui->setColormap(file->squareMeshGui->overlayColorChoice->value());
			iMAP->overlayAdjustment = true;
			break;
		case 1: // voronoi tessellation
			file->voronoiMeshGui->setColormap(file->voronoiMeshGui->overlayColorChoice->value());
			iMAP->overlayAdjustment = true;
			break;
		case 2: // quad-tree mesh
			file->treeMeshGui->setColormap(file->treeMeshGui->overlayColorChoice->value());
			iMAP->overlayAdjustment = true;
			break;
	}
}

void densityButtonCallback(Fl_Button*w,void*) {
	iMAP->startIntervalSlider->value(file->tMin);
	iMAP->endIntervalSlider->value(file->tMax);

	calculateDensity();
	file->plotType = 1;
	file->plotChanged = true;
	file->densityCalculated = true;
	iMAP->densityGui->densityButton->deactivate();
	file->drawLocalizations = true;
	iMAP->localizationsButton->value(1);


	file->animateTrajectories = false;
	iMAP->animateTrajectoriesButton->value(0);

	file->drawTrajectories = false;
	iMAP->drawTrajectoriesButton->value(0);
}

void neighbourRadiusCallback(Slider*w,void*) {
	file->densityCalculated = false;

	iMAP->densityGui->densityButton->activate();

	Fl::redraw();
}

DensityGui::DensityGui(int x, int y, int w, int h, const char* l) : Fl_Group(x,y,w,h,l) {

	this->begin();

	surroundingBox = new Fl_Box(x,y,w,h);
	surroundingBox->box(FL_BORDER_FRAME);
	surroundingBox->color(FL_WHITE);
	surroundingBox->show();

	densityButton = new Fl_Button(x+5,y+5,w/2-40,35,"Calculate\nDensity");
	densityButton->labelsize(12);
	densityButton->labelfont(iMAP->boldFont);
	densityButton->align(FL_ALIGN_CENTER);
	densityButton->callback((Fl_Callback*)densityButtonCallback);
	densityButton->show();

	neighbourRadiusSlider = new Slider(x+w/2-30,y+20,w/2+25,20,"Neighbor Radius [nm]");
	neighbourRadiusSlider->precision(0);
	neighbourRadiusSlider->labelsize(11);
	neighbourRadiusSlider->labelcolor(FL_WHITE);
	neighbourRadiusSlider->labelfont(iMAP->boldFont);
	neighbourRadiusSlider->value(100);
	neighbourRadiusSlider->callback((Fl_Callback*)neighbourRadiusCallback);
	neighbourRadiusSlider->bounds(20,1000);
	neighbourRadiusSlider->deactivate();
	neighbourRadiusSlider->show();

	this->end();

}

StereoVisionGui::StereoVisionGui() {
	const int w = 200;
	const int h = 40;

	window = new Fl_Window(40,50,w,h,"Stereo Vision");
	window->color(iMAP->bgColor);
	window->labelfont(iMAP->boldFont);
	window->begin();
	window->set_non_modal();

	eyeSeparationSlider = new Slider(5,15,w-10,20,"Separation [a.u.]");
	eyeSeparationSlider->bounds(0,2.0);
	eyeSeparationSlider->value(1.5);
	eyeSeparationSlider->precision(2);
	eyeSeparationSlider->labelcolor(FL_WHITE);
	eyeSeparationSlider->labelsize(10);
	eyeSeparationSlider->labelfont(iMAP->boldFont);
//			eyeSeparationSlider->callback((Fl_Callback*)eyeSeparationSliderCallback);
	eyeSeparationSlider->show();

	window->end();
//	window->show();
}
