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

#ifndef GUI_H_
#define GUI_H_

#include <stdio.h>			// Header File For Standard Input / Output
#include <stdarg.h>			// Header File For Variable Argument Routines
#include <stdlib.h>

#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Tabs.H>
#include <FL/fl_color_chooser.h>
#include <FL/Fl_Pixmap.H>

class Slider;
class FileInfo;
class SquareMeshGui;
class VoronoiMeshGui;
class TreeMeshGui;
class DistanceButton;
class ColorBar;
class SingleTrajectoryInferenceGui;
class CustomSelectionInferenceGui;
class StereoVisionGui;

#include "file.h"
extern File *file;

// mesh-based inference
void squareMeshWindowCallback(Fl_Widget*w,void*v);
void applySquareMeshCallback(Fl_Widget*w,void*v);
void minPointsSliderCallback(Slider*w,void*v);
void optimizationChoiceCallback(Fl_Choice*w,int*v);
void saveSquareMeshCallback(Fl_Button*w,int*v);
void cellsSliderCallback(Slider*w,void*);
void nullCallback(Fl_Widget*w,void*v);
void voronoiMeshWindowCallback(Fl_Widget *w,void*v);
void saveVoronoiMeshCallback(Fl_Button*w,int*v);
void clusteringChoiceCallback(Fl_Choice*w,int*v);
void colormapChoiceSquareCallback(Fl_Choice*w,int*v);
void colormapChoiceVoronoiCallback(Fl_Choice*w,int*v);
void colormapChoiceQuadTreeCallback(Fl_Choice*w,int*v);
void diffusionOverlayCallback(Fl_Check_Button*w,int*v);
void diffusionLogOverlayCallback(Fl_Check_Button*w,int*v);
void potentialOverlayCallback(Fl_Check_Button*w,int*v);
void forceMagnitudeOverlayCallback(Fl_Check_Button*w,int*v);
void pointNumberOverlayCallback(Fl_Check_Button*w,int*v);
void flipColormapCallback(Fl_Check_Button*w,int*v);
void treeMeshWindowCallback(Fl_Widget *w,void*v);
void updateTreeThresholdsCallback(Slider*w,void*);
void updateTreeThresholds(QuadTree *tree);
void generateTrajectoriesCallback(Fl_Widget*,int*v);
void saveTrajectoriesCallback(Fl_Widget*,int*v);
void selectionButtonCallback(Fl_Button*w, void*);
void gridColorButtonCallback(Fl_Button*w,int*v);
void simulationWidgetsCallback(Fl_Widget*w,int*v);
void optimizationFunctionChoiceCallback(Fl_Choice*w,int*v);
void smoothingPriorButtonCallback(Fl_Check_Button*w,int*v);
void jeffreysPriorButtonCallback(Fl_Check_Button*w,int*v);
void crossProduct(float *v, float v1x, float v1y, float v1z, float v2x, float v2y, float v2z);
void crossProduct(float*n, float*v1,float*v2);
void roButtonCallback(Fl_Check_Button*w,void*v);
void maximumNeighbourDistanceCallback(Slider*w,int*v);
void randomizedOptimizationToleranceCallback(Slider*w,int*v);
void axesCallback(Fl_Check_Button*w,void*);
void customSelectionNoiseSigmaSliderCallback(Fl_Widget*w,void*);
void samplePosteriorCallback(Fl_Button*w,int*v);
void savePosteriorCallback(Fl_Button*w,int*v);
void colorOverlayCallback(Fl_Widget*w,int*v);

void inferMeshCallback(Fl_Button*w,int*v);
void resetMeshCallback(Fl_Button*w,int*v);
void applyMeshCallback(Fl_Button*w,int*v);

void stopCalculationButtonCallback(Fl_Button*w,int*v);
void pauseCalculationButtonCallback(Fl_Button*w,int*v);

void voronoiLandscapeCallback(Fl_Check_Button*w,void*);
void generateVoronoiLandscape();
void updateVoronoiLandscapeCallback(Fl_Widget*w,void*v);
void alphaVoronoiLandscapeCallback(Fl_Widget*w,void*v);
void setLandscapeVoronoi();

void squareLandscapeCallback(Fl_Check_Button*w,void*);
void setLandscapeSquare();
void generateSquareLandscape();
void updateSquareLandscapeCallback(Fl_Widget*w,void*v);
void alphaSquareLandscapeCallback(Fl_Widget*w,void*v);

void quadTreeLandscapeCallback(Fl_Check_Button*w,void*);
void updateQuadTreeLandscapeCallback(Fl_Widget*w,void*);
void alphaQuadTreeLandscapeCallback(Fl_Widget*w,void*v);
void generateQuadTreeLandscape();
void countQuadTreeVertices(QuadTree *tree, int *count);
void setLandscapeQuadTree(QuadTree *tree);
void setCornerLandscapesQuadTree(double x, double y, double *value, int *num, QuadTree *tree);
void findCommonVertex(float *vertices, QuadTree *tree1, QuadTree *tree2);
void averageQuadTreeLandscapeNormals(QuadTree *tree, float x, float y, float *norms, int *num);
void defineQuadTreeLandscape(QuadTree *tree, int *i);
void countQuadTreeBorderVertices(QuadTree *tree, int *count);
void defineQuadTreeLandscapeBorders(QuadTree*tree,int*i);
void saveTreeMeshCallback(Fl_Button*w,int*v);
void densityButtonCallback(Fl_Button*w,void*);
void neighbourRadiusCallback(Slider*w,void*);

// custom selection inference
void customSelectionInferDFButtonCallback(Fl_Button*w,void*);
void customSelectionInferenceWindowCallback(Fl_Window*w,void*);
void customSelectionSaveCallback(Fl_Button*w,void*);
void customSelectionSampleDiffusionCallback(Fl_Button*w,void*);
void customSelectionSampleDiffusionSaveCallback(Fl_Button*w,void*);
void customSelectionSampleFxCallback(Fl_Button*w,void*);
void customSelectionSampleFxSaveCallback(Fl_Button*w,void*);
void customSelectionSampleFyCallback(Fl_Button*w,void*);
void customSelectionSampleFySaveCallback(Fl_Button*w,void*);

// posteriori tab related
void dPosteriorCallback(Fl_Button*w,void*v);
void fPosterioriCallback(Fl_Button*w,void*v);
void vPosterioriCallback(Fl_Button*w,void*v);
void potentialReferenceCallback(Fl_Button*w,void*v);

// trajectory inference
void minDiffusionSlider(Slider*w,void*);
void maxDiffusionSlider(Slider*w,void*);
void singleTrajectoryCallback(Fl_Button*w,void*);
void singleTrajectoryInferenceWindowCallback(Fl_Window*w,void*);
void singleTrajectoryStartCallback(Slider*w,void*);
void singleTrajectoryEndCallback(Slider*w,void*);
void singleTrajectoryInferenceCallback(Fl_Button*,void*);
void singleTrajectorySaveCallback(Fl_Button*,void*);
void singleTrajectorySampleDiffusionCallback(Fl_Button*,void*);
void singleTrajectorySampleDiffusionSaveCallback(Fl_Button*,void*);
void singleTrajectorySampleForceCallback(Fl_Button*,void*);
void singleTrajectorySampleForceSaveCallback(Fl_Button*,void*);

extern Fl_Menu_Item fileStatsItems[];
extern File *file;

class Slider : public Fl_Group {

    Fl_Float_Input *input;
    Fl_Slider *slider;
	int decimalPoints;
	float fieldWidth;

    // CALLBACK HANDLERS
    //    These 'attach' the input and slider's values together.
    //
    void Slider_CB2() {
        static int recurse = 0;
        if ( recurse ) {
            return;
        } else {
            recurse = 1;

            char s[80];
			switch (this->decimalPoints) {
				case 0:
					sprintf(s, "%.0f", slider->value());
					break;
				case 1:
					sprintf(s, "%.1f", slider->value());
					break;
				case 2:
					sprintf(s, "%.2f", slider->value());
					break;
				case 3:
					sprintf(s, "%.3f", slider->value());
					break;
				case 4:
					sprintf(s, "%.4f", slider->value());
					break;
			}

			// exceptional cases
			if (this->decimalPoints > 4) { sprintf(s, "%.4f", slider->value()); }
			if (this->decimalPoints < 0) { sprintf(s, "%.0f", slider->value()); }

            // fprintf(stderr, "SPRINTF(%d) -> '%s'\n", (int)(slider->value()+.5), s);
            input->value(s);          // pass slider's value to input
            recurse = 0;
        }
		// necessary to explicitly define callback condition of Fl_Group
		this->do_callback();
    }

    static void Slider_CB(Fl_Widget *w, void *data) {
        ((Slider*)data)->Slider_CB2();
        Fl::redraw();
    }

    void Input_CB2() {
        static int recurse = 0;
        if ( recurse ) {
            return;
        } else {
            recurse = 1;
            float val = 0;
            if ( sscanf(input->value(), "%f", &val) != 1 ) {
                val = 0;
            }
            // fprintf(stderr, "SCANF('%s') -> %d\n", input->value(), val);
            slider->value(val);         // pass input's value to slider
            recurse = 0;
        }
		// necessary to explicitly define callback condition of Fl_Group
		this->do_callback();
    }

    static void Input_CB(Fl_Widget *w, void *data) {
        ((Slider*)data)->Input_CB2();
        Fl::redraw();
    }

	public:
		// CTOR
		Slider(int x, int y, int w, int h, const char *l=0);

		// MINIMAL ACCESSORS --  Add your own as needed
		float value() const { return((float)slider->value()); }
		void value(float val) { slider->value(val); Slider_CB2(); }
		void minumum(float val) { slider->minimum(val); }
		float minumum() const { return((float)slider->minimum()); }
		void maximum(float val) { slider->maximum(val); }
		float maximum() const { return((float)slider->maximum()); }
		void bounds(float low, float high) { slider->bounds(low, high); }
		void type(int t) { slider->type(t); }
		void precision(int decimalPoints) {
			this->decimalPoints = decimalPoints;
			if (decimalPoints == 0) { slider->step(1); }
		}
		void deactivate() {
//			slider->labelcolor(FL_RED);
//			input->labelcolor(FL_RED);
//			input->color2(FL_RED);
//			slider->color2(FL_RED);
//			slider->selection_color(FL_RED);
//			input->selection_color(FL_RED);
			slider->color(FL_RED);
//			input->color(FL_RED);
//			slider->redraw_label();
//			input->redraw_label();
			slider->deactivate();
			input->deactivate();
		}
		void activate() {
//			slider->labelcolor(FL_GREEN);
//			input->labelcolor(FL_GREEN);
//			input->color2(FL_GREEN);
//			slider->color2(FL_GREEN);
//			slider->selection_color(FL_GREEN);
//			input->selection_color(FL_GREEN);
			slider->color(FL_GREEN);
//			input->color(FL_GREEN);
			slider->activate();
			input->activate();
		}

};

class FileInfo : public Fl_Group {

	// labels
	Fl_Box *surroundingBox;
	Fl_Box*filesLabelBox;
	Fl_Box *dimensionsLabelBox;
	Fl_Box *totalDetectionsLabelBox;
	Fl_Box *averageStepLabelBox;
	Fl_Box *durationLabelBox;
	Fl_Box *acquisitionTimeLabelBox;

	Fl_Box *filesBox;
	Fl_Box *dimensionsBox;
	Fl_Box *totalDetectionsBox;
	Fl_Box *averageStepBox;
	Fl_Box *durationBox;
	Fl_Box *acquisitionTimeBox;

	public:

	void update();
	void updatePoints();

	FileInfo(int x,int y,int w,int h);
};

class DensityGui : public Fl_Group {

	public:
		Fl_Box *surroundingBox;
		Slider *neighbourRadiusSlider;
		Fl_Button *densityButton;

		DensityGui(int x, int y, int w, int h, const char* l = 0);

		void fileOpen() {
			this->activate();
			neighbourRadiusSlider->activate();
			densityButton->activate();
		}

		double getRadius() { return neighbourRadiusSlider->value()/1000.0; }
};

class SquareMeshGui {

	private:
		bool meshApplied;
		Fl_Window *window;
		Slider *dxdySlider;
		int colormap;

		Fl_Tabs *tabs;
		Fl_Group *inferenceGroup;
		Fl_Group *annotationsGroup;
		Fl_Group *advancedGroup;
		Fl_Group *landscapeGroup;

	public:
		SquareMeshGui();

		Slider *minPointsSlider;
		Slider *polynomialOrderSlider;
		Slider *betaSlider;

		Fl_Button *saveButton;
		Fl_Button *resetButton;
		Fl_Button *applyButton;
		Fl_Button *inferButton;
		Fl_Button *stopButton;
		Fl_Button *pauseButton;
		Fl_Choice *inferenceModeChoice;
		Fl_Box *variablesBox;

		Fl_Check_Button *overlayButton;
		Fl_Check_Button *overlayPointNumberButton;
		Fl_Check_Button *overlayDiffusionButton;
		Fl_Check_Button *overlayDiffusionLogButton;
		Fl_Check_Button *overlayPotentialButton;
		Fl_Check_Button *overlayForceArrowsButton;
		Fl_Check_Button *overlayForceMagnitudeButton;
		Fl_Choice *overlayColorChoice;
		Fl_Check_Button *flipColormapButton;

		Fl_Check_Button *umUnitsButton;
		Fl_Check_Button *nmUnitsButton;
		Fl_Check_Button *localizationNumberLabelButton;
		Fl_Button *infoOverlayButton;
		Slider *overlayAlphaSlider;
		Slider *fontScaleSlider;
		Slider *gridAlphaSlider;
		Fl_Button *gridColorButton;
		Fl_Check_Button *displayGridButton;
		Fl_Check_Button *spotVisualizationButton;
		Slider *spotVisualizationScaleSlider;

		double gridRGB[3];

		// posteriori group
		Fl_Round_Button *fPosteriorButton;
		Fl_Round_Button *vPosteriorButton;
		Fl_Round_Button *dPosterioriButton;
		Slider *minPosteriorSlider;
		Slider *maxPosteriorSlider;
		Slider *posterioriSampleNumberSlider;
		Fl_Button *savePosteriorButton;
		Fl_Button *samplePosteriorButton;
		Fl_Button *potentialReferenceButton;

		// advanced group
		Slider *localizationPrecisionSlider;
		Fl_Check_Button *roButton;
		Slider *roRadiusSlider;
		Slider *roToleranceSlider;
		Slider *roMaximumIterationsSlider;
		Fl_Check_Button *neighbourDistanceViewButton;
		Slider *neighbourDistanceSlider;

		// Work calculator
		Fl_Button *workButton;
		Fl_Button *saveWorkButton;
		Fl_Box *workBox;

		// prior group
		Fl_Group *priorGroup;
		Slider *muSlider;
		Slider *lambdaSlider;
		Fl_Box *priorImageBox;
		Fl_Check_Button *smoothingPriorButton;
		Fl_Check_Button *jeffreysPriorButton;

		// simulation group
		Fl_Group *simulationGroup;
		Slider *numberOfTrajectoriesSlider;
		Slider *deltaTSlider;
		Slider *timeStepsSlider;
		Fl_Button *generateTrajectoriesButton;
		Fl_Button *saveTrajectoriesButton;

		// landscape view
		Fl_Check_Button *landscapeButton;
		Slider *scaleLandscapeSlider;
		Slider *landscapeAlphaSlider;
		Slider *xLightPositionSlider;
		Slider *yLightPositionSlider;
		Slider *zLightPositionSlider;
		Fl_Check_Button *fogButton;
		Slider *fogStartSlider;
		Slider *fogEndSlider;
		Fl_Check_Button *landscapeAxisButton;

		int visible() { return window->visible(); }
		void show() {
			window->show();
			polynomialOrderSlider->value(2);
			betaSlider->value(2.0f);
			inferenceModeChoice->value(0);
		}
		void hide() { window->hide(); }
		void close() { window->do_callback(); }

		void applySliders() {
			switch(inferenceModeChoice->value()) {
			case 0: // d
				roButton->deactivate();
				roRadiusSlider->activate();
				roToleranceSlider->activate();
				roMaximumIterationsSlider->activate();
				betaSlider->deactivate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 1: // df
				roButton->deactivate();
				roRadiusSlider->deactivate();
				roToleranceSlider->deactivate();
				roMaximumIterationsSlider->deactivate();
				betaSlider->activate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 2: // ddr
				roButton->deactivate();
				roRadiusSlider->activate();
				roToleranceSlider->activate();
				roMaximumIterationsSlider->activate();
				betaSlider->deactivate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 3: // dv
				roButton->activate();
				roRadiusSlider->activate();
				roToleranceSlider->activate();
				roMaximumIterationsSlider->activate();
				betaSlider->deactivate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->activate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 4: // poly
				roButton->deactivate();
				roRadiusSlider->deactivate();
				roToleranceSlider->deactivate();
				roMaximumIterationsSlider->deactivate();
				betaSlider->deactivate();
				polynomialOrderSlider->activate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->deactivate();
				smoothingPriorButton->deactivate();
				break;
			}
			annotationsGroup->activate();
			advancedGroup->activate();
			landscapeGroup->activate();

			inferenceModeChoice->activate();
			resetButton->activate();
			applyButton->deactivate();
			dxdySlider->deactivate();
			minPointsSlider->activate();
			inferButton->activate();
		}
		void resetSliders() {
			dxdySlider->activate();
			minPointsSlider->deactivate();
			resetButton->deactivate();
//			priorGroup->deactivate();
			annotationsGroup->deactivate();
			advancedGroup->deactivate();
			landscapeGroup->deactivate();
			polynomialOrderSlider->deactivate();
			inferButton->deactivate();
			simulationGroup->deactivate();
		}
		int overlay() { return overlayButton->value(); }

		int getPolynomialOrder() { return (int)polynomialOrderSlider->value(); }
		double getBeta() { return betaSlider->value(); }

		float getDx() {
			if (dxdySlider->value() > 1.0f) { return dxdySlider->value()/1000.0f;	}
			else { return 1.0f/1000.0f; }

		}

		int getPosteriorType() {
			if (this->dPosterioriButton->value()) { return 0; }
			if (this->fPosteriorButton->value()) { return 1; }
			if (this->vPosteriorButton->value()) { return 2; }
			fprintf(stderr,"invalid posterior selection\n");
			return -1;
		}

		void updateVariables(int vars) {
			char label [100];
			sprintf(label,"%i Zones",vars);
			variablesBox->copy_label(label);
		}

		int getColormap() { return colormap; }
		void setColormap(int c) { colormap = c; }

		// posteriori sampling
		Fl_Group *posterioriGroup;
		double *posterioriValues;
		int posterioriSamples;
		double posterioriRange;
		double posterioriIncrement;

		void drawPosteriorPlot(int type,int sup);
		void samplePosterior(int type);
		void savePosteriori(int type);

		void saveSquareMesh();
		void loadSquareMesh();
};

class VoronoiMeshGui {

	private:
		bool meshApplied;
		Fl_Window *window;
		int colormap;
		Fl_Group *inferenceGroup;
		Fl_Group *annotationsGroup;
		Fl_Group *advancedGroup;
		Fl_Tabs *tabs;
		Fl_Group *landscapeGroup;

	public:
		VoronoiMeshGui();

		// inference
		Fl_Button *saveButton;
		Fl_Button *resetButton;
		Fl_Button *applyButton;
		Fl_Button *inferButton;
		Fl_Button *stopButton;
		Fl_Button *pauseButton;
		Fl_Choice *inferenceModeChoice;
		Slider *cellsSlider;
		Slider *polynomialOrderSlider;
		Fl_Check_Button *overlayPointNumberButton;
		Fl_Check_Button *overlayDiffusionButton;
		Fl_Check_Button *overlayDiffusionLogButton;
		Fl_Check_Button *overlayPotentialButton;
		Fl_Check_Button *overlayForceArrowsButton;
		Fl_Check_Button *overlayForceMagnitudeButton;
		Fl_Box *variablesBox;

		// annotations
		Fl_Check_Button *overlayButton;
		Fl_Choice *overlayColorChoice;
		Fl_Check_Button *flipColormapButton;
		Fl_Check_Button *umUnitsButton;
		Fl_Check_Button *nmUnitsButton;
		Fl_Check_Button *localizationNumberLabelButton;
		Fl_Button *infoOverlayButton;
		Slider *fontScaleSlider;
		Fl_Check_Button *clickInformationButton;
		Slider *overlayAlphaSlider;
		Slider *gridAlphaSlider;
		Fl_Button *gridColorButton;
		Fl_Check_Button *displayGridButton;
		Fl_Check_Button *spotVisualizationButton;
		Slider *spotVisualizationScaleSlider;
		double gridRGB[3];

		// advanced
		Slider *betaSlider;
		Slider *maxIterationsSlider;
		Slider *minPointsSlider;
		Fl_Choice *clusteringChoice;
		DistanceButton *distanceButton;
		Slider *localizationPrecisionSlider;
		Fl_Check_Button *roButton;
		Slider *roRadiusSlider;
		Slider *roToleranceSlider;
		Slider *roMaximumIterationsSlider;
		Fl_Check_Button *neighbourDistanceViewButton;
		Slider *neighbourDistanceSlider;

		// Work calculator
		Fl_Button *workButton;
		Fl_Button *saveWorkButton;
		Fl_Box *workBox;

		// landscape view
		Fl_Check_Button *landscapeButton;
		Slider *scaleLandscapeSlider;
		Slider *landscapeAlphaSlider;
		Slider *xLightPositionSlider;
		Slider *yLightPositionSlider;
		Slider *zLightPositionSlider;
		Fl_Check_Button *fogButton;
		Slider *fogStartSlider;
		Slider *fogEndSlider;
		Fl_Check_Button *landscapeAxisButton;

		// posteriori
		Fl_Button *dPosterioriButton;
		Fl_Button *fPosteriorButton;
		Fl_Button *vPosteriorButton;
		Slider *minPosteriorSlider;
		Slider *maxPosteriorSlider;
		Slider *posterioriSampleNumberSlider;
		Fl_Button *samplePosteriorButton;
		Fl_Button *savePosteriorButton;
		Fl_Button *potentialReferenceButton;

		// prior group
		Fl_Group *priorGroup;
		Slider *muSlider;
		Slider *lambdaSlider;
		Fl_Box *priorImageBox;
		Fl_Check_Button *smoothingPriorButton;
		Fl_Check_Button *jeffreysPriorButton;

		// simulation group
		Fl_Group *simulationGroup;
		Slider *numberOfTrajectoriesSlider;
		Slider *deltaTSlider;
		Slider *timeStepsSlider;
		Fl_Button *generateTrajectoriesButton;
		Fl_Button *saveTrajectoriesButton;

		void show();
		int visible() {	return window->visible(); }
		void hide() { window->hide(); }
		void close() { window->do_callback(); }
		void applySliders() {
			switch(inferenceModeChoice->value()) {
			case 0: // d
				roButton->deactivate();
				roRadiusSlider->activate();
				roToleranceSlider->activate();
				roMaximumIterationsSlider->activate();
				betaSlider->deactivate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 1: // df
				roButton->deactivate();
				roRadiusSlider->deactivate();
				roToleranceSlider->deactivate();
				roMaximumIterationsSlider->deactivate();
				betaSlider->activate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 2: // ddr
				roButton->deactivate();
				roRadiusSlider->activate();
				roToleranceSlider->activate();
				roMaximumIterationsSlider->activate();
				betaSlider->deactivate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 3: // dv
				roButton->activate();
				roRadiusSlider->activate();
				roToleranceSlider->activate();
				roMaximumIterationsSlider->activate();
				betaSlider->deactivate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->activate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 4: // poly
				roButton->deactivate();
				roRadiusSlider->deactivate();
				roToleranceSlider->deactivate();
				roMaximumIterationsSlider->deactivate();
				betaSlider->deactivate();
				polynomialOrderSlider->activate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->deactivate();
				smoothingPriorButton->deactivate();
				break;
			}
			annotationsGroup->activate();
			advancedGroup->activate();
			landscapeGroup->activate();

			inferenceModeChoice->activate();
			resetButton->activate();
			applyButton->deactivate();
			cellsSlider->deactivate();
			minPointsSlider->activate();
			inferButton->activate();
			// optimizationChoice->deactivate();
		}
		void resetSliders() {
			cellsSlider->deactivate();
			minPointsSlider->deactivate();
			resetButton->deactivate();
			priorGroup->deactivate();
			annotationsGroup->deactivate();
			advancedGroup->deactivate();
			landscapeGroup->deactivate();
			polynomialOrderSlider->deactivate();
			inferButton->deactivate();
			simulationGroup->deactivate();
		}
		void updateVariables(int vars) {
			char label [100];
			sprintf(label,"%i Zones",vars);
			variablesBox->copy_label(label);
		}
		int overlay() { return overlayButton->value(); }

		int getCells() { return (int)cellsSlider->value(); }
		int getPolynomialOrder() { return (int)polynomialOrderSlider->value(); }
		int getMaxIterations() { return (int)maxIterationsSlider->value(); }
		double getBeta() { return betaSlider->value(); }

		int getColormap() { return colormap; }
		void setColormap(int c) { colormap = c; }

		// posteriori sampling
		Fl_Group *posterioriGroup;
		double *posterioriValues;
		int posterioriSamples;
		double posterioriRange;
		double posterioriIncrement;

		int getPosteriorType() {
			if (this->dPosterioriButton->value()) { return 0; }
			if (this->fPosteriorButton->value()) { return 1; }
			if (this->vPosteriorButton->value()) { return 2; }
			fprintf(stderr,"invalid posterior selection\n");
			return -1;
		}

		void drawPosteriorPlot(int type,int sup);
		void samplePosterior(int type);
		void savePosteriori(int type);

		void saveVoronoiMesh();
		void loadVoronoiMesh();
};

class TreeMeshGui {
	private:
		Fl_Tabs *tabs;
		Fl_Group *inferenceGroup;
		Fl_Group *annotationsGroup;
		Fl_Group *advancedGroup;
		bool meshApplied;
		int colormap;
		Fl_Group *landscapeGroup;

	public:
		TreeMeshGui();

		Fl_Window *window;
		Slider *minPointsSlider;
		Slider *minCapacitySlider;
		Slider *minLeafPowerSlider;
		Slider *polynomialOrderSlider;
		Slider *betaSlider;

		Fl_Button *saveButton;
		Fl_Button *resetButton;
		Fl_Button *applyButton;
		Fl_Button *inferButton;
		Fl_Button *stopButton;
		Fl_Button *pauseButton;
		Fl_Choice *inferenceModeChoice;
		Fl_Box *variablesBox;
		Slider *minSideSizeSlider;

		Fl_Check_Button *overlayButton;
		Fl_Check_Button *overlayDiffusionButton;
		Fl_Check_Button *overlayDiffusionLogButton;
		Fl_Check_Button *overlayPotentialButton;
		Fl_Check_Button *overlayForceArrowsButton;
		Fl_Check_Button *overlayForceMagnitudeButton;
		Fl_Check_Button *overlayPointNumberButton;
		Fl_Choice *overlayColorChoice;
		Fl_Check_Button *flipColormapButton;

		Fl_Check_Button *umUnitsButton;
		Fl_Check_Button *nmUnitsButton;
		Fl_Check_Button *localizationNumberLabelButton;
		Fl_Button *infoOverlayButton;
		Slider *overlayAlphaSlider;
		Slider *fontScaleSlider;
		Slider *gridAlphaSlider;
		Fl_Button *gridColorButton;
		Fl_Check_Button *displayGridButton;
		Fl_Check_Button *spotVisualizationButton;
		Slider *spotVisualizationScaleSlider;
		double gridRGB[3];

		Fl_Check_Button *roButton;
		Slider *roRadiusSlider;
		Slider *roToleranceSlider;
		Slider *roMaximumIterationsSlider;
		Fl_Check_Button *neighbourDistanceViewButton;
		Slider *neighbourDistanceSlider;

		// Work calculator
		Fl_Button *workButton;
		Fl_Button *saveWorkButton;
		Fl_Box *workBox;

		// landscape view
		Fl_Check_Button *landscapeButton;
		Slider *scaleLandscapeSlider;
		Slider *landscapeAlphaSlider;
		Slider *xLightPositionSlider;
		Slider *yLightPositionSlider;
		Slider *zLightPositionSlider;
		Fl_Check_Button *fogButton;
		Slider *fogStartSlider;
		Slider *fogEndSlider;
		Fl_Check_Button *landscapeAxisButton;

		Slider *localizationPrecisionSlider;

		// posteriori
		Fl_Group *posterioriGroup;
		Fl_Button *fPosteriorButton;
		Fl_Button *vPosteriorButton;
		Fl_Button *dPosteriorButton;
		Slider *minPosteriorSlider;
		Slider *maxPosteriorSlider;
		Slider *posterioriSampleNumberSlider;
		Fl_Button *savePosteriorButton;
		Fl_Button *samplePosteriorButton;
		Fl_Button *potentialReferenceButton;

		// prior group
		Fl_Group *priorGroup;
		Slider *muSlider;
		Slider *lambdaSlider;
		Fl_Box *priorImageBox;
		Fl_Check_Button *smoothingPriorButton;
		Fl_Check_Button *jeffreysPriorButton;

		// simulation group
		Fl_Group *simulationGroup;
		Slider *numberOfTrajectoriesSlider;
		Slider *deltaTSlider;
		Slider *timeStepsSlider;
		Fl_Button *generateTrajectoriesButton;
		Fl_Button *saveTrajectoriesButton;

		int visible() {
			return window->visible();
		}
		void show() {
			window->show();
			polynomialOrderSlider->value(2);
			betaSlider->value(2.0f);
			inferenceModeChoice->value(0);
		}
		void hide() { window->hide(); }
		void close() { window->do_callback(); }
		void applySliders() {
			switch(inferenceModeChoice->value()) {
			case 0: // d
				roButton->deactivate();
				roRadiusSlider->activate();
				roToleranceSlider->activate();
				roMaximumIterationsSlider->activate();
				betaSlider->deactivate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 1: // df
				roButton->deactivate();
				roRadiusSlider->deactivate();
				roToleranceSlider->deactivate();
				roMaximumIterationsSlider->deactivate();
				betaSlider->activate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 2: // ddr
				roButton->deactivate();
				roRadiusSlider->activate();
				roToleranceSlider->activate();
				roMaximumIterationsSlider->activate();
				betaSlider->deactivate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 3: // dv
				roButton->activate();
				roRadiusSlider->activate();
				roToleranceSlider->activate();
				roMaximumIterationsSlider->activate();
				betaSlider->deactivate();
				polynomialOrderSlider->deactivate();
				priorGroup->activate();
//				lambdaSlider->activate();
//				muSlider->activate();
				smoothingPriorButton->activate();
				break;
			case 4: // poly
				roButton->deactivate();
				roRadiusSlider->deactivate();
				roToleranceSlider->deactivate();
				roMaximumIterationsSlider->deactivate();
				betaSlider->deactivate();
				polynomialOrderSlider->activate();
				priorGroup->activate();
//				lambdaSlider->deactivate();
//				muSlider->deactivate();
				smoothingPriorButton->deactivate();
				break;
			}
			annotationsGroup->activate();
			advancedGroup->activate();
			landscapeGroup->activate();

			inferenceModeChoice->activate();
			minSideSizeSlider->deactivate();
			resetButton->activate();
			applyButton->deactivate();
			this->minCapacitySlider->deactivate();
			this->minLeafPowerSlider->activate();
			minPointsSlider->activate();
			inferButton->activate();
			// optimizationChoice->deactivate();
		}
		void resetSliders() {
			this->minCapacitySlider->activate();
			this->minLeafPowerSlider->deactivate();
			minPointsSlider->deactivate();
			resetButton->deactivate();
//			priorGroup->deactivate();
			annotationsGroup->deactivate();
			advancedGroup->deactivate();
			landscapeGroup->deactivate();
			polynomialOrderSlider->deactivate();
			inferButton->deactivate();
			simulationGroup->deactivate();
		}
		int getCount() {
			return (int)minCapacitySlider->value();
		}
		int overlay() { return overlayButton->value(); }

		int getPolynomialOrder() { return (int)polynomialOrderSlider->value(); }
		double getBeta() { return betaSlider->value(); }

		int getPosteriorType() {
			if (this->dPosteriorButton->value()) { return 0; }
			if (this->fPosteriorButton->value()) { return 1; }
			if (this->vPosteriorButton->value()) { return 2; }
			fprintf(stderr,"invalid posterior selection\n");
			return -1;
		}

		int getMinPoints() { return (int)minPointsSlider->value(); }
		int getMinLeafPower() { return (int)minLeafPowerSlider->value(); }

		void updateVariables(int vars) {
			char label [100];
			sprintf(label,"%i Zones",vars);
			variablesBox->copy_label(label);
		}

		int getColormap() { return colormap; }
		void setColormap(int c) { colormap = c; }

		// posteriori sampling
		double *posterioriValues;
		int posterioriSamples;
		double posterioriRange;
		double posterioriIncrement;

		void drawPosteriorPlot(int type,int sup);
		void samplePosterior(int type);
		void savePosteriori(int type);

		void saveQuadTreeMesh();

};

class DistanceButton : public Fl_Group {
	Fl_Box *labelBox;
	Fl_Round_Button *l1Button;
	Fl_Round_Button *l2Button;

	public:

	int radio;

	DistanceButton(int x, int y, int w, int h, const char* l = 0);

	int value() {
		if (l1Button->value()) { return 1; }
		else { return 2; }
	}
};

class ColorBar : public Fl_Box {
	void draw();
	public:
		// CONSTRUCTOR
		ColorBar(int x,int y,int w,int h,const char*l=0) : Fl_Box(x,y,w,h,l) {
//			box(FL_BORDER_BOX);
		}
		// DESTRUCTOR
//		~MyColorBar();
};

class SingleTrajectoryInferenceGui{
	public:
		SingleTrajectoryInferenceGui(char *filename);

		Fl_Window *window;
		Fl_Tabs *tabs;

		// single trajectory inference
		Fl_Group *inferenceGroup;
		Fl_Box *filenameBox;
		Fl_Box *pointsBox;
		Fl_Box *diffusionBox;
		Fl_Box *xForceBox;
		Fl_Box *yForceBox;
		Fl_Box *forceMagnitudeBox;
		Slider *trajectoryStartSlider;
		Slider *trajectoryEndSlider;
		Slider *noiseSigmaSlider;
		Fl_Button *inferenceButton;
		Fl_Button *saveButton;

		// diffusion group
		Slider *minDiffusionSlider;
		Slider *maxDiffusionSlider;
		Slider *diffusionSampleNumberSlider;
		Fl_Button *sampleDiffusionButton;
		Fl_Group *diffusionGroup;
		Fl_Button *saveDPosterioriButton;
		double minSampledDiffusion;
		double maxSampledDiffusion;
		double diffusionRange;
		double *diffusionValues;
		int diffusionSamples;
		double diffusionIncrement;

		// Fx group
		Slider *minForceSlider;
		Slider *maxForceSlider;
		Slider *forceSampleNumberSlider;
		Fl_Button *sampleForceButton;
		Fl_Group *forceGroup;
		Fl_Button *saveForcePosterioriButton;
		double minSampledForce;
		double maxSampledForce;
		double forceRange;
		double *forceValues;
		int forceSamples;
		double forceIncrement;

		// data
		int startIndex;
		int endIndex;
		bool valuesSampled;

		// functions
		void hide() { window->hide(); }
		void close() { window->do_callback(); }
		void show() { window->show(); }
		int visible() { return window->visible(); }
		void setStart(double start);
		void setEnd(double end);
		void draw();
		void infer();
		void save();
		void drawDiffusionPosterioriPlot();
		void sampleDiffusionPosteriori();
		void saveDiffusionPosteriori();
		void drawForcePosterioriPlot();
		void sampleForcePosteriori();
		void saveForcePosteriori();

	private:
		double optimizationArray[3];
};

class CustomSelectionInferenceGui : public Fl_Group {

	Fl_Box *surroundingBox;

	public:
		// custom selection inference
		Fl_Box *pointsBox;
		Fl_Box *diffusionBox;
		Fl_Box *forceBox;
		Fl_Box *areaBox;
		Fl_Box *forceMagnitudeBox;
		Fl_Box *pointsLabelBox;
		Fl_Box *diffusionLabelBox;
		Fl_Box *forceLabelBox;
		Fl_Box *areaLabelBox;
		Fl_Box *forceMagnitudeLabelBox;
		Slider *noiseSigmaSlider;
		Fl_Button *inferDFButton;
		Fl_Button *inferDDrButton;
//		Fl_Button *saveButton;
		Fl_Button *selectionButton;

		// diffusion group
		Slider *minDiffusionSlider;
		Slider *maxDiffusionSlider;
		Slider *diffusionSampleNumberSlider;
		Fl_Button *sampleDiffusionButton;
		Fl_Button *saveDPosterioriButton;
		double minSampledDiffusion;
		double maxSampledDiffusion;
		double diffusionRange;
		double *diffusionValues;
		int diffusionSamples;
		double diffusionIncrement;

		// Fx group
		Slider *minFxSlider;
		Slider *maxFxSlider;
		Slider *fxSampleNumberSlider;
		Fl_Button *sampleFxButton;
		Fl_Button *saveFxPosterioriButton;
		double minSampledFx;
		double maxSampledFx;
		double fxRange;
		double *fxValues;
		int fxSamples;
		double fxIncrement;

		// Fy group
		Slider *minFySlider;
		Slider *maxFySlider;
		Slider *fySampleNumberSlider;
		Fl_Button *saveFyPosterioriButton;
		double minSampledFy;
		double maxSampledFy;
		double fyRange;
		double *fyValues;
		int fySamples;
		double fyIncrement;

		// data
		bool valuesSampled;
		bool DFinferred;
		bool DDrinferred;

		CustomSelectionInferenceGui(int x, int y, int w, int h);

		// functions
		void save();
		void drawDiffusionPosterioriPlot();
		void sampleDiffusionPosteriori();
		void saveDiffusionPosteriori();
		void drawFxPosterioriPlot();
		void sampleFxPosteriori();
		void saveFxPosteriori();
		void drawFyPosterioriPlot();
		void sampleFyPosteriori();
		void saveFyPosteriori();
		void disableSelection();
		double getSigma() { return (double)noiseSigmaSlider->value()/1000.0; }

	private:
		double optimizationArray[3];
};

class StereoVisionGui {

	private:
		Fl_Window *window;
		Slider *eyeSeparationSlider;

	public:
		StereoVisionGui();
		void setEyeSeparation(float sep) { eyeSeparationSlider->value(sep); }
		float getEyeSeparation() { return eyeSeparationSlider->value(); }
		void show() { window->show(); }
		void hide() { window->hide(); }
		int visible() { return window->visible(); }
};

/* XPM of Smoothing Prior*/
//static char* prior_xpm[] = {
//"460 38 10 1",
//"  c #0a0a0a",
//"! c #d8d8d8",
//"# c #707070",
//"$ c #3e3e3e",
//"% c #989898",
//"& c #272727",
//"' c #898989",
//"( c #585858",
//") c #a0a0a0",
//"* c #181818",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%$'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!($%!!!!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'(#)!!!!!!%(#'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#(#!!!!!!!!'(#)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)($)!!!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'&#*#!!!!!)$($$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'&#*'!!!!!!%$(&'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)$(!!!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!($'&#!!!!!'&'$(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$('&'!!!!!!#&'$'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'&#!!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#&)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%*#!)!!!!!!$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!''')!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' '!)!!!!!!)&(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!''%!!!!!!!!!!!!!!!!!!!!!!!!!!(&%!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)&(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!# '!!!!!!!%*(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)(   &'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(*%!!!!!!!!'*#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#* &#!!!!!!!!!!!!!!!!!!!!!!!!!%&(!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#&%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$&)!!!!!!!# '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!# ('#*$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!&$!!!!!!!!!( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%&('$*'!!!!!!!!!!!!!!!!!!!!!!!!!(&%!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)&(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%&$!!!!!!!!(*'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$&)!!(*%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' (!!!!!!!!)$&%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#*%!' #!!!!!!!!!!!!!!!!!!!!!!!!!%&#!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'*'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' #!!!!!!!!&&!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'#!!)$&%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!# '!!!!!!!!%*$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'#!!# #!!!!!!!!!!!!!!!!!!!!!!!!!!$&)!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)$$)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!( '!!!!!!!' (!!!!!!!!!%()!##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)'!!!!!!!$)'$!!!!!'($$'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$*%!!!!!!!!# #!!!!!!!!!('!%$)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)'!!!!!!)$%!'(!!!!'$&#!!!!!!!!!!!!!!!!!!!!!!!!!!!'*#!",
//"!!'(((((('!!!)#%!!!%(((((('!!!!!#'!#####(&(!!!!!!!!!!!!!!#)!!!!!%#'!!!'#!)#####($!!!)#%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%&&)!!!!!!!( '!!!!!!!!!%$%!#(!!!!!!!!!!!!!!!!!!!#(((((#)!!!%#)!#####(&('!!!!!$%'$!!!!(*$'!!!!!!!!!!!!!)#'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' $!!!!!!!!)$ '!!!!!!!!!$#!%$)!!!!!!!!!!!!!!!!!'('!!!!!!##)!!'#!'####(&()!!!!)$%!#(!!!#&#)!!!!!!!!!!!!!!%('!!!!!!!!!!!)&$!",
//"!!$       (!!#&%!!!#       (!!!'&'!#####$&&'!!!!!!!!!!!!!&'!!!!!( #!!)$(!!#####$&'!!!$(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!# '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' $!!!!!!!!& '!!!!!!!!!%$%!#(!!%###########)!!!)&      (!!!$$!!#####$&&$)!!!!$%'$!!!# (!!!!!!!!!!!!!!!' #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)##!!!!!!!!!!!!!!!# $!!!!!!!!%*&%!!!!!!!!!$#!%$)!!###########'!!!# (!!!!!' $)!)$(!'####&$&(!!!!)$%!#(!!'&#!!!!!!!!!!!!!!!!# '!!!!!!!!!!!!(*%",
//"!)&*'%%'( &!%&#!!!!$ #%%'(* '!)$$!!!!!!!''$$)!!!!!!!!!!!!$#!!!!'*&)!!#&'!!!!!!!##$!!!'*'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(*%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!# (!!!!!!!' &)!!!!!!!!!%$%!#(!!' *&&&&&&&**%!!!' $%%''$&!!'&'!!!!!!!#''*'!!!!$%'$!!!$*((((#)!!!!!!!!!!# #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#&'!!!!!!!!!!!!!!( #!!!!!!!!' $!!!!!!!!!!$#!%$)!!$  &&&&&&& #!!!' $)!!!)$ '!!#&%!!!!!!#)#*'!!!)$%!#(!!# ((((%!!!!!!!!!!!!( '!!!!!!!!!!!!( '",
//"!' &)!!!%&&!#*%!!!!&*%!!!!( (!'*#!!!!!!))!'*'!!!!!!!!!!!!$(!!!!$ #!!!$$)!!!!!!))'*'!!)&(!!!!!!!!!!!!!!!!!!!!!!!!!!!!%)!!!!!!!!!!!!!!)%!!!!!!!!!!$&)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$ #!!!!!!!' $!!!!!!!!!!%$%!#(!!!&#%%%%%%%#(!!!!# (!!!!'*%!$$)!!!!!!))!!$$!!!!$%'$!!!$*****$%!!!!!!!!%!$*%!!!%)!!!!!!!!!!!!!!!!!!!!!!!!!)&$!!!!!!!!!!!!!!& '!!!!!!!!# (!!!!!!!!!!$#!%$)!!# $%%%%%%'$%!!!%*&%!!!# (!!)&(!!!!!!))!)&(!!!)$%!#(!!#****&'!!!!!!!!)%)!&&!!!))!!!!!!!!# #",
//"!# (!!!!'*&!&$!!!!' &!!!!!' $!(&%!!%(#(&#!)&(!!!!!!!!!!!!$$!!!' $!!!'*(!!!!#(#$$)&(!!!(*%!!!!!!!!!!!!!!!!!!!!!!!!!#$&&()!#$'!!##)(#(&&$'!!!!!!!)&&)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)& '!!!!!!!( (!!!!!!!!!!%$%!#(!!!($)!!!!!!$'!!!!( '!!!!)&''*(!!!!#(#$$)!#*'!!!$%'$!!!)%%%%%)!!!!!!!#$&($&)!#$&&(%!!!!!!!!!!!!!!!!!!!!!!!!(&%!!!!!!!!!!!!%**%!!!!!!!!$ #!!!!!!!!!!$#!%$)!!!& '!!!!!'#!!!!)&*%!!)&*'!!' #!!!((#&$!!(&)!!)$%!#(!!)%%%%%!!!!!!!)#$&$#&$!)#$&()!!!!!!' (",
//"!( #!!!'$ () #!!!!' $!!!!!' $)&$!!!#  *&'!!$&)!!!!!!!!!!!($!!)$ '!!!#*'!!!)&  &(!$&!!!' '!!!!!!!!)''%!!!%'%!!!!!!( *$& $)' $!!&'%*  $$ &'!!!!!!%*$!!!!!!!!!!!!!!!!!!!!!!)%!!))!!!!!!!!!!!' &)!!!!!!!& '!!!!!!!!!!%$%!#(!!!# '!!!!!#(!!!!!& '!!!!)&'# '!!!)&  &(!!' #!!!$%'$!!!!!!!!!!!!!!!!( &$* $!# &$& (!!!!!!!!!!!!!!!!!!!!!!!!'*#!!!!!!!!!!!!' &!!!!!!!!)& '!!!!!!!!!!$#!%$)!!!# (!!!!!#%!!!!!$ '!!# (!!!(*%!!)*  &(!!# '!!)$%!#(!!!!!!!!!!!!!!)$ &$*  #!( && (!!!!!!%*(",
//"!$ &$$$* $)% #!!!!( #!!!!!# (%*(!!!( &')!!!(*'!!!!!!!!!!!#&%!# (!!!!(*%!!!'* (%!!(*'!!' (!!!!!!!)($$&#!'$#'!!!!!'*&'!%$ '!( ##$!%* #!)# $!!!!!!% (!!!!!!!!!!!!!!!!!!!!!!#$!!##!!!!!!!!!!!# &!!!!!!!!& '!!!!!!!!!!%$%!#(!!!% $)!!!)$'!!!!%*&!!!!!)&'(*%!!!%* (%!!!%*(!!!$%'$!!!!!!!!!!!!!!!# $)!( #%&&'!'&$!!!!!!!!!!!!!!!!!!!!!!!!)$$)!!!!!!!!!!!# $!!!!!!!!' *%!!!!!!!!!!$#!%$)!!!!& '!!!!#!!!!!!( #!%&&%!!!$&)!!%  #)!!!' #!!)$%!#(!!!!!!!!!!!!!!( ()!'  '' $)%&(!!!!!!%*(",
//"%& &&&&&#)!' '!!!!$ '!!!!!( #' (!!!$ '!!!!!# #!!!!!!!!!!!# ')&&'!!!)$&%!!!' $!!!!# '!!)&$!!!!!!!##!)#*#()!!!!!!!# #!!!' (!)&  #!%*&)!!)&*%!!!!!' (!!!!!!!!!!!!!!!!!!!!!)$(!!('!!!!!!!!!!!( (!!!!!!!% &)!!!!!!!!!!%$%!#(!!!!$ '!!!#(!!!!!' $!!!!!'&)$&%!!!' $)!!!!)&&!!!$%'$!!!!!!!!!!!!!!)&*'!!# '' &')!))!!!!!!!!!!!!!!!!!!!!!!!!!(*'!!!!!!!!!!!$ #!!!!!!!!# $)!!!!!!!!!!$#!%$)!!!!' (!!!%#!!!!!!# (!( #!!!)&&)!!' (!!!!!%*$!!)$%!#(!!!!!!!!!!!!!%&&)!!!& '' &'!%)!!!!!!%*(",
//"' $)!!!!!!!' '!!!!&&)!!!!)&&%# #!!)&&)!!!!!' #!!!!!!!!!!!# '# (!!!!)&&)!!!( #!!!!' #!!!&&!!!!!!)(!!!!( #!!!!!!!!( &&&&& $!!'  %!%*$!!!!( '!!!!!' (!!!!!!!)(((((((((#)!!%&#!!$%!!!!!!!!!!!$ #!!!!!!!' $!!!!!!!!!!!%$%!#(!!!!'*$)!)$%!!!!!# #!!!!!#&)&&)!!!( #!!!!!!$&!!!$%'$!!!!!!!!!!!!!!' $!!!#*%!$  &$'!!!!!!#(((((((('!!!!!!!!!#$*(!!!!!!!!!!!& '!!!!!!!!( (!!!!!!!!!!!$'!'$)!!!!!$ '!!#'!!!!!!# ('&$)!!!%&$!!!# '!!!!!!&$!!)$%!#(!!!!!!!!!!!!!' #!!!!&&!)$  $'!!!!!!!%*(",
//"# (!!!!!!!!' '!!!' &!!!!!# (!' (!!) (!!!!!!' #!!!!!!!!!!!' #&&)!!!!)&&)!!!$*%!!!!' #!!!&$!!!!!!)(!!!!%&#!!!!!!!!( $$$$$$#!!#  '!%*$!!!!$ '!!!!!' (!!!!!!!)#########')!!#*'!!$!!!!!!!!!!!)& '!!!!!!!' (!!!!!!!!!!!%$%!#(!!!!!( #!#(!!!!!!$ '!!!!%&()$&%!!!$*%!!!!!!&&!!!$%'$!!!!!!!!!!!!!!# (!!!(&)!!'(&  '!!!!!'########'!!!!!!!!#$'$&)!!!!!!!!!%  '!!!!!!!!$ #!!!!!!!!!!!$'!%$)!!!!!' $!%$!!!!!!!% $$ #!!!!%&$!!!#&!!!!!!)&$!!)$%!#(!!!!!!!!!!!!!( #!!!' $!!)'(*&%!!!!!!%*(",
//"( #!!!!!!!!' '!!!# (!!!!# &%!' (!!' #!!!!!!# '!!!!!!!!!!!' * #!!!!!!$&%!!%&&)!!!!# '!!%*$!!!!!!)(!!!!'&$!!!!!!!!# #!!!%#%!' $&(!%*&%!!%&&%!!!!!% (!!!!!!!!!!!!!!!!!!!!!(&!!%$'%!!!!!!!!!' &)!!!!!!!# '!!!!!!!!!!!%$%!#(!!!!!'*&%$%!!!!!!&*%!!!'$ '!$*%!!)&&)!!!!!)&$!!!$%'$!!!!!!!!!!!!!!# (!!!$()(#!!!( #!!!!!!!!!!!!!!!!!!!!!!#$'!#*'!!!!!!!!!' &!!!!!!!!)& '!!!!!!!!!!!$#!%$)!!!!!!$ '(#!!!!!!!%*  $!!!!!)&&)!!($!!!!!!' (!!)$%!#(!!!!!!!!!!!!!( #!!!( #%('!!# '!!!!!!%*(",
//"$ '!!!!!!!!% #!!!( (##($ &'!!)&$!!' '!!!!!!(*%!'&(!!!!!!!%  $)!!!!!!(*'!!' (!!!!!$&)!!' #!!!!!!!#'!!'$#$#)!!!!!!%&&')%$ ')&*''&)%* (%%( (!!!!!!%*$!!!!!!!!!!!!!!!!!!!!)&$!!#(#%!!!!!!!!!' &!!!!!!!!( '!!!!!!!!!!!%$%!#(!!!!!!( $#!!!!!!' &###$* (!!( '!!' (!!!!!!' (!!!$%'$!!!!!!!!!!!!!!' &')# #' &')'& '!!!!!!!!!!!!!!!!!!!!!#&'!!%&(!!!!!!!!!# $!!!!!!!!%*&%!!!!!!!!!!!$#!%$)!!!!!!' &$)!!!!!!!)& *'!!!!!!$&)!!$#!!!!!!' #!!)$%!#(!!!!!!!!!!!!!# $%)#  '' $%)(&%!!!!!!%*(",
//"*&)!!!!!!!!)&(!!!$     &('!!!!$&)!#&%!!!!!!$$!!' &!!!!!!!)&&'!!!!!!!'*#!!# #!!!!!&$!!!# '!!!!!!!)(##(#!'&$'!!!!!!#* &* (!( #!!&#%* *** &'!!!!!!%&$)!!!!!!!!!!!!!!!!!!!'*&###(#!!!!!!!!!!( (!!!!!!!!$&)!!!!!!!!!!!%$%!#(!!!!!!'*&%!!!!!!'      &#!!!'*#!!# #!!!!!!' '!!!$%'$!!!!!!!!!!!!!!)$ *&**')$ *&* (!!!!!!!!!!!!!!!!!!!!!%$#!!!!(&%!!!!!!!!$ #!!!!!!!!# $!!!!!!!!!!!!$#!%$)!!!!!!!$ #!!!!!!!!!(*(!!!!!!!#*'!!$'!!!!!!# '!!)$%!#(!!!!!!!!!!!!!)& *&**&%%& *&*#!!!!!!!' (",
//"'%!!!!!!!!!!$&)!!%'''''%!!!!!!#*'!'#!!!!!!'*#!!!#$!!!!!!!!%%!!!!!!!!%&(!!'(%!!!!)*#!!!$&!!!!!!!!!)''%!!!%'%!!!!!!!'($('!!##!!!'#%*$'($(%!!!!!!!)$&)!!!!!!!!!!!!!!!!!!!#&''%!%)!!!!!!!!!!$ #!!!!!!!!$$!!!!!!!!!!!!%$%!#(!!!!!!!##!!!!!!!)'''''')!!!!)&(!!'(%!!!!!!(&)!!!$%'$!!!!!!!!!!!!!!!%($(''!!%($$('!!!!!!!!!!!!!!!!!!!!!!!%!!!!!)'!!!!!!!!!& '!!!!!!!!# (!!!!!!!!!!!!$#!%$)!!!!!!!'()!!!!!!!!!)%)!!!!!!!'*#!!#)!!!!!!$$!!!)$%!#(!!!!!!!!!!!!!!%($('''!!'($$#!!!!!!!!# #",
//"!!!!!!!!!!!!'*'!!!!!!!!!!!!!!!%&(!!!!!!!!!#&%!!!((!!!!!!!!!!!!!!!!!!!(&%!!!!!!!!'&%!!%*#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%*$!!!!!!!!!!!!!(*%!!!!!!!!!!!!!!!!!!!($!!!!!!!!!!!!!!!)&*'!!!!!!!!&#!!!!!!!!!!!!%$%!#(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(&%!!!!!!!!!%&(!!!!$%'$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%*&)!!!!!!!!$ #!!!!!!!!!!!!$#!%$)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$&)!!!!!!!!'*#!!!)$%!#(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!( '",
//"!!!!!!!!!!!!!$$!!!!!!!!!!!!!!!!(&%!!!!!!!)&(!!!%(%!!!!!!!!!!!!!!!!!!!'&#!!!!!!!!#(!!!#&)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%*$!!!!!!!!!!!!!(*'!!!!!!!!!!!!!!!!!!)&(!!!!!!!!!!!!!!!' &)!!!!!!!!&'!!!!!!!!!!!!%$%!#(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%&#!!!!!!!!!#&'!!!!$%'$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' $!!!!!!!!)&*'!!!!!!!!!!!!$#!%$)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'&'!!!!!!!!(&)!!!)$%!#(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(*%",
//"!!!!!!!!!!!!!'$'!!!!!!!!!!!!!!!)$#!!!!!!!'$%!!!!!!!!!!!!!!!!!!!!!!!!!!#$)!!!!!!!(%!!)$#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%&(!!!!!!!!!!!!!' #!!!!!!!!!!!!!!!!!!!#%!!!!!!!!!!!!!!!' $!!!!!!!!!&%!!!!!!!!!!!!%$%!#(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#$)!!!!!!!)$#!!!!!$%'$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!# (!!!!!!!!%*&)!!!!!!!!!!!!$#!%$)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!((!!!!!!!)$'!!!!)$%!#(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)$&)",
//"!!!!!!!!!!!!!!)!!!!!!!!!!!!!!!!!))!!!!!!!))!!!!!!!!!!!!!!!!!!!!!!!!!!!!)!!!!!!!!)!!!!)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!))!!!!!!!!!!!!!%*(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!# #!!!!!!!!%&!!!!!!!!!!!!!%$%!#(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)!!!!!!!!!)!!!!!!$%'$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!( '!!!!!!!!' (!!!!!!!!!!!!!$#!%$)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)!!!!!!!!)!!!!!)$%!#(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%*(!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(&%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$*'!!!!!!!!#$!!!!!!!!!!!!!)#!!''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#!)#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!&&)!!!!!!!!( '!!!!!!!!!!!!!'%!)#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#)!%'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#&'!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%&(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%*$!!!!!!!!!$(!!!)#$$()!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' (!!!!!!!!)&&)!!!%($$'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%&(!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(&)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' #!!!!!!!!!&#!!!(*$$ (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!# '!!!!!!!!%&(!!!)$*$&$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#*'!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'&#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(*'!!!!!!!!' '!!%*(!)(#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$&)!!!!!!!!#*'!!!'*#!%(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)$$!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(&)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$$!!!!!!!!!#*)!!'*&(')!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!&(!!!!!!!!!(&)!!!' &#'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#&'!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%&'!!!!!!!!!!!!!!!!!!!!!!!!!!!!''!'&#!!!!!)#%!$$!!!!'$* $%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!''%*'!!!!!)#))&(!!!!!#$* #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%&(!!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'(!!!!!!!!!!!!!!!!!!!!!!!!!!!!(&'($)!!!!!' #'&'!!!#'!%( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$&#$!!!!!!# '#&%!!!)#'!'$(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!($)!!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!('!!!!!!!!!!!!!!!!!!!!!!!!!!!#$($'!!!!!!%$#$#!!!)&$%'$*'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#$(%!!!!!!'$#$'!!!!'*$%'$(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#$'!!!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'$%!!!!!!!!!!!!!!!!!!!!!!!!!!!%')!!!!!!!!)'%!!!!!( &&&#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%!!!!!!!!%'%!!!!!)$ &&&'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#$'!!!!!!",
//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'##%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)'##%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#'!!!!!!!"};

//static Fl_Pixmap priorImage(prior_xpm);

#endif /* GUI_H_ */

