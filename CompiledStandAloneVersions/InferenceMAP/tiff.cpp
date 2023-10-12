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
#include "tiff.h"
#include <FL/fl_draw.H>
#include "file.h"

extern Globals *iMAP;
extern File *file;

#include "inference.h"

// TIFF DEFINITIONS
int writetiff(char *filename, char *description, int x, int y, int width, int height, int compression)
{
        TIFF *file;
        GLubyte *image, *p;
        int i;

        file = TIFFOpen(filename, "w");
        if (file == NULL) {
                return 1;
        }
        image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);

        glPixelStorei(GL_PACK_ALIGNMENT, 1);

        glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
        TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
        TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
        TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(file, TIFFTAG_COMPRESSION, compression);
        TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
        TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
        TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
        TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, description);
        p = image;
        for (i = height - 1; i >= 0; i--) {
        if (TIFFWriteScanline(file, p, i, 0) < 0) {
          free(image);
          TIFFClose(file);
          return 1;
        }
        p += width * sizeof(GLubyte) * 3;
        }
        TIFFClose(file);
        return 0;
}

void captureScreen(char* captureScreenName) {
        glViewport(0,0,iMAP->glWindow->w(),iMAP->glWindow->h());
        GLint viewport[4];                                      // Where The Viewport Values Will Be Stored
        glGetIntegerv(GL_VIEWPORT, viewport);                   // Retrieves The Viewport Values (X, Y, Width, Height)

        char filename [ FILENAME_MAX ];
        sprintf(filename, "%s%s", captureScreenName,".tif");

        // compression options: COMPRESSION_NONE, COMPRESSION_LZW, COMPRESSION_PACKBITS
        writetiff(filename, "OpenGL Screen Capture", 0, 0, viewport[2], viewport[3], COMPRESSION_NONE);
}

void writetiffOffscreen(char *offscreenFileName, Fl_Widget* widget, char* extension) {
        // storing recent widget size:
        int X = widget->x();
        int Y = widget->y();
        int W = widget->w();
        int H = widget->h();

        Fl_Offscreen o = fl_create_offscreen(W, H);
        fl_begin_offscreen(o);

        static unsigned int bufferSize=0;
        static unsigned char * buffer=0;

        if (bufferSize<3*(unsigned int)W*(unsigned int)H){
                bufferSize=3*W*H;
                if (buffer) delete[] buffer;
                buffer = new unsigned char[bufferSize];
        }

        iMAP->mainWindow->make_current();

        fl_read_image(buffer, X, Y, W, H, 0);

        fl_end_offscreen();
        fl_delete_offscreen(o); // cleanup

//       Define an image
        TIFF *image;

        char filename [ FILENAME_MAX ];
        sprintf(filename, "%s%s%s", offscreenFileName,extension,".tif");

        // Open the TIFF file
        if((image = TIFFOpen(filename, "w")) == NULL){
        printf("Could not open output.tif for writing\n");
        exit(42);
        }

        // We need to set some values for basic tags before we can add any data
        TIFFSetField(image, TIFFTAG_IMAGEWIDTH, W);
        TIFFSetField(image, TIFFTAG_IMAGELENGTH, H);
        TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(image, TIFFTAG_PLANARCONFIG,1);
        TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
        TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 3);
        TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, 1);
        TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, "");

    for(size_t i=0; i<(unsigned int)H; i++) {
       TIFFWriteScanline(image, buffer+i*(unsigned int)W*3, i);
    }

        // Close the file
        TIFFClose(image);
}

void saveScreenSequence() {
	char sequenceName[FILENAME_MAX];
	strcpy(sequenceName,iMAP->captureScreenName);
	strcat(sequenceName,"_");
	char screenNumber[10];
	sprintf(screenNumber,"%04i",iMAP->tiffCount);
	strcat(sequenceName,screenNumber);

	captureScreen(sequenceName);
	iMAP->tiffCount++;
}

OverlayTiffGui::OverlayTiffGui() {

	if (iMAP->fileNumber > 0) {

		const int w = 240;
		const int h = 290;

		count = 0;
		playEnable = false;

		window = new Fl_Window(w,h,"Overlay TIFF");
		window->callback((Fl_Callback*)overlayTiffWindowCallback);
		window->color(iMAP->bgColor);
		window->begin();
		window->set_non_modal();

		tabs = new Fl_Tabs(0,0,w,h);
		tabs->labelsize(12);
		tabs->box(FL_BORDER_FRAME);
		tabs->labelcolor(FL_WHITE);
		tabs->labelfont(iMAP->normalFont);
		tabs->begin();

		geometryGroup = new Fl_Group(0,25,w,h,"Geometry");
		geometryGroup->labelsize(12);
		geometryGroup->labelcolor(FL_DARK1);
		geometryGroup->labelfont(iMAP->normalFont);
		geometryGroup->box(FL_BORDER_FRAME);
		geometryGroup->begin();
		{
			resolutionSlider = new Slider(10,45,w-20,20,"Resolution [nm/px]");
			resolutionSlider->precision(0);
			resolutionSlider->labelsize(12);
			resolutionSlider->labelfont(iMAP->boldFont);
			resolutionSlider->labelcolor(FL_WHITE);
			resolutionSlider->value(160);
			resolutionSlider->bounds(100,300);
			resolutionSlider->show();

			xFlipButton = new Fl_Button(10,75,w/2-10,25,"x-Flip");
			xFlipButton->labelsize(12);
			xFlipButton->labelfont(iMAP->boldFont);
			xFlipButton->type(FL_TOGGLE_BUTTON);
			xFlipButton->callback((Fl_Callback*)xFlipButtonCallback);
			xFlipButton->show();

			yFlipButton = new Fl_Button(w/2,75,w/2-10,25,"y-Flip");
			yFlipButton->labelsize(12);
			yFlipButton->labelfont(iMAP->boldFont);
			yFlipButton->type(FL_TOGGLE_BUTTON);
			yFlipButton->callback((Fl_Callback*)yFlipButtonCallback);
			yFlipButton->show();

			xAlignSlider = new Slider(10,120,w/2-12,20,"x-Align [nm]");
			xAlignSlider->precision(0);
			xAlignSlider->labelsize(12);
			xAlignSlider->labelfont(iMAP->boldFont);
			xAlignSlider->labelcolor(FL_WHITE);
			xAlignSlider->value(0);
			xAlignSlider->bounds(-4000,4000);
			xAlignSlider->show();

			yAlignSlider = new Slider(w/2+2,120,w/2-12,20,"y-Align [nm]");
			yAlignSlider->precision(0);
			yAlignSlider->labelfont(iMAP->boldFont);
			yAlignSlider->labelcolor(FL_WHITE);
			yAlignSlider->labelsize(12);
			yAlignSlider->value(0);
			yAlignSlider->bounds(-4000,4000);
			yAlignSlider->show();

			transposeButton = new Fl_Button(10,150,w-20,25,"Transpose");
			transposeButton->labelsize(12);
			transposeButton->labelfont(iMAP->boldFont);
			transposeButton->type(FL_TOGGLE_BUTTON);
			transposeButton->callback(adjustImageCallback);
			transposeButton->show();
		}
		geometryGroup->end();

		imageGroup = new Fl_Group(0,25,w,h,"Image");
		imageGroup->labelsize(12);
		imageGroup->labelcolor(FL_DARK1);
		imageGroup->labelfont(iMAP->normalFont);
		imageGroup->box(FL_BORDER_FRAME);
		imageGroup->begin();
		{
			alphaSlider = new Slider(10,45,w-20,20,"Alpha");
			alphaSlider->precision(2);
			alphaSlider->labelsize(12);
			alphaSlider->labelfont(iMAP->boldFont);
			alphaSlider->labelcolor(FL_WHITE);
			alphaSlider->value(1.0);
			alphaSlider->bounds(0.0,1.0);
			alphaSlider->callback(adjustImageCallback);
			alphaSlider->show();

			contrastSlider = new Slider(10,85,w-20,20,"Contrast");
			contrastSlider->precision(2);
			contrastSlider->labelsize(12);
			contrastSlider->labelfont(iMAP->boldFont);
			contrastSlider->labelcolor(FL_WHITE);
			contrastSlider->value(0);
			contrastSlider->bounds(0.0,200);
			contrastSlider->callback(adjustImageCallback);
			contrastSlider->show();

			gammaSlider = new Slider(10,125,w-20,20,"Gamma");
			gammaSlider->precision(2);
			gammaSlider->labelsize(12);
			gammaSlider->labelfont(iMAP->boldFont);
			gammaSlider->labelcolor(FL_WHITE);
			gammaSlider->value(1.0);
			gammaSlider->bounds(0.001,5);
			gammaSlider->callback(adjustImageCallback);
			gammaSlider->show();

			brightnessSlider = new Slider(10,165,w-20,20,"Brightness");
			brightnessSlider->precision(0);
			brightnessSlider->labelsize(12);
			brightnessSlider->labelfont(iMAP->boldFont);
			brightnessSlider->labelcolor(FL_WHITE);
			brightnessSlider->value(0);
			brightnessSlider->bounds(1,iMAP->tiffSequence->saturation);
			brightnessSlider->callback(adjustImageCallback);
			brightnessSlider->show();

			minSlider = new Slider(10,205,w/2-10,20,"Min");
			minSlider->precision(0);
			minSlider->labelsize(12);
			minSlider->labelfont(iMAP->boldFont);
			minSlider->labelcolor(FL_WHITE);
			minSlider->value(1);
			minSlider->bounds(1,iMAP->tiffSequence->saturation);
			minSlider->callback(adjustImageCallback);
			minSlider->show();

			maxSlider = new Slider(w/2,205,w/2-10,20,"Max");
			maxSlider->precision(0);
			maxSlider->labelsize(12);
			maxSlider->labelfont(iMAP->boldFont);
			maxSlider->labelcolor(FL_WHITE);
			maxSlider->value(65384);
			maxSlider->bounds(1,iMAP->tiffSequence->saturation);
			maxSlider->callback(adjustImageCallback);
			maxSlider->show();

			normalizeButton = new Fl_Button(10,235,w-20,25,"Normalize");
			normalizeButton->labelsize(12);
			normalizeButton->labelfont(iMAP->boldFont);
			normalizeButton->type(FL_TOGGLE_BUTTON);
			normalizeButton->callback(adjustImageCallback);
			normalizeButton->show();

			clampButton = new Fl_Check_Button(10,265,w/2-60,20,"Clamp");
			clampButton->labelcolor(FL_WHITE);
			clampButton->labelfont(iMAP->boldFont);
			clampButton->labelsize(12);
			clampButton->show();
		}
		imageGroup->end();

		sequenceGroup = new Fl_Group(0,25,w,h,"Sequence");
		sequenceGroup->labelsize(12);
		sequenceGroup->labelcolor(FL_DARK1);
		sequenceGroup->labelfont(iMAP->normalFont);
		sequenceGroup->box(FL_BORDER_FRAME);
		sequenceGroup->begin();
		{
			sequenceSlider = new Slider(10,45,w-20,20,"Image Number");
			sequenceSlider->precision(0);
			sequenceSlider->labelsize(12);
			sequenceSlider->labelfont(iMAP->boldFont);
			sequenceSlider->labelcolor(FL_WHITE);
			sequenceSlider->bounds(1,iMAP->tiffSequence->stackSize);
			sequenceSlider->value(0);
			sequenceSlider->callback((Fl_Callback*)sequenceSliderCallback);
			sequenceSlider->show();

//			playButton = new Fl_Button(10,75,50,25,"Play");
//			playButton->labelcolor(iMAP->bgColor);
//			playButton->labelfont(iMAP->boldFont);
//			playButton->label("@>");
//			playButton->labelsize(12);
//			playButton->callback((Fl_Callback*)sequencePlayCallback, (void*)this);

//			stopButton = new Fl_Button(70,75,50,25,"Stop");
//			stopButton->labelcolor(iMAP->bgColor);
//			stopButton->labelfont(iMAP->boldFont);
//			stopButton->label("@square");
//			stopButton->labelsize(12);
//			stopButton->callback((Fl_Callback*)sequenceStopCallback, (void*)this);

			synchronizeButton = new Fl_Check_Button(10,75,80,25,"Synchronize");
			synchronizeButton->labelsize(12);
			synchronizeButton->labelfont(iMAP->boldFont);
			synchronizeButton->labelcolor(FL_WHITE);
			synchronizeButton->callback((Fl_Callback*)syncButtonCallback);
			synchronizeButton->show();

			// disable playback functions if only one TIFF opened
			if (count == 1) {
				sequenceSlider->deactivate();
				playButton->deactivate();
				stopButton->deactivate();
			}

		}
		sequenceGroup->end();
		tabs->end();
		window->end();
		window->show();
	}
}

void overlayTiffWindowCallback(Fl_Window*w,void*) {
	w->hide();
	iMAP->synchronizeTiff = false;
	iMAP->tiffOverlay = false;
	delete iMAP->overlayTiffGui;
	iMAP->overlayTiffGui = NULL;
}

void adjustImageCallback(Fl_Widget*,void*) {
	iMAP->tiffSequence->adjust();
}

void xFlipButtonCallback(Fl_Button*w, void*v) {
	if (w->value()) {
		iMAP->tiffSequence->xFlip = 1.0;
	} else {
		iMAP->tiffSequence->xFlip = 0.0;
	}
}

void yFlipButtonCallback(Fl_Button*w, void*v) {
	if (w->value()) {
		iMAP->tiffSequence->yFlip = 1.0;
	} else {
		iMAP->tiffSequence->yFlip = 0.0;
	}
}

void sequenceSliderCallback(Slider*w,void*v) {
	iMAP->tiffSequence->setSlice((int)(w->value()-1));
	iMAP->tiffSequence->adjust();
	Fl::redraw();
}

void clampButtonCallback(Fl_Check_Button*w,void*v) {

}

void getTiffAttributes(char *name, int *width, int *height, unsigned short *bitsPerSample, unsigned short *samplesPerPixel, int *n) {
		TIFFSetWarningHandler(NULL);

        TIFF* tif = TIFFOpen(name, "r");
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, height);
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, width);
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, bitsPerSample);
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, samplesPerPixel);
        if (!*samplesPerPixel) { *samplesPerPixel=1; } /* assume one sample per pixel if nothing is specified */

        // count images in stack
        int dircount = 0;
    if (tif) {
                do { dircount++; }
                while (TIFFReadDirectory(tif));
				if (iMAP->updateDisplay) {
					iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
					char fileName[300];
					sprintf(fileName,"%s\n",name);
					textDisplayUpdate(fileName);
					char progress[100];
					sprintf(progress,"%d-Image Stack\n",dircount);
					textDisplayUpdate(progress);
				}
                TIFFClose(tif);
        }
        *n = dircount;
}

void readTiffImage(char *name, int width, int height, unsigned short bitsPerSample, unsigned int samplesPerPixel, int slice, GLuint *texture) {

	TIFF *tif = NULL;
	unsigned short *buffer;

	buffer = new unsigned short [width*height];

	// set to desired directory of image stack
	TIFFSetWarningHandler(NULL);

	tif = TIFFOpen(name, "r");
	TIFFSetDirectory(tif, slice);

	// elimintate warnings

	glGenTextures(1, texture);

	for (uint32 r = 0; (int)r < height; r++) { TIFFReadScanline(tif, &buffer[r*width], r); }

	glBindTexture(GL_TEXTURE_2D, *texture);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
	glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_INTENSITY, width, height, 0, GL_LUMINANCE, GL_UNSIGNED_SHORT, buffer);

	delete [] buffer;
	TIFFClose(tif);
}

TiffSequence::TiffSequence(char *loadedFilename) {

	#ifdef __APPLE__
//	strcpy(filename,iMAP->cwd);
	strcat(filename,loadedFilename);
	#elif _WIN32
	strcpy(filename,loadedFilename);
	#endif

	fprintf(stderr,"%s\n",loadedFilename);

	slice = 0;
	xFlip = yFlip = 0.0;

	getTiffAttributes(filename,&width,&height,&bitsPerSample,&samplesPerPixel,&stackSize);
	readTiffImage(filename,width,height,bitsPerSample,samplesPerPixel,slice,&texture);

    saturation = 1;
    for (int b = 0; b < bitsPerSample; b++) { saturation *= 2; }

}

void TiffSequence::draw(int slice) {

    resolution = iMAP->overlayTiffGui->resolutionSlider->value()/1000.0;

    xMin = - (float)width/2.0*resolution;
    xMax = - xMin;
    yMin = - (float)height/2.0*resolution;
    yMax = - yMin;

    // apply alignment translation
    xMin += iMAP->overlayTiffGui->xAlignSlider->value()/1000.0;
    xMax += iMAP->overlayTiffGui->xAlignSlider->value()/1000.0;
    yMin += iMAP->overlayTiffGui->yAlignSlider->value()/1000.0;
    yMax += iMAP->overlayTiffGui->yAlignSlider->value()/1000.0;

    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);

	// Draw a textured quad
	glPushMatrix();
		updateView();

		glPushAttrib(GL_COLOR_BUFFER_BIT   |
				     GL_DEPTH_BUFFER_BIT   |
				     GL_ENABLE_BIT         |
				     GL_LIGHTING_BIT       |
				     GL_POLYGON_BIT        |
				     GL_TEXTURE_BIT);

		glDisable(GL_LIGHTING);
		glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		// enable alpha blending
		glEnable   (GL_BLEND);
		glDepthMask(GL_FALSE);
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

		if (iMAP->overlayTiffGui->clampButton->value()) { clip(); }

		glColor4f(1.0,1.0,1.0,iMAP->overlayTiffGui->alphaSlider->value());
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, texture);

        glBegin( GL_QUADS ) ;

			glTexCoord2f( xFlip, yFlip );
			glVertex2d( xMin, yMin );

			glTexCoord2f( 1.0-xFlip, yFlip );
			glVertex2d( xMax, yMin );

			glTexCoord2f( 1.0-xFlip, 1.0-yFlip );
			glVertex2d( xMax, yMax );

			glTexCoord2f( xFlip, 1.0-yFlip );
			glVertex2d( xMin, yMax );

        glEnd() ;

		unclip();

	glPopMatrix();

    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);

    glBindTexture(GL_TEXTURE_2D,0);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    glDisable (GL_TEXTURE_2D);
    glPopAttrib() ;
}

// sets the clipping planes (uses the 6 first)
void TiffSequence::clip() const {
	// clip the 6 faces of the cube
	GLdouble plane[4] ;

	plane[0] = +1. ;  plane[1] =  0. ;  plane[2] =  0. ;  plane[3] = -file->xMin ;
	glEnable( GL_CLIP_PLANE0 ) ;
	glClipPlane( GL_CLIP_PLANE0, plane ) ;

	plane[0] = -1. ;  plane[1] =  0. ;  plane[2] =  0. ;  plane[3] = file->xMax ;
	glEnable( GL_CLIP_PLANE1 ) ;
	glClipPlane( GL_CLIP_PLANE1, plane ) ;

	plane[0] =  0. ;  plane[1] = +1. ;  plane[2] =  0. ;  plane[3] = -file->yMin ;
	glEnable( GL_CLIP_PLANE2 ) ;
	glClipPlane( GL_CLIP_PLANE2, plane ) ;

	plane[0] =  0. ;  plane[1] = -1. ;  plane[2] =  0. ;  plane[3] = file->yMax ;
	glEnable( GL_CLIP_PLANE3 ) ;
	glClipPlane( GL_CLIP_PLANE3, plane ) ;
}

// unsets the clipping planes
void TiffSequence::unclip() const {
  // disable cube clip plane
  glDisable( GL_CLIP_PLANE0 ) ;
  glDisable( GL_CLIP_PLANE1 ) ;
  glDisable( GL_CLIP_PLANE2 ) ;
  glDisable( GL_CLIP_PLANE3 ) ;
}

void TiffSequence::adjust() {
	TIFF *tif = NULL;
	unsigned short *buffer;

	buffer = new unsigned short [width*height];

	// set to desired directory of image stack
	tif = TIFFOpen(filename, "r");
	TIFFSetDirectory(tif, slice);

	glGenTextures(1, &texture);

	for (uint32 r = 0; (int)r < height; r++) { TIFFReadScanline(tif, &buffer[r*width], r); }

//    int saturation = 1;
//    for (int b = 0; b < bitsPerSample; b++) { saturation *= 2; }

    // adjust brightness
    if (iMAP->overlayTiffGui->brightnessSlider->value() != 0) {
		const int brightness = (int)iMAP->overlayTiffGui->brightnessSlider->value();
		for (int v = 0; v < width*height; v++) {
				if ((int)(buffer[v])+(int)(brightness) >= saturation) { buffer[v] = saturation-1; }
				else if ((int)(buffer[v])+(int)(brightness) < 0) { buffer[v] = 0; }
				else { buffer[v] += brightness; }
		}
    }

    // adjust contrast
    if (iMAP->overlayTiffGui->contrastSlider->value() != 0) {
		float contrast = iMAP->overlayTiffGui->contrastSlider->value();
		float pixel = 0;
		contrast = (100.0+contrast)/100.0;
		contrast *= contrast;
		for (int v = 0; v < width*height; v++) {
				pixel = (float)buffer[v]/(float)saturation;
				pixel -= 0.5;
				pixel *= contrast;
				pixel += 0.5;
				pixel *= (float)saturation;
				if ((int)(pixel) > saturation) { buffer[v] = saturation-1; }
				else if ((int)(pixel) < 0) { buffer[v] = 0; }
				else { buffer[v] = (int)pixel; }
		}
    }

    // adjust min
    if (iMAP->overlayTiffGui->minSlider->value() != 0) {
		const int min = (int)iMAP->overlayTiffGui->minSlider->value();
		for (int v = 0; v < width*height; v++) {
				if (buffer[v] < min) { buffer[v] = 0; }
		}
    }

    // adjust max
    if (iMAP->overlayTiffGui->maxSlider->value() != saturation) {
		const int max = (int)iMAP->overlayTiffGui->maxSlider->value();
		for (int v = 0; v < width*height; v++) {
				if (buffer[v] > max) { buffer[v] = saturation-1; }
		}
    }

    // adjust gamma
    if (iMAP->overlayTiffGui->gammaSlider->value() != 1) {
		int *gamma;
		gamma = new int[saturation];
		for (int v = 0; v < saturation; v++) {
				gamma[v] = (int)( (saturation-1)*powf( (float)v/(float)(saturation-1), (float)1/(float)iMAP->overlayTiffGui->gammaSlider->value() ) );
		}
		for (int v = 0; v < width*height; v++) {
				buffer[v] = gamma[(int)buffer[v]];
				if (buffer[v] >= saturation) { buffer[v] = saturation-1; }
		}
		delete [] gamma;
    }

    // normalize
    if (iMAP->overlayTiffGui->normalizeButton->value()) {
		int max = -1;
		for (int v = 0; v < width*height; v++) {
				if ( max < (int)buffer[v] ) { max = buffer[v]; }
		}
		for (int v = 0; v < width*height; v++) {
				buffer[v] = (int)( (float)buffer[v]/(float)(max)*(float)(saturation-1) );
		}
    }

    // apply rotation
    if (iMAP->overlayTiffGui->transposeButton->value()) {
		unsigned short *bufferTemp = new unsigned short[width*height];

		for (int b = 0; b < width; b++) {
			for (int a = 0; a < height; a++) {
					bufferTemp[a + b*height] = buffer[a*width + b];
			}
		}

		for (int d = 0; d < width*height; d++) {
				buffer[d] = bufferTemp[d];
		}

		int temp = width;
		width = height;
		height = temp;

		delete [] bufferTemp;
    }

	glBindTexture(GL_TEXTURE_2D, texture);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
	glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_INTENSITY, width, height, 0, GL_LUMINANCE, GL_UNSIGNED_SHORT, buffer);

	delete [] buffer;
	TIFFClose(tif);
}

void syncButtonCallback(Fl_Check_Button*w,void*v) {
	if (w->value()) { 
		iMAP->synchronizeTiff = true; 
		iMAP->overlayTiffGui->sequenceSlider->deactivate();
	}
	else { 
		iMAP->synchronizeTiff = false;
		iMAP->overlayTiffGui->sequenceSlider->activate();
	}
	Fl::redraw();
}
