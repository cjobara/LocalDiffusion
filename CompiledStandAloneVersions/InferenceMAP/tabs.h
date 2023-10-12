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

#ifndef TABS_H_
#define TABS_H_

// FLTK LIBRARIES (GUI)
#include <FL/Fl.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Menu.H>
#include <FL/FL_Menu_Window.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Menu_Window.H>
#include <FL/Fl_Box.H>

class PopupWindow;
class Tabs;

void closeFile(Tabs*t);
void selectTab(Tabs*t);

class PopupWindow : public Fl_Menu_Window {
    Fl_Box *output;
    // Size window to just fit output's label text
    void SizeToText() {
        int W=0, H=0;
        fl_font(output->labelfont(), output->labelsize());
        fl_measure(output->label(), W, H, 0);
        resize(x(), y(), W+10, H+10);                           // +10: leaves +5 margin on all sides
        output->resize(0, 0, W+10, H+10);
    }
	public:
		PopupWindow() : Fl_Menu_Window(10,10) {
			output = new Fl_Box(0, 0, w(), h());                    // box will have the text of user's msg
			output->box(FL_UP_BOX);                                 // popup window will have an 'Up Box' border
			end();
			hide();
			border(0);                                              // popup will be borderless
			output->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);           // text should be left aligned
			output->label("No text defined");                       // (default msg if none defined)
			SizeToText();
		}
		// Change text in box
		void text(const char*s) {
			output->label(s);                                       // set message text
			SizeToText();                                           // resize window to size of text
		}
		// Pop up window at current mouse position
		void popup() {
			position(Fl::event_x_root(), Fl::event_y_root());       // position window at cursor
			show();
		}
};

class Tabs : public Fl_Button {

    static void closeTabCallback(Fl_Widget*, void *userdata) {
        Tabs *tab = (Tabs*)userdata;
		closeFile(tab);
    }
    static void selectTabCallback(Fl_Widget*, void *userdata) {
        Tabs *tab = (Tabs*)userdata;
		selectTab(tab);
    }

	int handle(int e) {
		switch (e) {
			case FL_PUSH:
				if ( Fl::event_button() == FL_RIGHT_MOUSE ) {
					Fl_Menu_Item rclick_menu[] = {
						{this->label(),0,selectTabCallback,(void*)this,FL_MENU_DIVIDER},
						{ "Close", 0, closeTabCallback, (void*)this },
						{ 0 }
					};
					const Fl_Menu_Item *m = rclick_menu->popup(Fl::event_x(), Fl::event_y(), 0, 0, 0);
					if ( m ) m->do_callback(0, m->user_data());
					return(1);          // (tells caller we handled this event)
				}
				break;
		}
		return(Fl_Button::handle(e));    // let Fl_Input handle all other events
	}

	public:
		Tabs(int x,int y,int w,int h,const char*l=0);

};

#endif /* TABS_H_ */
