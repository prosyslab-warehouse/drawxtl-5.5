// $Id: Draw_Fl_Input.h 900 2009-08-13 20:00:45Z larry $
//
#ifndef Draw_Fl_Input_h
#define Draw_Fl_Input_h

#include <FL/Fl_Input.H>

class Draw_Fl_Input:public Fl_Input
{
  public:
  Draw_Fl_Input (int x, int y, int w, int h, const char *l = 0):
    Fl_Input (x, y, w, h, l) {
    }
    int handle (int e);
};

#endif
