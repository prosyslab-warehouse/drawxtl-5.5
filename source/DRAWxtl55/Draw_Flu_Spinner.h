// $Id: Draw_Flu_Spinner.h 900 2009-08-13 20:00:45Z larry $
//
#ifndef Draw_Flu_Spinner_h
#define Draw_Flu_Spinner_h

#include "Flu_Spinner.h"
#include "Draw_Fl_Input.h"

class Draw_Flu_Spinner:public Flu_Spinner
{
  public:
  Draw_Flu_Spinner (int x, int y, int w, int h, const char *l = 0):
    Flu_Spinner (x, y, w, h, l) {
    }
    int handle (int e);
};

#endif
