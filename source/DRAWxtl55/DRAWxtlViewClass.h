// $Id: DRAWxtlViewClass.h 900 2009-08-13 20:00:45Z larry $
//
#ifndef DRAWxtlViewClass_h
#define DRAWxtlViewClass_h

class Draw_Fl_Input:public Fl_Input
{
  public:
  Draw_Fl_Input (int x, int y, int w, int h, const char *l = 0):
    Fl_Input (x, y, w, h, l) {
    }
  private:
    int handle (int e);
};

#endif
