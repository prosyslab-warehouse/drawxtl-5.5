#
# DRAWxtl V5.5 makefile using the Fast Light Tool Kit (FLTK).
#
# Copyright 2007-2011 by Larry Finger, Martin Kroeker, and Brian Toby
#

OPTIM = error
TARGET = ../../exe/DRAWxtl55
PREINC = 
RESFILE =

# post processing command - overwritten for Mac
POSTLINK = touch

OS=$(shell uname)

# Cygwin produces a string of the form CYGWIN_dd-x.y where dd is NT 95 98 or ME and
#    x and y are the major and minor versions, respectively. The next 'if' will effectively
#    truncate that to CYGWIN, no matter which version of Windows is being used. Other OS's
#    will be unaffected.

ifneq (,$(findstring CYGWIN,$(OS)))
    OS=CYGWIN
endif

# Map all BSD variants to Linux for now (seems to work for the Debian folks at least)
ifneq (,$(findstring BSD,$(OS)))
    OS=Linux
endif

ifeq ($(OS),Linux)

# section to define the various symbols needed to build DRAWxtl for Linux systems

  prefix          = /usr/local
  includedir      = ${prefix}/include
  libdir          = /usr/lib/

  BUILD_CPU=$(shell uname -m)
  ifneq (,$(findstring 64,$(BUILD_CPU)))
      libdir	  = /usr/lib64/
  endif

# compiler name:
  CXX             = g++

# flags for C++ compiler:
  OPTIM           = -g -Wall -Wunused  -fno-exceptions
# work around a bug in freeglut 2.4 - only enable this if you have to use a version 
# of fltk before 1.1.8 AND exactly version 2.4 of freeglut
#  GLUT		  = -DFREEGLUT24
  CXXFLAGS        = $(OPTIM)  -I/usr/X11R6/include $(GLUT)

# libraries to link with:
  GLDLIBS  =  -L/usr/X11R6/lib  -lGLU -lGL -lXinerama -lXft -lpthread -lm  -lXext -lX11  -lXpm \
		-lfltk \
# -lglut # only needed with old versions of fltk before 1.1.8
#       -lXi -lgdkglext-x11-1.0a  # Fedora Core 2 needs these libraries
  LINKFLTKGL  = $(libdir)libfltk.a $(libdir)libfltk_gl.a

# The extension to use for executables...
  EXEEXT          =

endif

ifeq ($(OS),CYGWIN)

# section to define the various symbols needed to build DRAWxtl for cygwin running on Windows

  prefix          = /usr/
  includedir      = 
  libdir          = /fltk/lib

# compiler name:
  CXX             = g++

# flags for C++ compiler:
  OPTIM           = -g -O2 -Wall -Wunused -pedantic -fno-exceptions -I /usr/i686-pc-cygwin/include/mingw
  CXXFLAGS        =  $(OPTIM) -I /fltk -I /FLU -mno-cygwin -mwindows -DWIN32
  PREINC 	  = -I/fltk -I/FLU

# libraries to link with:
  GLDLIBS  =  -mwindows -mno-cygwin -lglut32 -lglu32 -lopengl32 -lole32 -luuid -lcomctl32 -lwsock32 -lsupc++
  LINKFLTKGL  = -L/FLU/lib -lflu  -L/fltk/lib -lfltk_gl -lfltk

# The extension to use for executables...
  EXEEXT          = ".EXE"

  # Resource file (the icon)
  RESFILE	  = Drawv54.rc

endif

ifeq ($(OS),IRIX64)               # use GNU g++ compiler on SGI

# this uses the GNU g++ compiler
# libglut was loaded from the SGI freeware site: freeware.sgi.com
# but FLTK was built from source, since g++ and CC (MIPS c++) may be incompatible
# source patch needed to FLTK: change fltk-1.1.5/FL/Fl_Double_Window.H
# string  type(FL_DOUBLE_WINDOW) to type(0xF1) (two places & no, I don't understand)

  prefix          = /usr/user2/toby/fltk-1.1.5
  includedir      = ${prefix}/
  libdir          = ${prefix}/lib

  # compiler name:
  CXX             = g++

  # flags for C++ compiler:
  OPTIM           = -g -O2 -Wall -Wunused  -fno-exceptions
  CXXFLAGS        = $(OPTIM) 

  # libraries to link with:
  GLDLIBS  =  -L/usr/freeware/lib32 -lglut -L/usr/lib32/ -lGL -lGLU -lm -lXpm -lX11 
  LINKFLTKGL  = -L$(libdir) -lfltk_gl -lfltk

  # The extension to use for executables...
  EXEEXT          =

endif

ifeq ($(OS),Darwin)

# Mac stuff goes here
  #prefix          = /usr/local
  includedir      = ../../../fltk
  libdir          = ../../../fltk/lib

# compiler name:
  CXX             = g++

# flags for C++ compiler:
  OPTIM           = -g -O2 -Wall -Wunused  -fno-exceptions
  CXXFLAGS        = $(OPTIM)  
# in 10.2 I think I needed to add this to the above: 
# -framework AGL -framework GLUT -framework OpenGL -framework Carbon -framework ApplicationServices

# libraries to link with:
#  GLDLIBS  =  -framework AGL -framework GLUT -framework OpenGL  -framework Carbon -framework ApplicationServices -lsupc++ -lobjc
# libsupc++ not used with gcc v4.0
  GLDLIBS  =  -framework AGL -framework GLUT -framework OpenGL  -framework Carbon -framework ApplicationServices -lobjc
  LINKFLTKGL  = -L$(libdir) -static -lfltk_gl -lfltk
  PREINC 	  = -I$(includedir)
# The extension to use for executables...
  EXEEXT          =

# post processing command
  POSTLINK = /Developer/Tools/Rez -t APPL $(includedir)/FL/mac.r -o

endif

##################################################################################
#
# ALL Operating Specific material MUST be specified before this point
#
##################################################################################

ifeq ($(OPTIM),error)
  OPTIM=-O
  CXX=cc
  CXXFLAGS=$(OPTIM)
  GLDLIBS=-lfltk_gl -lfltk -lglut -lGLU -lGL
endif

# Be quiet when building...
.SILENT:

# Build commands and filename extensions...
.SUFFIXES:	.rc .cxx .h .o $(EXEEXT)

.rc.o:
	echo Building resource $< ...
	windres Drawv54.rc Drawv54.o

.cxx.o:
	echo Compiling $< ...
	$(CXX) -I$(includedir) $(CXXFLAGS) -c $<

CPPFILES =\
	CrystalView.cxx \
	cursor.cxx \
	drawmap.cxx \
	DRAWxtl1.cxx \
	DRAWxtl2.cxx \
	DRAWxtl3.cxx \
	DRAWxtlViewUI.cxx \
	EditView.cxx \
	Edit1.cxx \
	Edit2.cxx \
	Ellipsoids.cxx \
	Flu_Combo_Box.cxx \
	Flu_Combo_List.cxx \
	Flu_Spinner.cxx \
	Fourier.cxx \
	Generate_Drawing.cxx \
	gl2gif.cxx \
	Help.cxx \
	Import.cxx \
	MC.cxx \
	mpVector.cxx \
	Read_Input.cxx \
	symmtry.cxx \
	Tb_Window.cxx

OBJFILES = $(CPPFILES:.cxx=.o) $(RESFILE:.rc=.o)

all:	$(TARGET) 

clean:
	-@ rm -f $(TARGET) *.o core *~ .deps ../../examples/*.{out,cns,pov,wrl,tmp,tga,lst,fl,bmp} \
		*.orig *.rej

# Define the dependencies of the object files on the header files

.deps:  $(CPPFILES)
	echo "Building dependencies..."
	echo >.deps
	for i in $(CPPFILES); do gcc -MM $(PREINC) $$i >>.deps; done

-include .deps

# DRAWxtl 

$(TARGET) : .deps $(OBJFILES) 
	echo Linking $@...
	mkdir -p ../../exe
	$(CXX) $(CXXFLAGS) -o $@ $(OBJFILES) \
		$(LINKFLTKGL) $(GLDLIBS) 
	$(POSTLINK) $@

static: .deps $(OBJFILES) 
	echo Linking $@...
	mkdir -p ../../exe
	$(CXX) $(CXXFLAGS) -static -o DRAWxtl.$@ $(OBJFILES) \
		$(LINKFLTKGL) $(GLDLIBS) -lXxf86vm -ldl
	strip $@

