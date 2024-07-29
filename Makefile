
########################################################################
# Makefile for Linux
########################################################################

CXX = g++

ifneq (,$(wildcard libs))
    LIBDIR := libs
else
    ifneq (,$(wildcard ../libs))
        LIBDIR := ../libs
    else
        ifneq (,$(wildcard ../../libs))
            LIBDIR := ../../libs
        else
            LIBDIR := ../../../libs
        endif
    endif
endif

OPTIMIZE = -g -O4


CXXFLAGS = $(OPTIMIZE) -std=c++17 -I/home/gherron/projects/assimp/include -I. -I$(LIBDIR)/glfw/include -I$(LIBDIR)/bvh-v1  -I$(LIBDIR)/glm -I$(LIBDIR) -I/usr/include   -Wnarrowing -I.  -fopenmp -msse3 

LIBS = -L/home/gherron/projects/assimp/bin -L$(LIBDIR) -L/usr/lib -lassimp -lglbinding -lX11 -lGLU -lGL `pkg-config --static --libs glfw3`


target = raytrace.exe

headers = geom.h raytrace.h  realtime.h rgbe.h acceleration.h
src = geom.cpp main.cpp raytrace.cpp readAssimpFile.cpp realtime.cpp rgbe.cpp acceleration.cpp
extras = raytrace.vcxproj Makefile realtime.vert realtime.frag dwarf.jpg dwarf2.jpg axe.jpg

scenes = testscene.scn letterX.ply letterY.ply bunny.ply dwarf.x

pkgDir = /home/gherron/packages
pkgFiles = $(src) $(headers) $(extras) $(scenes)
pkgName = CS500-framework

objects = $(patsubst %.cpp,%.o,$(src))


$(target): $(objects)
	g++  $(CXXFLAGS) -o $@  $(objects) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(target)
	LD_LIBRARY_PATH="$(LIBDIR); $(LD_LIBRARY_PATH)" ./raytrace.exe testscene.scn

zip:
	rm -rf $(pkgDir)/$(pkgName) $(pkgDir)/$(pkgName).zip
	mkdir $(pkgDir)/$(pkgName)
	cp $(pkgFiles) $(pkgDir)/$(pkgName)
	cp -r ../libs $(pkgDir)/$(pkgName)
	rm -rf $(pkgDir)/$(pkgName)/libs/.hg 
	cd $(pkgDir);  zip -r $(pkgName).zip $(pkgName)


clean:
	rm -rf *.suo *.sdf *.orig Release Debug ipch *.o *~ raytrace dependencies *13*scn  *13*ppm 

dependencies: 
	g++ -MM $(CXXFLAGS)  $(src) > dependencies

include dependencies
