-include ../config

CFLAGS += -O2

all: libglutils.a libfb_osmesa.a libfb_glx.a

libglutils.a: trackball.o gl2ps.o utils.o polygonize.o \
		og_font.o og_stroke_mono_roman.o
	ar cr $@ $^

libfb_osmesa.a: fb_osmesa.o
	ar cr $@ $^

libfb_glx.a: fb_glx.o OffscreenContextGLX.o fbo.o
	ar cr $@ $^

gl2ps.o: gl2ps/gl2ps.c
	$(CC) $(CFLAGS) -Igl2ps -c gl2ps/gl2ps.c -o gl2ps.o

utils.o: utils.h

clean:
	rm -f *.o *.a
