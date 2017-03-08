export OPENGL=1
export LIBPNG=1
include ../../src/Makefile.defs
LIB+=-L/usr/X11/lib
OPT+=-I/usr/X11/include
ifeq ($(LIBPNG), 1)
PREDEF+= -DLIBPNG
LIB+= -lpng
endif

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I../../src/ -I./ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c spring.c m_output.c kepcart.c -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../../src/ 
	@-rm -f librebound.so
	@ln -s ../../src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
