platform=$(shell uname)

ifeq ($(platform),Darwin)
  ALL=bin/gl_server
  GL_OPTS=-framework OpenGL -framework GLUT -Wno-deprecated-declarations -ljpeg -L/opt/local/lib
  INCLUDE=-I/opt/local/include
else ifeq ($(platform),Linux)
  ALL=bin/gl_server
  GL_OPTS=-lGL -lglut -lGLU -lm -ljpeg
endif

all: $(ALL)

clean:
	rm -rf bin/*

bin/gl_server: src/gl_server.c src/opc_server.c src/cJSON.c
	mkdir -p bin
	gcc -o $@ $^ $(GL_OPTS) $(INCLUDE)
