all:nbody

CFLAGS = -O3 -lpthread -std=c++11

PLATFORM=OS_UNKNOWN
ifeq ($(shell uname), Linux)
	PLATFORM=linux
endif
ifeq ($(shell uname), linux)
	PLATFORM=linux
endif
ifeq ($(shell uname), Darwin)
	PLATFORM=darwin
endif

ifeq "$(OSTYPE)" "linux"
	PLATFORM=linux
endif
ifeq "$(OSTYPE)" "Linux"
	PLATFORM=linux
endif
ifeq "$(OSTYPE)" "linux-gnu"
	PLATFORM=linux
endif
ifeq "$(OSTYPE)" "darwin"
	PLATFORM=darwin
endif

ifeq "$(OS)" "Windows_NT"
	PLATFORM=win32
endif

# Linux
ifeq "$(PLATFORM)" "linux"
	LDFLAGS = -lGL -lGLU -lglut
endif

# OSX
ifeq "$(PLATFORM)" "darwin"
LDFLAGS = -framework OpenGL -framework GLUT -Wno-deprecated-declarations
endif

# Windows
ifeq "$(PLATFORM)" "win32"
	LDFLAGS = -lopengl32 -lglu32 -lglut32
endif

nbody: nbody.cpp
	g++ -o $@ $< $(CFLAGS) $(LDFLAGS)
