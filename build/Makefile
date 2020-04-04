SOURCES		= 	common.c \
				Smirnov.c \
				Davidson.c \
				SEM.c \
				Spectral.c
CC			=	gcc
CFLAGS		=	-Wall -O2 -fPIC -lm
LDFLAGS		=	-shared	
OBJECTS		=	$(SOURCES:.c=.o)
HEADERS		=	$(SOURCES:.c=.h)
TARGET		=	libstg.so


all: $(OBJECTS) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS)


$(OBJECTS): precompiled.h

%.o: %.c %.h
	$(CC) -c $(CFLAGS) $<


clean: 
	rm -f $(OBJECTS)
	rm -f $(TARGET)