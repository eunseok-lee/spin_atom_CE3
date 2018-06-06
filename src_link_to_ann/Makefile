# build cluster expansion associated codes.
CC = mpicc

#CFLAGS = -g -w
CFLAGS = -w

#LFLAGS = -L/opt/local/lib
LFLAGS = 
LIBS = -lm 

#INCLUDES = -I/opt/local/include
INCLUDES = 

SRCS = link_to_ann.c feed_forward.c dot_vectors.c sigmoidf.c
 
OBJS = $(SRCS:.c=.o)

TARGET = run_link_to_ann

default: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) $(TARGET) *.o
