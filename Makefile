
CC = g++

CPPFLAGS = -g -O3 -Wall -std=c++11

LIBS = -lpthread -ltbb

SRCS = read.cpp hash.cpp compress.cpp fqreader.cpp decompress.cpp main.cpp

OBJS = $(SRCS: .cpp = .o)

EXEC = Joseph

$(EXEC) : $(OBJS)
	$(CC) $(CPPFLAGS) $(LIBS) $^ -o $@

%.o : %.cpp
	$(CC) -c $(CPPFLAGS) $<

clean:
	rm -f *.o $(EXEC) *.out
