CC = g++
CFLAGS = -O2 -Wall -std=gnu++14
LDFLAGS =
INCLUDES =

SOURCES = main.cpp Household.cpp McmcObject.cpp
HEADERS = Household.h McmcObject.h
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = program.out

covid: $(EXECUTABLE) $(SOURCES) $(HEADERS)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(INCLUDES) $(OBJECTS) -o $@

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean

clean:
	rm -rf *.o *.out
