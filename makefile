EXECUTABLE=lab2-sequential
CC=g++
CXXFLAGS=-std=c++98 -Wall -Wextra -pedantic -Werror
RELEASEDIR=release
HFILES=$(wildcard *.h)
CPPFILES=$(wildcard *.cpp)
OBASENAMES=$(CPPFILES:.cpp=.o)
RELEASEOFILES=$(addprefix $(RELEASEDIR)/,$(OBASENAMES))

.PHONY: all
all: release

.PHONY: release
release: create_dir_release
release: $(RELEASEDIR)/$(EXECUTABLE)

.PHONY: create_dir_release
create_dir_release:
	test -d $(RELEASEDIR) || mkdir $(RELEASEDIR)

$(RELEASEDIR)/$(EXECUTABLE): $(RELEASEOFILES)
	$(CC) -o $@ $(RELEASEDIR)/*.o

$(RELEASEDIR)/%.o: %.cpp $(HFILES)
	$(CC) -c -o $@ $(CXXFLAGS) $<

.PHONY: clean
clean:
	rm -rf $(RELEASEDIR)/
