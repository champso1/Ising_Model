CXX = g++
CXXFLAGS = -Wall -Iinclude `pkg-config --cflags raylib` `pkg-config --cflags gsl` `root-config --cflags` -MMD
RM = rm -rf

LDFLAGS = -g
LDLIBS = `pkg-config --libs raylib` `pkg-config --libs gsl` `root-config --libs`

SRCDIR = src
BUILDDIR = build
TARGET = bin/isingmodel

SRCEXT = cpp
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
DEPS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.d))


$(TARGET): $(OBJS)
	@echo "Linking ..."
	@echo "$(CXX) $(LDFLAGS) $^ -o $(TARGET) $(LDLIBS)"; $(CXX) $(LDFLAGS) $^ -o $(TARGET) $(LDLIBS)


-include $(DEPS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo "$(CXX) $(CXXFLAGS) -c $< -o $@"; $(CXX) $(CXXFLAGS) -MMD -c $< -o $@
	


clean:
	@echo "Cleaning ..."
	@echo "$(RM) $(BUILDDIR) $(TARGET)"; $(RM) $(BUILDDIR) $(TARGET)

.PHONY: clean
