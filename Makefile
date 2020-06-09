

####
#### determine what platform: linux, osx or windows
####
ifeq ($(OS),Windows_NT)
	detected_OS := Windows
else
	detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

ifeq ($(detected_OS),Windows)
		PLATFORM = WINDOWS
endif
ifeq ($(detected_OS),Darwin)
		PLATFORM = OSX
endif
ifeq ($(detected_OS),Linux)
		PLATFORM = LINUX
endif


# create directory command
MKDIR = mkdir -p


# Set build directory
BUILDDIR = build


# Set source directory and source files extension
SRCDIR = src
SRCEXT = cpp


# Set executable name and directory
TARGETDIR = bin
TARGET = worm


# set library directory
LIBDIR = libs


# create compiler flags
CXX       = g++ -std=c++1z -Wall --pedantic
INCFLAGS  = -I $(LIBDIR)/json/include -I /usr/local/include
LIBFLAGS  = -L /usr/local/lib
DEF       = -D PLATFORM=$(PLATFORM)

####
#### options
####

####
#### debug mode
####
ifeq ($(DEBUG), 1)
	CXX += -g3 -O0
else
	CXX += -g0 -O3
endif





####
#### ------------- HERE FOLLOW TARGETS -------------
####

# Explore source tree
SRCS := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')


# Create build tree structure: src/../*.cpp -> build/src../*.o
OBJS = $(SRCS:%.cpp=$(BUILDDIR)/%.o)


# Define dependencies files for all objects
DEPS = $(OBJS:.o=.d)


# Indicate to make which targets are not files
.PHONY: all clean setup platform


# The default target triggered by the "make" command
# First the setup rule is run, then the build rule
all: setup $(TARGETDIR)/$(TARGET)


# The target creating all necessary directories for the build process
# Triggered by "make setup"
setup:
	@$(MKDIR) $(dir $(OBJS))
	@$(MKDIR) $(TARGETDIR)


# The build process target, which is a file
$(TARGETDIR)/$(TARGET): $(OBJS)
	@echo linking: $@
	@$(CXX) $(LIBFLAGS) -o $@ $(OBJS)


# Defines the function which will generate the rule for each object (.o file)
define generateRules =
$(1): $(subst .o,.cpp,$(subst $(BUILDDIR),./,$(1)))
	@echo compiling: $$<
	@$$(CXX) $(INCFLAGS) $$(DEF) -c -o $$@ $$< -MMD
endef


# Generates the rules for each .o file
$(foreach obj,$(OBJS),$(eval $(call generateRules,$(obj))))


# Removes build tree and the binary file
# Triggered by "make clean"
clean:
	@$(RM) -rf $(TARGETDIR)/$(TARGET)
	@$(RM) -fr $(BUILDDIR)

# rebuild rule
rebuild: clean all

run: all
	@echo ------------[$(TARGETDIR)/$(TARGET)]------------
	@$(TARGETDIR)/$(TARGET)

# outputs what platform
platform:
	@echo $(PLATFORM)

# Include dependencies
-include $(DEPS)