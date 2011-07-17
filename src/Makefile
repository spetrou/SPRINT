# Ensure the following line is uncommented for debug logging
#_DEBUG = -D_DEBUG
#_PROF = -D_PROF

# Pull in the system specifics which the configuration script cooked up

include Makevars

MAKEFILE = Makefile


# * ********************************************************************** *
# *  We have two top level targets. The library which will be linked to R  *
# *  and the computation node executable which will form the cluster.      *
# * ********************************************************************** *

SHLIB_OBJS = $(patsubst %.c,%.o,$(wildcard *.c))
SPRINT_OBJS = $(patsubst %.c,%.o,$(wildcard sprint/*.c))


# * **************************************************************** *
# *  Each algorithm provided by the compute cluster is stored in a   *
# *  seperate directory. This allows shared resources, in the local  *
# *  root, and both interface and implementation code in eponimus    *
# *  directories.                                                    *
# * **************************************************************** *

ALGORITHM_DIRS = $(wildcard algorithms/*)

INTERFACE_OBJS = $(patsubst %.c,%.o,$(foreach dir,$(ALGORITHM_DIRS),$(wildcard $(dir)/interface/*.c)))
IMPLEMENTATION_OBJS = $(patsubst %.c,%.o,$(foreach dir,$(ALGORITHM_DIRS),$(wildcard $(dir)/implementation/*.c)))
SHARED_ALGORITHM_OBJS = $(patsubst %.c,%.o,$(foreach dir,$(ALGORITHM_DIRS),$(wildcard $(dir)/*.c)))


# * **************************************************************** *
# *  Compiler flags additional to the ones R automatically provides  *
# * **************************************************************** *

LOCAL_CFLAGS = -pedantic -O3 -fPIC $(_DEBUG) $(_PROF)

# The universal target

.PHONY: all
all: sprint.so


# * **************************************** *
# *  This is the library which links with R  *
# * **************************************** *

sprint.so: $(INTERFACE_OBJS) $(SHARED_ALGORITHM_OBJS) $(SHLIB_OBJS) $(IMPLEMENTATION_OBJS)
ifeq (mpicc,$(findstring mpicc,$(MPICC)))
	$(MPICC) -shared -o $@ $^ $(PKG_LIBS)
else
	$(CC) -shared -o $@ $^ $(MPILIBS)
endif


# * ****************************************** *
# *  Create object files for each main target  *
# * ****************************************** *

$(SPRINT_OBJS) $(INTERFACE_OBJS) $(IMPLEMENTATION_OBJS) $(SHARED_ALGORITHM_OBJS) $(SHLIB_OBJS): %.o: %.c $(MAKEFILE)
	$(MPICC) $(LOCAL_CFLAGS) $(PKG_CFLAGS) -I$(R_INCLUDE_DIR) -c -o $@ $<

# * ****************************** *
# *  Clean up the build directory  *
# * ****************************** *

.PHONY: clean
clean:
	find . -name "*.o" -exec rm -f {} \;
	rm -f sprint.so
