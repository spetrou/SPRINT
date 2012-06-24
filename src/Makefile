# Ensure the following line is uncommented for debug logging
#_DEBUG = -D_DEBUG
#_PROF = -D_PROF

# Pull in the system specifics which the configuration script cooked up

include Makevars

MAKEFILE = Makefile

# * **************************************************************** *
# *  Each algorithm provided by the compute cluster is stored in a   *
# *  seperate directory. This allows shared resources, in the local  *
# *  root, and both interface and implementation code in eponimus    *
# *  directories.                                                    *
# * **************************************************************** *

SHLIB_OBJS = sprint.o
ALGORITHM_DIRS = algorithms/common algorithms/papply algorithms/pboot algorithms/pcor algorithms/phamming algorithms/pmaxT algorithms/ppam algorithms/prandomForest algorithms/pRP algorithms/ptest algorithms/rng
INTERFACE_OBJS = algorithms/papply/interface/papply.o algorithms/pboot/interface/pboot.o algorithms/pcor/interface/pcor.o algorithms/phamming/interface/phamming.o algorithms/pmaxT/interface/pmaxT.o algorithms/ppam/interface/ppam.o algorithms/prandomForest/interface/random_forest.o algorithms/pRP/interface/boot-rp.o algorithms/pRP/interface/rp.o algorithms/ptest/interface/ptest.o algorithms/rng/interface/sprint-rng.o algorithms/mpi/interface/mpi_calls.o algorithms/openmp/interface/openmp.o
IMPLEMENTATION_OBJS = algorithms/papply/implementation/apply.o algorithms/papply/implementation/comms.o algorithms/papply/implementation/ffapply.o algorithms/papply/implementation/lapply.o algorithms/papply/implementation/main.o algorithms/pboot/implementation/boot.o algorithms/pcor/implementation/correlation.o algorithms/pcor/implementation/kernel.o algorithms/phamming/implementation/hamming.o algorithms/phamming/implementation/hamming_kernel.o algorithms/pmaxT/implementation/block_sampling_fixed.o algorithms/pmaxT/implementation/bootloop.o algorithms/pmaxT/implementation/mt.o algorithms/pmaxT/implementation/pairt_sampling.o algorithms/pmaxT/implementation/pairt_sampling_fixed.o algorithms/pmaxT/implementation/permutation.o algorithms/pmaxT/implementation/random.o algorithms/pmaxT/implementation/Rpack.o algorithms/pmaxT/implementation/sampling.o algorithms/pmaxT/implementation/sampling_fixed.o algorithms/pmaxT/implementation/stat_func.o algorithms/pmaxT/implementation/stat_order.o algorithms/ppam/implementation/bswap.o algorithms/ppam/implementation/cstat.o algorithms/ppam/implementation/pamedoids.o algorithms/ppam/implementation/silhouette.o algorithms/prandomForest/implementation/random_forest_driver.o algorithms/pRP/implementation/boot-rank-product.o algorithms/pRP/implementation/rank-product.o algorithms/pRP/implementation/rp-generic.o algorithms/ptest/implementation/test.o algorithms/rng/implementation/sprint-rng-worker.o
SHARED_ALGORITHM_OBJS = algorithms/common/functions.o algorithms/common/mmap.o algorithms/common/serialize.o algorithms/common/utils.o algorithms/openmp/implementation/setThreads.o

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
