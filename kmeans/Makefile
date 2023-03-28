
TASK := kmeans
SRCDIRS := .
SRCEXT := cpp
SOURCES := $(wildcard $(addsuffix /*.${SRCEXT}, ${SRCDIRS}))
OBJECTS := $(SOURCES:%.${SRCEXT}=%.o)
DEPENDENCIES := $(OBJECTS:%.o=%.d)
EXECUTABLE := ${TASK}.exe
DIRNAME := $(notdir ${CURDIR})

# INPUT_DIR needs to be set according to your local file placement
INPUT_DIR := input

CXX ?= clang++
COMPILER = ${CXX}
FLAGS_OFFLOAD_OPENMP ?= -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target -march=sm_70
FLAGS = -std=c++14 -g ${FLAGS_OFFLOAD_OPENMP}
FLAGS_FAST = -O3
FLAGS_DEBUG = -O0 -Wall -Wextra
INCLUDES = $(addprefix -I, ${SRCDIRS})
LDLIBS =

NTHREADS ?= 1
GROUP ?= 2

# set default build target
build: release

# build for debugging
debug: FLAGS += ${FLAGS_DEBUG}
debug: ${EXECUTABLE}

# build for performance
release: FLAGS += ${FLAGS_FAST}
release: ${EXECUTABLE}

${EXECUTABLE}: ${OBJECTS}
	${COMPILER} ${FLAGS} -o $@ $^ ${LDLIBS}

%.o: %.${SRCEXT}
	${COMPILER} ${INCLUDES} -MMD -MP ${FLAGS} -c -o $@ $<

run-small: NPOINTS = 1000
run-small: NCENTERS ?= 5
run-small: NITERS ?= 20
run-small: release
	OMP_NUM_THREADS=${NTHREADS} ./${EXECUTABLE} ${INPUT_DIR}/small.in ${NPOINTS} ${NCENTERS} ${NITERS}

run-large: NPOINTS = 1000000
run-large: NCENTERS ?= 5000
run-large: NITERS ?= 50
run-large: release
	OMP_NUM_THREADS=${NTHREADS} ./${EXECUTABLE} ${INPUT_DIR}/large.in ${NPOINTS} ${NCENTERS} ${NITERS}

vis-small:
	./utils/make-gif.py ${INPUT_DIR}/small.in memory.out clusters.gif

vis-large:
	./utils/make-gif.py ${INPUT_DIR}/large.in memory.out clusters.gif

archive: clean
	find . -maxdepth 1 -type f -exec tar --transform 's|^|${DIRNAME}-group-${GROUP}/|g' -cvzf ${DIRNAME}-group-${GROUP}.tar.gz {} +

.PHONY: clean build debug release run-small run-large vis-small vis-large archive
clean:
	${RM} ${EXECUTABLE}
	${RM} ${OBJECTS}
	${RM} ${DEPENDENCIES}
	${RM} *.out *.gif

-include ${DEPENDENCIES}
