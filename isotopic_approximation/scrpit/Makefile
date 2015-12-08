SHELL = /bin/sh

DATA_DIR =  /home/wegatron/workspace/geometry/data
SIMP_OUTPUT_DIR = /home/wegatron/tmp/simp_tol

TEST_SIMP_TOL_EXEC = /home/wegatron/workspace/geometry/isotopic_approximation/release/bin/test_simptolerance
TEST_SIMP_TOL_EXEC_DEBUG = /home/wegatron/workspace/geometry/isotopic_approximation/build-debug/bin/test_simptolerance
DEBUG_SIMP_JPTS_EXEC = /home/wegatron/workspace/geometry/isotopic_approximation/release/bin/debug_simp_jpts

SPHERE = ${DATA_DIR}/sphere.stl
SPHERE_OUTPUT_DIR = ${SIMP_OUTPUT_DIR}/sphere/
SPHERE_OUTPUT_PREFIX = ${SIMP_OUTPUT_DIR}/sphere/sphere
SPHERE_THICK = 0.1
SPHERE_SAMPLE_R = 0.03

FANDISK = ${DATA_DIR}/fandisk.obj
FANDISK_OUTPUT_DIR = ${SIMP_OUTPUT_DIR}/fandisk/
FANDISK_OUTPUT_PREFIX = ${SIMP_OUTPUT_DIR}/fandisk/fandisk
FANDISK_THICK = 0.008
FANDISK_SAMPLE_R = 0.002

BUNNY = ${DATA_DIR}/bunny.obj
BUNNY_OUTPUT_DIR = ${SIMP_OUTPUT_DIR}/bunny/
BUNNY_OUTPUT_PREFIX = ${SIMP_OUTPUT_DIR}/bunny/bunny
BUNNY_THICK = 0.007
BUNNY_SAMPLE_R = 0.0035

CUBE = ${DATA_DIR}/cube.stl
CUBE_OUTPUT_DIR = ${SIMP_OUTPUT_DIR}/cube
CUBE_OUTPUT_PREFIX = ${SIMP_OUTPUT_DIR}/cube/cube
CUBE_THICK = 0.06
CUBE_SAMPLE_R = 0.006

CYLINDER = ${DATA_DIR}/cylinder.ply
CYLINDER_OUTPUT_DIR = ${SIMP_OUTPUT_DIR}/cylinder
CYLINDER_OUTPUT_PREFIX = ${SIMP_OUTPUT_DIR}/cylinder/cylinder
CYLINDER_THICK = 0.02
CYLINDER_SAMPLE_R =  0.006

FERTILITY =  ${DATA_DIR}/fertility.stl
FERTILITY_OUTPUT_DIR = ${SIMP_OUTPUT_DIR}/fertility
FERTILITY_OUTPUT_PREFIX = ${SIMP_OUTPUT_DIR}/fertility/fertility
FERTILITY_THICK = 0.5
FERTILITY_SAMPLE_R = 0.12

ARMADILLO = ${DATA_DIR}/Armadillo.obj
ARMADILLO_OUTPUT_DIR = ${SIMP_OUTPUT_DIR}/Armadillo
ARMADILLO_OUTPUT_PREFIX = ${SIMP_OUTPUT_DIR}/Armadillo/Armadillo
ARMADILLO_THICK = 0.6
ARMADILLO_SAMOKE_R = 0.12

.PHONY : default
default:
	echo "test simp tolerance, please specify a case to run!"

.PHONY : debug_simp_jpts
debug_simp_jpts: clean_simp_fandisk
	${DEBUG_SIMP_JPTS_EXEC} ${FANDISK} ${FANDISK_OUTPUT_DIR} ${FANDISK_THICK} ${FANDISK_SAMPLE_R} 35402 35403

.PHONY : test_simp_armadillo
test_simp_armadillo: 
	${TEST_SIMP_TOL_EXEC} ${ARMADILLO} ${ARMADILLO_OUTPUT_PREFIX} ${ARMADILLO_THICK} ${ARMADILLO_SAMOKE_R}
	
.PHONY : test_simp_cylinder
test_simp_cylinder: 
	${TEST_SIMP_TOL_EXEC} ${CYLINDER} ${CYLINDER_OUTPUT_PREFIX} ${CYLINDER_THICK} ${CYLINDER_SAMPLE_R}
	
.PHONY : test_simp_fertility
test_simp_fertility: 
	${TEST_SIMP_TOL_EXEC} ${FERTILITY} ${FERTILITY_OUTPUT_PREFIX} ${FERTILITY_THICK} ${FERTILITY_SAMPLE_R}

.PHONY : test_simp_sphere
test_simp_sphere: clean_simp_sphere
	${TEST_SIMP_TOL_EXEC} ${SPHERE} ${SPHERE_OUTPUT_PREFIX} ${SPHERE_THICK} ${SPHERE_SAMPLE_R}

.PHONY : test_simp_sphere_debug
test_simp_sphere_debug: clean_simp_sphere
	${TEST_SIMP_TOL_EXEC_DEBUG} ${SPHERE} ${SPHERE_OUTPUT_PREFIX} ${SPHERE_THICK} ${SPHERE_SAMPLE_R}

.PHONY : test_simp_fandisk
test_simp_fandisk: clean_simp_fandisk
	${TEST_SIMP_TOL_EXEC} ${FANDISK} ${FANDISK_OUTPUT_PREFIX} ${FANDISK_THICK} ${FANDISK_SAMPLE_R}
	
.PHONY : test_simp_fandisk
test_simp_fandisk_debug: clean_simp_fandisk
	${TEST_SIMP_TOL_EXEC_DEBUG} ${FANDISK} ${FANDISK_OUTPUT_PREFIX} ${FANDISK_THICK} ${FANDISK_SAMPLE_R}
	
.PHONY : test_simp_bunny
test_simp_bunny: clean_simp_bunny
	${TEST_SIMP_TOL_EXEC} ${BUNNY} ${BUNNY_OUTPUT_PREFIX} ${BUNNY_THICK} ${BUNNY_SAMPLE_R}

.PHONY : test_simp_cube
test_simp_cube: clean_simp_cube
	${TEST_SIMP_TOL_EXEC} ${CUBE} ${CUBE_OUTPUT_PREFIX} ${CUBE_THICK} ${CUBE_SAMPLE_R}

.PHONY : clean_simp_sphere
clean_simp_sphere: clean_log
	rm -rf ${SPHERE_OUTPUT_DIR}*; \
	mkdir ${SIMP_OUTPUT_DIR}/sphere/debug

.PHONY : clean_simp_fandisk
clean_simp_fandisk: clean_log
	rm -rf ${FANDISK_OUTPUT_DIR}*; \
	mkdir ${FANDISK_OUTPUT_DIR}jpts; \
	mkdir ${FANDISK_OUTPUT_DIR}adjv; \

.PHONY : clean_simp_bunny
clean_simp_bunny: clean_log
	rm -rf ${BUNNY_OUTPUT_DIR}*

.PHONY : clean_simp_cube
clean_simp_cube: clean_log
	rm -rf ${CUBE_OUTPUT_DIR}

.PHONY : clean_log
clean_log:
	echo > /tmp/zswlog.log

.PHONY : clean_all
clean_all: clean_simp_sphere clean_simp_fandisk clean_simp_bunny clean_simp_cube 
