SHELL = /bin/sh
TEST_SIMP_TOL_EXEC = /home/wegatron/workspace/geometry/isotopic_approximation/build-release/bin/test_simptolerance

SPHERE = /home/wegatron/workspace/geometry/data/sphere.stl
SPHERE_OUTPUT_DIR = /home/wegatron/tmp/simp_tol/sphere/
SPHERE_OUTPUT_PREFIX = /home/wegatron/tmp/simp_tol/sphere/sphere
SPHERE_THICK = 0.1
SPHERE_SAMPLE_R = 0.01

FANDISK = /home/wegatron/workspace/geometry/data/fandisk.obj
FANDISK_OUTPUT_DIR = /home/wegatron/tmp/simp_tol/fandisk/
FANDISK_OUTPUT_PREFIX = /home/wegatron/tmp/simp_tol/fandisk/fandisk
FANDISK_THICK = 0.008
FANDISK_SAMPLE_R = 0.0008

BUNNY = /home/wegatron/workspace/geometry/data/bunny.obj
BUNNY_OUTPUT_DIR = /home/wegatron/tmp/simp_tol/bunny/
BUNNY_OUTPUT_PREFIX = /home/wegatron/tmp/simp_tol/bunny/bunny
BUNNY_THICK = 0.007
BUNNY_SAMPLE_R = 0.0035

CUBE = /home/wegatron/workspace/geometry/data/cube.stl
CUBE_OUTPUT_DIR = /home/wegatron/tmp/simp_tol/cube
CUBE_OUTPUT_PREFIX = /home/wegatron/tmp/simp_tol/cube/cube
CUBE_THICK = 0.06
CUBE_SAMPLE_R = 0.006

.PHONY : default
default:
	echo "test simp tolerance, please specify a case to run!"

.PHONY : test_simp_sphere
test_simp_sphere:
	${TEST_SIMP_TOL_EXEC} ${SPHERE} ${SPHERE_OUTPUT_PREFIX} ${SPHERE_THICK} ${SPHERE_SAMPLE_R}

.PHONY : test_simp_fandisk
test_simp_fandisk:
	${TEST_SIMP_TOL_EXEC} ${FANDISK} ${FANDISK_OUTPUT_PREFIX} ${FANDISK_THICK} ${FANDISK_SAMPLE_R}
	
.PHONY : test_simp_bunny
test_simp_bunny:
	${TEST_SIMP_TOL_EXEC} ${BUNNY} ${BUNNY_OUTPUT_PREFIX} ${BUNNY_THICK} ${BUNNY_SAMPLE_R}

.PHONY : test_simp_cube
test_simp_cube:
	${TEST_SIMP_TOL_EXEC} ${CUBE} ${CUBE_OUTPUT_PREFIX} ${CUBE_THICK} ${CUBE_SAMPLE_R}

.PHONY : clean_simp_sphere
clean_simp_sphere:
	rm -rf ${SPHERE_OUTPUT_DIR}*

.PHONY : clean_simp_fandisk
clean_simp_fandisk:
	rm -rf ${FANDISK_OUTPUT_DIR}*

.PHONY : clean_simp_bunny
clean_simp_bunny:
	rm -rf ${BUNNY_OUTPUT_DIR}*

.PHONY : clean_simp_cube
clean_simp_cube:
	rm -rf ${CUBE_OUTPUT_DIR}
	
.PHONY : clean_all
clean_all: clean_simp_sphere clean_simp_fandisk clean_simp_bunny clean_simp_cube 