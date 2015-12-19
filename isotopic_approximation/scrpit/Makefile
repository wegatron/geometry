SHELL = /bin/bash

DATA_DIR =  /home/wegatron/workspace/geometry/data
APP_OUTPUT_DIR = /home/wegatron/tmp/approximate
APP_EXEC = /home/wegatron/workspace/geometry/isotopic_approximation/build-my-debug-o3/bin/approximate
EXPORT_LIB_DIR = /home/wegatron/workspace/geometry/isotopic_approximation/build-my-debug-o3/lib/
NORMAL_COND_SCALE = 0.7
FLAT_THRESHOLD = 0.25

SPHERE = ${DATA_DIR}/sphere.stl
SPHERE_OUTPUT_DIR = ${APP_OUTPUT_DIR}/sphere/
SPHERE_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/sphere/sphere
SPHERE_THICK = 0.1
SPHERE_SAMPLE_R = 0.03

FANDISK = ${DATA_DIR}/fandisk_smoothed.ply
FANDISK_OUTPUT_DIR = ${APP_OUTPUT_DIR}/fandisk/
FANDISK_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/fandisk/fandisk
FANDISK_THICK = 0.004
FANDISK_SAMPLE_R = 0.002

BUNNY = ${DATA_DIR}/bunny.obj
BUNNY_OUTPUT_DIR = ${APP_OUTPUT_DIR}/bunny/
BUNNY_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/bunny/bunny
BUNNY_THICK = 0.002
BUNNY_SAMPLE_R = 0.0008

CUBE = ${DATA_DIR}/cube.stl
CUBE_OUTPUT_DIR = ${APP_OUTPUT_DIR}/cube
CUBE_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/cube/cube
CUBE_THICK = 0.06
CUBE_SAMPLE_R = 0.03

CYLINDER = ${DATA_DIR}/cylinder_smoothed.ply
CYLINDER_OUTPUT_DIR = ${APP_OUTPUT_DIR}/cylinder/
CYLINDER_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/cylinder/cylinder
CYLINDER_THICK = 0.02
CYLINDER_SAMPLE_R =  0.01

FERTILITY =  ${DATA_DIR}/fertility.stl
FERTILITY_OUTPUT_DIR = ${APP_OUTPUT_DIR}/fertility/
FERTILITY_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/fertility/fertility
FERTILITY_THICK = 0.5
FERTILITY_SAMPLE_R = 0.2

ARMADILLO = ${DATA_DIR}/Armadillo.obj
ARMADILLO_OUTPUT_DIR = ${APP_OUTPUT_DIR}/Armadillo
ARMADILLO_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/Armadillo/Armadillo
ARMADILLO_THICK = 0.4
ARMADILLO_SAMOKE_R = 0.16

SPHERE = ${DATA_DIR}/sphere.stl
SPHERE_OUTPUT_DIR = ${APP_OUTPUT_DIR}/sphere/
SPHERE_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/sphere/sphere
SPHERE_THICK = 0.1
SPHERE_SAMPLE_R = 0.03

FANDISK = ${DATA_DIR}/fandisk_smoothed.ply
FANDISK_OUTPUT_DIR = ${APP_OUTPUT_DIR}/fandisk/
FANDISK_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/fandisk/fandisk
FANDISK_THICK = 0.004
FANDISK_SAMPLE_R = 0.002

BUNNY = ${DATA_DIR}/bunny.obj
BUNNY_OUTPUT_DIR = ${APP_OUTPUT_DIR}/bunny/
BUNNY_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/bunny/bunny
BUNNY_THICK = 0.002
BUNNY_SAMPLE_R = 0.0008

CUBE = ${DATA_DIR}/cube.stl
CUBE_OUTPUT_DIR = ${APP_OUTPUT_DIR}/cube
CUBE_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/cube/cube
CUBE_THICK = 0.06
CUBE_SAMPLE_R = 0.03

CYLINDER = ${DATA_DIR}/cylinder_smoothed.ply
CYLINDER_OUTPUT_DIR = ${APP_OUTPUT_DIR}/cylinder/
CYLINDER_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/cylinder/cylinder
CYLINDER_THICK = 0.02
CYLINDER_SAMPLE_R =  0.01

FERTILITY =  ${DATA_DIR}/fertility.stl
FERTILITY_OUTPUT_DIR = ${APP_OUTPUT_DIR}/fertility/
FERTILITY_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/fertility/fertility
FERTILITY_THICK = 0.5
FERTILITY_SAMPLE_R = 0.2

ARMADILLO = ${DATA_DIR}/Armadillo.obj
ARMADILLO_OUTPUT_DIR = ${APP_OUTPUT_DIR}/Armadillo
ARMADILLO_OUTPUT_PREFIX = ${APP_OUTPUT_DIR}/Armadillo/Armadillo
ARMADILLO_THICK = 0.4
ARMADILLO_SAMOKE_R = 0.16

.PHONY: export
export:
	export LD_LIBRARY_PATH=${EXPORT_LIB_DIR}

.PHONY : clean_log
clean_log:
	echo > /tmp/zswlog.log

.PHONY : clean_sphere
clean_sphere: clean_log
	rm -rf ${SPHERE_OUTPUT_DIR}*; mkdir ${SPHERE_OUTPUT_PREFIX}_tmp
.PHONY : simp_sphere
simp_sphere: clean_sphere export
	${APP_EXEC} ${SPHERE} ${SPHERE_OUTPUT_PREFIX} ${NORMAL_COND_SCALE} ${SPHERE_THICK} ${SPHERE_SAMPLE_R} ${FLAT_THRESHOLD}

.PHONY : clean_cylinder
clean_cylinder: clean_log
	rm -rf ${CYLINDER_OUTPUT_DIR}*; mkdir ${CYLINDER_OUTPUT_PREFIX}_tmp
.PHONY : simp_cylinder
simp_cylinder: clean_cylinder export
	${APP_EXEC} ${CYLINDER} ${CYLINDER_OUTPUT_PREFIX} ${NORMAL_COND_SCALE} ${CYLINDER_THICK} ${CYLINDER_SAMPLE_R} ${FLAT_THRESHOLD}

.PHONY : clean_fertility
clean_fertility: clean_log
	rm -rf ${FERTILITY_OUTPUT_DIR}*; mkdir ${FERTILITY_OUTPUT_PREFIX}_tmp
.PHONY :simp_fertility
simp_fertility: clean_fertility  export
	${APP_EXEC} ${FERTILITY} ${FERTILITY_OUTPUT_PREFIX} ${NORMAL_COND_SCALE} ${FERTILITY_THICK} ${FERTILITY_SAMPLE_R} ${FLAT_THRESHOLD}

.PHONY :clean_bunny
clean_bunny:clean_log
	rm -rf ${BUNNY_OUTPUT_DIR}*; mkdir ${BUNNY_OUTPUT_PREFIX}_tmp
.PHONY : simp_bunny
simp_bunny:clean_bunny  export
	${APP_EXEC} ${BUNNY} ${BUNNY_OUTPUT_PREFIX} ${NORMAL_COND_SCALE} ${BUNNY_THICK} ${BUNNY_SAMPLE_R} ${FLAT_THRESHOLD}

.PHONY:clean_fandisk
clean_fandisk : clean_log
	rm -rf ${FANDISK_OUTPUT_DIR}*; mkdir ${FANDISK_OUTPUT_PREFIX}_tmp
.PHONY:simp_fandisk
simp_fandisk:clean_fandisk export
	${APP_EXEC} ${FANDISK} ${FANDISK_OUTPUT_PREFIX} ${NORMAL_COND_SCALE} ${FANDISK_THICK} ${FANDISK_SAMPLE_R} ${FLAT_THRESHOLD}
