SHELL = /bin/bash

include pc_setting.in

SCRIPT_DIR = ${APP_ROOT_DIR}/isotopic_approximation2/script
APP_EXEC = ${APP_ROOT_DIR}/isotopic_approximation2/build-release/bin/appro
METRO = ${APP_ROOT_DIR}/3rd_tool/metro

.PHONY: export
export:
	export LD_LIBRARY_PATH=${EXPORT_LIB_DIR}

include ${SCRIPT_DIR}/suite.in