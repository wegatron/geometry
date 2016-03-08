SHELL = /bin/bash

include pc_setting.in

.PHONY: export
export:
	export LD_LIBRARY_PATH=${EXPORT_LIB_DIR}
.PHONY : clean_log
clean_log:
	echo > /tmp/zswlog.log

include data.in
include suite.in
