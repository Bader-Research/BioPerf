all: make_pom make_ce make_copy

make_ce:
	$(MAKE) -C ce

make_pom: 
	$(MAKE) -C pom 

make_copy:
	cp ce/ce  ./CE

use_exe_sun:
	cp ce/ce_sun ./CE
	cp pom/mkDB_sun pom/mkDB

use_exe_sgi:
	cp ce/ce_sgi ./CE
	cp pom/mkDB_sgi pom/mkDB
