AUTOMAKE_OPTIONS = foreign
SUBDIRS = lib bin test script

scan:
	script/scan.py lib > lib/Makefile.am
	script/scan.py bin > bin/Makefile.am
	echo > test/Makefile.am

scan-test:
	script/scan.py lib-test > lib/Makefile.am
	script/scan.py bin > bin/Makefile.am
	script/scan.py test > test/Makefile.am

scan-release:
	script/scan.py lib > lib/Makefile.am
	script/scan.py release > bin/Makefile.am
	echo > test/Makefile.am

unittest:
	script/run-unittest.py
