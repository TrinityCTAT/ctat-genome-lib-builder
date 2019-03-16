

ifeq ($(CTAT_GENOME_LIB), "")
	echo "You must set env var CTAT_GENOME_LIB to your CTAT genome lib dir before running tests."
	echo "Note, no building of software is required for general use according to documentation."
	exit 1
endif

test:
	cd testing && $(MAKE) test

clean:
	cd testing && $(MAKE) clean


