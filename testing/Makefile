DIRS = test_long_readthru_intron_removal \
       test_ctat_build \
       test_index_pfam

test:
	@for i in $(DIRS); do \
	echo "Running test in $$i..."; \
	(cd $$i; $(MAKE) test) || exit $$?; done

clean:
	@for i in $(DIRS); do \
	echo "Running test in $$i..."; \
	(cd $$i; $(MAKE) clean) || exit $$?; done

