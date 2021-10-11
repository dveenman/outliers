# If you are new to Makefiles: https://makefiletutorial.com

# This root Makefile only calls the Makefile in the `blog` subdirectory
# to prepare the visuals for our outlier blog post

SUBDIRS := blog

all: $(SUBDIRS)


clean:
	$(MAKE) -C blog clean

very-clean:
	$(MAKE) -C blog very-clean

dist-clean:
	$(MAKE) -C blog dist-clean


$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: all clean $(SUBDIRS)


