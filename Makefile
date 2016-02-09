host=${DSBUILD}
ifeq ($(host),)
host=$(shell hostname |  cut -d'.' -f 1)
endif
include Makefile.$(host)
include Makefile.main