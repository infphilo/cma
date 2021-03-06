#
# Copyright 2018, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of CMA.
#
# CMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Makefile for hisat2-align, hisat2-build, hisat2-inspect
#

INC =
GCC_PREFIX = $(shell dirname `which gcc`)
GCC_SUFFIX =
CC = $(GCC_PREFIX)/gcc$(GCC_SUFFIX)
CPP = $(GCC_PREFIX)/g++$(GCC_SUFFIX)
CXX = $(CPP)
HEADERS = $(wildcard *.h)
BOWTIE_MM = 1
BOWTIE_SHARED_MEM = 0

# Detect Cygwin or MinGW
WINDOWS = 0
CYGWIN = 0
MINGW = 0
ifneq (,$(findstring CYGWIN,$(shell uname)))
	WINDOWS = 1 
	CYGWIN = 1
	# POSIX memory-mapped files not currently supported on Windows
	BOWTIE_MM = 0
	BOWTIE_SHARED_MEM = 0
else
	ifneq (,$(findstring MINGW,$(shell uname)))
		WINDOWS = 1
		MINGW = 1
		# POSIX memory-mapped files not currently supported on Windows
		BOWTIE_MM = 0
		BOWTIE_SHARED_MEM = 0
	endif
endif

MACOS = 0
ifneq (,$(findstring Darwin,$(shell uname)))
	MACOS = 1
endif

EXTRA_FLAGS += -DPOPCNT_CAPABILITY
INC += -I third_party

MM_DEF = 

ifeq (1,$(BOWTIE_MM))
	MM_DEF = -DBOWTIE_MM
endif

SHMEM_DEF = 

ifeq (1,$(BOWTIE_SHARED_MEM))
	SHMEM_DEF = -DBOWTIE_SHARED_MEM
endif

PTHREAD_PKG =
PTHREAD_LIB = 

ifeq (1,$(MINGW))
	PTHREAD_LIB = 
else
	PTHREAD_LIB = -lpthread
endif

SEARCH_LIBS =
BUILD_LIBS = 
INSPECT_LIBS =

ifeq (1,$(MINGW))
	BUILD_LIBS = 
	INSPECT_LIBS = 
endif

SERACH_INC = 

LIBS = $(PTHREAD_LIB) -ltiff

SHARED_CPPS = 
SEARCH_CPPS = 

CMA_CPPS_MAIN = $(SEARCH_CPPS) hisat2_main.cpp
CMA_BUILD_CPPS_MAIN = $(BUILD_CPPS) hisat2_build_main.cpp

SEARCH_FRAGMENTS = $(wildcard search_*_phase*.c)
VERSION = $(shell cat VERSION)

ifeq (1,$(MACOS))
	LIBS += -framework Accelerate
	LIBS += -framework OpenCL
	EXTRA_FLAGS += -I/Systm/Library/Frameworks/Accelerate.framework/Versions/Current/Headers
else
	DYNAMIC_LIBS += -lOpenCL
endif

# Convert BITS=?? to a -m flag
BITS=32
ifeq (x86_64,$(shell uname -m))
BITS=64
endif
# msys will always be 32 bit so look at the cpu arch instead.
ifneq (,$(findstring AMD64,$(PROCESSOR_ARCHITEW6432)))
	ifeq (1,$(MINGW))
		BITS=64
	endif
endif
BITS_FLAG =

ifeq (32,$(BITS))
	BITS_FLAG = -m32
endif

ifeq (64,$(BITS))
	BITS_FLAG = -m64
endif
SSE_FLAG=-msse2

DEBUG_FLAGS    = -O0 -g3 $(BIToS_FLAG) $(SSE_FLAG)
DEBUG_DEFS     = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(EXTRA_FLAGS)\""
RELEASE_FLAGS  = -O3 $(BITS_FLAG) $(SSE_FLAG) -funroll-loops -g3
RELEASE_DEFS   = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
FILE_FLAGS     = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

ifeq (1,$(USE_SRA))
	ifeq (1, $(MACOS))
		DEBUG_FLAGS += -mmacosx-version-min=10.6
		RELEASE_FLAGS += -mmacosx-version-min=10.6
	endif
endif


CMA_BIN_LIST = cma-bin
CMA_BIN_LIST_AUX = cma-bin-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
	$(wildcard scripts/*.pl) \
	$(wildcard *.py) \
	$(wildcard hisatgenotype_modules/*.py) \
	$(wildcard hisatgenotype_scripts/*.py) \
	doc/manual.inc.html \
	doc/README \
	doc/style.css \
	$(wildcard example/index/*.ht2) \
	$(wildcard example/reads/*.fa) \
	example/reference/22_20-21M.fa \
	example/reference/22_20-21M.snp \
	$(PTHREAD_PKG) \
	cma \
	cma-bin \
	AUTHORS \
	LICENSE \
	NEWS \
	MANUAL \
	MANUAL.markdown \
	TUTORIAL \
	VERSION

ifeq (1,$(WINDOWS))
	CMA_BIN_LIST := $(CMA_BIN_LIST) hisat2.bat hisat2-build.bat hisat2-inspect.bat 
endif

# This is helpful on Windows under MinGW/MSYS, where Make might go for
# the Windows FIND tool instead.
FIND=$(shell which find)

SRC_PKG_LIST = $(wildcard *.h) \
	$(wildcard *.hh) \
	$(wildcard *.c) \
	$(wildcard *.cpp) \
	doc/strip_markdown.pl \
	Makefile \
	$(GENERAL_LIST)

BIN_PKG_LIST = $(GENERAL_LIST)

.PHONY: all allall both both-debug

all: $(CMA_BIN_LIST)

allall: $(CMA_BIN_LIST) $(CMA_BIN_LIST_AUX)

both: cma-bin

both-debug: cma-bin-debug

DEFS=-fno-strict-aliasing \
     -DCMA_VERSION="\"`cat VERSION`\"" \
     -DBUILD_HOST="\"`hostname`\"" \
     -DBUILD_TIME="\"`date`\"" \
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
     $(FILE_FLAGS) \
     $(PREF_DEF) \
     $(MM_DEF) \
     $(SHMEM_DEF)

cma-bin: cma.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SEARCH_LIBS)

cma-bin-debug: cma.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SEARCH_LIBS)


cma: ;

cma.bat:
	echo "@echo off" > hisat2.bat
	echo "perl %~dp0/cma %*" >> hisat2.bat

hisat2-build.bat:
	echo "@echo off" > hisat2-build.bat
	echo "python %~dp0/hisat2-build %*" >> hisat2-build.bat

hisat2-inspect.bat:
	echo "@echo off" > hisat2-inspect.bat
	echo "python %~dp0/hisat2-inspect %*" >> hisat2-inspect.bat


.PHONY: cma-src
cma-src-pkg: $(SRC_PKG_LIST)
	chmod a+x scripts/*.sh scripts/*.pl
	mkdir .src.tmp
	mkdir .src.tmp/hisat2-$(VERSION)
	zip tmp.zip $(SRC_PKG_LIST)
	mv tmp.zip .src.tmp/hisat2-$(VERSION)
	cd .src.tmp/hisat2-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .src.tmp ; zip -r hisat2-$(VERSION)-source.zip hisat2-$(VERSION)
	cp .src.tmp/hisat2-$(VERSION)-source.zip .
	rm -rf .src.tmp

.PHONY: cma-bin
cma-bin-pkg: $(BIN_PKG_LIST) $(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX)
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir .bin.tmp/hisat2-$(VERSION)
	if [ -f hisat2.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX) ; \
	fi
	mv tmp.zip .bin.tmp/hisat2-$(VERSION)
	cd .bin.tmp/hisat2-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r hisat2-$(VERSION)-$(BITS).zip hisat2-$(VERSION)
	cp .bin.tmp/hisat2-$(VERSION)-$(BITS).zip .
	rm -rf .bin.tmp

.PHONY: doc
doc: doc/manual.inc.html MANUAL

doc/manual.inc.html: MANUAL.markdown
	pandoc -T "HISAT2 Manual" -o $@ \
	 --from markdown --to HTML --toc $^
	perl -i -ne \
	 '$$w=0 if m|^</body>|;print if $$w;$$w=1 if m|^<body>|;' $@

MANUAL: MANUAL.markdown
	perl doc/strip_markdown.pl < $^ > $@

.PHONY: clean
clean:
	rm -f $(CMA_BIN_LIST) $(CMA_BIN_LIST_AUX) \
	$(addsuffix .exe,$(CMA_BIN_LIST) $(CMA_BIN_LIST_AUX)) \
	cma-src.zip cma-bin.zip
	rm -f core.* .tmp.head
	rm -rf *.dSYM

.PHONY: push-doc
push-doc: doc/manual.inc.html
	scp doc/*.*html doc/indexes.txt salz-dmz:/ccb/salz7-data/www/ccb.jhu.edu/html/software/hisat2/
