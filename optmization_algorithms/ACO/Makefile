#
# 'make dist' : Package source code in tar.gz. Files included in the package
# 		must be listed explicitly below.
#
# 'make rsync' : Copy the whole source directory to a remote location. Useful
# 		 for copying the code to a computation cluster.
#
# 'make upload' : Package the source code and upload it to the web. It should
# 		  only be used when releasing a new version. If you only want
# 		  to edit the web, use the Makefile in web/ ('cd web/; make').
#
VERSION=0.99
DEBUG=0
TRACE=0

ifneq ($(DEBUG),0)
DEBUGFLAGS=-g3 -DMALLOC_CHECK_=3 -DDODEBUG=$(DEBUG)
else
DEBUGFLAGS=
endif

ifneq ($(TRACE),0)
TRACEFLAGS=-DTRACE=$(TRACE)
endif

OPT_FLAGS=-O2 -ffast-math
WARN_FLAGS=-Wall -Wextra -std=c99 -pedantic
CFLAGS=$(WARN_FLAGS) $(OPT_FLAGS) $(DEBUGFLAGS) $(TRACEFLAGS)
LDLIBS+= $(OPT_FLAGS) -lm
# OS:=$(shell uname)
# ifeq ($(OS),Linux)
# TIMER=unix
# else
TIMER=dos
# endif

.PHONY: default all clean rsync dist upload

default: main

main: main.o acoqap.o aco.o QAP.o utilities.o ants.o qap-ants.o InOut.o $(TIMER)_timer.o qap-ls.o parse.o adaptation.o

acoqap: acoqap.o aco.o QAP.o utilities.o ants.o qap-ants.o InOut.o $(TIMER)_timer.o qap-ls.o parse.o adaptation.o

acotsp: acotsp.o aco.o TSP.o utilities.o ants.o tsp-ants.o InOut.o $(TIMER)_timer.o ls.o parse.o adaptation.o


all: clean acoqap acotsp

clean:
	@$(RM) *.o acotsp acoqap

runner.o: runner.c

aco.o: aco.c

acotsp.o: acotsp.c

acoqap.o: acoqap.c

TSP.o: TSP.c TSP.h

QAP.o: QAP.c QAP.h

ants.o: ants.c ants.h aco-parameters.def

qap-ants.o: qap-ants.c ants.h aco-parameters.def

InOut.o: InOut.c InOut.h

utilities.o: utilities.c utilities.h

ls.o: ls.c ls.h

parse.o: parse.c parse.h aco-parameters.def

$(TIMER)_timer.o: $(TIMER)_timer.c timer.h

adaptation.o : adaptation.h

TAGS: *.c *.h
	etags *.c *.h

rsync : RNODE=iridiacluster
rsync : all
ifdef RNODE
	rsync -rlp -CIzc -L --delete --copy-unsafe-links \
	--exclude=.svn --progress --relative --exclude=acotsp --exclude=acoqap	\
	.                                       \
	$(RNODE):~/src/acoqap/
else
	@echo "ERROR: You must specify a remote node (e.g., RNODE=b21)"
	@exit 1
endif

dist : DIST_SRC_FILES=*.c *.h README Makefile gpl.txt *.def tsp/*.tsp qap/*.dat
dist : DISTFILENAME:=ACOTSPQAP-$(VERSION)
dist : DISTDIR:=../
dist : all
	@(mkdir -p $(DISTDIR)/$(DISTFILENAME)			\
	&& rsync --relative -rlpC --exclude=.svn $(DIST_SRC_FILES) $(DISTDIR)/$(DISTFILENAME)/ \
        && cd  $(DISTDIR) 	\
	&& tar cf - $(DISTFILENAME) | gzip -f9 > $(DISTFILENAME).tar.gz \
	&& rm -rf ./$(DISTFILENAME)					\
	&& echo " $(DISTDIR)/$(DISTFILENAME).tar.gz created." && cd $(CWD) )

upload:
	$(MAKE) DISTDIR=web dist
	$(MAKE) -C web update
	@echo "****"
	@echo "**** DO NOT forget to update the link to $(DISTFILENAME).tar.gz in web/index.html! (TODO: do this automatically using 'sed') ***"
	@echo "****"
