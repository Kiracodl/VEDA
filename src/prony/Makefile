#compile a matlab script.  modifed from the fermi example on rappture's wiki
#this is probably only going to work on nanohub itself, don't expect it to work
#on your local machine
APP = prony_fit
matlabversion = 2018a # recommended per Steven Clark
DELETE_FILES =
BIN = ../../bin
RP_DIR = $(shell . /etc/environ.sh; use -e -r rappture; echo $${RAPPTURE_PATH})

all:
	if [ -e  /etc/environ.sh ] ; \
	then \
	  . /etc/environ.sh; \
	  use -e -r matlab-$(matlabversion); \
	  make $(APP) ; \
	fi  #otherwise not running on nanohub so do nothing

$(APP): $(APP).m
	mcc -v -C -m -R -nojvm -I $(RP_DIR)/lib/matlab $(APP).m 
	patchrunscript.sh run_$(APP).sh run_$(APP).sh 1

install: all
	if [ -e  /etc/environ.sh ] ; \
	then \
	  install --mode 0755 -D run_$(APP).sh $(BIN) ;\
	  install --mode 0755 -D $(APP) $(BIN) ;\
	  install --mode 0644 -D readme.txt $(BIN)/$(APP).readme ;\
	  if [ -e $(APP).ctf ] ;then mv $(APP).ctf $(BIN); fi ;\
	fi

clean:
	rm -f $(APP) $(APP).ctf \
              $(APP).prj run_$(APP).sh \
              $(APP)_main.c $(APP)_mcc_component_data.c \
              mccExcludedFiles.log readme.txt \
              requiredMCRProducts.txt $(DELETE_FILES)

distclean: clean
	rm -rf $(BIN)/$(APP) \
               $(BIN)/$(APP).ctf \
               $(BIN)/$(APP).readme \
               $(BIN)/$(APP)_mcr \
               $(BIN)/run_$(APP).sh

