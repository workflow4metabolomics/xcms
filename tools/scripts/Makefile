# USAGE: make [install|clean]

# -------- VARIABLE --------

OBJ=camera_script-
MAIN_SCRIPT=../camera_w4m_script/CAMERA.r
DEP=$(MAIN_SCRIPT) ../camera_w4m_script/lib.r

HOST=brazil
HOST_PATH=/sbr2/download/w4m/


# ------------------------

all: $(DEP)
	chmod +x $(MAIN_SCRIPT)
	tar --exclude=".svn" -zchf $(OBJ)$$(grep version $(MAIN_SCRIPT) | sed -r 's/^.*=.(.*).$$/\1/').tgz $^

# ------------------------

install:
	@ssh $(HOST) stat $(HOST_PATH)/$(OBJ)$$(grep version $(MAIN_SCRIPT) | sed -r 's/^.*=.(.*).$$/\1/').tgz 2> /dev/null; if [ $$? -eq 1 ]; then scp *.tgz $(HOST):$(HOST_PATH); else echo "The file in this version already exist. Please, change the version"; fi

installforce:
	scp *.tgz $(HOST):$(HOST_PATH)
    
clean:
	rm *.tgz

