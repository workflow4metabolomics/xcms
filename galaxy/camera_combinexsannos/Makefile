# USAGE: make [install|clean]

# -------- VARIABLE --------

OBJ=camera_combinexsAnnos.tgz
DEP=abims_CAMERA_combinexsAnnos.xml tool_dependencies.xml repository_dependencies.xml static test-data


# ------------------------

all: $(OBJ)

$(OBJ): $(DEP)
	tar --exclude=".svn" -zchf $@ $^

# ------------------------

install: $(OBJ)
	mv *.tgz ~

clean:
	rm *.tgz

