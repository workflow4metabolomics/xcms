# USAGE: make [install|clean]

# -------- VARIABLE --------

OBJ=xcms_set.tgz
DEP=abims_xcms_xcmsSet.xml tool_dependencies.xml repository_dependencies.xml static test-data


# ------------------------

all: $(OBJ)

$(OBJ): $(DEP)
	tar --exclude=".svn" -zchf $@ $^

# ------------------------

install: $(OBJ)
	mv *.tgz ~

clean:
	rm *.tgz

