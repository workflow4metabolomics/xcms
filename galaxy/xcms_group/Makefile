# USAGE: make [install|clean]

# -------- VARIABLE --------

OBJ=xcms_group.tgz
DEP=abims_xcms_group.xml tool_dependencies.xml repository_dependencies.xml static test-data


# ------------------------

all: $(OBJ)

$(OBJ): $(DEP)
	tar --exclude=".svn" -zchf $@ $^

# ------------------------

install: $(OBJ)
	mv *.tgz ~

clean:
	rm *.tgz

