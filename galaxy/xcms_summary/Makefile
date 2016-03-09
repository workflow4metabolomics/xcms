# USAGE: make [install|clean]

# -------- VARIABLE --------

OBJ=xcms_summary.tgz
DEP=abims_xcms_summary.xml xcms_summary.r tool_dependencies.xml repository_dependencies.xml test-data


# ------------------------

all: $(OBJ)

$(OBJ): $(DEP)
	tar --exclude=".svn" -zchf $@ $^

# ------------------------

install: $(OBJ)
	mv *.tgz ~

clean:
	rm *.tgz

