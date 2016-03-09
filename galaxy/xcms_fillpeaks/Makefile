# USAGE: make [install|clean]

# -------- VARIABLE --------

OBJ=xcms_fillpeaks.tgz
DEP=abims_xcms_fillPeaks.xml tool_dependencies.xml repository_dependencies.xml static test-data


# ------------------------

all: $(OBJ)

$(OBJ): $(DEP)
	tar --exclude=".svn" -zchf $@ $^

# ------------------------

install: $(OBJ)
	mv *.tgz ~

clean:
	rm *.tgz

