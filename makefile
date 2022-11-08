.PHONY: all build check document test

all: document build check

build: doc
	R CMD build .

#check: build
#	R CMD check spatialcluster*tar.gz

clean:
	-rm -f spatialcluster*tar.gz
	-rm -fr spatialcluster.Rcheck
	-rm -fr src/*.{o,so}

doc: clean
	Rscript -e 'devtools::document()'
	Rscript -e 'rmarkdown::render("README.Rmd")'

test:
	Rscript -e 'devtools::test()'

check:
	Rscript -e 'library(pkgcheck); checks <- pkgcheck(); print(checks); summary (checks)'

install: clean
	R CMD INSTALL .
