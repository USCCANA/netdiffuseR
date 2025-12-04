help:
	@echo "Makefile for netdiffuseR package"
	@echo "Usage:"
	@echo "  make install     - Install the package"
	@echo "  make build       - Build the package source tarball"
	@echo "  make check       - Check the package with R CMD check"
	@echo "  make checkv      - Check the package with R CMD check using Valgrind"
	@echo "  make clean       - Clean up the build artifacts"
	@echo "  make docs        - Generate documentation"

install: 
	Rscript --vanilla -e 'devtools::install()'

build:
	R CMD build . 

README.md: README.qmd
	quarto render README.qmd

check:
	Rscript --vanilla -e 'devtools::check()'

checkv: netdiffuseR_$(VERSION).tar.gz
	R CMD check --as-cran --use-valgrind netdiffuseR_$(VERSION).tar.gz

clean:
	rm -rf netdiffuseR.Rcheck src/*.so src/*.o

docs:
	Rscript --vanilla -e 'devtools::document()'

.PHONY: check checkv clean install docs
