VERSION:=$(shell Rscript -e 'x<-readLines("DESCRIPTION");cat(gsub(".+[:]\\\s*", "", x[grepl("^Vers", x)]))')

help:
	@echo "Makefile for netdiffuseR package"
	@echo "Usage:"
	@echo "  make install     - Install the package"
	@echo "  make check       - Check the package with R CMD check"
	@echo "  make checkv      - Check the package with R CMD check using Valgrind"
	@echo "  make clean       - Clean up the build artifacts"
	@echo "  make docs        - Generate documentation"

install: netdiffuseR_$(VERSION).tar.gz
	R CMD INSTALL netdiffuseR_$(VERSION).tar.gz

netdiffuseR_$(VERSION).tar.gz: */*.R 
	R CMD build . 

README.md: README.Rmd
	Rscript -e 'rmarkdown::render("README.Rmd")'

check: netdiffuseR_$(VERSION).tar.gz
	R CMD check --as-cran netdiffuseR_$(VERSION).tar.gz

checkv: netdiffuseR_$(VERSION).tar.gz
	R CMD check --as-cran --use-valgrind netdiffuseR_$(VERSION).tar.gz

clean:
	rm -rf netdiffuseR.Rcheck

docs:
	Rscript --vanilla -e 'devtools::document()'

.PHONY: check checkv clean install docs