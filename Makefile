VERSION:=$(shell Rscript -e 'x<-readLines("DESCRIPTION");cat(gsub(".+[:]\\s*", "", x[grepl("^Vers", x)]))')

install: netdiffuseR_$(VERSION).tar.gz
	R CMD INSTALL netdiffuseR_$(VERSION).tar.gz

netdiffuseR_$(VERSION).tar.gz: */*.R inst/NEWS README.md
	R CMD build . 

inst/NEWS: NEWS.md
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
	head -n 80 inst/NEWS

README.md: README.Rmd
	Rscript -e 'rmarkdown::render("README.Rmd")'

.PHONY: check checkv clean

check: netdiffuseR_$(VERSION).tar.gz
	R CMD check --as-cran netdiffuseR_$(VERSION).tar.gz

checkv: netdiffuseR_$(VERSION).tar.gz
	R CMD check --as-cran --use-valgrind netdiffuseR_$(VERSION).tar.gz

clean:
	rm -rf netdiffuseR.Rcheck

man/moran.Rd: R/* src/*.cpp src/*.h
	Rscript --vanilla -e 'roxygen2::roxygenize()'

