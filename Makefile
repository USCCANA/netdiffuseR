netdiffuseR.tar.gz: */*.R
	$(MAKE) clean ; \
	R CMD build . && \
		mv netdiffuseR*.tar.gz netdiffuseR.tar.gz

inst/NEWS: NEWS.md
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
	head -n 80 inst/NEWS

README.md: README.Rmd
	Rscript -e 'rmarkdown::render("README.Rmd")'

.PHONY: check checkv clean

check: netdiffuseR.tar.gz
	R CMD check --as-cran netdiffuseR.tar.gz

checkv: netdiffuseR.tar.gz
	R CMD check --as-cran --use-valgrind netdiffuseR.tar.gz

clean:
	rm -rf netdiffuseR.Rcheck ; rm -f netdiffuseR.tar.gz
