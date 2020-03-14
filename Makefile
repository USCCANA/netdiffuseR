netdiffuseR.tar.gz: */*.R
	rm netdiffuseR.tar.gz ; \
	R CMD build . && \
		mv netdiffuseR*.tar.gz netdiffuseR.tar.gz

inst/NEWS: NEWS.md
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
	head -n 80 inst/NEWS
check: netdiffuseR.tar.gz
	R CMD check --as-cran netdiffuseR.tar.gz

checkv: netdiffuseR.tar.gz
	R CMD check --as-cran --use-valgrind netdiffuseR.tar.gz

