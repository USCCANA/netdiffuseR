all:
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
	head -n 80 inst/NEWS
check:
	cd ..&&R CMD build netdiffuseR/ && \
		R CMD check --as-cran netdiffuseR*.tar.gz

checkv:
	cd ..&&R CMD build netdiffuseR/ && \
		R CMD check --as-cran --use-valgrind netdiffuseR*.tar.gz

