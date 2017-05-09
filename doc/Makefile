NEMSDOC=markdown.md BUILD.md NEWTEST.md OLDTEST.md
#MODELDOC=modeldoc.md  README.GFS.md README.NMM.md

ifneq ($(wildcard ../../doc/*.md),)
APPDOC=../../doc/*.md
else
APPDOC=
endif

ALLDOC=$(APPDOC) $(NEMSDOC) # $(MODELDOC)


all: README.html

clean:
	rm -f README.html

README.html: $(ALLDOC) Makefile
	./md2html.py $(ALLDOC) README.html
	head README.html

deliver: README.html README.css
	scp README.{css,html} \
	   samuel.trahan@dmzgw.ncep.noaa.gov:/home/www/emc/htdocs/projects/rt/doc/.