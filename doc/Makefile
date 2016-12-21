NEMSDOC=markdown.md BUILD.md NEWTEST.md OLDTEST.md
APPDOC=../../doc/README.md
#MODELDOC=modeldoc.md  README.GFS.md README.NMM.md

ALLDOC=$(APPDOC) $(NEMSDOC) # $(MODELDOC)


all: README.html

clean:
	rm -f README.html

README.html: $(ALLDOC)
	./md2html.py $(ALLDOC) README.html
	head README.html

deliver: README.html README.css
	scp README.{css,html} \
	   samuel.trahan@dmzgw.ncep.noaa.gov:/home/www/emc/htdocs/projects/rt/doc/.