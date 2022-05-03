all: package

guis:
	cd rappture; make

force_look:
	true

bin: force_look
	cd src;  make clean install

package: veda.zip

veda.zip: bin guis 
	zip veda$(VERSION).zip bin/*.pl bin/*.m bin/az_ddaskr bin/force_viewer rappture/*/tool.xml rappture/examples/*/*xml Release_Notes.txt LICENSE.txt README-bin.txt