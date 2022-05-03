
img_objs := $(patsubst %.png,%.png.txt,$(wildcard *.png))

tool_dirs := amp_red_avg pmscan scan_contact forceViewer ssac fmscan fmac fzcurves_new adac dac  freqsweep freqadv scanning scanbasic fixedpoint forcemod forcevolume peakforcebasic

tool_objs := $(patsubst %, %/tool.xml, $(tool_dirs))

all_xmls := $(wildcard *.xml)
template_xmls := $(wildcard *.template.xml)
depend_xmls := $(filter-out $(template_xmls), $(all_xmls))

all : $(img_objs) $(tool_objs)

install : all

clean :
	rm -f *.png.txt
	rm -f tmp

distclean : clean
	rm -f */tool.xml

%/tool.xml : %.template.xml $(img_objs) $(depend_xmls)
	cpp -nostdinc -C $<  2> /dev/null | sed '/\#.*/d' > $@ 

%.png.txt : %.png
	base64 $< > tmp
	echo -n "<current>" | cat - tmp > $@
	echo "</current>" >> $@
	rm tmp	


