4/30/09: the PNGs must be mimencoded in order to be loaded into the gui.
The generate_guis file will automatically do this for you.  to add a new 
just #include the resulting .txt file in the gui.

4/9/09: each tool now has it's own separate directory ('dac', 'freqsweep', etc).
within each directory there is a tool.xml.  generate_guis.sh automatically
updates the files in these directories. To test each tool, change to that
directory and run "rappture"

3/2/09: I'm introducing a new concept for the GUI files.  As far
as I can tell (and somebody correct me if I'm wrong) rappture does not
have an "include" statement.  But managing the multiple different versions
of the GUIs is getting cumbersome.  There is a lot of duplicated material
and it's geting annoying trying to keep it all synchronized.  So here's what
we do:

Do not edit scanning.xml, freqsweep.xml, or tool.xml by hand any more.
Instead, edit scanning.template.xml, freqsweep.template.xml and 
tool.template.xml.  You'll notice that I've made use of C style macros
in the template file. For example, because all three files should always
have the same tip-sample interaction menu, I've moved it out to a separate file
and then just said 

#include "tip_sample_interaction.xml"


  This way, you only have to edit the one tip-sample file to effect all three
GUIs.

  When you are done editing the files, run the generate_guis.sh file.  It does
not need any arguments or anything.  When you are done, check everything 
into SVN (the templates and the generated xml files).

Daniel