Map of the directory structure

bin - executables get installed here.  Technically you can run the executables directly without
      the GUI, but it's not for beginners.  See the rappture directory for the GUIs

doc - documentation lives here.  use 'make' to generate a PDF (assuming you have 
      make and the necessary latex software installed)

src - source code, self explanatory.

rappture - this is where the GUI definitions live.  Each tool has a template xml file named
           xxx.template.xml and a directory named xxx/, where xxx is the name of the tool.
           the script generate_guis.sh will create the files named xxx/tool.xml (really this
           should be a proper Makefile). To run the GUI, cd into the directory of the tool you
           wish to run and type 'rappture'.  Note: this assumes that you have rappture installed
           on your machine (https://nanohub.org/infrastructure/rappture/)  Currently rappture is
           only supported on Linux and Mac OS so Windows users cannot run the GUI locally.
           If somebody wanted to port to Windows, that would be appreciated (hint, hint ;) 

middleware - this is stuff that we need to run the code on nanohub.  If you are running on your
             local machine you can pretty much ignore this directory

examples - empty
data - empty
testcases - empty
