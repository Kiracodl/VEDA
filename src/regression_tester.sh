# automated regression tester for VEDA
# Daniel Kiracofe, 2008

# how to use this file
# generate known test case input.  These will be the driverXXXXX.xml
# that rappture generates.  rename them test1.xml test2.xml etc
# and place them in the testcases directory.
# run those input files on a known good version (i.e. "az_ddaskr test1.xml")
# and store the output files (the runXXXXX.xml files) as output1.xml
# output2.xml etc and store them in the testcases directory as well.
# then, this script will run the current build of az_ddaskr and
# compare the output to the previous output.  If your changes were good
# there should be no output (other than a list of all of the curves
# in the file). if your changes are bad, then there will be a 
# tremendous amount of output showing all the places where the 
# outputs differ.

#note: changing compiler optimizations (such as -O3) can cause the
#results to be different in 7th or 8th decimal place.  If you find
#differences, try to make sure you have not changed the Makefile

#usage: if you specify the argument "diffonly" then the comparison
# is re-run on previously generated output file

# bugs:currently hardcoded for 9 test cases. Need to make that variable
#      tolerance for comparison hardcoded.

# todo list: 


make
rm run*.xml

for I in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
  if [ $# -lt 1 ]; then 
      echo "Running case $I"
      ./az_ddaskr ../testcases/test${I}.xml 2>&1 > /dev/null
  else echo "Diff only case $I"
  fi


  mv run*.xml output${I}.xml 2> /dev/null

  echo "Comparing results"
  gzip -df ../testcases/output${I}.xml.gz 2> /dev/null

#      diff -y --suppress-common-lines output${I}.xml ../testcases/output${I}.xml | grep -v time | wc
  diff -y output${I}.xml ../testcases/output${I}.xml | grep -v time | ./compare.pl
  
done
