all: smoothSDE.so 

smoothSDE.so: smoothSDE.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('smoothSDE.cpp','-O0 -g -Wno-ignored-attributes -I ../inst/include/', safebounds=FALSE, framework='CppAD')"
	
clean:
	rm -rf *o