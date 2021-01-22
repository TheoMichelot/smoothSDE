all: smooth_sde.so 

smooth_sde.so: smooth_sde.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('smooth_sde.cpp','-O0 -g -Wno-ignored-attributes -I ../inst/include/', safebounds=FALSE)"
	
clean:
	rm -rf *o