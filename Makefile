clean_read_TFlux:
	-rm read_TFlux	
read_TFlux:clean_read_TFlux
	g++ -o  read_TFlux read_TFlux.cc -O -std=c++11 `root-config --libs` -I. -I/home/xji/data0/software/root_build/include -L/home/xji/data0/software/root_build/lib
	#g++ -o  read_TFlux read_TFlux.cc -O -std=c++17 `root-config --libs` -I. -I/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+3.10-2.17-e17-prof/include -L/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+3.10-2.17-e17-prof/lib
clean: clean_read_TFlux
	rm -f canv* *.d *.so *.pcm *~
