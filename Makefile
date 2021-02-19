clean_read_TFlux:
	-rm read_TFlux	
read_TFlux:clean_read_TFlux
	g++ -o  read_TFlux read_TFlux.cc -O -std=c++11 `root-config --libs` -I. -I/home/xji/data0/software/root_build/include -L/home/xji/data0/software/root_build/lib

clean: clean_read_TFlux
	rm -f canv* *.d *.so *.pcm *~
