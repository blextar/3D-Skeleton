GPP = g++ -O3 -g -frounding-math -fno-strict-aliasing

hierskel: calc_skeleton.o hull.o field_comp.o distTrans.o octree.o misc.o skel_mesh.o
		$(GPP) -lgsl -lgslcblas -I/opt/local/include/ -L/opt/local/lib /opt/local/lib/libCGAL_Core.dylib /opt/local/lib/libCGAL.dylib /opt/local/lib/libgmpxx.dylib /opt/local/lib/libmpfr.dylib /opt/local/lib/libgmp.dylib /opt/local/lib/libboost_thread-mt.dylib calc_skeleton.o field_comp.o distTrans.o octree.o misc.o skel_mesh.o hull.o -o hierskel

hull.o: *.cc *.h
	$(GPP) -I. -I/opt/local/include/ -c hull.cc

skel_mesh.o: *.cc *.h
	$(GPP) -I. -I/opt/local/include/ -c skel_mesh.cc
	
misc.o: *.cc *.h
	$(GPP) -c misc.cc

field_comp.o:	*.cc *.h
		$(GPP) -I. -I/opt/local/include/ -c field_comp.cc
	
distTrans.o:	*.cc *.h
	$(GPP) -I. -I/opt/local/include/ -c distTrans.cc

octree.o: *.cc *.h
	$(GPP) -I. -I/opt/local/include/ -c octree.cc

calc_skeleton.o: *.cc *.h
		$(GPP) -I. -Wno-non-template-friend -I/opt/local/include/ -c calc_skeleton.cc

clean: 
	rm *.o
