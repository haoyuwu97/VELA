target = Gt_Gw

SOURCES = correlator_mode0.cc correlator_mode1.cc correlator_Gw.cc transform_methods.cc lammps_mua_backend.cc main.cc
OBJECTS = $(SOURCES:.cc=.o)

CXXFLAGS = -O3
LDLIBS = -ldl

default: $(target)

$(target): $(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	-rm -f *.o $(target)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

depend:
	makedepend -Y. -o.o -- $(SOURCES)
# DO NOT DELETE
