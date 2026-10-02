CXXFLAGS = -Wall -Wextra -std=c++23 -O
LINK.o = $(CXX) $(LDFLAGS) # This does the same a the line below
LDLIBS = -lstdc++ -lm #L is for link

all: Out.txt data.svg intensities.svg

Out.txt : main
	./main > Out.txt 


main : main.o Bunto.o hf.o

main.o : main.cc Bunto.h hf.h

hf : hf.o

hf.o : hf.h hf.cc

Bunto : Bunto.o hf.o

Bunto.o : Bunto.h Bunto.cc hf.h

data.svg : data.txt Makefile
	echo '\
	set terminal svg ;\
	set output "$@" ;\
	set xrange[0:25] ;\
	set key outside right ;\
	set multiplot layout 3,1 ;\
	set ylabel "Amplitude" ;\
	set format y "%.0t" ;\
	set label "10^{-15}" at graph 0.01,1.05 ;\
	plot \
		"$<" using 1:2 index 0 with linespoints title "{/Symbol z} = 0" ,\
		"$<" using 1:2 index 1 with linespoints title "{/Symbol z} = r" ,\
		"$<" using 1:2 index 2 with linespoints title "{/Symbol z} = ir" ,\
		"$<" using 1:2 index 3 with linespoints title "{/Symbol z} = -r" ;\
	unset label ;\
	set ylabel "Contrast C" ;\
	set yrange[-0.2:1.3] ;\
	set ytics 1 ;\
	set mytics 5 ;\
	plot \
		"$<" using 1:3 index 0 with linespoints title "{/Symbol z} = 0" ,\
		"$<" using 1:3 index 1 with linespoints title "{/Symbol z} = r",\
		"$<" using 1:3 index 2 with linespoints title "{/Symbol z} = ir",\
		"$<" using 1:3 index 3 with linespoints title "{/Symbol z} = -r";\
	set ylabel "{/Symbol q}" ;\
	set xlabel "Energy" ;\
	set yrange[-1:1] ;\
	set ytics 1 ;\
	set mytics 5 ;\
	plot \
		"$<" using 1:4 index 0 with linespoints title "{/Symbol z} = 0" ,\
		"$<" using 1:4 index 1 with linespoints title "{/Symbol z} = r" ,\
		"$<" using 1:4 index 2 with linespoints title "{/Symbol z} = ir" ,\
		"$<" using 1:4 index 3 with linespoints title "{/Symbol z} = -r" ;\
	unset multiplot'|tee log.data|gnuplot

intensities.svg : intensities.txt Makefile
	echo '\
	set terminal svg ;\
	set output "$@" ;\
	set xrange[0:2*pi] ;\
	set key outside right ;\
	set format y "%.0t" ;\
	set multiplot layout 4,1 ;\
	set ylabel "Intensity spectrum" ;\
	set label "10^{-x}" at graph 0.01,1.05 ;\
	set label "a) {/Symbol z} = 0" at graph 0.1,0.8 ;\
	plot \
		"$<" using 1:2 index 0 with linespoints title "" ,\
		"$<" using 1:3 index 0 with linespoints title "" ,\
		"$<" using 1:4 index 0 with linespoints title "" ,\
		"$<" using 1:5 index 0 with linespoints title "" ;\
	unset label ;\
	set label "a) {/Symbol z} = r" at graph 0.1,0.8 ;\
	plot \
		"$<" using 1:2 index 1 with linespoints title "" ,\
		"$<" using 1:3 index 1 with linespoints title "" ,\
		"$<" using 1:4 index 1 with linespoints title "" ,\
		"$<" using 1:5 index 1 with linespoints title "" ;\
	unset label ;\
	set label "a) {/Symbol z} = ir" at graph 0.1,0.8 ;\
	plot \
		"$<" using 1:2 index 2 with linespoints title "" ,\
		"$<" using 1:3 index 2 with linespoints title "" ,\
		"$<" using 1:4 index 2 with linespoints title "" ,\
		"$<" using 1:5 index 2 with linespoints title "" ;\
	unset label ;\
	set label "a) {/Symbol z} = -r" at graph 0.1,0.8 ;\
	plot \
		"$<" using 1:2 index 3 with linespoints title "E = 1.892 eV" ,\
		"$<" using 1:3 index 3 with linespoints title "E = 8.088 eV" ,\
		"$<" using 1:4 index 3 with linespoints title "E = 14.29 eV" ,\
		"$<" using 1:5 index 3 with linespoints title "E = 20.48 eV" ;\
	unset multiplot'|tee log.data|gnuplot


.PHONY:clean 
clean: # Forces to remove the stuff listed
	$(RM) main main.o Out.txt *.o out* err *.svg *.log log.* *.Identifier times.txt data.txt

test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM) 
