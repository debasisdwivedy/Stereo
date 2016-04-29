all : segment stereo render

segment: CImg.h segment.cpp
	g++ -std=c++11 -pthread -Dcimg_display=0 segment.cpp -o segment -I. -O3

stereo: CImg.h stereo.cpp
	g++ -std=c++11 -pthread -Dcimg_display=0 stereo.cpp -o stereo -I. -O3

render: CImg.h render.cpp
	g++ -std=c++11 -pthread -Dcimg_display=0 render.cpp -o render -I. -O3

clean:
	rm segment stereo render
