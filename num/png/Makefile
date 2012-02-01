targets=mesh0 mesh1 mesh2  sol_u sol_q sol_t

all: ${targets:=.png}

%.png:../ps/%.ps
	gs -dNOPAUSE -sDEVICE=pngmono -sOutputFile=$@\
	     -g1600x1600 -r222 -q -dBATCH $<
	gs -dNOPAUSE -sDEVICE=pngmono -sOutputFile=-\
	     -g2400x2400 -r4440 -q -dBATCH $< |\
	     convert PNG:- -crop 1600x1600+580+220 -repage 0\
	     PNG8:-> ${@:.png=s.png}

