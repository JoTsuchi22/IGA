all: NURBS_input_s

clean:
	rm -f *.o NURBS_input_s.x

NURBS_input_s: constant.h NURBS_calc.c NURBS_calc.h NURBS_input_for_s-IGA.c
	gcc -c -o NURBS_calc.o NURBS_calc.c -Wall -mcmodel=large
	gcc -c -o NURBS_input_for_s-IGA.o NURBS_input_for_s-IGA.c -Wall -mcmodel=large
	gcc NURBS_calc.o NURBS_input_for_s-IGA.o -o NURBS_input_s.x -Wall -mcmodel=large
