CC				= g++
CFLAGS			= -O3 -lm -Wall -Wextra

MC3D:	./src/s_IGA_MC3D.cpp
	cd ./src && $(CC) -o s_IGA_MC3D.x s_IGA_MC3D.cpp $(CFLAGS)

# for debug
run_a:	MC3D
	cd ./textfile && ../src/s_IGA_MC3D.x ./test_boundary.txt ./test_MC3D_00_edited.txt ./test_MC3D_01_edited.txt

run_b:	MC3D
	cd ./textfile && ../src/s_IGA_MC3D.x ./test_boundary.txt ./test_MC3D_00_edited.txt ./test_MC3D_02_edited.txt

run_c:	MC3D
	cd ./textfile && ../src/s_IGA_MC3D.x ./test_boundary.txt ./test_MC3D_03_edited.txt ./test_MC3D_04_edited.txt

run_d:	MC3D
	cd ./textfile && ../src/s_IGA_MC3D.x ./test_2D_boundary.txt ./test_MC3D_05_edited.txt ./test_MC3D_06_edited.txt ./test_MC3D_07_edited.txt ./test_MC3D_08_edited.txt

run_e:
	cd ./textfile && ../src/s_IGA_MI3D.x ./test_MC3D_00.txt ./test_MC3D_01.txt ./test_MC3D_02.txt ./test_MC3D_03.txt ./test_MC3D_04.txt

run_f:
	cd ./textfile && ../src/s_IGA_MI3D.x ./test_MC3D_05.txt ./test_MC3D_06.txt ./test_MC3D_07.txt ./test_MC3D_08.txt