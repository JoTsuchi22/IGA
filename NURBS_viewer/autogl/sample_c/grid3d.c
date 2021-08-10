/*  AutoGL $B!c$*<j7Z!d(B $B$N%5%s%W%k%W%m%0%i%`Nc(B */

/* $B%f!<%F%#%j%F%#$N2D;k2=5!G=$NMxMQ(B */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL $B!I$*$F$,$k!I$r;H$&$K$O(B,$B$3$N%X%C%@%U%!%$%k$r%$%s%/%k!<%I$9$k!#(B*/
#include <autogl.h>

/* $B;0<!85:9J,K!2r@O7k2L$N2D;k2=(B */
/* $BEyCMLL$HCGLL$rI=<($9$k!#(B */
/* $B$3$l$O!"E57?E*$J;0<!85%"%W%j%1!<%7%g%s$NNc$G$"$k!#(B */

/* $B%f!<%F%#%j%F%#$N2D;k2=!JEyCMLL(B,$BCGLL!K$NMxMQ(B */
/* $B0z?t$J$I>\$7$$>pJs$O(B, autogl_visualization.h$B$r;2>H!#(B*/
/* $B$J$*!"$3$l$O(Bautogl.h$B$+$i(Binclude$B$5$l$F$$$k!#(B*/



/* $B%S%e!<>e$KI=<($5$l$k%b%G%k$rI=8=$9$k$?$a$NJQ?t72(B */

/* $BD>8r:9J,3J;R$K$h$k;0<!85%9%+%i!<>l(B */

/* $B3J;R$NBg$-$5!JH>J,!K(B */
#define HALF_GRID_SIZE 30.0

/* $B3J;R$NJ,3d?t(B */
#define CELLS 40

/* $B3J;R>e$N%9%+%i!<CM(B */
static double GridValues[CELLS + 1][CELLS + 1][CELLS + 1];

/* $B%3%s%?(B-$B$N:GBg!&:G>.CM(B */
static double MinRange = -1.0;
static double MaxRange = 5.0;

/* $B??$sCf$NEyCMLL$NCM(B */
static double LevelValue = 2.0;
/* $BI}!JB>$NFsKg$NEyCMLL$O!"(BLevelValue +- LevelWidth $B$H$J$k!#!K(B */
static double LevelWidth = 2.0;

/* $BCGLL$NJ}Dx<0(B Ax+By+Cz+D=0 $B$N78?t(B */
static double A = 1.0;
static double B = 1.0;
static double C = 1.0;
static double D = 0.0;

/* $BCGLL$rI=<($9$k$+(B */
static int SectionIsShown = 0;

/* $BEyCMLL$rI=<($9$k$+(B */
static int IsosurfaceIsShown = 0;



/* $B;0<!85%9%+%i!<>l$NDj5A(B */
static double Function (double x, double y, double z)
{
  return sin (x * M_PI) + y * y + z * z;
}

/* $B;0<!853J;R$r=i4|2=$9$k!#(B */
static void InitializeGrid ()
{
  int i, j, k;

  for (i = 0; i <= CELLS; i++) {
    for (j = 0; j <= CELLS; j++) {
      for (k = 0; k <= CELLS; k++) {
	double scale = 4.0 / CELLS;
	
	GridValues[i][j][k] 
	  = Function ((i - CELLS / 2) * scale, 
		      (j - CELLS / 2) * scale,
		      (k - CELLS / 2) * scale);
      }
    }
  }
}

/* $B3J;R$NCGLL%3%s%?(B-$B$rIA2h$9$k!#(B */
static void PlotSection (void)
{
  double cellSize = HALF_GRID_SIZE * 2 / CELLS;
  int gx, gy, gz;

  for (gx = 0; gx < CELLS; gx++) {
    for (gy = 0; gy < CELLS; gy++) {
      for (gz = 0; gz < CELLS; gz++) {
	double x000 = (gx - CELLS / 2) * cellSize;
	double y000 = (gy - CELLS / 2) * cellSize;
	double z000 = (gz - CELLS / 2) * cellSize;
	double value000 = GridValues[gx][gy][gz];
	
	double x100 = x000 + cellSize;
	double y100 = y000;
	double z100 = z000;
	double value100 = GridValues[gx + 1][gy][gz];

	double x110 = x000 + cellSize;
	double y110 = y000 + cellSize;
	double z110 = z000;
	double value110 = GridValues[gx + 1][gy + 1][gz];

	double x010 = x000;
	double y010 = y000 + cellSize;
	double z010 = z000;
	double value010 = GridValues[gx][gy + 1][gz];
	
	double x001 = x000;
	double y001 = y000;
	double z001 = z000 + cellSize;
	double value001 = GridValues[gx][gy][gz + 1];

	double x101 = x000 + cellSize;
	double y101 = y000;
	double z101 = z000 + cellSize;
	double value101 = GridValues[gx + 1][gy][gz + 1];

	double x111 = x000 + cellSize;
	double y111 = y000 + cellSize;
	double z111 = z000 + cellSize;
	double value111 = GridValues[gx + 1][gy + 1][gz + 1];
	
	double x011 = x000;
	double y011 = y000 + cellSize;
	double z011 = z000 + cellSize;
	double value011 = GridValues[gx][gy + 1][gz + 1];

	/* $BO;LLBN$NCGLL%3%s%?(B-$B$rIA2h$9$k!#(B*/
	AutoGL_DrawSectionSmoothContourHexahedron 
	  (A, B, C, D, MinRange, MaxRange,
	   x000, y000, z000, value000,
	   x100, y100, z100, value100,
	   x110, y110, z110, value110,
	   x010, y010, z010, value010,
	   x001, y001, z001, value001,
	   x101, y101, z101, value101,
	   x111, y111, z111, value111,
	   x011, y011, z011, value011);
      }
    }
  }
}

/* $B3J;R$NEyCMLL$rIA2h$9$k!#(B */
static void PlotIsosurfaceOfLevel (double levelValue)
{
  double cellSize = HALF_GRID_SIZE * 2 / CELLS;
  double grade;
  double red, green, blue;
  int gx, gy, gz;

  /* 0.0 - 1.0 $B$N%0%l!<%ICM$KBP1~$9$k?'$rF@$k!#(B */
  grade = (levelValue - MinRange) / (MaxRange - MinRange);
  AutoGL_GetContourColor (&red, &green, &blue,
			  grade);
  AutoGL_SetColor (red, green, blue);

  for (gx = 0; gx < CELLS; gx++) {
    for (gy = 0; gy < CELLS; gy++) {
      for (gz = 0; gz < CELLS; gz++) {
	double x000 = (gx - CELLS / 2) * cellSize;
	double y000 = (gy - CELLS / 2) * cellSize;
	double z000 = (gz - CELLS / 2) * cellSize;
	double value000 = GridValues[gx][gy][gz];
	
	double x100 = x000 + cellSize;
	double y100 = y000;
	double z100 = z000;
	double value100 = GridValues[gx + 1][gy][gz];

	double x110 = x000 + cellSize;
	double y110 = y000 + cellSize;
	double z110 = z000;
	double value110 = GridValues[gx + 1][gy + 1][gz];

	double x010 = x000;
	double y010 = y000 + cellSize;
	double z010 = z000;
	double value010 = GridValues[gx][gy + 1][gz];
	
	double x001 = x000;
	double y001 = y000;
	double z001 = z000 + cellSize;
	double value001 = GridValues[gx][gy][gz + 1];

	double x101 = x000 + cellSize;
	double y101 = y000;
	double z101 = z000 + cellSize;
	double value101 = GridValues[gx + 1][gy][gz + 1];

	double x111 = x000 + cellSize;
	double y111 = y000 + cellSize;
	double z111 = z000 + cellSize;
	double value111 = GridValues[gx + 1][gy + 1][gz + 1];
	
	double x011 = x000;
	double y011 = y000 + cellSize;
	double z011 = z000 + cellSize;
	double value011 = GridValues[gx][gy + 1][gz + 1];

	/* $BO;LLBN$NEyCMLL$rIA2h$9$k!#(B*/
	AutoGL_DrawIsosurfaceHexahedron 
	  (levelValue, levelValue, 0,
	   x000, y000, z000, value000,
	   x100, y100, z100, value100,
	   x110, y110, z110, value110,
	   x010, y010, z010, value010,
	   x001, y001, z001, value001,
	   x101, y101, z101, value101,
	   x111, y111, z111, value111,
	   x011, y011, z011, value011);
      }
    }
  }
}

/* $BEyCMLL$rIA2h$9$k!#(B */
static void PlotIsosurface (void)
{
  /* $B;0Kg$NEyCMLL$rIA2h$9$k!#(B */
  PlotIsosurfaceOfLevel (LevelValue - LevelWidth);
  PlotIsosurfaceOfLevel (LevelValue);
  PlotIsosurfaceOfLevel (LevelValue + LevelWidth);
}



/* $B%S%e!<$N:FIA2h$N$?$a$N%3!<%k%P%C%/4X?t(B */
/* $B%S%e!<$,:FI=<($5$l$k$4$H$K8F$P$l$k!#(B */
static void RedrawView (void) 
{
  /* $B2hLL:82<$K;0<!85:BI87O$rIA2h$9$k!#(B */
  AutoGL_DrawCoordinateSystem3D (20);

  if (SectionIsShown) {
    /* $B2hLL1&B&$K%3%s%?(B-$B$N%P!<$rIA2h$9$k!#(B */
    AutoGL_DrawContourMap (80, MinRange, MaxRange);
  }

  /* $B3J;R$NOH$rIA2h$9$k!#(B */
  {
    double size = HALF_GRID_SIZE;

    AutoGL_SetColor (1, 1, 1);
    AutoGL_DrawBox3D (-size, -size, -size,
		      size, size, size);
  }

  /* $BH>F)L@I=<($rM-8z$K$9$k!#(B */
  AutoGL_TurnOnTransparency ();

  /* $B%G%#%9%W%l%$%j%9%H$KF~$l$?IA2hL?Na$r<B:]$KIA2h$9$k!#(B */
  AutoGL_DrawDisplayList ();
}

/* $BIA2h%"%C%W%G!<%H%3%^%s%I$N$?$a$N%3!<%k%P%C%/4X?t(B */
/* "UpdateModel"$B%\%?%s$rA*$V$H8F$P$l$k!#(B */
static void UpdateModel (void) 
{
  /* $B%G%#%9%W%l%$%j%9%H$r3+$-!"$3$l$r6u$K$9$k!#(B */
  AutoGL_OpenDisplayList ();

  /* $B0J2<$NIA2hL?Na$r%G%#%9%W%l%$%j%9%H$K=g<!F~$l$F$$$/!#(B */
  /* $B<B:]$NIA2h$O$3$3$G$O$^$@9T$o$l$J$$!#(B */

  if (SectionIsShown) {
    PlotSection ();
  }

  if (IsosurfaceIsShown) {
    PlotIsosurface ();
  }

  /* $B%G%#%9%W%l%$%j%9%H$rJD$8$k!#(B */
  AutoGL_CloseDisplayList ();

  /* $B%S%e!<$r:FIA2h$9$k!#(B */
  AutoGL_DrawView (); 
  /* $BEPO?$5$l$?%S%e!<$N:FIA2h4X?t!J$3$N>l9g$O(BRedrawView ()$B!K$,8F$P$l$k!#(B */
}



/* $B4X?t(BAutoGL_SetUp$B$O%f!<%6!<B&%W%m%0%i%`$4$H$KI,$:0l$DMQ0U$9$k$3$H!#(B*/
/* $B$3$3$G!"%3!<%k%P%C%/4X?t$d@)8fJQ?t$J$I$rEPO?$9$k!#(B*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  InitializeGrid ();
  
  /* $B%S%e!<$r=`Hw$9$k!#(B */

  AutoGL_SetViewSize (80); 

  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetDefaultCallbackInMode3D (NULL);   /* $B;0<!85%"%W%jMQ%$%Y%s%H=hM}(B */

  /* $B$b$7IA2h$,=E$1$l$P%I%i%C%0Cf$OIA2h$5$;$J$$$[$&$,$h$$!#(B */
  AutoGL_EnableDragInMode3D ();

#if 0
  /* $B%Q%M%k(Bviewing$B$r=`Hw$9$k!#(B */
  /* $B%f!<%6!<$,IaCJMQ$$$J$$Hf3SE*%^%$%J!<$J5!G=!#(B */

  AutoGL_AddGroup ("viewing");    
  AutoGL_SetPanelInView3D ();        /* $B;0<!85%"%W%jMQ%S%e!<%$%s%0JQ?t(B */
  AutoGL_SetPanelForInterpreter ();  /* $B%3%^%s%I%i%$%s=hM}(B */
  AutoGL_SetPanelForSave ();         /* $B2hLL%;!<%V5!G=(B */
#endif

  /* $B%Q%M%k(Bmodel$B$r=`Hw$9$k!#(B */
  /* $B%f!<%6!<$,IQHK$KMQ$$$k5!G=!#(B */
  /* $B$*$h$S!"3F%W%m%0%i%`$4$H$KFCM-$J5!G=!#(B */

  AutoGL_AddGroup ("model");    
  AutoGL_SetPanelInMode3D ();
  AutoGL_AddComment ();
  AutoGL_SetLabel ("Model Parameters");

  /* $BCGLL%3%s%?(B- */
  AutoGL_AddBoolean (&SectionIsShown, "SectionIsShown"); 
  AutoGL_SetLabel ("section");
  AutoGL_AddComment ();
  AutoGL_SetLabel ("equation Ax+By+Cz+D=0");
  AutoGL_AddReal (&A, "A");      
  AutoGL_AddReal (&B, "B");      
  AutoGL_AddReal (&C, "C");      
  AutoGL_AddReal (&D, "D");      
  AutoGL_AddReal (&MinRange, "MinRange");      
  AutoGL_AddReal (&MaxRange, "MaxRange");      

  /* $BEyCMLL(B */
  AutoGL_AddBoolean (&IsosurfaceIsShown, "IsosurfaceIsShown"); 
  AutoGL_SetLabel ("iso-surface");
  AutoGL_AddReal (&LevelValue, "LevelValue");      
  AutoGL_AddReal (&LevelWidth, "LevelWidth");      

  AutoGL_AddCallback (UpdateModel, "UpdateModel");
}