/*******************************************************************************
 * Macro Name: Circle
 * Creator  : Joseph Suttle
 *
 *		Generates a circle, dude.
 *
 * Revision History:
 * Fo Tweny 2012
 *******************************************************************************/
#define Scale 1000 
#define PI 3.141592653589796
 
module Circle_macro_module
{
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <malloc.h>
#include "ldata.h"
/* Begin -- Remove this block if you are not using L-Comp. */
#include "lcomp.h"
/* End */

	void Circle(void)
	{
		LCell	pCell	=	LCell_GetVisible();
		LFile	pFile	=	LCell_GetFile(pCell);
		
		
		
		int idx, jdx, divisions,i;
		float width,theta, r, Xc, Yc;
		
		LCoord Xi, Yi;
		LLayer pLayer;
		
		
		
		
		LDialogItem Dialog_Items[5] = {{"Wire Width", "5"}, {"X Center", "0"}, {"Y Center", "0"},{"Layer", "Boundary"},{"Radius","50800"}};
	
		if ( LDialog_MultiLineInputBox ( "Curve Parameters", Dialog_Items, 5 ) )
      {
          width = atof(Dialog_Items[0].value)*Scale;
          Xc = atof(Dialog_Items[1].value)*Scale;
          Yc = atof(Dialog_Items[2].value)*Scale;
		  pLayer = LLayer_Find (pFile, Dialog_Items[3].value);
		  r = atof(Dialog_Items[4].value)*Scale;
          
      }  
      else
      {
      	return;
      }
		static const int points = (8+4*1000)*1000;
		LPoint ptArray[48];
		LWireConfig derp = {width,0,0,0};
		divisions=48;
		i=0;
		ptArray[i] = LPoint_Set(Xc, Yc+r);
		i++;
		
		for (idx=0; idx<(divisions); idx++)
		{	
			theta=idx*PI/24;
			Xi = Xc+r*sin(theta); Yi = Yc+r*cos(theta);
			ptArray[i] = LPoint_Set(Xi, Yi);
			i++;
	 	}
		
	 	LWire_New(pCell, pLayer, &derp , 1, ptArray, i);
		
	}

	void Circle_register(void)
	{
		LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "Circle", "Circle", NULL /*hotkey category*/);
	}

}
Circle_register();
