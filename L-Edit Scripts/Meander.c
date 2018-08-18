/*******************************************************************************
 * Macro Name: Meander
 * Creator  : Joseph Suttle
 *
 *		Generates a set of parallel transmission lines according to user input for line #, spacing, and thickness
 *
 * Revision History:
 * April 18th 2012
 *******************************************************************************/
#define Scale 1000 
 
module Meander_macro_module
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

	void Meander(void)
	{
		LCell	pCell	=	LCell_GetVisible();
		LFile	pFile	=	LCell_GetFile(pCell);
		
		
		
		int nLines, idx;
		float width, spacing, l, Xc, Yc;
		
		LCoord Xi, Yi;
		LLayer pLayer;
		LPoint ptArray[5];
		
		
		
		LDialogItem Dialog_Items[7] = {{"Number of Lines", "1000"}, {"Wire Width", "5"}, {"Wire Spacing", "2"}, {"X Center", "0"}, {"Y Center", "0"},{"Layer", "base"},{"Transmission Line Length","10000"}};
	
		if ( LDialog_MultiLineInputBox ( "Curve Parameters", Dialog_Items, 7 ) )
      {
          nLines = atoi(Dialog_Items[0].value);
          width = atof(Dialog_Items[1].value)*Scale;
          spacing = atof(Dialog_Items[2].value)*Scale;
          Xc = atof(Dialog_Items[3].value)*Scale;
          Yc = atof(Dialog_Items[4].value)*Scale;
		  pLayer = LLayer_Find (pFile, Dialog_Items[5].value);
          l = atof(Dialog_Items[6].value)*Scale;
          
      }  
      else
      {
      	return;
      }

		/*LObject archetype = LWire_New(pCell, pLayer,NULL, 0, ptArray, 0);
		LWire_SetWidth( pCell, archetype, width);
		
		LCoord w = LWire_GetWidth(archetype);
		LJoinType j = LWire_GetJoinType(archetype);
		LCapType  cap = LWire_GetCapType(archetype);
		LCoord   miter = LWire_GetMiterAngle(archetype);*/
		
		LWireConfig derp = {width,0,0,0};
		
		for (idx=0; idx<nLines/2; idx++)
		{	
                        //if(idx%2==0)
			//{
				Xi = Xc; Yi = Yc-(width+spacing)*2*idx;
				ptArray[0] = LPoint_Set(Xi, Yi);
				Xi = Xc+l; Yi = Yc-(spacing+width)*2*idx;
				ptArray[1] = LPoint_Set(Xi, Yi);				
				Xi = Xc+l; Yi = Yc-(spacing+width)*(2*idx+1);
				ptArray[2] = LPoint_Set(Xi, Yi);
				Xi = Xc; Yi = Yc-(spacing+width)*(2*idx+1);
				ptArray[3] = LPoint_Set(Xi, Yi);
				Xi = Xc; Yi = Yc-(spacing+width)*(2*(idx+1));
				ptArray[4] = LPoint_Set(Xi, Yi);
				LWire_New(pCell, pLayer, &derp , 1, ptArray, 5);
			/*} else
			{
				Xi = 2*d+l; Yi = width*idx;
				ptArray[0] = LPoint_Set(Xi, Yi);
				Xi = d+l; Yi = (spacing+width)*idx;
				ptArray[1] = LPoint_Set(Xi, Yi);	
				Xi = d; Yi = (spacing+width)*idx;
				ptArray[2] = LPoint_Set(Xi, Yi);
				Xi = 0; Yi = width*idx;
				ptArray[3] = LPoint_Set(Xi, Yi);
				LWire_New(pCell, pLayer,&derp, 1, ptArray, 4);
			}*/
		   
	 	}
	 	
	 
	 	
	 	//LObject wire = LWire_New(pCell, pLayer,NULL, 0, ptArray, 4*idx);
		//LWire_SetWidth( pCell, wire, width);	

	}

	void Meander_register(void)
	{
		LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "Meander", "Meander", NULL /*hotkey category*/);
	}

}
Meander_register();
