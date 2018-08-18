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
#define PI 3.141592653589796
 
module CPWParaMeander_macro_module
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

	void CPWParaMeander(void)
	{
		LCell	pCell	=	LCell_GetVisible();
		LFile	pFile	=	LCell_GetFile(pCell);
		
		
		
		int nLines, idx,jdx,kdx,nMe;
		float width, spacing, Xc, Yc, l,d, g, r1, r2, CPx, CPy, ls;
		
		LCoord Xi, Yi;
		LLayer pLayer;
		LLayer pLayer2;
		LPoint ptArray[20];
		
		
		LDialogItem Dialog_Items[11] = {{"Number of Lines (Odd for sym)", "21"}, {"Wire Width", "5"}, {"Wire Spacing", "2"}, {"X Center", "0"}, {"Y Center", "0"},{"Layer", "base"},{"Number of Meanders","100"},{"Length of Each","100"},{"Diagonal Length","20"},{"CPW Gap","5"},{"Temp Gap Layer","temp"}};
		
		if ( LDialog_MultiLineInputBox ( "Curve Parameters", Dialog_Items, 11 ) )
      {
          nLines = atoi(Dialog_Items[0].value);
          width = atof(Dialog_Items[1].value)*Scale;
          spacing = atof(Dialog_Items[2].value)*Scale;
          Xc = atof(Dialog_Items[3].value)*Scale;
          Yc = atof(Dialog_Items[4].value)*Scale;
		  pLayer = LLayer_Find (pFile, Dialog_Items[5].value);
          nMe = atof(Dialog_Items[6].value);
		  l = atof(Dialog_Items[7].value)*Scale;
		  d = atof(Dialog_Items[8].value)*Scale;
		  g = atof(Dialog_Items[9].value)*Scale;
		  pLayer2 = LLayer_Find (pFile, Dialog_Items[10].value);
          
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
		ls = (width+spacing)*nLines-width; //ls is the line spacing between meanders less one parallel spacing
		
		LWireConfig derp = {width,0,0,0};
		LWireConfig herp = {width*nLines,0,0,0};
		LWireConfig berp = {width+2*g,0,0,0};
		
		Xi = Xc+l/2+d;
		Yi = Yc + l/6+d;
		r1 = (Xc+3*l/4+d-width*nLines/2)-Xi;
		r2=  Yc - (nLines/2)*spacing-width*(nLines-1)/2-Yi;
		CPx = Xi+r1;
		CPy =  Yi;
		ptArray[0] = LPoint_Set(Xi, Yi+width*nLines*1.5);
		ptArray[1] = LPoint_Set(Xi, Yi);
		
		
		for(jdx=0;jdx<=16;jdx++)
		{
			Xi = CPx-r1*cos(PI*jdx/32); Yi = CPy+r2*sin(PI*jdx/32);
			ptArray[jdx+2] = LPoint_Set(Xi, Yi);
		}
		
		Xi = Xc+3*l/4+d-width*nLines/2;
		Yi = Yc - (nLines/2)*spacing-width*(nLines-1)/2;
		
		ptArray[19] = LPoint_Set(Xi, Yi);
			LWire_New(pCell, pLayer, &herp , 1, ptArray, 20);
			LWire_New(pCell, pLayer2, &berp, 1, ptArray, 20);
		
		for(jdx=0;jdx<nLines;jdx++)
			{
			Xi = Xc+3*l/4+d; Yi = Yc - (nLines/2)*spacing-jdx*width;
			ptArray[0] = LPoint_Set(Xi, Yi);
			Xi = Xc+3*l/4+2*d; Yi = Yc-(width+spacing)*jdx;
			ptArray[1] = LPoint_Set(Xi, Yi);
			Xi = Xc+l-width/2; Yi = Yc-(width+spacing)*jdx;
			ptArray[2] = LPoint_Set(Xi, Yi);
			
			LWire_New(pCell, pLayer, &derp , 1, ptArray, 3);
			LWire_New(pCell, pLayer2, &berp , 1, ptArray, 3);
			}
		
		
		for (idx=0; idx<nMe/2; idx++)
		{	
		    for(jdx=0;jdx<nLines;jdx++)
			{
                
				r1 = (width+spacing)*(2*nLines-2*jdx-1)/2+ls/2;
				r2 = (width+spacing)*(2*jdx+1)/2+ls/2;
				Xi = Xc+width; Yi = Yc-(width+spacing)*jdx-(width+spacing)*nLines*2*idx-2*(ls)*idx;
				if(idx==0)
				{
					Xi = Xi+3*l/4+2*d;
				}
				ptArray[0] = LPoint_Set(Xi, Yi);
				Xi = Xc+l; Yi = Yc-(width+spacing)*jdx-(width+spacing)*nLines*2*idx-2*(ls)*idx;
				ptArray[1] = LPoint_Set(Xi, Yi);
				//insert a nice smooth curve
				CPx = Xi; CPy = Yi - r1; //set centerpoint of circle
				for(kdx=1;kdx<=8;kdx++)
				{
				
					Xi = CPx+r1*sin(PI*kdx/8); Yi = CPy+r1*cos(PI*kdx/8);
					ptArray[1+kdx] = LPoint_Set(Xi, Yi);
					
				
				}
				Xi = Xc+l; Yi = Yc-(width+spacing)*(2*nLines-jdx-1)-ls-(width+spacing)*nLines*2*idx-2*(ls)*idx;
				ptArray[10] = LPoint_Set(Xi, Yi);
				Xi = Xc+width/2; Yi = Yi = Yc-(width+spacing)*(2*nLines-jdx-1)-ls-(width+spacing)*nLines*2*idx-2*(ls)*idx;
				ptArray[11] = LPoint_Set(Xi, Yi);
				//some more nice smooth curves
				CPx = Xi; CPy = Yi - r2; //set centerpoint of circle
				for(kdx=1;kdx<=8;kdx++)
				{
				
					Xi = CPx-r2*sin(PI*kdx/8); Yi = CPy+r2*cos(PI*kdx/8);
					ptArray[11+kdx] = LPoint_Set(Xi, Yi);
					
				
				}
				Xi = Xi+d;
				//Xi = Xc+(width+spacing)*(nLines-jdx); Yi = Yc-(width+spacing)*(2*nLines+jdx)-(width+spacing)*nLines*2*idx;
				ptArray[20] = LPoint_Set(Xi, Yi);
				LWire_New(pCell, pLayer, &derp , 1, ptArray, 21);
				LWire_New(pCell, pLayer2, &berp , 1, ptArray, 21);
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
	 	}
		
		for(jdx=0;jdx<nLines;jdx++)
			{
			Xi = Xc+width; Yi = Yc-(width+spacing)*(nMe*nLines+jdx)-ls*2*idx;
			ptArray[0] = LPoint_Set(Xi, Yi);
			Xi = Xc+l/4; Yi =Yc-(width+spacing)*(nMe*nLines+jdx)-ls*2*idx;
			ptArray[1] = LPoint_Set(Xi, Yi);
			Xi = Xc+l/4+d; Yi =Yc-(width+spacing)*(nMe*nLines)-nLines/2*spacing-jdx*width-ls*2*idx;
			ptArray[2] = LPoint_Set(Xi, Yi);
			
			LWire_New(pCell, pLayer, &derp , 1, ptArray, 3);
			LWire_New(pCell, pLayer2, &berp , 1, ptArray, 3);
			}
			
		Xi = Xc+l/4+d+width*nLines/2;
		Yi = Yc-(width+spacing)*(nMe*nLines)-ls*2*idx-(nLines-1)*width/2-(nLines-1)*spacing/2;
		ptArray[0] = LPoint_Set(Xi, Yi);
		r1 = Xc+l/2+d-Xi;
		r2 = Yc-(width+spacing)*(nMe*nLines)-ls*2*idx-(nLines-1)*width/2-(nLines-1)*spacing/2-l/6-d-Yi;
		CPx = Xi;
		CPy = Yi+r2;
		
		
			for(jdx=0;jdx<=16;jdx++)
		{
			Xi = CPx+r1*sin(PI*jdx/32); Yi = CPy - r2*cos(PI*jdx/32);
			ptArray[jdx+1] = LPoint_Set(Xi, Yi);
		}
		
		Xi = Xc+l/2+d;
		Yi = Yc-(width+spacing)*(nMe*nLines)-ls*2*idx-(nLines-1)*width/2-(nLines-1)*spacing/2-l/6-d;
		ptArray[17] = LPoint_Set(Xi, Yi);
		Yi = Yi - (width)*nLines*1.5;
		ptArray[18] = LPoint_Set(Xi, Yi);
	    LWire_New(pCell, pLayer, &herp , 1, ptArray, 19);
	    LWire_New(pCell, pLayer2, &berp , 1, ptArray, 19);
		
		

	}

	void CPWParaMeander_register(void)
	{
		LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "CPWParaMeander", "CPWParaMeander", NULL /*hotkey category*/);
	}

}
CPWParaMeander_register();
