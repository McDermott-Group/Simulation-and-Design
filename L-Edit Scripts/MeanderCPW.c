/*******************************************************************************
 * Macro Name: MeanderCPW
 * Creator  : David Hover
 *
 * Revision History:
 * 28 April 2010	Generated by L-Edit
 *******************************************************************************/

#define PI 3.141592653589796
#define Scale 1000



module Macro1_macro_module
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
	
	LPoint *ArcSegment(float rRadius, float rX_0, float rY_0, float theta_0, int nSegments, LFile pFile)
	{
		LPoint *points = (LPoint*)calloc(nSegments, sizeof(LPoint));
		float dTheta;
		float xt, yt;
		int idx;
		
		dTheta = (PI/4.)/(float)nSegments;
		
		for (idx=0; idx<(nSegments+1); idx++)
		{
			xt = rX_0 + rRadius*cos(theta_0 + idx*dTheta);
			yt = rY_0 + rRadius*cos(theta_0 + idx*dTheta);
			
			points[idx].x = LFile_DispUtoIntU(pFile, xt); 
			
		}
		
		return points;
	}
		

	void MeanderCPW(void)
	{
		LCell	pCell	=	LCell_GetVisible();
		LFile	pFile	=	LCell_GetFile(pCell);
	
		LLayer baseLayer = LLayer_Find(pFile, "Base");
		LLayer tempLayer = LLayer_Find(pFile, "Temp");
		
		LDialogItem Dialog_Items[5] = { {"Number of Meanders", "1"}, {"Length of Straight Segments (um)", "1000"}, {"Turn Radius (um)", "125"}, {"Trace Width (um)", "10"}, {"Gap Width (um)", "6"} } ;
		int nMeander;
		float tRadius;
      float tLength;
      float tWidth;
      float tGap;
      
      LCoord X0;
      LCoord Y0;
      LLayer pLayer;
      
      int arcSegments = 8;
      
      int idx;
      
      //Scale = 1000;
      
      if ( LDialog_MultiLineInputBox ( "Resonator Parameters", Dialog_Items, 6 ) )
      {
      	 Radius = atoi(Dialog_Items[0].value)*Scale;
          numSegments = atoi(Dialog_Items[1].value);
          traceWidth = atof(Dialog_Items[2].value)*Scale;
          X0 = atof(Dialog_Items[3].value);
          Y0 = atof(Dialog_Items[4].value);
          pLayer = LLayer_Find (pFile, Dialog_Items[5].value);
      }   
      else
      	return;
      	
		LPoint pt_Array[1000];   
		LPoint *temp = (LPoint*)calloc(numSegments, sizeof(LPoint));
		 
		 
		 
		for (idx=0; idx < (nMeander+1)
	}


	void MeanderCPW_register(void)
	{
		LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "MeanderCPW", "MeanderCPW", NULL /*hotkey category*/);
	}

}
MeanderCPW_register();