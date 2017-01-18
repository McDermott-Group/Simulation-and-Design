/*******************************************************************************
 * Macro Name: CPW2
 * Creator  : Matthew Beck
 *
 *		Generates a meandering CPW based upon user input
 *
 * Revision History:
 * May 17, 2013
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

void CPW2(void)
{
	/*The following are L-Edit calls to determine what cell/file to draw in */
	LCell	pCell	=	LCell_GetVisible();
	LFile	pFile	=	LCell_GetFile(pCell);
	
	float cpwLength; 		// (user defined) The OAL of the cpw
	float meanderWidth; 	// (user defined) defines the length of the straight portion of the meander line
	float meanderSpace; 	//	(user defined) defines the spacing between the meander triplet lines.
	float gpSpacing;		// (user defined) defines the spacing between the ground plane and center conductor; 
	float wireWidth;		// (user defined) defines the width of the ground plane and center conductor wires;
	float meanderRadius; // (calculated)	  defines the radius needed to turn the meander
	float xCenter;			// (user defined) defines the cartesian x center point for the meander;
	float yCenter;		   // (user defined) defines the cartesian y center point for the meander;
	float numMeander;		// length of the CPW / length of meander plus radius.
	float radiusCenterX;	// This determines the center of the radius need to draw the connection between wires.
	float radiusCenterY;
	float meanderLength;
	//The following are L-Edit calls do define what layer needs to be drawn 
	LCoord centerX, centerY, xOffSet, yOffSet, leftCornerX, leftCornerY, rightCornerX, rightCornerY;
	LCoord centerX1, centerY1, leftCornerX1, leftCornerY1, rightCornerX1, rightCornerY1;
	LCoord centerX2, centerY2, leftCornerX2, leftCornerY2, rightCornerX2, rightCornerY2;
	LLayer pLayer;
	LLayer pLayer2;
	LPoint ptArray[20];
		
		
	LDialogItem Dialog_Items[8] = {
										{"Overall CPW length (um)", "21"}, 
										{"meander length (um)", "5"}, 
										{"Spacing between meanders","50"},
										{"Ground Plane Spacing", "2"}, 
										{"Wire Width", "10"},
										{"X Center", "0"}, 
										{"Y Center", "0"},
										{"Layer", "Base"},
									};
		
	if ( LDialog_MultiLineInputBox ( "Curve Parameters", Dialog_Items, 8 ) )
	{
		cpwLength		= atoi(Dialog_Items[0].value);
		meanderWidth 	= atof(Dialog_Items[1].value)*Scale;
		meanderSpace 	= atof(Dialog_Items[2].value)*Scale;
		gpSpacing 		= atof(Dialog_Items[3].value)*Scale;
		wireWidth 		= atof(Dialog_Items[4].value)*Scale;
		xCenter 			= atof(Dialog_Items[5].value)*Scale;
		yCenter 			= atof(Dialog_Items[6].value)*Scale;
		pLayer2 			= LLayer_Find (pFile, Dialog_Items[7].value);
			
	}  
	else
	{
		return;
	}
	  
	//	float meanderSpacing = (wireWidth + gpSpacing)*3 - wireWidth;
	  
	//	meanderRadius = 2*(3*wireWidth/2 + gpSpacing) + meanderSpacing;
	//	LWireConfig centerTrace = {wireWidth,0,0,0};
	//	currentX = xCenter;
	//	currentY = yCenter + radius;
	//	ptArray[0] = LPoint_Set(currentX, currentY);
	//	ptArray[1] = LPoint_Set(xCenter,yCenter);
		float IR;
		float OR;
		float AR;
		float Pi = 3.14159;
		float xCenterRadius;
		IR = (2*gpSpacing+2*wireWidth + meanderSpace)/2;
		OR = IR + wireWidth;
		AR = (OR + IR)/2 / Scale;
	  	meanderLength = meanderWidth/Scale + 2*Pi*AR;
		numMeander = cpwLength / meanderLength;
	  	float newWireOffSetX = 3*wireWidth + 2*gpSpacing + meanderSpace;
	  	float startAngle, endAngle;
	  	int n=0;
	  	LCoord groundWireOffSet;
	  	groundWireOffSet = wireWidth/2 + gpSpacing/2;
	  	//Place the wires for the center trace
	  	/***********************************************************************************************/
	leftCornerX1 = xCenter - cpwLength/2;
	leftCornerY1 = yCenter + wireWidth/2;
	rightCornerX1= xCenter + cpwLength/2;
	rightCornerY1= yCenter + wireWidth/2 + gpSpacing;
	
	LBox_New(pCell, pLayer2, leftCornerX, leftCornerY, rightCornerX, rightCornerY);
	
	leftCornerX2 = xCenter - cpwLength/2;
	leftCornerY2 = yCenter - wireWidth/2;
	rightCornerX2= xCenter + cpwLength/2;
	rightCornerY2= yCenter - wireWidth/2 - gpSpacing;
	
	LBox_New(pCell, pLayer2, leftCornerX, leftCornerY, rightCornerX, rightCornerY);  
		/**************************************************************************************************/	
		xCenterRadius = 3*wireWidth/2 + gpSpacing + meanderSpace/ 2;
	//	float radiusCenterY1 = (meanderWidth / 2 + yCenter);
	//	float radiusCenterY2 = (meanderWidth / 2 - yCenter);
		float radiusCenterY1 = yCenter + meanderWidth/2 ;
		float radiusCenterY2 = yCenter - meanderWidth/2;
		float smallRadiusIR, smallRadiusOR, largeRadiusIR, largeRadiusOR;
		
		smallRadiusIR = meanderSpace/2 + wireWidth;
		smallRadiusOR = smallRadiusIR + gpSpacing;
		largeRadiusIR = 2*wireWidth + gpSpacing + meanderSpace/2;
		largeRadiusOR = largeRadiusIR + gpSpacing;
		//place the curves for the center trace
		/*****************************************************************************************************/
		for(n = 0; n < numMeander-1; n++)
		{
		//	radiusCenterX = (2*n+1)*xCenterRadius + xCenter;
			radiusCenterX = xCenter + (2*n+1)*xCenterRadius;
		//	radiusCenterY = -1*radiusCenterY;
		/*	if(radiusCenterY > 0){
				startAngle = 0;
				endAngle = 180;
			}
			else{
				startAngle = 180;
				endAngle = 360;
			}*/
			if(n%2 ==0){
			
				radiusCenterY = radiusCenterY1;
				startAngle = 0;
				endAngle = 180;
			}
			else{
			
				radiusCenterY = radiusCenterY2;
				startAngle = 180;
				endAngle = 360;
			}
				
			//Center Trace
			LTorusParams tParams;
			tParams.ptCenter = LPoint_Set(radiusCenterX,radiusCenterY);
			tParams.nInnerRadius = IR;
			tParams.nOuterRadius = OR;
			tParams.dStartAngle = startAngle;
			tParams.dStopAngle = endAngle;
			//Ground Wire 1
			LTorusParams tParams1;
			tParams1.ptCenter = LPoint_Set(radiusCenterX,radiusCenterY);
			tParams1.nInnerRadius = smallRadiusIR;
			tParams1.nOuterRadius = smallRadiusOR;
			tParams1.dStartAngle = startAngle;
			tParams1.dStopAngle = endAngle;
			//Ground Wire 2
			LTorusParams tParams2;
			tParams2.ptCenter = LPoint_Set(radiusCenterX,radiusCenterY);
			tParams2.nInnerRadius = largeRadiusIR;
			tParams2.nOuterRadius = largeRadiusOR;
			tParams2.dStartAngle = startAngle;
			tParams2.dStopAngle = endAngle;
			
			//Draw the curves
			//LTorus_CreateNew(pCell, pLayer2, &tParams);
			 LTorus_CreateNew(pCell, pLayer2, &tParams1);
			 LTorus_CreateNew(pCell, pLayer2, &tParams2);
			//LObject_ConvertToPolygon(pCell , myTorus1, 1);
		}
		/***********************************************************************************************************/
		
		// place in the ground wires
		//float i = truncf(numMeander);
		int i = (int)numMeander;
		float remain = numMeander - i;
		float traceLength = numMeander*meanderLength - remain*meanderLength;
		remain = cpwLength - traceLength;
		
	   	LDialog_MsgBox(LFormat("Avg Radius = %f \n meander length %f. \n Trace Length %f \n remaining wire = %f",AR, meanderLength,traceLength,remain));
}
	  
	  void CPW_register(void)
	  {
		LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "CPW2", "CPW2", NULL /*hotkey category*/);
	  }
	  
}
CPW_register();
	  
	  
	  