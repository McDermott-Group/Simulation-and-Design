/*******************************************************************************
* Macro Name: CPW_DrawStraight
* Creator   : Edward Leonard
*
*		Generates straight CPW based on user input
*
* Revision History:
* June 25, 2016
*******************************************************************************/
#define Scale 1000
#define PI 3.141592653589796

module CPWDrawStraight_macro_module
{
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <malloc.h>
#include "ldata.h"
#include "lcomp.h"


	/* ---------------------
	/  Primary draw function
	/  --------------------- */
	void CPW_DrawStraight(void)
	{
		/*The following are L-Edit calls to determine what cell/file to draw in */
		LCell	pCell	=	LCell_GetVisible();
		LFile	pFile	=	LCell_GetFile(pCell);

		float cpwLength; 		// (user defined) The overall length of the cpw
		float gpSpacing;		// (user defined) defines the spacing between the ground plane and center conductor
		float wireWidth;		// (user defined) defines the width of the center trace
		float xCenter;			// (user defined) defines the cartesian x center point
		float yCenter;		  // (user defined) defines the cartesian y center point

		//The following are L-Edit objects used to define geometries to be drawn
		LCoord centerX, centerY, xOffSet, yOffSet, leftCornerX, leftCornerY, rightCornerX, rightCornerY;
		LLayer pLayer;

		// Pop-up dialog box asking user for input on CPW paramters
		LDialogItem Dialog_Items[6] =
		{
			{"CPW length (um)", "1000"},
			{"Ground Plane Gap (um)", "2"},
			{"Wire Width (um)", "10"},
			{"X Center (um)", "0"},
			{"Y Center (um)", "0"},
			{"Layer", "M0"},
		};

		if ( LDialog_MultiLineInputBox ( "CPW Parameters", Dialog_Items, 6 ) )
		{
			cpwLength		= atof(Dialog_Items[0].value)*Scale;
			gpSpacing 		= atof(Dialog_Items[1].value)*Scale;
			wireWidth 		= atof(Dialog_Items[2].value)*Scale;
			xCenter 		= atof(Dialog_Items[3].value)*Scale;
			yCenter 		= atof(Dialog_Items[4].value)*Scale;
			pLayer	 		= LLayer_Find (pFile, Dialog_Items[5].value);
		}
		else
		{
			return;
		}

		// Define and draw the GND cut on the +Y side of the center trace
		leftCornerX = xCenter - cpwLength/2;
		leftCornerY = yCenter + wireWidth/2;
		rightCornerX= xCenter + cpwLength/2;
		rightCornerY= yCenter + wireWidth/2 + gpSpacing;

		LBox_New(pCell, pLayer, leftCornerX, leftCornerY, rightCornerX, rightCornerY);

		// Define and draw the GND cut on the -Y side of the center trace
		leftCornerX = xCenter - cpwLength/2;
		leftCornerY = yCenter - wireWidth/2;
		rightCornerX= xCenter + cpwLength/2;
		rightCornerY= yCenter - wireWidth/2 - gpSpacing;

		LBox_New(pCell, pLayer, leftCornerX, leftCornerY, rightCornerX, rightCornerY);
	}

	/* ------------------
	/  Macro Registration
	/  ------------------ */
	void CPW_register(void)
	{
		LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "CPW_Straight", "CPW_DrawStraight", NULL /*hotkey category*/);
	}

}
CPW_register();
