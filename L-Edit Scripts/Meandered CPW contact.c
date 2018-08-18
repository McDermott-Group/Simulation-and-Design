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
    
    void fuck(void)
    {
        /*The following are L-Edit calls to determine what cell/file to draw in */
       
         LCell	pCell	=	LCell_GetVisible();
         LFile	pFile	=	LCell_GetFile(pCell);
         
        //int w[] = {25,30,35,40,45,50,55,60,65,70};
        //int s[] = {3,4,5,6,7,8,9,10};
        //char buffer[32];
        //int ss = s[1];
        //snprintf(buffer,sizeof(char)*32,"30 um trace %i um gap",ss);
       // puts(filename);
        //pCell = LCell_New(pFile,buffer);
        
        
        
        
        //float cpwLength; 		// (user defined) The OAL of the cpw
        //float meanderWidth; 	// (user defined) defines the length of the straight portion of the meander line
        //float meanderSpace; 	//	(user defined) defines the spacing between the meander triplet lines.
        //float gpSpacing;		// (user defined) defines the spacing between the ground plane and center conductor;
        //float wireWidth;		// (user defined) defines the width of the ground plane and center conductor wires;
        //float meanderRadius; // (calculated)	  defines the radius needed to turn the meander
        //float xCenter;			// (user defined) defines the cartesian x center point for the meander;
        //float yCenter;		   // (user defined) defines the cartesian y center point for the meander;
        //float numMeander;		// length of the CPW / length of meander plus radius.
        //float radiusCenterX;	// This determines the center of the radius need to draw the connection between wires.
        //float radiusCenterY;
        //float meanderLength;
        //The following are L-Edit calls do define what layer needs to be drawn
        LCoord centerX, centerY, xOffSet, yOffSet, leftCornerX, leftCornerY, rightCornerX, rightCornerY;
        LCoord centerX1, centerY1, leftCornerX1, leftCornerY1, rightCornerX1, rightCornerY1;
        LCoord centerX2, centerY2, leftCornerX2, leftCornerY2, rightCornerX2, rightCornerY2;
        LCoord centerX3, centerY3, leftCornerX3, leftCornerY3, rightCornerX3, rightCornerY3;
        LLayer pLayer;
        LLayer pLayer2;
        LPoint ptArray[20];
        
        
        LDialogItem Dialog_Items[3] = {
           // {"Overall CPW length (um)", "21"},
            //{"meander length (um)", "5"},
            //{"Spacing between meanders","50"},
            {"Ground Plane Spacing", "2"},
            {"Wire Width", "10"},
            {"Total Length (<6600)","6600"},
            //{"X Center", "0"},
            //{"Y Center", "0"},
           // {"Layer", "Base"},
        };
        
        float width = 50*Scale;
        int length = 6600;
        int gap = 10*Scale;
        int l0 = 6600;
        
        if ( LDialog_MultiLineInputBox ( "Curve Parameters", Dialog_Items, 3 ) )
        {
            //cpwLength		= atoi(Dialog_Items[0].value);
            //meanderWidth 	= atof(Dialog_Items[1].value)*Scale;
           // meanderSpace 	= atof(Dialog_Items[2].value)*Scale;
            gap 		= atof(Dialog_Items[0].value)*Scale;
            width       = atof(Dialog_Items[1].value)*Scale;
            length      = atof(Dialog_Items[2].value);
           // xCenter 			= atof(Dialog_Items[5].value)*Scale;
           // yCenter 			= atof(Dialog_Items[6].value)*Scale;
            //pLayer2 			= LLayer_Find (pFile, Dialog_Items[7].value);
            
        }
        else
        {
            return;
        }
        
     
        pLayer = LLayer_Find(pFile,"Dielectric");
        pLayer2 = LLayer_Find (pFile, "Metal");
        int n=0;
        
        
        int dL = (l0 - length)*Scale;
        int outRad = 87.5*Scale + width/2;
        int inRad = 87.5*Scale - width/2;
        for(n = 0; n < 7; n++)
        {
         
            //Center Trace
            centerX = 175*Scale;
            centerY = (262.5+n*350.)*Scale;
           // LTorusParams tParams;
         //   tParams.ptCenter = LPoint_Set(xCenter,262.5);
         //   tParams.nInnerRadius = IR;
         //   tParams.nOuterRadius = OR;
         //   tParams.dStartAngle = startAngle;
         //   tParams.dStopAngle = endAngle;
            //Ground Wire 1
           /*******************************************/
            /******** Write The 1st curve *****************/
            LTorusParams tParams1;
            tParams1.ptCenter = LPoint_Set(centerX,centerY);
            tParams1.nInnerRadius = inRad;
            tParams1.nOuterRadius = outRad;
            tParams1.dStartAngle = 90.;
            tParams1.dStopAngle = 270.;
            
            
            LTorus_CreateNew(pCell, pLayer2, &tParams1);
            /*******************************************/
           
            /*******************************************/
           
            /********* Write the Gap (1/2) ************/
            
            tParams1.ptCenter = LPoint_Set(centerX,centerY);
            tParams1.nInnerRadius = outRad;
            tParams1.nOuterRadius = outRad + gap;
            tParams1.dStartAngle = 90.;
            tParams1.dStopAngle = 270.;
            
             LTorus_CreateNew(pCell, pLayer, &tParams1);
           
            /********************************************/
            /********* Write the Gap (2/2) ************/
           
            tParams1.ptCenter = LPoint_Set(centerX,centerY);
            tParams1.nInnerRadius = inRad;
            tParams1.nOuterRadius = inRad - gap;
            tParams1.dStartAngle = 90.;
            tParams1.dStopAngle = 270.;
          
            
          
            LTorus_CreateNew(pCell, pLayer, &tParams1);
            /********************************************/
            
            
            /********************************************/
            /************* Write the 2nd curve **********/
            
            centerX = 250*Scale;
            centerY = (87.5 + n*350.)*Scale;
            LTorusParams tParams2;
            tParams2.ptCenter = LPoint_Set(centerX,centerY);
            tParams2.nInnerRadius = inRad;
            tParams2.nOuterRadius = outRad;
            tParams2.dStartAngle = 270.;
            tParams2.dStopAngle = 90.;
            
            LTorus_CreateNew(pCell, pLayer2, &tParams2);
            /*********************************************/
            
            /********* Write the Gap (1/2) ************/
            
            tParams1.ptCenter = LPoint_Set(centerX,centerY);
            tParams1.nInnerRadius = outRad;
            tParams1.nOuterRadius = outRad + gap;
            tParams1.dStartAngle = 270;
            tParams1.dStopAngle = 90;
            
            LTorus_CreateNew(pCell, pLayer, &tParams1);
            /************************************************
             /********* Write the Gap (2/2) ************/
            
            tParams1.ptCenter = LPoint_Set(centerX,centerY);
            tParams1.nInnerRadius = inRad;
            tParams1.nOuterRadius = inRad - gap;
            tParams1.dStartAngle = 270;
            tParams1.dStopAngle = 90;
            
            
            
            LTorus_CreateNew(pCell, pLayer, &tParams1);
            /********************************************/
             
        }
        
        leftCornerX1 = -250*Scale;
        rightCornerX1 = 250*Scale;
        leftCornerY1 = -width/2;
        rightCornerY1 = width/2;
        
         LBox_New(pCell, pLayer2, leftCornerX1, leftCornerY1, rightCornerX1, rightCornerY1);
        
        leftCornerX1 = -250*Scale;
        rightCornerX1 = 250*Scale;
        leftCornerY1 = -width/2;
        rightCornerY1 = -width/2 - gap;
        
        LBox_New(pCell, pLayer, leftCornerX1, leftCornerY1, rightCornerX1, rightCornerY1);
        
        leftCornerX1 = -250*Scale;
        rightCornerX1 = 250*Scale;
        leftCornerY1 = width/2;
        rightCornerY1 = width/2 + gap;
        
        LBox_New(pCell, pLayer, leftCornerX1, leftCornerY1, rightCornerX1, rightCornerY1);
        
        leftCornerX1 = -250*Scale-gap;
        rightCornerX1 = -250*Scale;
        leftCornerY1 = -width/2 -gap;
        rightCornerY1 = width/2 + gap;
        
         LBox_New(pCell, pLayer, leftCornerX1, leftCornerY1, rightCornerX1, rightCornerY1);
        /*
        leftCornerX2 = 175*Scale;
        rightCornerX2 = 425*Scale;
        leftCornerY2 = 2437.5*Scale;
        rightCornerY2 = 2462.5*Scale;
        */
        
        leftCornerX2 = 175*Scale;
        rightCornerX2 = 425*Scale;
        leftCornerY2 = (2450*Scale-width/2);
        rightCornerY2 = (2450*Scale + width/2);
        
        LBox_New(pCell, pLayer2, leftCornerX2, leftCornerY2, rightCornerX2, rightCornerY2);
        
        leftCornerX2 = 175*Scale;
        rightCornerX2 = 425*Scale;
        leftCornerY2 = (2450*Scale-width/2) - gap;
        rightCornerY2 = (2450*Scale - width/2);
        
        LBox_New(pCell, pLayer, leftCornerX2, leftCornerY2, rightCornerX2, rightCornerY2);
        
        leftCornerX2 = 175*Scale;
        rightCornerX2 = 425*Scale;
        leftCornerY2 = (2450*Scale+width/2);
        rightCornerY2 = (2450*Scale + width/2) + gap;
        
        LBox_New(pCell, pLayer, leftCornerX2, leftCornerY2, rightCornerX2, rightCornerY2);
        
        
        /*
        leftCornerX3 = 500*Scale;
        rightCornerX3 = 525*Scale;
        leftCornerY3 = 1475*Scale;
        rightCornerY3 = 2362.5*Scale;
        */
        
        leftCornerX3 = 512.5*Scale - width/2;
        rightCornerX3 = 512.5*Scale + width/2;
        leftCornerY3 = 1475*Scale + dL;
        rightCornerY3 = 2362.5*Scale;
        
        LBox_New(pCell, pLayer2, leftCornerX3, leftCornerY3, rightCornerX3, rightCornerY3);
        
        leftCornerX3 = 512.5*Scale - width/2;
        rightCornerX3 = 512.5*Scale - width/2 - gap;
        leftCornerY3 = 1475*Scale + dL;
        rightCornerY3 = 2362.5*Scale;
        
        LBox_New(pCell, pLayer, leftCornerX3, leftCornerY3, rightCornerX3, rightCornerY3);
        
        leftCornerX3 = 512.5*Scale + width/2;
        rightCornerX3 = 512.5*Scale + width/2 + gap;
        leftCornerY3 = 1475*Scale+dL;
        rightCornerY3 = 2362.5*Scale;
        
        LBox_New(pCell, pLayer, leftCornerX3, leftCornerY3, rightCornerX3, rightCornerY3);
        
        leftCornerX3 = 512.5*Scale - width/2 - gap;
        rightCornerX3 = 512.5*Scale + width/2 + gap;
        leftCornerY3 = 1475*Scale+dL;
        rightCornerY3 = 1475*Scale+dL - gap;
        
         LBox_New(pCell, pLayer, leftCornerX3, leftCornerY3, rightCornerX3, rightCornerY3);
        
        centerX = 425*Scale;
        centerY = 2362.5*Scale;
        
        LTorusParams tParams3;
        tParams3.ptCenter = LPoint_Set(centerX,centerY);
        tParams3.nInnerRadius = 87.5*Scale - width/2;
        tParams3.nOuterRadius = 87.5*Scale + width/2;
        tParams3.dStartAngle = 0;
        tParams3.dStopAngle = 90.;
         LTorus_CreateNew(pCell, pLayer2, &tParams3);
        
        tParams3.ptCenter = LPoint_Set(centerX,centerY);
        tParams3.nInnerRadius = inRad - gap;
        tParams3.nOuterRadius = inRad;
        tParams3.dStartAngle = 0;
        tParams3.dStopAngle = 90.;
        LTorus_CreateNew(pCell, pLayer, &tParams3);
        
        tParams3.ptCenter = LPoint_Set(centerX,centerY);
        tParams3.nInnerRadius = outRad;
        tParams3.nOuterRadius = outRad + gap;
        tParams3.dStartAngle = 0;
        tParams3.dStopAngle = 90.;
        LTorus_CreateNew(pCell, pLayer, &tParams3);
        
        n = 0;
        float yBottom = 0;
        float yTop = 0;
        for(n = 0; n<13; n++)
        {
            yBottom = 175*Scale-width/2+n*175*Scale;
            yTop = 175*Scale+width/2+n*175*Scale;
            
            /* Write the Trace */
            leftCornerX = 175*Scale;
            rightCornerX = 250*Scale;
            leftCornerY = yBottom;
            rightCornerY = yTop;
            LBox_New(pCell, pLayer2, leftCornerX, leftCornerY, rightCornerX, rightCornerY);
            
            /* Write the top gap */
            leftCornerY = yTop;
            rightCornerY = yTop+gap;
            LBox_New(pCell, pLayer, leftCornerX, leftCornerY, rightCornerX, rightCornerY);
            
            /* Write the Bottom Gap */
            leftCornerY = yBottom;
            rightCornerY = yBottom - gap;
            LBox_New(pCell, pLayer, leftCornerX, leftCornerY, rightCornerX, rightCornerY);
        }
        /***********************************************************************************************************/
        
        // place in the ground wires
        //float i = truncf(numMeander);
        //int i = (int)numMeander;
       // float remain = numMeander - i;
       // float traceLength = numMeander*meanderLength - remain*meanderLength;
       // remain = cpwLength - traceLength;
        
     //   LDialog_MsgBox(LFormat("Avg Radius = %f \n meander length %f. \n Trace Length %f \n remaining wire = %f",AR, meanderLength,traceLength,remain));
    }
    
    void CPW_register(void)
    {
        LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "fuck", "fuck", NULL /*hotkey category*/);
    }
    
}
CPW_register();


