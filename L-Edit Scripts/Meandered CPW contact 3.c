/*******************************************************************************
 * Macro Name: CPW2
 * Creator  : Matthew Beck
 *
 *		Generates several multiplexed CPW chips of differing trace widths and gaps
 *
 * Revision History:
 * Aug 24, 2015
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
    
    
    
    void buildChip(void) 
    {
    	int gapArray[5][6] = {	{1,2,3,4,5,6},
    							{7,8,9,10,11,12},
    							{13,14,15,16,17,18},
    							{19,20,21,22,23,24},
    							{25,26,27,28,29,30} 
    						};
    					//		{13,14,15,16,17,18},
    					//		{19,20,21,22,23,24},
    					//		{25,26,27,28,29,30} 
    						 							
    	// gapArray[1][6] = {7,8,9,10,11,12};
    	// gapArray[2][6] = {13,14,15,16,17,18};
    	// gapArray[3][6] = {19,20,21,22,23,24};
    	// gapArray[4][6] = {25,26,27,28,29,30};
    	 
    	 
    	int offSet = -2023.25*Scale;	
    	int traceStep = 780*Scale;		
    	int yOffSet = -43.5*Scale;		
    	int lengthArray[6] = {6300,6170,6080,6010,5940,5880}; 
    	int index = 0;							
    	int traceIndex = 0;
    	int chipIndex = 0;							
    	int traceArray[10] = {5,10,15,20,25,30,35,40,45,50};
    	int i = 0;
    	for(traceIndex = 0; traceIndex < 10; traceIndex++) 
    	{
    		LCell	pCell	=	LCell_GetVisible();			
       		LFile	pFile	=	LCell_GetFile(pCell);		
    		//gapIndex = 0;
    		int traceWidth = traceArray[traceIndex]*Scale;		
    		char cellName[] = "5 um chip, 1 - 6 um gaps";
    		char chipLabel[] = "5_um_chip_1_6_um_gaps";	
    		for(chipIndex = 0; chipIndex < 6; chipIndex++)
    		{
       			sprintf(cellName,"%d um chip, %d - %d um gaps",traceArray[traceIndex],gapArray[chipIndex][0],gapArray[chipIndex][5]);
       			
    			//LCell pCell3 = LCell_New(pFile, cellName ); 	
    			//LCell_Delete(pCell3);
    			LCell pCell2 = LCell_New(pFile, cellName );
    			sprintf(chipLabel,"%d_um_chip_%d_%d_um_gaps",traceArray[traceIndex],gapArray[chipIndex][0],gapArray[chipIndex][5]);	
    			
    			for(index= 0; index < 6; index++) 
    			{
    				buildWires(traceWidth,gapArray[chipIndex][index]*Scale,lengthArray[index],offSet + traceStep*index, yOffSet, index); 
    			}
    			
    			for(i = 0; i< 23; i++)
    			{
    				nameChip(chipLabel[i],i);
    			}
    			LWindow pWindow = LWindow_GetVisible();
    			LWindow_Close(pWindow);
    		}
    		
    	}
    }
    
    void buildWires(int width, int gap, int length, int offSet, int yOff, int flip)
    {
        /*The following are L-Edit calls to determine what cell/file to draw in */
       
        LCell	pCell	=	LCell_GetVisible();
        LFile	pFile	=	LCell_GetFile(pCell);
        //LPoint instCenter = LPoint_Set(centerX,centerY);
        //LCoord xCoord = 0;
        //LCoord yCoord = 0;
        struct LMagnification magnify;
        magnify.num = 1;
        magnify.denom = 1;
        LCell baseCell = LCell_Find(pFile,"Chip Outline Mask Version 2");
        LTransform trans = LTransform_Set(0,0,LNormalOrientation,magnify);
        LPoint point = LPoint_Set(1,1);
        LPoint delta = LPoint_Set(0,0);
        LInstance_New(pCell,baseCell,trans,point,delta);
       // LCell pCell = LCell_New(pFile,sprint("%ium chip,
        int l0 = 6300;
      	//int yOffSet = 43.5*Scale;
        //The following are L-Edit calls do define what layer needs to be drawn
        LCoord centerX, centerY, xOffSet, yOffSet, leftCornerX, leftCornerY, rightCornerX, rightCornerY;
        LCoord centerX1, centerY1, leftCornerX1, leftCornerY1, rightCornerX1, rightCornerY1;
        LCoord centerX2, centerY2, leftCornerX2, leftCornerY2, rightCornerX2, rightCornerY2;
        LCoord centerX3, centerY3, leftCornerX3, leftCornerY3, rightCornerX3, rightCornerY3;
        LLayer pLayer;
        LLayer pLayer2;
        LPoint ptArray[20];
        
        
     //   LDialogItem Dialog_Items[3] = {
           // {"Overall CPW length (um)", "21"},
            //{"meander length (um)", "5"},
            //{"Spacing between meanders","50"},
      //      {"Ground Plane Spacing", "2"},
      //      {"Wire Width", "10"},
       //     {"Total Length (<6300)","6300"},
            //{"X Center", "0"},
            //{"Y Center", "0"},
           // {"Layer", "Base"},
      //  };
        
       // float width = 50*Scale;
       // int length = 6300;
       // int gap = 10*Scale;
        
        
        //if ( LDialog_MultiLineInputBox ( "Curve Parameters", Dialog_Items, 3 ) )
       // {
            //cpwLength		= atoi(Dialog_Items[0].value);
            //meanderWidth 	= atof(Dialog_Items[1].value)*Scale;
           // meanderSpace 	= atof(Dialog_Items[2].value)*Scale;
         //   gap 		= atof(Dialog_Items[0].value)*Scale;
          //  width       = atof(Dialog_Items[1].value)*Scale;
           // length      = atof(Dialog_Items[2].value);
           // xCenter 			= atof(Dialog_Items[5].value)*Scale;
           // yCenter 			= atof(Dialog_Items[6].value)*Scale;
            //pLayer2 			= LLayer_Find (pFile, Dialog_Items[7].value);
            
       // }
       // else
       // {
        //    return;
       // }
        
     
        pLayer = LLayer_Find(pFile,"Dielectric");
        pLayer2 = LLayer_Find (pFile, "Metal");
        int n=0;
        
        
        int dL = (l0 - length)*Scale;
        int outRad = 87.5*Scale + width/2;
        int inRad = 87.5*Scale - width/2;
        for(n = 0; n < 7; n++)
        {
         
            //Center Trace
            centerX = 175*Scale + offSet;
            centerY = ((262.5 + n*350.)*Scale + yOff)*pow(-1,flip);
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
            }
            
            /********************************************/
            /************* Write the 2nd curve **********/
    	for(n = 0; n < 7; n++)
    	{
    		int startAng = 0;
    		int stopAng  = 90;
    		if (n==0)
    		{
    			if (pow(-1,flip)<1)
    			{
    				startAng = 270;
    				stopAng = 0;
    			}
    			else
    			{
    			
    				startAng = 0;
    				stopAng = 90;
    			}
    		}
    		else
    		{
    			startAng = 270.;
    			stopAng = 90;
    		}
            centerX = 250*Scale + offSet;
            centerY = ((87.5 + n*350.)*Scale + yOff)*pow(-1,flip);
            
            if(n==0)
            {
            	float radius = (inRad + outRad)/2;
            	leftCornerX1 = centerX + inRad - gap;
            	rightCornerX1 = centerX + inRad + width/2 - 3*Scale;
            	leftCornerY1 = centerY - 2*Scale*pow(-1,flip);
            	rightCornerY1 = centerY;
            	LBox_New(pCell, pLayer, leftCornerX1, leftCornerY1, rightCornerX1, rightCornerY1);
            
            	leftCornerX1 = centerX + outRad - width/2 + 3*Scale;
            	rightCornerX1 = centerX + outRad + gap;
            	leftCornerY1 = centerY - 2*Scale*pow(-1,flip);
            	rightCornerY1 = centerY;
            	LBox_New(pCell, pLayer, leftCornerX1, leftCornerY1, rightCornerX1, rightCornerY1);
            
            }
            LTorusParams tParams2;
            tParams2.ptCenter = LPoint_Set(centerX,centerY);
            tParams2.nInnerRadius = inRad;
            tParams2.nOuterRadius = outRad;
            tParams2.dStartAngle = startAng;
            tParams2.dStopAngle = stopAng;
            
            LTorus_CreateNew(pCell, pLayer2, &tParams2);
            /*********************************************/
            
            /********* Write the Gap (1/2) ************/
            
            tParams2.ptCenter = LPoint_Set(centerX,centerY);
            tParams2.nInnerRadius = outRad;
            tParams2.nOuterRadius = outRad + gap;
            tParams2.dStartAngle = startAng;
            tParams2.dStopAngle = stopAng;
            
            LTorus_CreateNew(pCell, pLayer, &tParams2);
            /************************************************
             /********* Write the Gap (2/2) ************/
            
            tParams2.ptCenter = LPoint_Set(centerX,centerY);
            tParams2.nInnerRadius = inRad;
            tParams2.nOuterRadius = inRad - gap;
            tParams2.dStartAngle = startAng;
            tParams2.dStopAngle = stopAng;
            
            
            
            LTorus_CreateNew(pCell, pLayer, &tParams2);
            /********************************************/
             
        }
        
     
        
        
        leftCornerX2 = 175*Scale + offSet;
        rightCornerX2 = 425*Scale + offSet;
        leftCornerY2 = ((2450*Scale-width/2) + yOff)*pow(-1,flip);
        rightCornerY2 = ((2450*Scale + width/2) + yOff)*pow(-1,flip);
        
        LBox_New(pCell, pLayer2, leftCornerX2, leftCornerY2, rightCornerX2, rightCornerY2);
        
        leftCornerX2 = 175*Scale + offSet;
        rightCornerX2 = 425*Scale + offSet;
        leftCornerY2 = ((2450*Scale-width/2) - gap + yOff)*pow(-1,flip);
        rightCornerY2 = ((2450*Scale - width/2) + yOff)*pow(-1,flip);
        
        LBox_New(pCell, pLayer, leftCornerX2, leftCornerY2, rightCornerX2, rightCornerY2);
        
        leftCornerX2 = 175*Scale + offSet;
        rightCornerX2 = 425*Scale + offSet;
        leftCornerY2 = ((2450*Scale+width/2) + yOff)*pow(-1,flip);
        rightCornerY2 = ((2450*Scale + width/2) + gap + yOff)*pow(-1,flip);
        
        LBox_New(pCell, pLayer, leftCornerX2, leftCornerY2, rightCornerX2, rightCornerY2);
        
        
       
        
        leftCornerX3 = 512.5*Scale - width/2 + offSet;
        rightCornerX3 = 512.5*Scale + width/2 + offSet;
        leftCornerY3 = (850*Scale + dL + yOff)*pow(-1,flip);
        rightCornerY3 = (2362.5*Scale + yOff)*pow(-1,flip);
        
        LBox_New(pCell, pLayer2, leftCornerX3, leftCornerY3, rightCornerX3, rightCornerY3);
        
        leftCornerX3 = 512.5*Scale - width/2 + offSet;
        rightCornerX3 = 512.5*Scale - width/2 - gap + offSet;
        leftCornerY3 = (850*Scale + dL + yOff)*pow(-1,flip);
        rightCornerY3 = (2362.5*Scale + yOff)*pow(-1,flip);
        
        LBox_New(pCell, pLayer, leftCornerX3, leftCornerY3, rightCornerX3, rightCornerY3);
        
        leftCornerX3 = 512.5*Scale + width/2 + offSet;
        rightCornerX3 = 512.5*Scale + width/2 + gap + offSet;
        leftCornerY3 = (850*Scale+dL + yOff)*pow(-1,flip);
        rightCornerY3 = (2362.5*Scale + yOff)*pow(-1,flip);
        
        LBox_New(pCell, pLayer, leftCornerX3, leftCornerY3, rightCornerX3, rightCornerY3);
        
     
        
        centerX = 425*Scale + offSet;
        centerY = (2362.5*Scale + yOff)*pow(-1,flip);
        
        int startAng = 0;
        int stopAng = 0;
        if (pow(-1,flip)<1)
        {
        	startAng = 270;
        	stopAng = 0;
        }
        else
        {
        	startAng = 0;
        	stopAng = 90;
        
        
        }
        LTorusParams tParams3;
        tParams3.ptCenter = LPoint_Set(centerX,centerY);
        tParams3.nInnerRadius = 87.5*Scale - width/2;
        tParams3.nOuterRadius = 87.5*Scale + width/2;
        tParams3.dStartAngle = startAng;
        tParams3.dStopAngle = stopAng;
         LTorus_CreateNew(pCell, pLayer2, &tParams3);
        
        tParams3.ptCenter = LPoint_Set(centerX,centerY);
        tParams3.nInnerRadius = inRad - gap;
        tParams3.nOuterRadius = inRad;
        tParams3.dStartAngle = startAng;
        tParams3.dStopAngle = stopAng;
        LTorus_CreateNew(pCell, pLayer, &tParams3);
        
        tParams3.ptCenter = LPoint_Set(centerX,centerY);
        tParams3.nInnerRadius = outRad;
        tParams3.nOuterRadius = outRad + gap;
        tParams3.dStartAngle = startAng;
        tParams3.dStopAngle = stopAng;
        LTorus_CreateNew(pCell, pLayer, &tParams3);
        
        n = 0;
        float yBottom = 0;
        float yTop = 0;
        for(n = 0; n<13; n++)
        {
            yBottom = 175*Scale-width/2+n*175*Scale;
            yTop = 175*Scale+width/2+n*175*Scale;
            
            /* Write the Trace */
            leftCornerX = 175*Scale + offSet;
            rightCornerX = 250*Scale + offSet;
            leftCornerY = (yBottom + yOff)*pow(-1,flip);
            rightCornerY = (yTop + yOff)*pow(-1,flip);
            LBox_New(pCell, pLayer2, leftCornerX, leftCornerY, rightCornerX, rightCornerY);
            
            /* Write the top gap */
            leftCornerY = (yTop + yOff)*pow(-1,flip);
            rightCornerY = (yTop+gap + yOff)*pow(-1,flip);
            LBox_New(pCell, pLayer, leftCornerX, leftCornerY, rightCornerX, rightCornerY);
            
            /* Write the Bottom Gap */
            leftCornerY = (yBottom + yOff)*pow(-1,flip);
            rightCornerY = (yBottom - gap + yOff)*pow(-1,flip);
            LBox_New(pCell, pLayer, leftCornerX, leftCornerY, rightCornerX, rightCornerY);
        }
       
        
    
    }
    
    void nameChip(char letter, int i)
    {
    
    
    	LCell	pCell	=	LCell_GetVisible();
       	LFile	pFile	=	LCell_GetFile(pCell);
       	
    	char letterCellName[] = "alpha_A";
       
        const char* instancePointer = (const char*)malloc(2);
        
        sprintf(letterCellName,"alpha_%c",toupper(letter));
       
        
        strcpy(instancePointer,letterCellName);
      
        
        LCell letterCell = LCell_Find(pFile,instancePointer);
       	
       	struct LMagnification magnify;
        magnify.num = 4;
        magnify.denom = 1;
        
        LTransform trans = LTransform_Set(-1000*Scale + 75*i*Scale,-2700*Scale,LNormalOrientation,magnify);
        LPoint point = LPoint_Set(1,1);
        LPoint delta = LPoint_Set(0,0);
        LInstance_New(pCell,letterCell,trans,point,delta);
    
    
    
    }
    
    void CPW_register(void)
    {
        LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "CPW Meander Script", "buildChip", NULL /*hotkey category*/);
    }
    
}
CPW_register();


