#define Scale 1000 
#define PI 3.141592653589796
#define cSpeed 3.0e8
#define Scale 1000 

module CuCPW_macro_module
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

void CuCPW(void)
{
	LCell	pCell	=	LCell_GetVisible();
	LFile	pFile	=	LCell_GetFile(pCell);
	
	float startAngle;
	float endAngle;
	float freq;
	float perm;
	float OAL;
	float OAW;
	float cpWidth;
	float gpWidth;
	float lambda;
	float meanderLength;
	float resPoint;
	float arcLength;
	float straight;
	LCoord leftCornerX, leftCornerY, rightCornerX, rightCornerY;
	float endLength = 500;
	float IR, OR;
	float permEff;
	float resWidth = 3000;
	int numMeander1;
	int numMeander2;
	float centralRadius;
	float c = 300000000;
	LLayer pLayer2;
	
	
									
	LDialogItem Dialog_Items[9] = {
										{"Frequency (GHz)", "5"}, 
										{"Relative Permitivity", "9.6"}, 
										{"Overall Length (um)","8000"},
										{"Overall Width (um)", "6000"}, 
										{"Center Pin Width (um)", "20"},
										{"Ground Plane Spacing (um)","100"},
										{"num Meanders B4 straight","5"},
										{"num meander after straight","1"},
										{"Layer", "Copper"},
									};
	
	if ( LDialog_MultiLineInputBox ( "CPW Parameters", Dialog_Items, 9 ) )
	{
		freq				= atof(Dialog_Items[0].value);
		perm			 	= atof(Dialog_Items[1].value);
		OAL			 	= atof(Dialog_Items[2].value);
		OAW		 		= atof(Dialog_Items[3].value);
		cpWidth			= atof(Dialog_Items[4].value);
		gpWidth 			= atof(Dialog_Items[5].value);
		numMeander1		= atof(Dialog_Items[6].value);
		numMeander2 	= atof(Dialog_Items[7].value);
		pLayer2 			= LLayer_Find (pFile, Dialog_Items[8].value);
	}  
	else
	{
		return;
	}
	freq = freq*1000000000;
	centralRadius = (OAL - resWidth - 2*endLength) / ((numMeander1+numMeander2 + 4)*2);
//	centralRadius /= Scale;
	permEff = (1 + perm)/2;
	lambda = c / (freq*sqrt(permEff))*1e6;
	meanderLength = 3*lambda/4;
	resPoint = 2 * meanderLength / 3;
	arcLength = PI*centralRadius;
	straight = (meanderLength - resWidth - 2*endLength -(numMeander1 + numMeander2 + 4)*PI*centralRadius + 4*centralRadius)/(2 + numMeander1 + numMeander2);
//	straight /= Scale;
	IR = centralRadius - cpWidth/2;
	OR = centralRadius + cpWidth/2;
	float xOffSet, yOffSet;
	float newMeanderOffSet = 2*centralRadius;
	xOffSet = cpWidth / 2;
	yOffSet = straight / 2;
	
	//	LDialog_MsgBox(LFormat("IR = %f \n OR = %f",IR,OR));
	//Draw the beginning straight piece;
	LBox_New(pCell, pLayer2, -OAL*Scale/2, -cpWidth*Scale/2, endLength*Scale-OAL*Scale/2, cpWidth*Scale/2);
	//Draw its upturn to the first meander
	LTorusParams tParams;
			tParams.ptCenter = LPoint_Set(endLength*Scale-OAL*Scale/2,centralRadius*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 270;
			tParams.dStopAngle = 360;
	LTorus_CreateNew(pCell, pLayer2, &tParams);
	//Draw lin segment up to first meander
	LBox_New(pCell, pLayer2, (endLength+centralRadius-xOffSet)*Scale-OAL*Scale/2, centralRadius*Scale, (endLength+centralRadius+xOffSet)*Scale-OAL*Scale/2, straight*Scale/2);
	
	//Draw first set of meanders
	int i = 0;
	float xCenter, yCenter,xCenterTorus,yCenterTorus;
	LDialog_MsgBox(LFormat("radius = %f \n straight = %f",centralRadius,straight));
	LDialog_MsgBox(LFormat("freq = %f \n Wave Length = %f \n Meander Length = %f \n straight section = %f",freq, lambda,meanderLength,straight));
	
	for(i = 0; i < numMeander1; i++){
	
			xCenter = endLength + centralRadius + (i+1)*newMeanderOffSet;
	  		yCenter = 0;	
	  		leftCornerX = xCenter - xOffSet;
	  		leftCornerY = yCenter - yOffSet;
	  		rightCornerX = xCenter + xOffSet;
	  		rightCornerY = yCenter + yOffSet;
	  		LBox_New(pCell, pLayer2, leftCornerX*Scale-OAL*Scale/2, leftCornerY*Scale, rightCornerX*Scale-OAL*Scale/2, rightCornerY*Scale);
	  		
//	}
//	float torusStart = xCenter + 2*centralRadius*Scale;
//	for(i = 0; i < numMeander1+1; i++){
//			xCenter = torusStart - (2*i + 1)*centralRadius*Scale;
			xCenterTorus = xCenter - centralRadius;
			if(i%2 ==0){
				yCenterTorus = straight/2;
				startAngle = 0;
				endAngle = 180;
			}
			else{
				yCenterTorus = -straight/2;
				startAngle = 180;
				endAngle = 360;
			}
			tParams.ptCenter = LPoint_Set(xCenterTorus*Scale-OAL*Scale/2,yCenterTorus*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = startAngle;
			tParams.dStopAngle = endAngle;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	}

	//draw the stem coming down from last meander
	//xCenter = xCenter + newMeanderOffSet;
	xCenterTorus = xCenter + centralRadius;
	//yCenterTorus = -yCenterTorus;
	if (yCenterTorus > 0){
	
		tParams.ptCenter = LPoint_Set(xCenterTorus*Scale-OAL*Scale/2,-yCenterTorus*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 180;
			tParams.dStopAngle = 360;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	//Draw connecting Stub
	//xCenter = xCenterTorus+3*centralRadius;
	xCenter = xCenter + newMeanderOffSet;
	LBox_New(pCell,pLayer2,(xCenter - xOffSet)*Scale-OAL*Scale/2,-straight*Scale/2,xCenter*Scale+xOffSet*Scale-OAL*Scale/2,-centralRadius*Scale);
	//Draw 90 degree turnaround to final section
	
		tParams.ptCenter = LPoint_Set((xCenter+centralRadius)*Scale-OAL*Scale/2,-centralRadius*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 90;
			tParams.dStopAngle = 180;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	}
	else{
		tParams.ptCenter = LPoint_Set(xCenterTorus*Scale-OAL*Scale/2,-yCenterTorus*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 0;
			tParams.dStopAngle = 180;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
			xCenter = xCenter + newMeanderOffSet;
		//	xCenter = xCenterTorus+centralRadius;
	LBox_New(pCell, pLayer2, (xCenter - xOffSet)*Scale-OAL*Scale/2, straight*Scale/2, (xCenter+xOffSet)*Scale-OAL*Scale/2, centralRadius*Scale);
	tParams.ptCenter = LPoint_Set((xCenter)*Scale + centralRadius*Scale-OAL*Scale/2,centralRadius*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 180;
			tParams.dStopAngle = 270;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	
	
	}
	//Calculate how far along we are...
	
	float distTrav = endLength+(numMeander1+2)*PI*centralRadius + numMeander1*straight+straight-2*centralRadius;
	
	float diff = resPoint - distTrav;
	LDialog_MsgBox(LFormat("distance = %f \n resPoint = %f \n diff = %f",distTrav,resPoint,diff));
	//Draw resonant section
	LBox_New(pCell,pLayer2,xCenter*Scale + centralRadius*Scale-OAL*Scale/2,-cpWidth*Scale/2,xCenter*Scale + centralRadius*Scale+(diff-50)*Scale-OAL*Scale/2,cpWidth*Scale/2);
	pLayer2 = LLayer_Find (pFile, "Resonant spot");
	LBox_New(pCell,pLayer2,xCenter*Scale + centralRadius*Scale+(diff-50)*Scale-OAL*Scale/2,-cpWidth*Scale/2,xCenter*Scale + centralRadius*Scale+(diff+50)*Scale-OAL*Scale/2,cpWidth*Scale/2);
	pLayer2 = LLayer_Find (pFile, Dialog_Items[8].value);
	LBox_New(pCell,pLayer2,xCenter*Scale + centralRadius*Scale+(diff+50)*Scale-OAL*Scale/2,-cpWidth*Scale/2,xCenter*Scale + centralRadius*Scale+resWidth*Scale-OAL*Scale/2,cpWidth*Scale/2);
	
	//Draw the upturn to the next set of meanders
	xCenter = xCenter + centralRadius+resWidth;
	tParams.ptCenter = LPoint_Set(xCenter*Scale-OAL*Scale/2,centralRadius*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 270;
			tParams.dStopAngle = 360;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	//Draw stem leading to firt meander		
	LBox_New(pCell,pLayer2,xCenter*Scale+centralRadius*Scale-xOffSet*Scale-OAL*Scale/2,centralRadius*Scale,xCenter*Scale+centralRadius*Scale+xOffSet*Scale-OAL*Scale/2,straight*Scale/2);

	//Draw the 2nd set of meanders
	float newCenter = xCenter+centralRadius;
		for(i = 0; i < numMeander2; i++){
	
			xCenter = newCenter + (i+1)*newMeanderOffSet;
	  		yCenter = 0;	
	  		leftCornerX = xCenter - xOffSet;
	  		leftCornerY = yCenter - yOffSet;
	  		rightCornerX = xCenter + xOffSet;
	  		rightCornerY = yCenter + yOffSet;
	  		LBox_New(pCell, pLayer2, leftCornerX*Scale-OAL*Scale/2, leftCornerY*Scale, rightCornerX*Scale-OAL*Scale/2, rightCornerY*Scale);
	  		
//	}
//	float torusStart = xCenter + 2*centralRadius*Scale;
//	for(i = 0; i < numMeander1+1; i++){
//			xCenter = torusStart - (2*i + 1)*centralRadius*Scale;
			xCenterTorus = xCenter - centralRadius;
			if(i%2 ==0){
				yCenterTorus = straight/2;
				startAngle = 0;
				endAngle = 180;
			}
			else{
				yCenterTorus = -straight/2;
				startAngle = 180;
				endAngle = 360;
			}
			tParams.ptCenter = LPoint_Set(xCenterTorus*Scale-OAL*Scale/2,yCenterTorus*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = startAngle;
			tParams.dStopAngle = endAngle;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	}
	
		if (yCenterTorus > 0){
	
		tParams.ptCenter = LPoint_Set(xCenterTorus*Scale+2*centralRadius*Scale-OAL*Scale/2,-yCenterTorus*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 180;
			tParams.dStopAngle = 360;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	//Draw connecting Stub
	//xCenter = xCenterTorus+3*centralRadius;
	xCenter = xCenter + newMeanderOffSet;
	LBox_New(pCell,pLayer2,(xCenter - xOffSet)*Scale-OAL*Scale/2,-straight*Scale/2,xCenter*Scale+xOffSet*Scale-OAL*Scale/2,-centralRadius*Scale);
	//Draw 90 degree turnaround to final section
	
		tParams.ptCenter = LPoint_Set((xCenter+centralRadius)*Scale-OAL*Scale/2,-centralRadius*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 90;
			tParams.dStopAngle = 180;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	}
	else{
		tParams.ptCenter =LPoint_Set(xCenterTorus*Scale+2*centralRadius*Scale-OAL*Scale/2,-yCenterTorus*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 0;
			tParams.dStopAngle = 180;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
			xCenter = xCenter + newMeanderOffSet;
		//	xCenter = xCenterTorus+centralRadius;
	LBox_New(pCell, pLayer2, (xCenter - xOffSet)*Scale-OAL*Scale/2, straight*Scale/2, (xCenter+xOffSet)*Scale-OAL*Scale/2, centralRadius*Scale);
	tParams.ptCenter = LPoint_Set((xCenter)*Scale + centralRadius*Scale-OAL*Scale/2,centralRadius*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 180;
			tParams.dStopAngle = 270;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	
	
	}
	
	
/*	//Draw upTurn to final section
		tParams.ptCenter = LPoint_Set(xCenterTorus*Scale+2*centralRadius*Scale,-yCenterTorus*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 180;
			tParams.dStopAngle = 360;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	//Draw connecting Stub
	xCenter = xCenterTorus+3*centralRadius;
	LBox_New(pCell,pLayer2,xCenter*Scale - xOffSet*Scale,-straight*Scale/2,xCenter*Scale+xOffSet*Scale,-centralRadius*Scale);
	//Draw 90 degree turnaround to final section
	xCenter = xCenter + centralRadius;
		tParams.ptCenter = LPoint_Set(xCenter*Scale,-centralRadius*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 90;
			tParams.dStopAngle = 180;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
			
			*/
	//Draw final Section
	LBox_New(pCell,pLayer2,(xCenter+centralRadius)*Scale-OAL*Scale/2,-cpWidth*Scale/2,(xCenter+centralRadius)*Scale+endLength*Scale-OAL*Scale/2,cpWidth*Scale/2);

}
/*
double drawUpTurn(LCell pCell, LLayer pLayer2,float IR,float OR,float radius, float length, float yCenter, float xCenter,float xOffSet,float Scale)
{
	float xCenter;
	LTorusParams tParams;
	//Draw upTurn 
		tParams.ptCenter = LPoint_Set(xCenter*Scale+2*radius*Scale,-yCenter*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 180;
			tParams.dStopAngle = 360;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	//Draw connecting Stub
	xCenter = xCenterTorus+3*radius;
	LBox_New(pCell,pLayer2,xCenter*Scale - xOffSet*Scale,-length*Scale/2,xCenter*Scale+xOffSet*Scale,-radius*Scale);
	//Draw 90 degree turnaround to final section
	xCenter = xCenter + radius;
		tParams.ptCenter = LPoint_Set(xCenter*Scale,-radius*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 90;
			tParams.dStopAngle = 180;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
}

void drawDownTurn(LCell pCell, LLayer pLayer2, float IR, float OR,float radius, float length, float yCenter, float xCenter,float xOffSet,float Scale){
	
	LTorusParams tParams;
	//Draw Downturn
	tParams.ptCenter = LPoint_Set(xCenter*Scale,yCenter*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 0;
			tParams.dStopAngle = 180;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
	LBox_New(pCell, pLayer2, (xCenter - xOffSet)*Scale, length*Scale/2, (xCenter+xOffSet)*Scale, radius*Scale);
	tParams.ptCenter = LPoint_Set(xCenter*Scale + radius*Scale,radius*Scale);
			tParams.nInnerRadius = IR*Scale;
			tParams.nOuterRadius = OR*Scale;
			tParams.dStartAngle = 180;
			tParams.dStopAngle = 270;
			LTorus_CreateNew(pCell, pLayer2, &tParams);
			
}
*/
void CuCPW_register(void){
	LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "CuCPW", "CuCPW", NULL /*hotkey category*/);
}
	  
}

CuCPW_register();