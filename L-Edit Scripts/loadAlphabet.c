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
    
    void fuckOff(void)
    {
    	LCell	pCell	=	LCell_GetVisible();
       	LFile	pFile	=	LCell_GetFile(pCell);
       	//char *a = malloc(256);
       	//char *stringa;
       	//sprintf(stringa,"myCell %d",15);
		//strcpy(a, stringa);
       	char stringa[] = "Help 1";
       	sprintf(stringa,"Help %d",15);
       //char *cellName;
       //sprintf(cellName,"%d um chip,1 - 6 um gaps",15);
       //	printf(cellName);
       	LCell 	pCell2 = LCell_New(pFile, stringa );
       	char F[] = "F";
       	write(F[0]);
       	//LOGO_Main("Fuck");
       	//LFile_OpenCell(pFile,a);
       	//free(a);
       //	LCell 	pCell2 = LCell_New(pFile, "15 um chip, 1 - 6 um gaps" );
    
    
    }
    
    void write(char letter)
    {
    	LCell	pCell	=	LCell_GetVisible();
       	LFile	pFile	=	LCell_GetFile(pCell);
       	LFile 	pFile2	=	LFile_Open("C:\Program Files (x86)\Tanner EDA\Tanner Tools v13.0\AddIns\alphabet.tdb",LTdbFile);
    	struct LMagnification magnify;
        magnify.num = 4;
        magnify.denom = 1;
        char alphabet[] = "_0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        int i = 0;
        //char letter[0] = "";
        char sourceCellName[] = "_A";
        char destCellName[] = "alpha_A";
        for(i = 0; i < 37; i++)
        {
        	//letter = alphabet[i];
        	sprintf(sourceCellName,"_%c",toupper(alphabet[i]));
        	sprintf(destCellName,"alpha_%c",toupper(alphabet[i]));
        	LCell sourceCell = LCell_Find(pFile2,sourceCellName);
        	LCell_Copy(pFile2,sourceCell,pFile,destCellName);
        
        }
        
        char letterCellName[] = "alpha_A";
       
        const char* instancePointer = (const char*)malloc(2);
        
        sprintf(letterCellName,"alpha_%c",toupper(letter));
       
        
        strcpy(instancePointer,letterCellName);
      
        
        LCell letterCell = LCell_Find(pFile,instancePointer);
       
        LTransform trans = LTransform_Set(0,0,LNormalOrientation,magnify);
        LPoint point = LPoint_Set(1,1);
        LPoint delta = LPoint_Set(0,0);
        LInstance_New(pCell,letterCell,trans,point,delta);
    }
    
    
    
    
    
    void CPW_register(void)
    {
        LMacro_BindToMenuAndHotKey_v9_30("Tools", NULL /*hotkey*/, "loadAlphabet", "fuckOff", NULL /*hotkey category*/);
    }
    
}
CPW_register();