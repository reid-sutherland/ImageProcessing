#include "mex.h"
#include "matrix.h"
#include "levelset.h"
#include <string.h>


/******************************  image size  ******************************/
int img_size (char* fichier, int ordre)
{
int size;
double* out;

mxArray *input[1];
mxArray *output[1];

   input[0] = mxCreateString( fichier );

   if( mexCallMATLAB( 1, output, 1, input, "imgsize" ) )  
       	mexErrMsgTxt("Unable to perform imread");

   mxDestroyArray(input[0]);

   out=mxGetPr(output[0]);
   size= (int)out[ordre];

	return (size);
}
/******************************************************************************************/
/************************ Energy computation  ***************************/
double cal_energy(List_LevelSet *p_lls, Image *p_image, unsigned char vect)
{
int v, h, k, i, v_limit, h_limit, offset;
double Energy;
float* position;

h_limit = p_image->xSize;
v_limit = p_image->ySize;
Energy = 0;

if (vect==IMG_FALSE)
{
for(v=0;v<v_limit;v++) 
{
    for(h=0;h<h_limit;h++)
    {
      offset = v*h_limit+h;	
	k=0;
	i=0;
	while (k==0)
	{
	  position = p_lls->l_levelset[i].U;
	  if (position[offset] >=0)
	  {
	    Energy += (p_image->I0[offset] - p_lls->l_levelset[i].p_param[2]) * (p_image->I0[offset] - p_lls->l_levelset[i].p_param[2]);
	    k=1;
	  }
	    else
	    {
	      if (i==p_lls->Num-1)
	      {
		Energy += (p_image->I0[offset] - p_lls->p_param_out[2]) * (p_image->I0[offset] - p_lls->p_param_out[2]); 
		k=1;
	      }
		else
		  i=i+1;
	      
	    }
	}
     }
}
}
return Energy;
 }  
/**************************************************************************************/
/************************  The mean of each region  ***************************/
double cal_mean(List_LevelSet *p_lls, Image *p_image, int offset, unsigned char vect)
{
int k, i;
double mean;
float* position;
if (vect==IMG_FALSE)
{
	k=0;
	i=0;
	while (k==0)
	{
	  position = p_lls->l_levelset[i].U;
	  if (position[offset] >=0)
	  {
	    mean = p_lls->l_levelset[i].p_param[2];
	    k=1;
	  }
	    else
	    {
	      if (i==p_lls->Num-1)
	      {
		mean = p_lls->p_param_out[2];
		k=1;
	      }
		else
		  i=i+1;
	      
	    }
	}
}
return mean;
}
/**************************************************************************************/




