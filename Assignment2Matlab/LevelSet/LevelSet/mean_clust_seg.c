
   

/* This function evolves multiple curves so that the curves define a partition.
The method embeds directly a simple partition constraint in the evolution equations. 
Starting from an arbitrary initial partition, the constraints implements the rule 
that if a point leaves a region, it is claimed by a single other region. The scheme 
is fast, and results in a significant reduction in the computational load 
(the complexity grows linearly as a function of the number of regions). 
The method is used in conjunction with the piecewise constant model.

References:

I. Ben Ayed, A. Mitiche, and Z. Belhadj, Polarimetric image segmentation via maximum 
likelihood approximation and efficient multiphase level sets,” IEEE Transactions on 
Pattern Analysis and Machine Intelligence, vol. 28, no. 9, pp. 1493–1500, 2006. 

I. Ben Ayed, A. Mitiche: A Partition Constrained Minimization Scheme for Efficient 
Multiphase Level Set Image Segmentation. IEEE ICIP: 1641-1644, 2006

I. Ben Ayed and A. Mitiche, “A region merging prior for variational level set image 
segmentation,” IEEE  Transactions on Image Processing, vol. 17, no. 12,   
pp. 2301–2313, 2008. */ 


#include "levelset.h"

/**********************************************************************************
                    Compute the data-driven velocity
************************************************************************************/

float *Seg_Clust_Mean(List_LevelSet *p_lls, Image *p_image, int p, float *F)
{

  int k, l, h_limit;
  double uu,vv, back; 
  double I;
  int uu_nb, vv_nb;

h_limit=p_image->xSize;
I = p_image->I0[p];

uu_nb=p_lls->Num;
vv_nb=p_lls->Num;

for(k=0;k<p_lls->Num;k++){
if (p_lls->l_levelset[k].U[p] > 0)
{ uu_nb=k;
  uu = p_lls->l_levelset[k].p_param[2];
 break;
}}

if (uu_nb==p_lls->Num)
{
uu = p_lls->p_param_out[2];
vv = 2000000000;
for(l=0;l<p_lls->Num;l++){
back = p_lls->l_levelset[l].p_param[2];
if(fabs(I-back)<fabs(I-vv))
{ vv = back;
  vv_nb=l;}
}
}
else
{
vv=p_lls->p_param_out[2];
for(l=0;l<p_lls->Num;l++){
back = p_lls->l_levelset[l].p_param[2];
if(fabs(I-back)<fabs(I-vv))
{ vv = back;
  vv_nb=l;}
}}
F[0]=uu_nb;
F[1]=vv_nb;
F[2]=(I-uu)*(I-uu)-(I-vv)*(I-vv); 
F[3]=-(I-uu)*(I-uu)+(I-vv)*(I-vv);
return(F);
}
/*************************************************************************************/
/* Calculate the mean and the area of domains inside the curves and in the background */
/**************************************************************************************/


int Init_Clust_Mean(List_LevelSet *p_lls, Image *p_image)
{
  const int Num_Param = 3;
  int v,h,offset,i;
  int flag;
  
  for(i=0;i<p_lls->Num;i++){
    p_lls->l_levelset[i].Par_Num = Num_Param;
    if((p_lls->l_levelset[i].p_param = 
	(double*)calloc(Num_Param,sizeof(double)))==NULL)
      return(IMG_FATAL);
  }
  p_lls->Par_Num = Num_Param;
  if((p_lls->p_param_out = (double*)calloc(Num_Param,sizeof(double)))==NULL)
      return(IMG_FATAL);

  flag = IMG_FALSE;
  for(v=0;v<p_image->ySize;v++) {
    for(h=0;h<p_image->xSize;h++) {
      offset = v*p_image->xSize+h;
      flag = IMG_FALSE;
      for(i=0;i<p_lls->Num;i++){
	if (p_lls->l_levelset[i].U[offset] > 0.0) {
	  p_lls->l_levelset[i].p_param[0] += 1.0;
	  p_lls->l_levelset[i].p_param[1] += (double)p_image->I0[offset];
	 
	  flag = IMG_TRUE;
	}
      }
      if(!flag){
	p_lls->p_param_out[0] += 1.0;
	p_lls->p_param_out[1] += (double)p_image->I0[offset];

      }      
    }
  }
  for(i=0;i<p_lls->Num;i++){
    p_lls->l_levelset[i].p_param[2] = 
      p_lls->l_levelset[i].p_param[1]/p_lls->l_levelset[i].p_param[0];
  
  }
  p_lls->p_param_out[2] = p_lls->p_param_out[1]/p_lls->p_param_out[0];

  return(IMG_NO_ERROR);
  
}
