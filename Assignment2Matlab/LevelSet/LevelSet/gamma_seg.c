/**********************************************************************************/
/* This function implements level set image segmentation by minimizing a data term 
   which measures the conformity of image data within each region to the Gamma model

   Reference:

   I. Ben Ayed, A. Mitiche, and Z. Belhadj, “Multiregion level set partitioning of 
   synthetic aperture radar images,” IEEE Transactions on Pattern Analysis and Machine
   Intelligence, vol. 27, no. 5, pp. 793–800, 2005.
************************************************************************************/

#include "levelset.h"

 /**********************************************************************************
                    Compute the data-driven velocity
************************************************************************************/

float *Seg_Gamma(List_LevelSet *p_lls, Image *p_image, int p, float *F)
{

  int k,l;
  float uu,vv = 2000000000,I,back,back_m;
  float uu_m;
  
  I = p_image->I0[p];

  for(k=0;k<p_lls->Num;k++){
    
    vv = 2000000000;
    uu_m = p_lls->l_levelset[k].p_param[2];
    uu   = log(uu_m)+I/(uu_m);

      for(l=0;l<p_lls->Num;l++){
	if(l==k)
	  continue;
	back_m = p_lls->l_levelset[l].p_param[2];
	back = log(back_m)+I/back_m;
	
	if(back<vv)
	  vv = back;
      }
      
       back_m = p_lls->p_param_out[2];
       back = log(back_m)+I/back_m; 
     
       if(back<vv)
 	vv = back;
      
     F[k] = 1000*(uu-vv);
     
  }

  return(F);
  
}


/***************************************************************************/
/* Calculate the mean and the area of domains inside and outside the curve */
/***************************************************************************/

int Init_Gamma(List_LevelSet *p_lls, Image *p_image)
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

/*****************************************/
/*        Update parametres              */
/*****************************************/
int Update_Gamma(List_LevelSet *p_lls, Image *p_image,
		   int p){

  int k;
  int flag_change = IMG_FALSE;
  int flag_before = IMG_TRUE;
  int flag_after = IMG_TRUE;

  for(k=0;k<p_lls->Num;k++){
    
    if(p_lls->l_levelset[k].U[p]>0.0)
      flag_before = IMG_FALSE;
    if(p_lls->l_levelset[k].Un[p]>0.0)
      flag_after = IMG_FALSE;   
    
    if(p_lls->l_levelset[k].U[p]*p_lls->l_levelset[k].Un[p] < 0.0){
    
      flag_change = IMG_TRUE;
      if(p_lls->l_levelset[k].Un[p]>0.0){
	p_lls->l_levelset[k].p_param[0] += 1.0;
	p_lls->l_levelset[k].p_param[1] += (double)p_image->I0[p];
      }
      else{
	p_lls->l_levelset[k].p_param[0] -= 1.0;
	p_lls->l_levelset[k].p_param[1] -= (double)p_image->I0[p];
      }
      if((int)p_lls->l_levelset[k].p_param[0] == 0){
	p_lls->l_levelset[k].p_param[2] = 0.0;
      }
      else{
	p_lls->l_levelset[k].p_param[2] = 
	  p_lls->l_levelset[k].p_param[1] / p_lls->l_levelset[k].p_param[0];
      }
    }
  }
  
  if(flag_change){
    if(flag_before){
      p_lls->p_param_out[0] -= 1.0;
      p_lls->p_param_out[1] -= (double)p_image->I0[p];
    }
    if(flag_after){
      p_lls->p_param_out[0] += 1.0;
      p_lls->p_param_out[1] += (double)p_image->I0[p];
    }
    if((int)p_lls->p_param_out[0] == 0){
      p_lls->p_param_out[2] = 0.0;
    }
    else{
      p_lls->p_param_out[2] =
	p_lls->p_param_out[1] / p_lls->p_param_out[0];
    }
  }

  return(flag_change);
}
