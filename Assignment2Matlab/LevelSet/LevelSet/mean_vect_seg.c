
/****************************************************************** 
This function implements level set image segmentation by minimizing 
the classical Chan-Vese data term in the case of vector-valued images. 
*******************************************************************/

#include "levelset.h"



/**********************************************************************************
                    Compute the data-driven velocity
************************************************************************************/

float *Seg_Vec_Mean(List_LevelSet *p_lls, Image *p_image, int p, float *F)
{

  int k,l;
  float uu,vv = 2000000000,I,back;
  
  I = p_image->I0[p];

  for(k=0;k<p_lls->Num;k++){
    
    vv = 2000000000;
    uu = p_lls->l_levelset[k].p_param[2];
    
      for(l=0;l<p_lls->Num;l++){
	back = p_lls->l_levelset[l].p_param[2];
	if(l==k)
	  continue;
	if(fabs(I-back)<fabs(I-vv))
	  vv = back;
      }
      
      
      back = p_lls->p_param_out[2];
      if(fabs(I-back)<fabs(I-vv))
	vv = back;      
      
      F[k] = -2*(uu-vv)*(I - (uu + vv)/2);
    
  }
  return(F);
  
}


/***************************************************************************/
/* Calculate the mean and the area of domains inside and outside the curve */
/***************************************************************************/

int Init_Vec_Mean(List_LevelSet *p_lls, Image *p_image)
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

int Update_Vec_Mean(List_LevelSet *p_lls, Image *p_image,
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
	p_lls->l_levelset[k].p_param[1] += p_image->I0[p];
      }
      else{
	p_lls->l_levelset[k].p_param[0] -= 1.0;
	p_lls->l_levelset[k].p_param[1] -= p_image->I0[p];	
      }
      if(p_lls->l_levelset[k].p_param[0]==0)
	p_lls->l_levelset[k].p_param[2] = 0;
      else
	p_lls->l_levelset[k].p_param[2] =
	  p_lls->l_levelset[k].p_param[1]/p_lls->l_levelset[k].p_param[0];
    }
  }
  
  if(flag_change){
    if(flag_before){
      p_lls->p_param_out[0] -= 1.0;
      p_lls->p_param_out[1] -= p_image->I0[p];
    }
    if(flag_after){
      p_lls->p_param_out[0] += 1.0;
      p_lls->p_param_out[1] += p_image->I0[p];
    }
    if(p_lls->p_param_out[0]==0)
      p_lls->p_param_out[2] = 0;
    else
      p_lls->p_param_out[2] = p_lls->p_param_out[1]/p_lls->p_param_out[0];
  }

  return(flag_change);
  
}
