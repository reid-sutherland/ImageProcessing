
/**********************************************************************************/
/* This function implements multiphase level set image segmentation by enforcing an
   explicit partition constraint at all times during curve evolution, and uses a 
   Gaussian image model in the case of vector-valued images. 

   Reference:

    A. Mansouri, A. Mitiche, and C. Vazquez, “Multiregion competition: A level set 
	extension of region competition to multiple region partioning,” Computer Vision 
	and Image Understanding, vol. 101, no. 3, pp. 137–150, 2006.


/**********************************************************************************/


#include "levelset.h"
#include "levelutil.h"
#define C 0.0


/**********************************************************************************
                    Compute the data-driven velocity
************************************************************************************/

float *Seg_Vec_Part_Gauss(List_LevelSet *p_lls, Image *p_image, int p, float *F)
{

  int i,k,l,N;
  double vv;
  double *I;
  double uu;
  
  N = p_image->tSize;
  if((I = (double*)malloc(N*sizeof(double)))==NULL)
    mexErrMsgTxt("Gauss_Vect_Part , Memory allocation error");
  for(i=0;i<N;i++)
    I[i] = (double)p_image->Il[i][p];

  for(k=0;k<p_lls->Num;k++){
    
    vv = cal_uu(I,p_lls->p_param_out,N);
    uu = cal_uu(I,p_lls->l_levelset[k].p_param,N);
    F[k] = (uu-vv);
    
    for(l=0;l<p_lls->Num;l++){
      if(/*(p_lls->l_levelset[l].U[p-p_image->xSize-1]>0.0)||
	   (p_lls->l_levelset[l].U[p-p_image->xSize]>0.0)||
	   (p_lls->l_levelset[l].U[p-p_image->xSize+1]>0.0)||*/
	 (p_lls->l_levelset[l].U[p]>0.0)/*||
	   (p_lls->l_levelset[l].U[p-1]>0.0)||
	   (p_lls->l_levelset[l].U[p+1]>0.0)||
	   (p_lls->l_levelset[l].U[p+p_image->xSize-1]>0.0)||
	   (p_lls->l_levelset[l].U[p+p_image->xSize]>0.0)||
	   (p_lls->l_levelset[l].U[p+p_image->xSize+1]>0.0)*/){
	if(l<k){
	  /*	  F[k] = C;
		  break;*/
	  continue;
	}
	if(l==k)
	  continue;
	vv = cal_uu(I,p_lls->l_levelset[l].p_param,N);
	F[k] = (uu-vv);
	break;
      }
    }
  }
  free(I);
  return(F);
}


/***************************************************************************/
/* Calculate the mean and the area of domains inside and outside the curve */
/***************************************************************************/
int Init_Vec_Part_Gauss(List_LevelSet *p_lls, Image *p_image)
{
  int Num_Param;    /* Area, Int. of I, Mean, Int. of I^2, Var */
  int v,h,offset,i,j,l,N;
  int flag;
  
  

  N = p_image->tSize;
  Num_Param = 2 * (N *(N + 1) + 1);

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
      for(i=0;i<p_lls->Num;i++) {
	if (p_lls->l_levelset[i].U[offset] > 0.0) {
	  p_lls->l_levelset[i].p_param[0] += 1.0;
	  for(j=0;j<N;j++){
	    p_lls->l_levelset[i].p_param[j+1] += (double)p_image->Il[j][offset];
	    for(l=0;l<N;l++){
	      p_lls->l_levelset[i].p_param[j*N+l + (N+1)] += 
		(double)p_image->Il[j][offset]*(double)p_image->Il[l][offset];
	    }
	  }
	  flag = IMG_TRUE;
	  break;
	}
      }
      if(!flag){
	p_lls->p_param_out[0] += 1.0;
	for(j=0;j<N;j++){
	  p_lls->p_param_out[j+1] += (double)p_image->Il[j][offset];
	  for(l=0;l<N;l++){
	    p_lls->p_param_out[j*N+l + (N+1)] += 
	      (double)p_image->Il[j][offset]*(double)p_image->Il[l][offset];
	  }
	}
      }      
    }
  }

  /* Compute means and covariance matrix */
  for(i=0;i<p_lls->Num;i++) {
    update_param(p_lls->l_levelset[i].p_param,N);
  }
  /* Background */
  update_param(p_lls->p_param_out,N);
  
  return(IMG_NO_ERROR);
}

/*****************************************/
/*        Update parametres              */
/*****************************************/
int Update_Vec_Part_Gauss(List_LevelSet *p_lls, Image *p_image,
		   int p){

  int N;
  int k,i,l;
  double *I;
  unsigned char flag_change = IMG_FALSE;
  unsigned char flag_finish = IMG_FALSE;

  N = p_image->tSize;
  if((I = (double*)malloc(N*sizeof(double)))==NULL)
    mexErrMsgTxt("Gauss_Vect_Part , Memory allocation error");
  for(i=0;i<N;i++)
    I[i] = (double)p_image->Il[i][p];


  /* Initialize to the first level-set */
  k = 0;
  while((!flag_finish)&&(k<p_lls->Num)){ /* Not finished yet */
    
    /* Is there a change in sign? */
    if(p_lls->l_levelset[k].U[p]*p_lls->l_levelset[k].Un[p] < 0.0){
      
      flag_change = IMG_TRUE;
      if(p_lls->l_levelset[k].Un[p]>0.0){

	/* The point enter the level-set, so it must be added and mean
	 recalculated */
	add_point(I,p_lls->l_levelset[k].p_param,N,1);
	update_param(p_lls->l_levelset[k].p_param,N);   

	/* Search from where the point is entering starting from the
	   next level-set */
	l = k+1;
	while((!flag_finish)&&(l<p_lls->Num)){
	  if(p_lls->l_levelset[l].U[p]>0.0){

	    /* Was inside level-set l? Must be substracted and updated
	     and all is done */
	    add_point(I,p_lls->l_levelset[l].p_param,N,-1);
	    update_param(p_lls->l_levelset[l].p_param,N);   
	    flag_finish = IMG_TRUE;
	  }
	  l++;
	}
	if(!flag_finish){

	  /* It wasn't inside any level-set, so it was in the
	     background, must be substracted and updated */
	  add_point(I,p_lls->p_param_out,N,-1);
	  update_param(p_lls->p_param_out,N);   
	  flag_finish = IMG_TRUE;
	}
      }
      else{
	/* The same but the point is leaving the level-set */
	add_point(I,p_lls->l_levelset[k].p_param,N,-1);
	update_param(p_lls->l_levelset[k].p_param,N);   
	l = k+1;
	while((!flag_finish)&&(l<p_lls->Num)){
	  if(p_lls->l_levelset[l].Un[p]>0.0){
	    add_point(I,p_lls->l_levelset[l].p_param,N,1);
	    update_param(p_lls->l_levelset[l].p_param,N);   
	    flag_finish = IMG_TRUE;
	  }
	  l++;
	}
	if(!flag_finish){
	  add_point(I,p_lls->p_param_out,N,1);
	  update_param(p_lls->p_param_out,N);
	  flag_finish = IMG_TRUE;
	}
      }
    }
    else{
      /* there was no change but the point was and still is inside the
	 level-set, there is nothing to do */
      if((p_lls->l_levelset[k].U[p] > 0.0)&&(p_lls->l_levelset[k].Un[p] > 0.0)){
	break;
      }
    }
    k++;
  }
  
  free(I);
  return(flag_change);
}



