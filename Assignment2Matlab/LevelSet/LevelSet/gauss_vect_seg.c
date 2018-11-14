/*************************************************************************/
/* This function implements a Gaussian generalization of the classical  
   piecewise constant model in the case of vector-valued images. */
/*************************************************************************/


#include "levelset.h"
#include "levelutil.h"



/**********************************************************************************
                    Compute the data-driven velocity
************************************************************************************/

float *Seg_Vec_Gauss(List_LevelSet *p_lls, Image *p_image, int p, float *F)
{

  int i,j,k,l,N;
  double vv = 1e30;
  double I, Ii;
  double uu, uu_m, uu_v, uu_mi, uu_sum, uu_sumi, uu_d;
  double back, back_m, back_v, back_mi, back_sum, back_sumi, back_d;

  N = p_image->tSize;

  for(k=0;k<p_lls->Num;k++){
    
    vv = 1e30;

    uu_sum = 0.0;
    for(i=0;i<N;i++){
      I = (double)p_image->Il[i][p];
      uu_m = p_lls->l_levelset[k].p_param[N*N+N+1+i];
      uu_sumi = 0.0;
      for(j=0;j<N;j++){
	Ii = (double)p_image->Il[j][p];
	uu_mi = p_lls->l_levelset[k].p_param[N*N+N+1+j];
	uu_v = p_lls->l_levelset[k].p_param[N*N+2*N+1+i*N+j];
	uu_sumi += (Ii-uu_mi)*uu_v;
      }
      uu_sum += (I-uu_m)*uu_sumi;
    }
    uu_d =  p_lls->l_levelset[k].p_param[2*N*N+2*N+1];
    
    uu = (log(uu_d) + uu_sum)/2.0;



    for(l=0;l<p_lls->Num;l++){

      if(l==k)
	continue;
      
      
      back_sum = 0.0;
      for(i=0;i<N;i++){
	I = (double)p_image->Il[i][p];
	back_m = p_lls->l_levelset[l].p_param[N*N+N+1+i];
	back_sumi = 0.0;
	for(j=0;j<N;j++){
	  Ii = (double)p_image->Il[j][p];
	  back_mi = p_lls->l_levelset[l].p_param[N*N+N+1+j];
	  back_v = p_lls->l_levelset[l].p_param[N*N+2*N+1+i*N+j];
	  back_sumi += (Ii-back_mi)*back_v;
	}
	back_sum += (I-back_m)*back_sumi;
      }
      back_d =  p_lls->l_levelset[l].p_param[2*N*N+2*N+1];
      
      back = (log(back_d) + back_sum)/2.0;
    
      
      if(back<vv)
	vv = back;
    }
      
     back_sum = 0.0; 
     for(i=0;i<N;i++){ 
       I = (double)p_image->Il[i][p]; 
       back_m = p_lls->p_param_out[N*N+N+1+i]; 
       back_sumi = 0.0; 
       for(j=0;j<N;j++){ 
 	Ii = (double)p_image->Il[j][p]; 
 	back_mi = p_lls->p_param_out[N*N+N+1+j]; 
 	back_v = p_lls->p_param_out[N*N+2*N+1+i*N+j]; 
 	back_sumi += (Ii-back_mi)*back_v; 
       } 
       back_sum += (I-back_m)*back_sumi; 
     } 
     back_d =  p_lls->p_param_out[2*N*N+2*N+1]; 
  
     back = (log(back_d) + back_sum)/2.0; 
    
     if(back<vv) 
       vv = back; 
    
    F[k] = (uu-vv);
    
  }
  return(F);
  
}


/***************************************************************************/
/* Calculate the mean and the area of domains inside and outside the curve */
/***************************************************************************/

int Init_Vec_Gauss(List_LevelSet *p_lls, Image *p_image)
{
  int Num_Param;    
  int v,h,offset,i,j,l,N;
  int flag;
  double **temp;
  int count = 0;
  

  N = p_image->tSize;
  Num_Param = 2 * (N *(N + 1) + 1);
  if((temp=(double**)malloc(N*sizeof(double*)))==NULL)
    mexErrMsgTxt("Gauss Vect , Memory allocation error");

  for(i=0;i<p_lls->Num;i++){
    p_lls->l_levelset[i].Par_Num = Num_Param;
    p_lls->l_levelset[i].Num_pos = 0;
    if((p_lls->l_levelset[i].p_param = 
	(double*)calloc(Num_Param,sizeof(double)))==NULL)
      return(IMG_FATAL);
  }
  p_lls->Par_Num = Num_Param;
  p_lls->Num_out = 0;
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
	  p_lls->l_levelset[i].Num_pos++;
	  for(j=0;j<N;j++){
	    p_lls->l_levelset[i].p_param[j+1] += (double)p_image->Il[j][offset];
	    for(l=0;l<N;l++){
	      p_lls->l_levelset[i].p_param[j*N+l + (N+1)] += 
		(double)p_image->Il[j][offset]*(double)p_image->Il[l][offset];
	    }
	  }
	  flag = IMG_TRUE;
	}
      }
      if(!flag){
	p_lls->p_param_out[0] += 1.0;
	p_lls->Num_out++;
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
  for(i=0;i<p_lls->Num;i++) {
    /* means */
    for(j=0;j<N;j++){
      p_lls->l_levelset[i].p_param[N*(N+1)+1 + j] = 
	p_lls->l_levelset[i].p_param[j+1]/p_lls->l_levelset[i].p_param[0];
    }
    /* covariance = Int. II^T/A - mean*mean^T */
    for(j=0;j<N;j++){
      for(l=0;l<N;l++){
	p_lls->l_levelset[i].p_param[N*(N+2)+1 + j*N+l] = /*(j==l)?1.0:0.0;*/
	  p_lls->l_levelset[i].p_param[(N+1) + j*N+l]/
	  p_lls->l_levelset[i].p_param[0] -
	  p_lls->l_levelset[i].p_param[N*(N+1)+1 + j]*
	  p_lls->l_levelset[i].p_param[N*(N+1)+1 + l];
      }
    } 

    /* Compute determinant */
    /* Set temporal matrix */
    for(j=0;j<N;j++)
      temp[j] = p_lls->l_levelset[i].p_param + N*N + 2*N + 1 + N*j;
    
    /* Compute inverse covariance */
    p_lls->l_levelset[i].p_param[2*N*N + 2*N + 1] = detinv(temp,N);
    
    
  }

  /* Background */
  /* means */
  for(j=0;j<N;j++){
    p_lls->p_param_out[N*(N+1)+1 + j] = 
      p_lls->p_param_out[j+1]/p_lls->p_param_out[0];
  }
  /* covariance = Int. II^T/A - mean*mean^T */
  for(j=0;j<N;j++){
    for(l=0;l<N;l++){
      p_lls->p_param_out[N*(N+2)+1 + j*N+l] = /*(l==j)?1.0:0.0;*/
      	p_lls->p_param_out[(N+1) + j*N+l]/
	p_lls->p_param_out[0] -
	p_lls->p_param_out[N*(N+1)+1 + j]*
	p_lls->p_param_out[N*(N+1)+1 + l];
    }
  } 
  /* Compute determinant */
  /* Set temporal matrix */
  for(i=0;i<N;i++)
    temp[i] = p_lls->p_param_out + N*N + 2*N + 1 + N*i;

  /* Compute inverse covariance */
  p_lls->p_param_out[2*N*N + 2*N + 1] = detinv(temp,N);;  
  
  
  free(temp);
  return(IMG_NO_ERROR);
  
}

/*****************************************/
/*        Update parametres              */
/*****************************************/
int Update_Vec_Gauss(List_LevelSet *p_lls, Image *p_image,
		   int p){

  int N;
  double **temp;
  int k,i,j,l;
  int flag_change = IMG_FALSE;
  int flag_before = IMG_TRUE;
  int flag_after = IMG_TRUE;

  N = p_image->tSize;
  if((temp=(double**)malloc(N*sizeof(double*)))==NULL)
    mexErrMsgTxt("Gauss Vect , Memory allocation error");


  for(k=0;k<p_lls->Num;k++){
    
    if(p_lls->l_levelset[k].U[p]>0.0)
      flag_before = IMG_FALSE;
    if(p_lls->l_levelset[k].Un[p]>0.0)
      flag_after = IMG_FALSE;   
    
    if(p_lls->l_levelset[k].U[p]*p_lls->l_levelset[k].Un[p] < 0.0){
    
      flag_change = IMG_TRUE;
      if(p_lls->l_levelset[k].Un[p]>0.0){
	p_lls->l_levelset[k].p_param[0] += 1.0;
	p_lls->l_levelset[k].Num_pos++;
	for(j=0;j<N;j++){
	  p_lls->l_levelset[k].p_param[j+1] += (double)p_image->Il[j][p];
	  for(l=0;l<N;l++){
	    p_lls->l_levelset[k].p_param[j*N+l + (N+1)] += 
	      (double)p_image->Il[j][p]*(double)p_image->Il[l][p];
	  }
	}
      }
      else{
	p_lls->l_levelset[k].p_param[0] -= 1.0;
	p_lls->l_levelset[k].Num_pos--;
	for(j=0;j<N;j++){
	  p_lls->l_levelset[k].p_param[j+1] -= (double)p_image->Il[j][p];
	  for(l=0;l<N;l++){
	    p_lls->l_levelset[k].p_param[j*N+l + (N+1)] -= 
	      (double)p_image->Il[j][p]*(double)p_image->Il[l][p];
	  }
	}
      }
      /* means */
      for(j=0;j<N;j++){
	if(p_lls->l_levelset[k].p_param[0]<EPSILON)
	  p_lls->l_levelset[k].p_param[N*(N+1)+1 + j] = 0.0;
	else
	  p_lls->l_levelset[k].p_param[N*(N+1)+1 + j] = 
	    p_lls->l_levelset[k].p_param[j+1]/p_lls->l_levelset[k].p_param[0];
      }
      /* covariance = Int. II^T/A - mean*mean^T */
      for(j=0;j<N;j++){
	for(l=0;l<N;l++){
	  if(p_lls->l_levelset[k].p_param[0]<EPSILON)
	    p_lls->l_levelset[k].p_param[N*(N+2)+1 + j*N+l] = 0.0;
	  else
	    p_lls->l_levelset[k].p_param[N*(N+2)+1 + j*N+l] = 
	      p_lls->l_levelset[k].p_param[(N+1) + j*N+l]/
	      p_lls->l_levelset[k].p_param[0] -
	      p_lls->l_levelset[k].p_param[N*(N+1)+1 + j]*
	      p_lls->l_levelset[k].p_param[N*(N+1)+1 + l];
	}
      } 
      /* Compute determinant */
      /* Set temporal matrix */
      for(i=0;i<N;i++)
	temp[i] = p_lls->l_levelset[k].p_param + N*N + 2*N + 1 + N*i;
      
      /* Compute inverse covariance */
      if(p_lls->l_levelset[k].p_param[0]<EPSILON)
	p_lls->l_levelset[k].p_param[2*N*N + 2*N + 1] = 1.0;
      else
	p_lls->l_levelset[k].p_param[2*N*N + 2*N + 1] = detinv(temp,N);
    }
  }
  
  if(flag_change){
    if(flag_before){
      p_lls->p_param_out[0] -= 1.0;
      p_lls->Num_out--;
      for(j=0;j<N;j++){
	p_lls->p_param_out[j+1] -= (double)p_image->Il[j][p];
	for(l=0;l<N;l++){
	  p_lls->p_param_out[j*N+l + (N+1)] -= 
	    (double)p_image->Il[j][p]*(double)p_image->Il[l][p];
	}
      }
    }
    if(flag_after){
      p_lls->p_param_out[0] += 1.0;
      p_lls->Num_out++;
      for(j=0;j<N;j++){
	p_lls->p_param_out[j+1] += (double)p_image->Il[j][p];
	for(l=0;l<N;l++){
	  p_lls->p_param_out[j*N+l + (N+1)] += 
	    (double)p_image->Il[j][p]*(double)p_image->Il[l][p];
	}
      }
    }

    /* Background */
    /* means */
    for(j=0;j<N;j++){
      p_lls->p_param_out[N*(N+1)+1 + j] = 
	p_lls->p_param_out[j+1]/p_lls->p_param_out[0];
    }
    /* covariance = Int. II^T/A - mean*mean^T */
    for(j=0;j<N;j++){
      for(l=0;l<N;l++){
	p_lls->p_param_out[N*(N+2)+1 + j*N+l] = 
	  p_lls->p_param_out[(N+1) + j*N+l]/
	  p_lls->p_param_out[0] -
	  p_lls->p_param_out[N*(N+1)+1 + j]*
	  p_lls->p_param_out[N*(N+1)+1 + l];
      }
    } 
    /* Compute determinant */
    /* Set temporal matrix */
    for(i=0;i<N;i++)
      temp[i] = p_lls->p_param_out + N*N + 2*N + 1 + N*i;
    
    /* Compute inverse covariance */
    p_lls->p_param_out[2*N*N + 2*N + 1] = detinv(temp,N);

    
  }
  
  free(temp);
  return(flag_change);
}

