



  

/**********************DESCRIPTION: 

Basic level set evolution routine which implements evolution equations in the form 
of the sum of two velocities, one is a function of the curvature of the evolving curve
and the other is a scalar function which depends on position and time. Formal and 
complete details on this implementation can be found in the book "Variational and 
Level Set Methods in Image Segmentation", by Mitiche and Ben Ayed 2010, Springer, 
1st edition, pages 22-25 and pages 50-51.

AUTHORS:
   Abdol-Reza Mansouri
   Carlos Vazquez
   Ismail Ben Ayed ************************/


#include "levelset.h"


int LS_LevelSetPDE (List_LevelSet *p_lls,Image *p_image, Data data)
{
  int i, j, upd;
  int k_pos, l_max;
  int p ;
  int k ;
  float *F;
  Nabla nabla;
  float curvature_term;  
  double delta,sum1,sum2,error;
  unsigned char change = IMG_FALSE;
  int flag_pos_change;
  int flag_max_change;

  /* The case of multiregion (multiphase) segmentation embeds directly a simple
     partition constraint in the equations. Starting from an arbitrary initial  
	 partition, the constraints implements the rule that if a point leaves a 
	 region, it is claimed by a single other region. */

if (data.Seg_Fun==Seg_Clust_Mean){
  F = (float*)malloc(sizeof(float)*4);
  if(F==NULL)
    return(IMG_FATAL);

  
  for (i = 1; i < p_image->ySize-1; i++) {
    for (j = 1; j < p_image->xSize-1; j++) {
      p = j + i * p_image->xSize ;
      flag_pos_change=IMG_FALSE;
      flag_max_change=IMG_FALSE;
      F = (*data.Seg_Fun)(p_lls,p_image,p,F); 

k_pos = F[0];
l_max = F[1]; 
if (k_pos==l_max)
continue;
if ((k_pos<p_lls->Num)) {
nabla = Nabla_cal(p_lls->l_levelset+k_pos,p);
curvature_term = Curv_cal(p_lls->l_levelset+k_pos,p);
delta = data.dt * ( -(imax(F[2], 0) * nabla.p  +
			      imin(F[2], 0) * nabla.m))+ data.K*curvature_term;
delta = (delta > 2 * MAX_VAL)? 2 * MAX_VAL : 
  (delta < 2*MIN_VAL )? 2* MIN_VAL : delta;
p_lls->l_levelset[k_pos].Un[p] =  p_lls->l_levelset[k_pos].U[p] +
	  delta;

 	p_lls->l_levelset[k_pos].Un[p] =  
 	  (p_lls->l_levelset[k_pos].Un[p] > 2000000000*MAX_VAL) ? 2000000000*MAX_VAL : 
 	  ((p_lls->l_levelset[k_pos].Un[p] < 2000000000*MIN_VAL) ? 2000000000*MIN_VAL :  
	   p_lls->l_levelset[k_pos].Un[p]); 
if(p_lls->l_levelset[k_pos].U[p]*p_lls->l_levelset[k_pos].Un[p] < 0.0) {
flag_pos_change=IMG_TRUE;
p_lls->l_levelset[k_pos].p_param[0] -= 1.0;
p_lls->l_levelset[k_pos].p_param[1] -= (double)p_image->I0[p];
p_lls->l_levelset[k_pos].p_param[2] = 
	  p_lls->l_levelset[k_pos].p_param[1] / p_lls->l_levelset[k_pos].p_param[0];
	}

p_lls->l_levelset[k_pos].U[p] = p_lls->l_levelset[k_pos].Un[p] ;
}


if ((l_max<p_lls->Num)) {
nabla = Nabla_cal(p_lls->l_levelset+l_max,p);
curvature_term = Curv_cal(p_lls->l_levelset+l_max,p);
delta = data.dt * ( -(imax(F[3], 0) * nabla.p  +
			      imin(F[3], 0) * nabla.m)) +data.K*curvature_term;

delta = (delta > 2 * MAX_VAL)? 2 * MAX_VAL : 
  (delta < 2*MIN_VAL )? 2* MIN_VAL : delta;
p_lls->l_levelset[l_max].Un[p] =  p_lls->l_levelset[l_max].U[p] +
	  delta;

 	p_lls->l_levelset[l_max].Un[p] =  
 	  (p_lls->l_levelset[l_max].Un[p] > 2000000000*MAX_VAL) ? 2000000000*MAX_VAL : 
 	  ((p_lls->l_levelset[l_max].Un[p] < 2000000000*MIN_VAL) ? 2000000000*MIN_VAL :  
	   p_lls->l_levelset[l_max].Un[p]);
if(p_lls->l_levelset[l_max].U[p]*p_lls->l_levelset[l_max].Un[p] < 0.0) {
flag_max_change=IMG_TRUE;

p_lls->l_levelset[l_max].p_param[0] += 1.0;
p_lls->l_levelset[l_max].p_param[1] += (double)p_image->I0[p];
p_lls->l_levelset[l_max].p_param[2] = 
	  p_lls->l_levelset[l_max].p_param[1] / p_lls->l_levelset[l_max].p_param[0];
 }

p_lls->l_levelset[l_max].U[p] = p_lls->l_levelset[l_max].Un[p] ;
 }

if ((flag_pos_change==IMG_TRUE) & (l_max==p_lls->Num))
 {
p_lls->p_param_out[0] += 1.0;
p_lls->p_param_out[1] += (double)p_image->I0[p];
 p_lls->p_param_out[2] =
	p_lls->p_param_out[1] / p_lls->p_param_out[0];
}   

if ((flag_pos_change==IMG_TRUE) & (flag_max_change==IMG_FALSE) & (l_max<p_lls->Num)) 
 {
p_lls->p_param_out[0] += 1.0;
p_lls->p_param_out[1] += (double)p_image->I0[p];
 p_lls->p_param_out[2] =
	p_lls->p_param_out[1] / p_lls->p_param_out[0];
}  

if ((k_pos==p_lls->Num)&(flag_max_change==IMG_TRUE))
{
p_lls->p_param_out[0] -= 1.0;
p_lls->p_param_out[1] -= (double)p_image->I0[p];
 p_lls->p_param_out[2] =
	p_lls->p_param_out[1] / p_lls->p_param_out[0];
}

 
    }
    }
    
}

/***************** The general case *********************************/

else {
  
  F = (float*)malloc(sizeof(float)*p_lls->Num);
  if(F==NULL)
    return(IMG_FATAL);

  
  for (i = 1; i < p_image->ySize-1; i++) {
    for (j = 1; j < p_image->xSize-1; j++) {
      p = j + i * p_image->xSize ;
      
      F = (*data.Seg_Fun)(p_lls,p_image,p,F); 
      
      for(k=0;k<p_lls->Num;k++){
	nabla = Nabla_cal(p_lls->l_levelset+k,p);
	curvature_term = Curv_cal(p_lls->l_levelset+k,p);
	delta = data.dt * ( -(imax(F[k], 0) * nabla.p  +
			      imin(F[k], 0) * nabla.m))+ data.K  * curvature_term;
	
	delta = (delta > 2 * MAX_VAL)? 2 * MAX_VAL : 
  (delta < 2*MIN_VAL )? 2* MIN_VAL : delta;

      p_lls->l_levelset[k].Un[p] =  p_lls->l_levelset[k].U[p] +
	  delta;

 	p_lls->l_levelset[k].Un[p] =  
 	  (p_lls->l_levelset[k].Un[p] > 2000000000*MAX_VAL) ? 2000000000*MAX_VAL : 
 	  ((p_lls->l_levelset[k].Un[p] < 2000000000*MIN_VAL) ? 2000000000*MIN_VAL :  
	   p_lls->l_levelset[k].Un[p]); 
      }
    }
    }
  
  for (i = 1; i < p_image->ySize-1; i++) {
    for (j = 1; j < p_image->xSize-1; j++) {
      p = j + i * p_image->xSize ;
      change |= (*data.Upd_Fun)(p_lls,p_image,p);
    }
  }

  for(k=0;k<p_lls->Num;k++) 
    for (p = 0; p < p_image->size; p++)  
      p_lls->l_levelset[k].U[p] = p_lls->l_levelset[k].Un[p] ; 
}
  
  free(F);
  return(IMG_TRUE);

}


/*****************************************/
/*        Curvature term                 */
/*****************************************/

float Curv_cal(LevelSet *p_levelset, int p){

  int pu, pl, pr, pd ;
  int pul, pur, pdl, pdr ;
  float u_x, u_y ;
  float uxx, uxy, uyy ;
  float norm_sq ;
  float curvature_term ;
  
  pl = p-1 ;
  pr = p+1 ;
  pu = p - p_levelset->size_x ;
  pd = p + p_levelset->size_x ;
  pul = pu-1 ;
  pur = pu+1 ;
  pdl = pd-1 ;
  pdr = pd+1 ;
  
  u_x = .5 * (p_levelset->U[pr] - p_levelset->U[pl]) ;
  u_y = .5 * (p_levelset->U[pd] - p_levelset->U[pu]) ;
  uxx = (p_levelset->U[pr] - 2.0 * p_levelset->U[p] +
	 p_levelset->U[pl]) ;
  uxy = (p_levelset->U[pdr] - p_levelset->U[pur] -
	 p_levelset->U[pdl] + p_levelset->U[pul]) / 4.0 ;
  uyy = (p_levelset->U[pd] - 2.0 * p_levelset->U[p] +
	 p_levelset->U[pu]) ;
  
  norm_sq = u_x * u_x + u_y * u_y ;
  
  if(fabs(norm_sq) < EPSILON)
    return(0.0);
  
  curvature_term = (uxx * ( u_y * u_y) - 2.0 * u_x * u_y * uxy
		    + uyy * ( u_x * u_x)) / norm_sq; 

  return(curvature_term);

}

/*****************************************/
/*                 Nabla                 */
/*****************************************/

Nabla Nabla_cal(LevelSet *p_levelset, int p){

  int pu, pl, pr, pd ;
  int pul, pur, pdl, pdr ;
  float Dpx, Dpy, Dmx, Dmy ;
  float tmx, tmy, tpx, tpy ;
  float Nabla_p, Nabla_m ;
  Nabla nabla;

  pl = p-1 ;
  pr = p+1 ;
  pu = p - p_levelset->size_x ;
  pd = p + p_levelset->size_x ;
  pul = pu-1 ;
  pur = pu+1 ;
  pdl = pd-1 ;
  pdr = pd+1 ;
  
  Dpx = p_levelset->U[pr] - p_levelset->U[p] ;
  Dmx = p_levelset->U[p] - p_levelset->U[pl] ;
  Dpy = p_levelset->U[pd] - p_levelset->U[p] ;
  Dmy = p_levelset->U[p] - p_levelset->U[pu] ;
  
  tmx = imax (Dmx, 0) ;
  tpx = imin (Dpx, 0) ;
  tmy = imax (Dmy, 0) ;
  tpy = imin (Dpy, 0) ;
  Nabla_p = tmx * tmx + tpx * tpx + tmy * tmy + tpy * tpy ;
  Nabla_p = (float) sqrt ((double) Nabla_p) ;
  
  tmx = imin (Dmx, 0) ;
  tpx = imax (Dpx, 0) ;
  tmy = imin (Dmy, 0) ;
  tpy = imax (Dpy, 0) ;
  Nabla_m = tmx * tmx + tpx * tpx + tmy * tmy + tpy * tpy ;
  Nabla_m = (float) sqrt ((double) Nabla_m) ;

  nabla.p = Nabla_p;
  nabla.m = Nabla_m;

  return(nabla);
}


/*********************************************/
/*        Initialize list of LevelSets       */
/*********************************************/

List_LevelSet *Init_ListLevelSet(Image *p_image,int num_ls, int initial)
{

  List_LevelSet *list;
  int v, h, i;
  int j, k, d, pas, x, y, pos, l1, l2;
  int *final;
  float min;

  float v_c, h_c, offset;
  int count_pos, count_neg, h_limit, v_limit;
  float coeff, value_min, value_max;
  float distance;
  
  printf ("initial  %i ",initial);
  
  h_limit=p_image->xSize;
  v_limit=p_image->ySize;
  

  /* Get memory for the list */

  if((list = (List_LevelSet *)malloc(sizeof(List_LevelSet)))==NULL)
    return((List_LevelSet *)NULL);  
  list->Num = num_ls;
  if((list->l_levelset = (LevelSet *)malloc(num_ls*sizeof(LevelSet)))==NULL)
    return((List_LevelSet *)NULL);  

  /* Get memory for each levelset function */

  for(i=0;i<num_ls;i++){

    /* Initialize values */

    list->l_levelset[i].size_x = h_limit;
    list->l_levelset[i].size_y = v_limit;

    /* Get memory for arrays */

    if((list->l_levelset[i].U=(float *)malloc(h_limit*v_limit*sizeof(float)))==NULL)
      return((List_LevelSet *)NULL);  
    if((list->l_levelset[i].Un=(float *)malloc(h_limit*v_limit*sizeof(float)))==NULL)
      return((List_LevelSet *)NULL);     
  }    

  /* Initialize arrays */

if (initial == 1)
{
  count_pos = 0;
  count_neg = 0;
  value_min = -100;
  value_max = 100;
  coeff = (value_max-value_min);
  for(v=0;v<v_limit;v++) {
    for(h=0;h<h_limit;h++) {
      for(i=0;i<num_ls;i++){
	offset = i * 2 * M_PI /num_ls;
	v_c = v_limit/2*(1-sin(offset)/2.0);
	h_c = h_limit/2*(1-cos(offset)/2.0);
	
	distance = dist(v, h, v_c, h_c)/
	  ((float) sqrt((double) ((h_limit*h_limit/4)+(v_limit*v_limit/4))));
	list->l_levelset[i].U[v*h_limit+h] = 
	  coeff * exp(-distance*distance/.1)+value_min;
      }
    }
  }
}


/********** Initialize arrays with tiny circles all over the image *******************/

else
{
pas=1;
 x= 5;
 d=(int)((h_limit-pas*x)/x);
 pos=(int)(pas + (d/2));
 y=(int)((v_limit)/(d+pas));
 
final=(int *)malloc(x*y*sizeof(int));

h=0;
for(j=0;j<y;j++) {
 for(i=0;i<x;i++){

final[j*x+i]=h;
h+=1; 
if (h==num_ls) h=0; }}

 for(k=0;k<num_ls;k++){list->l_levelset[k].Num_pos=0;};
  list->Num_out=0;

for(k=0;k<num_ls;k++){
for(v=0;v<v_limit;v++) {
    for(h=0;h<h_limit;h++) {

min=9999;
for(j=0;j<y;j++) {
  for(i=0;i<x;i++){if
   ( ((v-(pos+j*(pas+d)))*(v-(pos+j*(pas+d)))+(h-(pos+i*(pas+d)))*(h-(pos+i*(pas+d)))<min)&(final[j*x+i]==k)){
min=(v-(pos+j*(pas+d)))*(v-(pos+j*(pas+d)))+(h-(pos+i*(pas+d)))*(h-(pos+i*(pas+d)));
l1=i;l2=j;}}};

list->l_levelset[k].U[v*h_limit+h]
	  = ((d/2)-dist(v,h,(pos+l2*(pas+d)),(pos+l1*(pas+d))))/sqrt(v_limit*v_limit+h_limit*h_limit);

if (list->l_levelset[k].U[v*h_limit+h]>0)
list->l_levelset[k].Num_pos+=1;
}}};

list->Num_out+=h_limit*v_limit;
for(k=0;k<num_ls;k++){list->Num_out-=list->l_levelset[k].Num_pos;}


free(final); }

  return(list);
}


/**********************************/
/*** Free LevelSet List         ***/
/**********************************/
void Free_ListLevelSet(List_LevelSet *p_lls){
 
  int i;

  for(i=0;i<p_lls->Num;i++){
    
    /* Free parameters and functions */
    free(p_lls->l_levelset[i].U);
    free(p_lls->l_levelset[i].Un);
    free(p_lls->l_levelset[i].p_param);
       
  }

  free(p_lls->l_levelset);
  free(p_lls->p_param_out);
  
  free(p_lls);
  
}

