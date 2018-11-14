

/* Contributors: C. Vazquez, I. Ben Ayed, M. Ben Salah, A. Mansouri, and S. Mseddi*/

#include <ctype.h>
#include <string.h>
#include "levelset.h"
#include "mex.h"
#include "matrix.h"



/* Function prototypes */

void fill_border_table(float* tab, int h_limit, int v_limit);
Data read_args(int argc, char *argv[], Data data);
void help();
Data init_params();
void print_args(Data data);
int  Seg2_Stim(Data data);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


/*--------------------------------------------------------------------------------*/
/*  fonction      :  void fill_border_table(float*, int, int)                     */
/*  description   :  fill the edges values of an array use the nearest            */
/*                   value in the array to fill each compartment                  */
/*  @return       :  void                                                         */
/*  @param        :  float* tab the array                                         */
/*                :  int h_limit the horizontal size of tab                       */
/*                :  int v_limit the vertical size of tab                         */
/*--------------------------------------------------------------------------------*/

void fill_border_table(float* tab, int h_limit, int v_limit) {

  int h, v;

  tab[0] = tab[h_limit+1];
  
  for(h=1;h<h_limit-1;h++) {
    tab[h] = tab[h+h_limit];
  }
  tab[h_limit-1] = tab[2*h_limit-2];
  
  for(v=1;v<v_limit-1;v++) {
    tab[v*h_limit] = tab[v*h_limit+1];
    tab[v*h_limit+h_limit-1] = tab[v*h_limit+h_limit-2];
  }
  
  tab[(v_limit-1)*h_limit] = tab[(v_limit-1)*h_limit-h_limit+1];
  for(h=1;h<h_limit-1;h++) {
    tab[(v_limit-1)*h_limit+h] = tab[(v_limit-1)*h_limit+h-h_limit];
  }
  tab[h_limit*v_limit-1] = tab[h_limit*v_limit-2-h_limit];
}

/*------------------------------------------------------------------------------*/
/*  fonction      :  void Seg2_Stim()                                           */
/*  description   :  level set image segmentation                               */
/*  @return       :  void                                                       */
/*  @param        :  void                                                       */
/*------------------------------------------------------------------------------*/
             
int Seg2_Stim(Data data) {    
           
  /*******************************************************************************/
  /* Variables declaration                                                       */
  /*******************************************************************************/

  char prg[]= "Seg2_Stim";
  char *nom_couleur[NB_COLOR]; 
    
  float* position;
  double *img_var;
  double *Energy;
  double *mean;
  double *large;


  int h_limit, v_limit, t_limit, offset, mean_pos; 
  int comp_pcs_in[IMG_NB_COMP],comp_pcs_out[IMG_NB_COMP];
  int i,j,l,t,p,k,h,v,iter,length,energ;


  /* Noise */
  float v1,v2,rsq,g1;

  unsigned int flg_in, no_change;
  char flag;
  int compt;

/* MEX files */

mxArray *ptrput[1];
mxArray *image_input[1];
mxArray *image_output[2];
mxArray *img_moy[1];
mxArray *img_energy[1];
mxArray *img_large1[1];
mxArray *img_large2[3];


  List_LevelSet *p_lls;
  Image *p_image;

  p_lls = (List_LevelSet *)malloc(sizeof(List_LevelSet));
  p_image = (Image *)malloc(sizeof(Image));

  h_limit = img_size(data.in_fname, 1);
  v_limit = img_size(data.in_fname, 0);
  t_limit = img_size(data.in_fname, 2);

  if((p_image->I0 = (float*)malloc(h_limit*v_limit*sizeof(float)))==NULL)
    mexErrMsgTxt("pas bon 1");
  p_image->xSize = h_limit;
  p_image->ySize = v_limit;
  p_image->size = h_limit*v_limit;
  p_image->tSize = t_limit;

  if((mean = (double*)malloc(h_limit*v_limit*sizeof(double)))==NULL)
    mexErrMsgTxt("error mean");

  if((Energy = (double*)malloc(100*sizeof(double)))==NULL)
    mexErrMsgTxt("error Energy");

 if((data.vect==IMG_TRUE))
{
      if((p_image->Il = (float**)malloc(t_limit*sizeof(float*)))==NULL)
	 mexErrMsgTxt("Error img vect");
      for(i=0;i<t_limit;i++)
	{
	if((p_image->Il[i] = (float*)malloc(h_limit*v_limit*sizeof(float)))==NULL)
	   mexErrMsgTxt("Error img vect");
        }
    
  }

/*****************************************   Reading the image   ***********************************************/

ptrput[0] = mxCreateString( data.in_fname );

mexCallMATLAB(1, image_input, 1, ptrput,"reading");
mexCallMATLAB(1, img_moy, 1, ptrput,"reading"); 

img_var = mxGetPr(image_input[0]);
mean = mxGetPr(img_moy[0]); 

if((data.vect==IMG_TRUE))
{
  for (k = 0; k < p_image->tSize; k++)
  { 
    l= k*p_image->xSize*p_image->ySize;

    for (j = 0; j < p_image->ySize; j++) 
    {
       for (i = 0; i < p_image->xSize; i++) 
       {   
	  p = i + j * p_image->xSize; 
	  p_image->Il[k][p]=(float)img_var[p+l];     
       }
    }
  }
}

else
{
 for (j = 0; j < p_image->ySize; j++) 
 {
   for (i = 0; i < p_image->xSize; i++) 
   {   
	p = i + j * p_image->xSize;      
	p_image->I0[p]= (float)img_var[p];

   }
 }
}


  /* Add noise to the image */
  if(data.var >0){
    for(v=0;v<v_limit;v++){
      for(h=0;h<h_limit;h++){
	offset = v * h_limit + h;
	do{
	  v1 = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
	  v2 = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
	  rsq = v1*v1 + v2*v2;
	}while(rsq>=1.0 || rsq == 0.0);
	g1 = v1*sqrt(-2.0*log(rsq)/rsq);
	p_image->I0[offset] += sqrt(data.var)*g1*255;
      }
    }
  }

  /* Initialize List of levelsets */

  p_lls = Init_ListLevelSet(p_image,data.N-1, data.start);
  (*data.Init_Fun)(p_lls,p_image);

  p_lls->data = &data;

  /*********************  Initialize the energy  ****************************/

compt=0;

	  mexCallMATLAB(1, img_energy, 0, NULL,"initEnergy");	
	  Energy = mxGetPr(img_energy[0]);
	  energ=0;	
	  Energy[energ]=cal_energy(p_lls, p_image, data.vect);

mexCallMATLAB(1, img_large1, 0, NULL,"initlarge");


  /*******************************************************************************/
  /* Draw contours on the original image                                         */
  /*******************************************************************************/

for(v=0;v<v_limit;v++) {
    for(h=0;h<h_limit;h++){
      offset = v*h_limit+h;	
      for(i=p_lls->Num-1;i>=0;i--){
	position = p_lls->l_levelset[i].U;
	if ((position[offset] >= 0.0) &&
	    ((position[offset+1] < 0.0) ||
	     (position[offset-1] < 0.0) ||
	     (position[offset+h_limit] < 0.0) ||
	     (position[offset-h_limit] < 0.0))) {
	  img_var [offset] = (i+1)*1000;
	  length+=1;
	  
	}

      }
		
    }
  }

if((data.vect==IMG_FALSE))
{
Energy[0]=Energy[0]+data.K*length;
}

if( mexCallMATLAB( 0, NULL, 1, image_input, "img_aff" ) )  
        mexErrMsgTxt("Unable to perform img_aff");

if( mexCallMATLAB( 0, NULL, 0, NULL, "stopping" ) )  
        mexErrMsgTxt("Unable to perform arret");
  

  /****************************************************************************/
  /*                       The LevelSet algorithm                             */
  /****************************************************************************/
 
  no_change = 0;
  
  if(
     (data.Init_Fun == Init_Vec_Gauss)||
     (data.Init_Fun == Init_Vec_Part_Gauss))
    mean_pos = ((t_limit+1)*t_limit+1);
  else
    mean_pos = (t_limit+1);

  
  printf("\n");
  printf("Iter: %d",0);
  for(i=p_lls->Num-1;i>=0;i--){
    printf("  u%d: %6.2f",i,p_lls->l_levelset[i].p_param[mean_pos]);
    if(p_lls->l_levelset[i].Par_Num == 5)
      printf("  v%d: %6.2f",i,(p_lls->l_levelset[i].p_param[3]));
  }
  printf("  u%d: %6.2f",p_lls->Num,p_lls->p_param_out[mean_pos]);
  if(p_lls->Par_Num == 5)
    printf("  v%d: %6.2f",p_lls->Num,(p_lls->p_param_out[3]));
  printf("\n");
  fflush(stdout);

  for(iter=1;(iter<data.Nls)&(no_change<20*MAX_NO_CHANGE);iter++){
  
    if(LS_LevelSetPDE(p_lls, p_image, data)==IMG_FALSE)
      no_change++;
    else
      no_change = 0;
    for(i=0;i<p_lls->Num;i++)
      fill_border_table(p_lls->l_levelset[i].U,h_limit,v_limit);

    if(!((iter+1)%data.nb_frames)){

for (j = 0; j < p_image->ySize; j++) {
 for (i = 0; i < p_image->xSize; i++) {   
	p = i + j * p_image->xSize;   

	if (data.vect==IMG_FALSE)  
	img_var[p] = (double)p_image->I0[p];
	else
	{
	  for(l=0;l<p_image->tSize;l++)
	  {
	    img_var[p+l*v_limit*h_limit]= (double)p_image->Il[l][p];
	  }
	}


}
}
length=0;
energ+=1;
Energy[energ]=cal_energy(p_lls, p_image, data.vect);;
compt+=1;

      /************************************************************/
      /* Draw contours on the original image                      */
      /************************************************************/

         for(v=0;v<v_limit;v++) {
	for(h=0;h<h_limit;h++) {
	  offset = v*h_limit+h;

	  if (data.vect==IMG_FALSE)
	      mean[offset]= cal_mean(p_lls, p_image, offset, data.vect);
	  else
	  {
	      k=0;
	      i=0;
	      while (k==0)
	      {
	        position = p_lls->l_levelset[i].U;
	        if (position[offset] >=0)
	        {
		  for(l=0;l<p_image->tSize;l++)
		  {
	          mean[offset+l*v_limit*h_limit] = p_lls->l_levelset[i].p_param[p_image->tSize*(p_image->tSize+1)+1 + l];
		  }
	          k=1;
	        }
	          else
	          {
	            if (i==p_lls->Num-1)
	            {
		      for(l=0;l<p_image->tSize;l++)
		      {
	              mean[offset+l*v_limit*h_limit] = p_lls->p_param_out[p_image->tSize*(p_image->tSize+1)+1 + l];
		      }
		      k=1;
	            }
		      else
		        i=i+1;
	      
	          }
	      }


	  }

	  for(i=p_lls->Num-1;i>=0;i--){
	    position = p_lls->l_levelset[i].U;
	    if ((position[offset] >= 0.0) &&
		((position[offset+1] < 0.0) ||
		 (position[offset-1] < 0.0) ||
		 (position[offset+h_limit] < 0.0) ||
		 (position[offset-h_limit] < 0.0))) 
	    {		  	
		img_var [offset] = (i+1)*1000;
		length+=1;

	    }

	}

      }
}
Energy[energ]=Energy[energ]+data.K*length;

if( mexCallMATLAB( 0, NULL, 1, image_input, "img_aff" ) )  
        mexErrMsgTxt("Unable to perform imshow");

    } 
if( mexCallMATLAB( 0, NULL, 0, NULL, "stopping" ) )  
        mexErrMsgTxt("Unable to perform imshow"); 


 /********************************************************************************/
 /*                        Write the image                                       */
 /********************************************************************************/


     printf("Iter: %d",iter); 
     for(i=p_lls->Num-1;i>=0;i--){
       printf("  u0%d: %6.2f",i,p_lls->l_levelset[i].p_param[mean_pos]);
       if(p_lls->l_levelset[i].Par_Num == 5)
	 printf("  VV0%d: %6.2f",i,(p_lls->l_levelset[i].p_param[4]));
     }
     printf("  u0%d: %6.2f",p_lls->Num,p_lls->p_param_out[mean_pos]);
     if(p_lls->Par_Num == 5)
       printf("  VV0%d: %6.2f",p_lls->Num,(p_lls->p_param_out[4]));
     printf("\r"); 
     fflush(stdout); 

  }

if((data.vect==IMG_FALSE))
{
if( mexCallMATLAB( 0, NULL, 1, img_energy, "Energy") )  
        mexErrMsgTxt("Unable to perform Energy");
}

if( mexCallMATLAB( 0, NULL, 1, img_moy, "img_mean") )  
        mexErrMsgTxt("Unable to perform Energy");

 image_output[0] = mxDuplicateArray(img_moy[0]);

image_output[1] = mxCreateString(data.out_fname);

if( mexCallMATLAB(0, NULL, 2, image_output, "writing" ) )  
        mexErrMsgTxt("Unable to perform size");

mxDestroyArray(image_input);
mxDestroyArray(image_output);

for(j=1;j<3;j++)
{
  mexCallMATLAB(1, img_large1, 0, NULL,"initlarge");
  large = mxGetPr(img_large1[0]);
    for(i=p_lls->Num-1;i>=0;i--)
    {
	large[i]=p_lls->l_levelset[i].p_param[(j-1)*2];
    }
  large[p_lls->Num]=p_lls->p_param_out[(j-1)*2];
img_large2[j] = mxDuplicateArray(img_large1[0]);
}

  img_large2[0] = mxDuplicateArray(img_moy[0]);
  mexCallMATLAB(0, NULL, 3, img_large2,"thelargest");

  Free_ListLevelSet(p_lls);
  free(p_image->I0);
  if(data.vect==IMG_TRUE){
    for(i=0;i<p_image->tSize;i++)
      free(p_image->Il[i]);
    free(p_image->Il);
  }

  free(p_image);

  return(IMG_NO_ERROR);
  
}                  
   
/*------------------------------------------------------------------------------*/
/*  fonction      :  void help()                                                */
/*  description   :  print a help for seg2D users                               */
/*  @return       :  void                                                       */
/*  @param        :  void                                                       */
/*------------------------------------------------------------------------------*/
  
void help() {
  printf("User Commands\n\n");
  printf("NAME\n");
  printf("\tseg - segment a ViDS image using its intensity components\n\n");
  printf("SYNOPSIS\n");
  printf("\t./Seg2D_Stim -in cosinus -field 2 -vec yes -iter 2800 -fra 10 -dt 10 -l 0.025 -r 2\n\n");
  printf("OPTIONS\n");
  printf("\t[-in <ViDS_input_file>]    name of ViDS input file (without extension)\n");
  printf("\t[-invect <ViDS_input_file>]    name of ViDS vector file (without extension)\n");
  printf("\t[-field <ViDS_input_field> field number of input image\n");
  printf("\t[-out <ViDS_output_file>]  name of ViDS output file (without extension)\n");
  printf("\t[-vec <yes|no>]            vector valued image?\n");
  printf("\t[-iter <nb_iter>]          number of iterations\n");
  printf("\t[-f <nb_frames>]           period of written frames (one out of nb_frames)\n");
  printf("\t[-dt <delta_t>]            temporal discretisation term\n");
  printf("\t[-l <lambda>]              curvature factor\n");
  printf("\t[-r <regions>]             number of regions\n");
  printf("\t[-sf <strategy(mean|gauss)>]     segmentation strategy\n");
  printf("\t[-v <variance>]            variance of noise\n"); 
}




/*------------------------------------------------------------------------------*/
/*  fonction      :  Data init_params()                                         */
/*  description   :  initialize default parmeters                               */
/*  @return       :  Data                                                       */
/*  @param        :  void                                                       */
/*------------------------------------------------------------------------------*/

Data init_params(){
  
  Data data;

  data.dt = .1;
  data.K = 2;
  data.Nls = 200;
  strcpy(data.in_fname, "brain.tif");
  strcpy(data.vect_fname,"");
  data.in_field = 0;
  strcpy(data.out_fname, "test.tif");
  data.nb_frames = 20;
  data.N = 3;
  data.vect = IMG_FALSE;
  data.Seg_Fun = Seg_Mean;
  data.Upd_Fun = Update_Mean;
  data.Init_Fun = Init_Mean;
  data.var = 0.0;
  data.start = 1;
  return(data);
}


/*------------------------------------------------------------------------------*/
/*  fonction      :  data read_args(int, char**, Data)                          */
/*  description   :  read the arguments for the Seg2D_Stim fonction             */
/*  @return       :  Data                                                       */
/*  @param        :  int argc, the number of arguments                          */
/*                :  char** argv, the list of arguments                         */
/*                :  Data data, The structure to keep the parametres            */
/*------------------------------------------------------------------------------*/

Data read_args(int argc, char *argv[], Data data) {

  /* Variable declarations */
 int   i;          
 char temp[IMG_M_LEN];
 char strat[10];

 /* Function body */
 for (i = 1; i < argc; i++) {
   if (argv[i][0] == '-') {
     if (strcmp(argv[i],"-h") == 0) {
       help();    
       
     } 
     if (strcmp(argv[i],"-in") == 0)  strcpy(data.in_fname,argv[++i]);
     if (strcmp(argv[i],"-invect") == 0)  strcpy(data.vect_fname,argv[++i]);
     if (strcmp(argv[i],"-out") == 0) strcpy(data.out_fname,argv[++i]);
     if (strcmp(argv[i],"-vec") == 0) data.vect =
					(tolower(argv[++i][0])=='y')?IMG_TRUE:IMG_FALSE;

     if (strcmp(argv[i],"-field") == 0)  data.in_field = atoi(argv[++i]);
     if (strcmp(argv[i],"-iter") == 0)  data.Nls = atoi(argv[++i]);
     if (strcmp(argv[i],"-f") == 0)     data.nb_frames = atoi(argv[++i]);
     if (strcmp(argv[i],"-dt") == 0)    data.dt = atof(argv[++i]);
     if (strcmp(argv[i],"-l") == 0)     data.K = atof(argv[++i]);
     if (strcmp(argv[i],"-r") == 0)     data.N = atoi(argv[++i]);  
     if (strcmp(argv[i],"-st") == 0)     data.start = atoi(argv[++i]); 
     if (strcmp(argv[i],"-v") == 0)     data.var = atof(argv[++i]);    
     if (strcmp(argv[i],"-sf") == 0){    
       if (strcmp(argv[++i],"mean") == 0){
	 strcpy(strat,"mean");
       }
       else{
	 if (strcmp(argv[i],"gauss") == 0){
	   strcpy(strat,"gauss");
	 }
	 else{
	   if (strcmp(argv[i],"part") == 0){
	     strcpy(strat,"part");
	   }
         else{
	   if (strcmp(argv[i],"clust") == 0){
	     strcpy(strat,"clust");
	   }
         else{
	   if (strcmp(argv[i],"ker") == 0){
	     strcpy(strat,"ker");
	   }
         else{
	   if (strcmp(argv[i],"gamma") == 0){
	     strcpy(strat,"gamma");
	   }
         else{
	   if (strcmp(argv[i],"ker_part") == 0){
	     strcpy(strat,"ker_part");
	   }
	   else{	   
	     printf("Segmentation function not yet implemented\n");
	     mexErrMsgTxt("methods error");
	   }
	 }
       }
      }
     }
    }
   }}    
 }}

 /* Define the functions to use */
 if(data.vect==IMG_FALSE){
   if (strcmp(strat,"mean") == 0){
     data.Seg_Fun = Seg_Mean;
     data.Upd_Fun = Update_Mean;
     data.Init_Fun = Init_Mean;

   }
   else{
     if (strcmp(strat,"gauss") == 0){
       data.Seg_Fun = Seg_Gauss;
       data.Upd_Fun = Update_Gauss;
       data.Init_Fun = Init_Gauss;
     }
     else{
       if (strcmp(strat,"part") == 0){
	 data.Seg_Fun = Seg_Part_Mean;
	 data.Upd_Fun = Update_Part_Mean;
	 data.Init_Fun = Init_Part_Mean;
       }
    else{
       if (strcmp(strat,"clust") == 0){
	 data.Seg_Fun = Seg_Clust_Mean;
	 data.Init_Fun = Init_Clust_Mean;
       }
    else{
     if (strcmp(strat,"ker") == 0){
       data.Seg_Fun = Seg_Ker;
       data.Upd_Fun = Update_Ker;
       data.Init_Fun = Init_Ker;
     }
    else{
     if (strcmp(strat,"gamma") == 0){
       data.Seg_Fun = Seg_Gamma;
       data.Upd_Fun = Update_Gamma;
       data.Init_Fun = Init_Gamma;
     }
    else{
     if (strcmp(strat,"ker_part") == 0){
       data.Seg_Fun = Seg_Part_Ker;
       data.Upd_Fun = Update_Part_Ker;
       data.Init_Fun = Init_Part_Ker;
     }
       else{
	 printf("Segmentation function not yet implemented\n"); 
	 mexErrMsgTxt("methods error");
       } 
      }
     }
    }
  }
 }}
}

 else{
   if (strcmp(strat,"mean") == 0){
     data.Seg_Fun = Seg_Vec_Mean;
     data.Upd_Fun = Update_Vec_Mean;
     data.Init_Fun = Init_Vec_Mean;
   }
   else{
     if (strcmp(strat,"gauss") == 0){
       data.Seg_Fun = Seg_Vec_Gauss;
       data.Upd_Fun = Update_Vec_Gauss;
       data.Init_Fun = Init_Vec_Gauss;
     }
     else{
       if (strcmp(strat,"part") == 0){
	 data.Seg_Fun = Seg_Vec_Part_Gauss;
	 data.Upd_Fun = Update_Vec_Part_Gauss;
	 data.Init_Fun = Init_Vec_Part_Gauss;
       }
       else{
	 printf("Segmentation function not yet implemented\n");
	mexErrMsgTxt("error");
       }
     }
   }
 }

 if (strcmp(data.out_fname,"") == 0) {
   for (i=strlen(data.in_fname); i>=0; i--)
     if (data.in_fname[i] == '/')
       break;
   strcpy(data.out_fname, data.in_fname+i+1);
   sprintf(temp,"_out_i%u_f%u_dt%05.0f_l%05.0f_r%u_v%04.0f_%s",
	   data.Nls,data.nb_frames,100000*data.dt,100000*data.K,data.N,data.var,strat);
   strcat(data.out_fname,temp);
 }

 return(data);
} 
  

/*------------------------------------------------------------------------------*/
/*  fonction      :  void print_args(Data)                                      */
/*  description   :  print the arguments of the seg2D fonction                  */
/*  @return       :  void                                                       */
/*  @param        :  Data data, Structure of parameters                         */
/*------------------------------------------------------------------------------*/

void print_args(Data data){

  printf("\n Input file: %s, field: %d\n", data.in_fname, data.in_field);
  printf("Output file : %s\n", data.out_fname);
  printf("Vector valued : %s\n", (data.vect==IMG_TRUE)?"Yes":"No");
  if(data.vect==IMG_TRUE)
    printf("Vector file: %s\n", data.vect_fname);
  printf("Number of iterations : %i\n", data.Nls);
  printf("Display frequency : %i\n", data.nb_frames);
  printf("Smoothness weight: %f\n", data.K);
  printf("Time step : %f\n", data.dt);
  printf("Number of regions : %d\n", data.N);
  if(data.Seg_Fun == Seg_Mean)
    printf("Segmentation strategy: %s\n","Mean");
  if((data.Seg_Fun == Seg_Gauss)||(data.Seg_Fun == Seg_Vec_Gauss))
    printf("Segmentation strategy: %s\n","Gauss"); 
  printf("Noise variance : %f\n", data.var); 

}


/*------------------------------------------------------------------------------*/
/*  fonction      :  mexFunction()                                              */
/*  description   :  mex fonction                                               */
/*------------------------------------------------------------------------------*/

	
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	int i;
	int argc;
	char** argv;
	double* p;
	
	Data data;
	p = mxGetPr(prhs[0]);	
	argc = p[0];
	argv = (char**)malloc((argc+1)*sizeof(char*));
						
  	for (i=0;i<argc+1;i++) 
	{
	   argv[i] = (char *) malloc (sizeof(char) * IMG_S_LEN);
	}

strcpy(argv[1], "-in");
mxGetString ( mxGetField(prhs[1], 0, "in"), argv[2], IMG_S_LEN );

strcpy(argv[3], "-out");
mxGetString ( mxGetField(prhs[1], 0, "out"), argv[4], IMG_S_LEN );

strcpy(argv[5], "-sf");
mxGetString ( mxGetField(prhs[1], 0, "sf"), argv[6], IMG_S_LEN );

strcpy(argv[7], "-iter");
mxGetString ( mxGetField(prhs[1], 0, "Nls"), argv[8], IMG_S_LEN );

strcpy(argv[9], "-r");
mxGetString ( mxGetField(prhs[1], 0, "N"), argv[10], IMG_S_LEN );

strcpy(argv[11], "-f");
mxGetString ( mxGetField(prhs[1], 0, "nb"), argv[12], IMG_S_LEN );

strcpy(argv[13], "-v");
mxGetString ( mxGetField(prhs[1], 0, "var"), argv[14], IMG_S_LEN );

strcpy(argv[15], "-l");
mxGetString ( mxGetField(prhs[1], 0, "K"), argv[16], IMG_S_LEN );

strcpy(argv[17], "-dt");
mxGetString ( mxGetField(prhs[1], 0, "dt"), argv[18], IMG_S_LEN );

strcpy(argv[19], "-st");
mxGetString ( mxGetField(prhs[1], 0, "start"), argv[20], IMG_S_LEN );

strcpy(argv[21], "-vec");
mxGetString ( mxGetField(prhs[1], 0, "vect"), argv[22], IMG_S_LEN );

  	data = init_params();

  	data = read_args(argc, argv, data);
  	print_args(data);

	
		   	for (i=0;i<argc+1;i++) {
				free(argv[i]);				
			}
			free(argv);
			//free(p);
	
  	if(!Seg2_Stim(data))	 
		return;
	else
		 mexErrMsgTxt("mexfunction error");
		 

			

}
