#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#include <stdlib.h>

/* Constants definitions */

#define NB_COLOR 6
#define DEFAULT_COLOR 0
#define COLOR_M_LEN 10

#define EPSILON 1e-10
#define MAX_VAL  100000000.0
#define MIN_VAL -100000000.0

#define MAX_NO_CHANGE 500
#define sigma  10000


/* Macro definitions */
#define imin(x,y) ((x < y) ? x : y)
#define imax(x,y) ((x < y) ? y : x)
#define max(x,y)  ((x > y) ? x : y)
#define min(x,y)  ((x < y) ? x : y)
#define dist(x_a, y_a, x_b, y_b) ((float) sqrt((double) ((x_a-x_b)*(x_a-x_b)+(y_a-y_b)*(y_a-y_b))))

#define IMG_NB_COMP        3    /* Number of components in the image (max)  */
#define IMG_S_LEN         80    /* Length of a small string                 */
#define IMG_M_LEN        200    /* Length of a medium string                */
#define IMG_L_LEN        400    /* Length of a long string                  */


/* ENGLISH */ 
#define IMG_TRUE           1    /* Returned by functions value when true    */ 
#define IMG_FALSE          0    /* Returned by functions value when false   */ 

/* General flags */
#define IMG_NO_ERROR       0    /* Execution ok                             */
#define IMG_FATAL          1    /* Fatal error occured                      */


/* Structure definitions */


/* LevelSet structure */
typedef struct LevelSet
{
  float *U ;
  float *Un ;
  int size_x;
  int size_y;
  int Num_pos;         /* Number of positive points */
  int Par_Num;         /* Number of parameters */
  double *p_param;
} LevelSet;

/* List of LevelSet for multi-region segmentation */
typedef struct List_LevelSet
{
  int Num;
  struct LevelSet *l_levelset;
  int Num_out;        /* Number of points that are negatives for all
			             functions */
  int Par_Num;        /* Number of parameters for the exterior */
  double *p_param_out; /* List of parameters for the exterior */
  struct Data *data;
} List_LevelSet;

/* Image structure */
typedef struct Image
{
  int xSize;      /* X Size */
  int ySize;      /* Y Size */
  int size;       /* xSize*ySize */
  int tSize;      /* Number of components in vector valued image */
  float *I0;      /* Intensity image, to apear in output image */
  float **Il;     /* List of images for vector-valued */
} Image ;

/* LevelSet Algorithm parameters */
typedef struct Data
{
  int Nls;                   /* Number of iterations */
  float dt;                  /* delta time */
  float K;                   /* Curbature factor */
  int N;                     /* Number of regions */
  unsigned char vect;        /* Vector valued image? */
  char in_fname[IMG_M_LEN];  /* Input file name */
  int in_field;              /* field number of input image */
  char vect_fname[IMG_M_LEN];/* Input file name with vector images */
  char out_fname[IMG_M_LEN]; /* Output file name */
  int nb_frames;             /* Number of frames to skip before displaying */
  int (*Init_Fun)(List_LevelSet *, Image *); /*Init function*/
  float *(*Seg_Fun)(List_LevelSet *, Image *, int, float *); /*Segmentation function*/
  int (*Upd_Fun)(List_LevelSet *, Image *, int); /*Update function*/
  float var;                 /* variance of noise */
  int start;			     /* type of initalisation */
} Data;

/* Nabla structure */
typedef struct
{
  float p;
  float m;
} Nabla;

/* Force structure */
typedef struct
{
  float Force;          /* positive value unless background */
  int levelset;         /* position of the function to increase or 
			   -1 for the background */
} Force;


/****************************************************************************/
/*		  Image files		    */
/****************************************************************************/

int img_size (char* fichier, int ordre);
double cal_energy(List_LevelSet *p_lls, Image *p_image, unsigned char vect);
double cal_mean(List_LevelSet *p_lls, Image *p_image, int offset, unsigned char vect);
double* reception(Image *p_image);


/* Function prototypes */

/* Mean */
float *Seg_Mean(List_LevelSet *p_lls, Image *p_image, int p,float *F);
int    Init_Mean(List_LevelSet *p_lls, Image *p_image);
int    Update_Mean(List_LevelSet *p_lls, Image *p_image, int p);

float *Seg_Vec_Mean(List_LevelSet *p_lls, Image *p_image, int p,float *F);
int    Init_Vec_Mean(List_LevelSet *p_lls, Image *p_image);
int    Update_Vec_Mean(List_LevelSet *p_lls, Image *p_image, int p);

float *Seg_Part_Mean(List_LevelSet *p_lls, Image *p_image, int p,float *F);
int    Init_Part_Mean(List_LevelSet *p_lls, Image *p_image);
int    Update_Part_Mean(List_LevelSet *p_lls, Image *p_image, int p);

float *Seg_Clust_Mean(List_LevelSet *p_lls, Image *p_image, int p,float *F);
int    Init_Clust_Mean(List_LevelSet *p_lls, Image *p_image);

/* Gauss */
float *Seg_Gauss(List_LevelSet *p_lls, Image *p_image, int p,float *F);
int    Init_Gauss(List_LevelSet *p_lls, Image *p_image);
int    Update_Gauss(List_LevelSet *p_lls, Image *p_image, int p);

float *Seg_Vec_Gauss(List_LevelSet *p_lls, Image *p_image, int p,float *F);
int    Init_Vec_Gauss(List_LevelSet *p_lls, Image *p_image);
int    Update_Vec_Gauss(List_LevelSet *p_lls, Image *p_image, int p);

float *Seg_Vec_Part_Gauss(List_LevelSet *p_lls, Image *p_image, int p,float *F);
int    Init_Vec_Part_Gauss(List_LevelSet *p_lls, Image *p_image);
int    Update_Vec_Part_Gauss(List_LevelSet *p_lls, Image *p_image, int p);

/* Kernel */
float *Seg_Ker(List_LevelSet *p_lls, Image *p_image, int p, float *F);
int    Init_Ker(List_LevelSet *p_lls, Image *p_image);
int    Update_Ker(List_LevelSet *p_lls, Image *p_image, int p);

float *Seg_Part_Ker(List_LevelSet *p_lls, Image *p_image, int p, float *F);
int    Init_Part_Ker(List_LevelSet *p_lls, Image *p_image);
int    Update_Part_Ker(List_LevelSet *p_lls, Image *p_image, int p);

/* Gamma */
float *Seg_Gamma(List_LevelSet *p_lls, Image *p_image, int p, float *F);
int    Init_Gamma(List_LevelSet *p_lls, Image *p_image);
int    Update_Gamma(List_LevelSet *p_lls, Image *p_image, int p);

/* Irregular reconstruction */

float *Seg_TV(List_LevelSet *p_lls, Image *p_image, int p,float *F);
int    Init_TV(List_LevelSet *p_lls, Image *p_image);
int    Update_TV(List_LevelSet *p_lls, Image *p_image, int p);

/* LevelSet functions */
float          Curv_cal(LevelSet *p_levelset, int p);
float          Grad_cal(LevelSet *p_levelset, int p);
Nabla          Nabla_cal(LevelSet *p_levelset, int p);
List_LevelSet *Init_ListLevelSet(Image *p_image,int num_ls, int initial);
int            LS_LevelSetPDE (List_LevelSet *p_lls,Image *p_image, Data data);
void           Free_ListLevelSet(List_LevelSet *p_lls);
