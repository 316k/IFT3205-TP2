/*------------------------------------------------------*/
/* Prog    : Tp2_IFT3205-2-4.c                          */
/* Auteur  :                                            */
/* Date    : --/--/2010                                 */
/* version :                                            */ 
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo2.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/   
/*------------------------------------------------*/
#define NAME_VISUALISER "display "
#define NAME_IMG_IN1  "UdM_1"
#define NAME_IMG_IN2  "UdM_2"
#define NAME_IMG_OUT1 "image-TpIFT3205-2-1a"
#define NAME_IMG_OUT2 "image-TpIFT3205-2-1b"

/*------------------------------------------------*/
/* PROTOTYPE DE FONCTIONS  -----------------------*/   
/*------------------------------------------------*/
void PreFFT_Translation(float** Matrix, int length, int width) {
    int x,y;
    for(x=0; x<length; x++)
        for(y=0; y<width; y++)
            if((x + y) % 2 == 1)
                Matrix[x][y] = -Matrix[x][y];
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/
int main(int argc,char **argv)
 {
  int i,j,k;
  int length,width;
  float Theta0;
  int x0,y0;
  char BufSystVisuImg[100];

  //Constante
  length=512;
  width=512;
  
  //Allocation Memoire 
  float** MatriceImgI1=fmatrix_allocate_2d(length,width);
  float** MatriceImgM1=fmatrix_allocate_2d(length,width);
  float** MatriceImgI2=fmatrix_allocate_2d(length,width);
  float** MatriceImgM2=fmatrix_allocate_2d(length,width);
  float** MatriceImgI3=fmatrix_allocate_2d(length,width);
  float** MatriceImgM3=fmatrix_allocate_2d(length,width);
  float** MatriceImgR3=fmatrix_allocate_2d(length,width);
  float** MatriceImg3=fmatrix_allocate_2d(length,width);

  //Lecture Image 
  float** MatriceImg1=LoadImagePgm(NAME_IMG_IN1,&length,&width);
  float** MatriceImg2=LoadImagePgm(NAME_IMG_IN2,&length,&width);

 
  // .... .... .... .... .... .... .... .... .... .... .... ....
  PreFFT_Translation(MatriceImg1,length,width);
  PreFFT_Translation(MatriceImg2,length,width);
  // ... ... ... ... ... ... ... ... ... ... ... ... ... ... ...
  FFTDD(MatriceImg1,MatriceImgI1,length,width);
  FFTDD(MatriceImg2,MatriceImgI2,length,width);
  // .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..
  Mod(MatriceImgM1,MatriceImg1,MatriceImgI1,length,width);
  Mod(MatriceImgM2,MatriceImg2,MatriceImgI2,length,width);
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  /*Pour visu*/
    for(i=0;i<length;i++) 
        for(j=0;j<width;j++) 
        {
            MatriceImgM1[i][j] = log(1 + MatriceImgM1[i][j]);
            MatriceImgM2[i][j] = log(1 + MatriceImgM2[i][j]);
        }
  Recal(MatriceImgM1,length,width);
  Recal(MatriceImgM2,length,width);

  //Sauvegarde
  SaveImagePgm(NAME_IMG_OUT1,MatriceImgM1,length,width);
  SaveImagePgm(NAME_IMG_OUT2,MatriceImgM2,length,width);

  //Commande systeme: VISU
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_OUT1);
  strcat(BufSystVisuImg,".pgm&");
  printf(" %s",BufSystVisuImg);
  system(BufSystVisuImg);
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_OUT2);
  strcat(BufSystVisuImg,".pgm&");
  printf(" %s",BufSystVisuImg);
  system(BufSystVisuImg);


  //==End=========================================================

  //Liberation memoire 
  free_fmatrix_2d(MatriceImgI1);
  free_fmatrix_2d(MatriceImgM1);
  free_fmatrix_2d(MatriceImgI2);
  free_fmatrix_2d(MatriceImgM2);
  free_fmatrix_2d(MatriceImgR3);
  free_fmatrix_2d(MatriceImgI3);
  free_fmatrix_2d(MatriceImgM3);
  free_fmatrix_2d(MatriceImg1);
  free_fmatrix_2d(MatriceImg2);  
  free_fmatrix_2d(MatriceImg3);

  //retour sans probleme
  printf("\n C'est fini ... \n\n");
  return 0;
}
