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
#define NAME_IMG_OUT1 "image-TpIFT3205-2-5"

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

float matrix_val(float** m, int j, int i, int length, int width) {

    if(i < 0 || i >= length)
        return 0;
    
    if(j < 0 || j >= width)
        return 0;

    return m[i][j];
}

float matrix_valw(float** m, int j, int i, int length, int width) {

    if(i < 0 || i >= length)
        i = (i + length) % length;
    
    if(j < 0 || j >= width)
        j = (j + width) % width;

    return m[i][j];
}

void matrix_rotation_ppv(float** src, float** dest, float angle, int l, int w) {
    int x1,y1,x2,y2;
    float xp,yp;

    for(y2=0; y2<l; y2++) {
        for(x2=0; x2<w; x2++) {

            xp =  (x2 - w/2) * cos(-angle) + (y2 - l/2) * sin(-angle) + w/2;
            yp = -(x2 - w/2) * sin(-angle) + (y2 - l/2) * cos(-angle) + l/2;

            x1 = floor(xp);
            y1 = floor(yp);

            dest[y2][x2] = matrix_val(src, x1, y1, l, w);
        }
    }
}

void matrix_rotation_billy(float** src, float** dest, float angle, int l, int w) {
    int x1,y1,x2,y2;
    float xp,yp;
    
    float dx, dy;
    float f_xy, f_x1y, f_xy1, f_x1y1, f_xpy, f_xpy1;
    
    for(y2=0; y2<l; y2++) {
        for(x2=0; x2<w; x2++) {

            xp =  (x2 - w/2) * cos(-angle) + (y2 - l/2) * sin(-angle) + w/2;
            yp = -(x2 - w/2) * sin(-angle) + (y2 - l/2) * cos(-angle) + l/2;

            x1 = (int) xp;
            y1 = (int) yp;

            
            if(y1 < 0 || y1 >= l || x1 < 0 || x1 >= w) {
                dest[y2][x2] = 0;
                continue;
            }

    
            dx = xp - x1;
            dy = yp - y1;
            
            f_xy = matrix_valw(src, x1, y1, l, w);

            f_x1y = matrix_valw(src, x1 + 1, y1, l, w);

            f_xy1 = matrix_valw(src, x1, y1 + 1, l, w);

            f_x1y1 = matrix_valw(src, x1 + 1, y1 + 1, l, w);

            f_xpy = f_xy + dx * (f_x1y - f_xy);

            f_xpy1 = f_xy1 + dx * (f_x1y1 - f_xy1);

            dest[y2][x2] = f_xpy + dy * (f_xpy1 - f_xpy);
        }
    }
}

/*------------------------------------------------*/
/* ET MAINTENANT, NOTRE PROGRAMME PRINCIPAL   ----*/
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
  float** MatriceImgC=fmatrix_allocate_2d(length,width);
  float** MatriceImgIC=fmatrix_allocate_2d(length,width);
  float** MatriceImg4 = fmatrix_allocate_2d(length,width);
  float** MatriceImgI4 = fmatrix_allocate_2d(length,width);
  float** MatriceImgM4 = fmatrix_allocate_2d(length,width);
  
  //Lecture Image 
  float** MatriceImg1=LoadImagePgm(NAME_IMG_IN1,&length,&width);
  float** MatriceImg2=LoadImagePgm(NAME_IMG_IN2,&length,&width);

 
  // .... .... .... .... .... .... .... .... .... .... .... ....
  PreFFT_Translation(MatriceImg1,length,width);
  PreFFT_Translation(MatriceImg2,length,width);
  // ... ... ... ... ... ... ... ... ... ... ... ... ... ... ...
  FFTDD(MatriceImg1,MatriceImgI1,length,width);
  FFTDD(MatriceImg2,MatriceImgI2,length,width);
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  Mod(MatriceImgM1,MatriceImg1,MatriceImgI1,length,width);
  Mod(MatriceImgM2,MatriceImg2,MatriceImgI2,length,width);
  // .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..

  float angle, best_angle, min_error = INFINITY, error;

  for(angle = -PI/16; angle < PI/16; angle += 0.005) {

      // Rotation
      matrix_rotation_billy(MatriceImgM2, MatriceImgM3, angle, length, width);

      // Calcul de l'erreur
      error = 0;
      
      for(i=0;i<length;i++) 
          for(j=0;j<width;j++) 
          {
              error += fabs(MatriceImgM3[i][j] - MatriceImgM1[i][j]);
          }

      // Keep min
      if(error < min_error) {
        min_error = error;
        best_angle = angle;
      }
  }

  IFFTDD(MatriceImg2, MatriceImgI2, length, width);
  PreFFT_Translation(MatriceImg2, length, width);
  matrix_rotation_billy(MatriceImg2, MatriceImg3, best_angle, length, width);

  
  PreFFT_Translation(MatriceImg3, length, width);
  for(i=0; i<length; i++) {
    for(j=0; j<width; j++) { 
      MatriceImg4[i][j] = MatriceImg3[i][j];
      MatriceImgI4[i][j] = MatriceImgI3[i][j];
    }
  }
  for(i=0; i<length/2; i++) {
    for(j=0; j<width/2; j++) {
         float tmp = MatriceImg3[i][j];
         MatriceImg3[i][j] = MatriceImg3[length - i-1][length - j-1];
         MatriceImg3[length - i-1][length - j-1] = tmp;
    }
  }
  FFTDD(MatriceImg4, MatriceImgI4, length, width);
  Mod(MatriceImgM4,MatriceImg4,MatriceImgI4,length,width);
  
  FFTDD(MatriceImg3, MatriceImgI3, length, width);

  for(i=0; i<length; i++) {
    for(j=0; j<width; j++) {
        //Multiplication complexe
        MatriceImgC[i][j] = MatriceImg1[i][j] * MatriceImg3[i][j] - MatriceImgI1[i][j] * MatriceImgI3[i][j];
        MatriceImgIC[i][j] = MatriceImg1[i][j] * MatriceImgI3[i][j] - MatriceImgI1[i][j] * MatriceImg3[i][j];

        MatriceImg2[i][j] = MatriceImgC[i][j] / (MatriceImgM4[i][j]*MatriceImgM4[i][j]);
        MatriceImgI2[i][j] = MatriceImgIC[i][j] / (MatriceImgM4[i][j]*MatriceImgM4[i][j]);
    }
  }

  IFFTDD(MatriceImg2, MatriceImgI2, length, width);
  PreFFT_Translation(MatriceImg2, length, width);
  float max = -INFINITY;
  int col =0 ;
  int row=0;
  for(i=0; i<length; i++) {
    for(j=0; j<width; j++) {
        if(max < MatriceImg2[i][j]) {
            max = MatriceImg2[i][j];
            col = j - width/2;
            row = i - length/2;
        }
    }
  }
  printf("Ligne = %d; Colonne = %d\n", row, col);
  // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  Recal(MatriceImg2,length,width);
  //Recal(MatriceImg3,length,width);

  //Sauvegarde
  SaveImagePgm(NAME_IMG_OUT1,MatriceImgM3,length,width);

  //Commande systeme: VISU
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_OUT1);
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
  free_fmatrix_2d(MatriceImgIC);
  free_fmatrix_2d(MatriceImgC);
  //retour sans probleme
  printf("\n C'est fini ... \n\n");
  return 0;
}
