
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "spring.h"

#ifdef OPENGL
#include "display.h"
#ifdef LIBPNG
#include <png.h>
#endif // LIBPNG
#ifdef _APPLE
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif  // _APPLE
#endif  // OPENGL



extern int NS;

// write out springs to a file
void write_springs(struct reb_simulation* const r,char *fileroot){
   FILE *fpo;
   char filename[100];
   strcpy(filename,fileroot);
   strcat(filename,"_springs.txt");
   fpo = fopen(filename,"w");
   for(int i=0;i<NS;i++){
      double dr = spring_length(r,springs[i]);
      fprintf(fpo,"%d %d %.6f %.6f %.6f %.6f ",
        springs[i].i, springs[i].j, 
        springs[i].ks, springs[i].rs0, 
        springs[i].gamma, springs[i].smax);
      fprintf(fpo,"%.6f %.6f\n",
         strain(r,springs[i]), dr);
   }
   fclose(fpo);
   printf("\n write_springs: NS=%d %s\n",NS,filename);
}


void read_springs(struct reb_simulation* const r,char *fileroot){
   char filename[100];
   strcpy(filename,fileroot);
   strcat(filename,"_springs.txt");
   printf("\n reading in springs %s\n",filename);
   FILE *fpi;
   fpi = fopen(filename,"r");
   char string[300];
   struct spring spr;
   int  i,j; 
   double ks,rs0,gamma,smax;
   while(fgets(string,300,fpi) != NULL){
      sscanf(string,"%d %d %lf %lf %lf %lf",
        &i,&j,&ks,&rs0,&gamma,&smax);
      spr.i=i; spr.j=j;
      spr.ks=ks; spr.rs0=rs0; spr.gamma=gamma; spr.smax=smax; 
      springs_add(r,spr);
   }
   fclose(fpi);
   printf("read_springs: NS=%d\n",NS);
}


void read_particles(struct reb_simulation* const r,char *fileroot){
   char filename[100];
   strcpy(filename,fileroot);
   strcat(filename,"_particles.txt");
   printf("\n reading in particles %s\n",filename);
   FILE *fpi;
   fpi = fopen(filename,"r");
   char string[300];
   struct reb_particle pt;
   pt.ax = 0.0; pt.ay = 0.0; pt.az = 0.0;
   double x,y,z,vx,vy,vz,m,rad;
   while(fgets(string,300,fpi) != NULL){
      sscanf(string,"%lf %lf %lf %lf %lf %lf %lf %lf\n",
        &x, &y, &z, &vx, &vy, &vz, &rad, &m);
      pt.x = x;   pt.y = y;   pt.z = z;
      pt.vx = vx; pt.vy = vy; pt.vz = vz;
      pt.r = rad; pt.m = m; 
      reb_add(r,pt);
   }
   fclose(fpi);
   printf("read_particles: N=%d\n",r->N);
}

void write_particles(struct reb_simulation* const r,char *fileroot){
   FILE *fpo;
   char filename[100];
   strcpy(filename,fileroot);
   strcat(filename,"_particles.txt");
   fpo = fopen(filename,"w");
   for(int i=0;i< r->N;i++){
      fprintf(fpo,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6e\n",
        r->particles[i].x, r->particles[i].y, r->particles[i].z,
        r->particles[i].vx, r->particles[i].vy, r->particles[i].vz,
        r->particles[i].r, r->particles[i].m); 
   }
   fclose(fpo);
   printf("\n write_particles: N=%d %s\n",r->N,filename);
}

#ifdef LIBPNG
unsigned char* 	imgdata = NULL;
int output_png_num = 0;
void output_png(struct reb_simulation* const r,char* dirname){
	char filename[1024];
	sprintf(filename,"%s%09d.png",dirname,output_png_num);
	output_png_num++;
	output_png_single(r,filename);
}

void output_png_single(struct reb_simulation* const r, char* filename){
	// Read Image
	// if (display_init_done==0) return;
        if (r->t < 3.0*r->dt) return; 
        printf("got here \n"); // got here
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport); // bombs here!
	int width = viewport[2];
	int height = viewport[3];
        printf("png: %d %d\n",width,height);
	glReadBuffer(GL_BACK);
	//glReadBuffer(GL_FRONT);
	if (imgdata==NULL){
		imgdata = calloc(width*height*3,sizeof(unsigned char));
	}
	png_byte* row_pointers[height];
	for (int h = 0; h < height; h++) {
		row_pointers[height-h-1] = (png_bytep) &imgdata[width*3*h];
	}

	glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,imgdata);

	/* open the file */
	FILE *fp;
	png_structp png_ptr;
	png_infop info_ptr;
	fp = fopen(filename, "wb");
	if (fp==NULL){
		printf("\n\nError while opening file '%s'.\n",filename);
		return;
	}

	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (png_ptr == NULL) {
		fclose(fp);
		return;
	}

	/* Allocate/initialize the image information data.  REQUIRED */
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fclose(fp);
		png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
		return;
	}

	/* Set error handling.  REQUIRED if you aren't supplying your own
	* error hadnling functions in the png_create_write_struct() call.
	*/
	/*
	if (setjmp(png_ptr->jmpbuf))
	{
	fclose(fp);
	png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
	return;
	}
	*/

	/* I/O initialization functions is REQUIRED */
	/* set up the output control if you are using standard C streams */
	png_init_io(png_ptr, fp);

	/* Set the image information here.  Width and height are up to 2^31,
	* bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
	* the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
	* PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
	* or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
	* PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
	* currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
	*/
	png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB,
	PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	/* Write the file  */
	png_write_info(png_ptr, info_ptr);
	png_write_image(png_ptr, row_pointers);
	png_write_end(png_ptr, info_ptr);

	/* if you allocated any text comments, free them here */

	/* clean up after the write, and free any memory allocated */
	png_destroy_write_struct(&png_ptr, (png_infopp)NULL);

	/* close the file */
	fclose(fp);
}
#endif

