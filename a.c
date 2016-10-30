#include <stdio.h>
#include <setjmp.h>
#include <string.h>
#include <stdlib.h>
#include "jpeglib.h"

void binarizate_image(int height, int width, unsigned char* img) {
  unsigned char *p = malloc(sizeof(unsigned char) * width);
  double *g0 = malloc(sizeof(double) * width);
  double *g1 = malloc(sizeof(double) * width);
  double *tmp;

  int i, j;
  int c;
  double t;
  unsigned char *data = img;
  int s = width / 8;
  
  g1[width - 1] = 127 * s;
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      p[(i & 1) ? width - j - 1: j] = data[j];
    }
    
    for (j = 0; j < width; j++) {
      g0[j] = (j ? g0[j - 1] : g1[width - 1]) * ((double) (s - 1) / s) + p[j];
      t = i ? (g0[j] + g1[width - j - 1]) / 2 : g0[i];
      p[j] = (unsigned char) (p[j] < t / s * .63);  // .8 is a parameter
    }
    
    for (j = 0; j < width; j++) {
      data[j] = p[(i & 1) ? width - j - 1 : j] ? 0 : 255;
    }

    tmp = g0;
    g0 = g1;
    g1 = tmp;
    data += width;
  }
  
  free(p);
  free(g0);
  free(g1);
}

#define PUT_2B(array,offset,value)  \
        (array[offset] = (char) ((value) & 0xFF), \
         array[offset+1] = (char) (((value) >> 8) & 0xFF))
#define PUT_4B(array,offset,value)  \
        (array[offset] = (char) ((value) & 0xFF), \
         array[offset+1] = (char) (((value) >> 8) & 0xFF), \
         array[offset+2] = (char) (((value) >> 16) & 0xFF), \
         array[offset+3] = (char) (((value) >> 24) & 0xFF))
#define adjust4(x) x = ((x + 3) & ~3)

void write_bmp_header(j_decompress_ptr cinfo, FILE *output_file)
{
  char bmpfileheader[14];
  char bmpinfoheader[40];
  long headersize, bfSize;
  int bits_per_pixel, cmap_entries;

  int step;

  bits_per_pixel = 1;
  cmap_entries = 0;

  step = (cinfo->output_width + 7) / 8;

  adjust4(step);

  /* File size */
  headersize = 14 + 40 + 8; /* Header and colormap */

  bfSize = headersize + (long) step * (long) cinfo->output_height;

  /* Set unused fields of header to 0 */
  memset(bmpfileheader, 0, sizeof(bmpfileheader));
  memset(bmpinfoheader, 0 ,sizeof(bmpinfoheader));

  /* Fill the file header */
  bmpfileheader[0] = 0x42;/* first 2 bytes are ASCII 'B', 'M' */
  bmpfileheader[1] = 0x4D;
  PUT_4B(bmpfileheader, 2, bfSize); /* bfSize */
  /* we leave bfReserved1 & bfReserved2 = 0 */
  PUT_4B(bmpfileheader, 10, headersize); /* bfOffBits */

  /* Fill the info header (Microsoft calls this a BITMAPINFOHEADER) */
  PUT_2B(bmpinfoheader, 0, 40);   /* biSize */
  PUT_4B(bmpinfoheader, 4, cinfo->output_width); /* biWidth */
  PUT_4B(bmpinfoheader, 8, cinfo->output_height); /* biHeight */
  PUT_2B(bmpinfoheader, 12, 1);   /* biPlanes - must be 1 */
  PUT_2B(bmpinfoheader, 14, bits_per_pixel); /* biBitCount */
  /* we leave biCompression = 0, for none */
  /* we leave biSizeImage = 0; this is correct for uncompressed data */
  if (cinfo->density_unit == 2) { /* if have density in dots/cm, then */
    PUT_4B(bmpinfoheader, 24, (INT32) (cinfo->X_density*100)); /* XPels/M */
    PUT_4B(bmpinfoheader, 28, (INT32) (cinfo->Y_density*100)); /* XPels/M */
  }
  PUT_2B(bmpinfoheader, 32, cmap_entries); /* biClrUsed */
  /* we leave biClrImportant = 0 */

  if (fwrite(bmpfileheader, 1, 14, output_file) != (size_t) 14) {
    printf("write bmpfileheader error\n");
  }
  if (fwrite(bmpinfoheader, 1, 40, output_file) != (size_t) 40) {
    printf("write bmpinfoheader error\n");
  }
  
  char colors[] = {0, 0, 0, 0, 255, 255, 255, 0};
  if (fwrite(colors, 1, 8, output_file) != (size_t) 8) {
    printf("write colormap error\n");
  }
}

void write_pixel_data(j_decompress_ptr cinfo, unsigned char *output_buffer, FILE *output_file)
{
  int rows, cols;
  int row_width;
  int step;
  unsigned char *tmp = NULL;

  unsigned char *pdata;

  row_width = cinfo->output_width * cinfo->output_components;
  step = row_width;
  adjust4(step);

  int counter = 0;
  unsigned char *gray = malloc(cinfo->output_height * cinfo->output_width * (sizeof (unsigned char)));
  tmp = output_buffer + row_width * (cinfo->output_height - 1);
  for (rows = 0; rows < cinfo->output_height; rows++) {
    for (cols = 0; cols < row_width; cols += 3) {
      unsigned char 
        r = tmp[cols + 0],
        g = tmp[cols + 1],
        b = tmp[cols + 2];
      gray[counter] = (unsigned char) ((0.3 * r) + (0.59 * g) + (0.11 * b));
      counter++;
    }
    tmp -= row_width;
  }

  binarizate_image(cinfo->output_height, cinfo->output_width, gray);

  // Save image

  int step2 = (cinfo->output_width + 7) / 8;
  adjust4(step2);

  pdata = (unsigned char *)malloc(step2);

  counter = 0;
  for (rows = 0; rows < cinfo->output_height; rows++) {
    memset(pdata, 0, step2);
    for (cols = 0; cols < cinfo->output_width; cols++) {
      int c = gray[counter];
      counter++;

      if (c > 0) {
        pdata[cols / 8] |= 1 << (7 - cols % 8);
      }
    }
    fwrite(pdata, 1, step2, output_file);
  }

  free(gray);
  free(pdata);
}

int read_jpeg_file(const char *input_filename, const char *output_filename)
{
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  FILE *input_file;
  FILE *output_file;
  JSAMPARRAY buffer;
  int row_width;

  unsigned char *output_buffer;
  unsigned char *tmp = NULL;

  cinfo.err = jpeg_std_error(&jerr);

  if ((input_file = fopen(input_filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", input_filename);
    return -1;
  }

  if ((output_file = fopen(output_filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s\n", output_filename);
    return -1;
  }

  jpeg_create_decompress(&cinfo);

  /* Specify data source for decompression */
  jpeg_stdio_src(&cinfo, input_file);

  /* Read file header, set default decompression parameters */
  (void) jpeg_read_header(&cinfo, TRUE);

  /* Start decompressor */
  (void) jpeg_start_decompress(&cinfo);

  row_width = cinfo.output_width * cinfo.output_components;

  buffer = (*cinfo.mem->alloc_sarray)
          ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_width, 1);

  write_bmp_header(&cinfo, output_file);

  output_buffer = (unsigned char *)malloc(row_width * cinfo.output_height);
  memset(output_buffer, 0, row_width * cinfo.output_height);
  tmp = output_buffer;

  /* Process data */
  while (cinfo.output_scanline < cinfo.output_height) {
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);

    memcpy(tmp, *buffer, row_width);
    tmp += row_width;
  }

  write_pixel_data(&cinfo, output_buffer, output_file);

  free(output_buffer);

  (void) jpeg_finish_decompress(&cinfo);

  jpeg_destroy_decompress(&cinfo);

  /* Close files, if we opened them */
  fclose(input_file);
  fclose(output_file);
  return 0;
}

int main(int argc, char *argv[])
{
  if (argc < 3) {
    read_jpeg_file("1.jpg", "1.bmp");
  } else {
    read_jpeg_file(argv[1], argv[2]);
  }
  return 0;
}

