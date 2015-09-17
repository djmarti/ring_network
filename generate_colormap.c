#include <stdio.h>
#include <stdlib.h>
#include "colormap.h"

const int width = 20;
const int height = 256;
const int MAXVAL = 255;         /* maximum value for color codes */

void save_array_ppm(double *a, const char *fname)
{
    unsigned int rgb_bone[MAXVAL + 1][3];
    fill_cm_lookuptable_bone(rgb_bone, MAXVAL);

    FILE *fin = fopen(fname, "w");
    fprintf(fin, "P3\n%d %d\n%d\n", width, height, MAXVAL);

    int pos;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            pos = (int) rint((double) MAXVAL * a[width*i + j]);
            fprintf(fin, "% 3d % 3d % 3d   ", rgb_bone[pos][0], rgb_bone[pos][1],
                    rgb_bone[pos][2]);
        }
        fprintf(fin, "\n");
    }
    fclose(fin);
}

int main(int argc, char *argv[])
{
    argc = 0;
    argv = NULL;
    int N = width * height;
    double array[N];
    int status = 0;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            array[width * i + j] = 1.0 - (double) i / (double) MAXVAL;
        }
    }
    save_array_ppm(array, "colormap_bone.ppm");
    status = system("convert colormap_bone.ppm colormap_bone.jpg");
    return status;
}
