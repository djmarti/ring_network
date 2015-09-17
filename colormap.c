#include "colormap.h"

/* How the dicts are organized (colors.py in matploblib)
 *
 *        Example: suppose you want red to increase from 0 to 1 over
        the bottom half, green to do the same over the middle half,
        and blue over the top half.  Then you would use::

            cdict = {'red':   [(0.0,  0.0, 0.0),
                               (0.5,  1.0, 1.0),
                               (1.0,  1.0, 1.0)],

                     'green': [(0.0,  0.0, 0.0),
                               (0.25, 0.0, 0.0),
                               (0.75, 1.0, 1.0),
                               (1.0,  1.0, 1.0)],

                     'blue':  [(0.0,  0.0, 0.0),
                               (0.5,  0.0, 0.0),
                               (1.0,  1.0, 1.0)]}

        Each row in the table for a given color is a sequence of
        *x*, *y0*, *y1* tuples.  In each sequence, *x* must increase
        monotonically from 0 to 1.  For any input value *z* falling
        between *x[i]* and *x[i+1]*, the output value of a given color
        will be linearly interpolated between *y1[i]* and *y0[i+1]*::

            row i:   x  y0  y1
                           /
                          /
            row i+1: x  y0  y1

        Hence y0 in the first row and y1 in the last row are never used.

*/


void fill_cm_lookuptable_bone(unsigned int table[][3], const int N)
{
    /* original bone specs. From _cm.py (matplotlib)
     *
     * {'red':   ((0., 0., 0.),(0.746032, 0.652778, 0.652778),(1.0, 1.0, 1.0)),
     *  'green': ((0., 0., 0.),(0.365079, 0.319444, 0.319444),
     *            (0.746032, 0.777778, 0.777778),(1.0, 1.0, 1.0)),
     *  'blue':  ((0., 0., 0.),(0.365079, 0.444444, 0.444444),(1.0, 1.0, 1.0))} 
     *
     */
    
    /* bone clippled at 15% and 90% */
    struct intval_node red_desc[] = {
        {0.0, 0.15}, /* {0.0, 0.0} */
        {0.746032, 0.652778}, /* {0.746032, 0.652778} */
        {1.0, 0.9}
    };

    /* Green */
    struct intval_node green_desc[] = {
        {0.0, 0.15},
        {0.365079, 0.319444},
        {0.746032, 0.777778},
        {1.0, 0.9}
    };

    /* Blue */
    struct intval_node blue_desc[] = {
        {0.0, 0.15},
        {0.365079, 0.444444},
        {1.0, 0.9},
    };

    struct color_desc red = { 3, red_desc };
    struct color_desc green = { 4, green_desc };
    struct color_desc blue = { 3, blue_desc };

    interpolate_col(&red, table, N, 0);
    interpolate_col(&green, table, N, 1);
    interpolate_col(&blue, table, N, 2);
}

void fill_cm_lookuptable_jet(unsigned int table[][3], const int N)
{
    /* Red */
    struct intval_node red_desc[] = {
        {0.00, 0.0},
        {0.35, 0.0},
        {0.66, 1.0},
        {0.89, 1.0},
        {1.00, 0.5}
    };

    /* Green */
    struct intval_node green_desc[] = {
        {0.00, 0.0},
        {0.125, 0.0},
        {0.375, 1.0},
        {0.64, 1.0},
        {0.91, 0.0},
        {1.00, 0.0}
    };

    /* Blue */
    struct intval_node blue_desc[] = {
        {0.00, 0.5},
        {0.11, 1.0},
        {0.34, 1.0},
        {0.65, 0.0},
        {1.00, 0.0}
    };

    struct color_desc red = { 5, red_desc };
    struct color_desc green = { 6, green_desc };
    struct color_desc blue = { 5, blue_desc };

    interpolate_col(&red, table, N, 0);
    interpolate_col(&green, table, N, 1);
    interpolate_col(&blue, table, N, 2);
}


void interpolate_col(const struct color_desc *color, unsigned int table[][3],
                     const int N, const int col_index)
{
    double x_l, x_r;            /* endpoints position interval */
    int ind_lo, ind_hi;         /* indices for subintervals */
    double val_l, val_r;        /* endpoints value interval */


    for (int l = 0; l < color->n_subintervals; l++) {
        x_l = color->nodes[l].position;
        x_r = color->nodes[l + 1].position;
        ind_lo = (int) rint(N * x_l);
        ind_hi = (int) (rint(N * x_r) + 1.);
        val_l = color->nodes[l].value;
        val_r = color->nodes[l + 1].value;
        if (val_l == val_r) {
            for (int i = ind_lo; i < ind_hi; i++)
                table[i][col_index] = (unsigned int) rint(N * val_l);
        } else {
            /* Interpolate linearly */
            for (int i = ind_lo; i < ind_hi; i++) {
                table[i][col_index] =
                    (unsigned int) rint(N * (val_l + (val_r - val_l)
                                             * (i - (double) ind_lo) / (ind_hi -
                                                                        (double)
                                                                        ind_lo)));
            }
        }
    }
}
