#include "CONST.hpp"

void bit_reverse_short(short *x_r, short *x_i);
void fft_short(short *x_r, short *x_i);
void bit_reverse_short_pointer(short *x_r, short *x_i);
void fft_short_pointer(short *x_r, short *x_i);
void create_table_short();
short add_sin_short(int n);
short add_cos_short(int n);

double add_sin(int n);
double add_cos(int n);
void create_table();
void fft(double x_r[N], double x_i[N]);
void bit_reverse(double x_r, double x_i);
void fft_pointer(double *x_r, double *x_i);
void bit_reverse_pointer(double *x_r, double *x_i);
