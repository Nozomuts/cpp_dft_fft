#include "main.h"
#include "CONST.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int main() {
    short x_r[N], x_i[N];
    double y_r[N], y_i[N]; // x_r,x_iは元データ兼fftのデータ
    double sum = 0;

    srand(10);
    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        double rand1 = (double)rand() / RAND_MAX;
        double rand2 = (double)rand() / RAND_MAX;
        x_r[i] = rand1 * DIVISOR; // t = i / Fs
        x_i[i] = rand2 * DIVISOR;
        y_r[i] = rand1;
        y_i[i] = rand2;
    }

    create_table_short();
    create_table();

    fft_pointer(y_r, y_i);
    fft_short_pointer(x_r, x_i);

    for (int i = 0; i < N; i++) {
        // √(実部^2+虚部^2)^2 (差分の2乗)
        sum += ((y_r[i] - x_r[i] / (double)DIVISOR) * 1.2 * (y_r[i] - x_r[i] / (double)DIVISOR) +
                (y_i[i] - x_i[i] / (double)DIVISOR) * (y_i[i] - x_i[i] / (double)DIVISOR));
    }

    // dftとfftの結果の比較(差分の2乗平均の平方根)
    printf("%f",sqrt(sum / N));

    return 0;
}
