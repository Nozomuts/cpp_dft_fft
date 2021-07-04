#include "CONST.hpp"
#include <iostream>
#include <math.h>

using namespace std;

short sin_table_short[N / 4 / M + 1];
short table_short[2];
short sini_short, cosi_short, tmp_short;

short add_sin_short(int n) {
    n %= N;
    tmp_short = 0;

    if (n == 0) {
        return 0;
    } else if (n == 1 || n == N / 2 - 1) {
        return table_short[0];
    } else if (n == N / 4 - 1 || n == N / 2 - 1) {
        return table_short[1];
    } else if (n == N / 2 + 1 || n == N - 1) {
        return -table_short[0];
    } else if (n == 3 * N / 4 - 1) {
        return -table_short[1];
    } else if (n == N / 4) {
        return 1;
    } else if (n == 3 * N / 4) {
        return -1;
    } else if (n <= N / 4) {
        if (n % M == 0) {
            return sin_table_short[n / M];
        }
        sini_short = sin_table_short[(n - n % M) / M];
        cosi_short = sin_table_short[(N / 2 - (n + N / 4) + n % M) / M];
        for (int i = 0; i < n % M; i++) {
            tmp_short = sini_short * table_short[1] / DIVISOR + cosi_short * table_short[0] / DIVISOR;
            cosi_short = cosi_short * table_short[1] / DIVISOR - sini_short * table_short[0] / DIVISOR;
            sini_short = tmp_short;
        }
        return sini_short;
    } else if (n <= N / 2) {
        if (n % M == 0) {
            return sin_table_short[(N / 2 - n) / M];
        }
        sini_short = sin_table_short[(N / 2 - n + n % M) / M];
        cosi_short = -sin_table_short[((n + N / 4) - N / 2 - n % M) / M];
        for (int i = 0; i < n % M; i++) {
            tmp_short = sini_short * table_short[1] / DIVISOR + cosi_short * table_short[0] / DIVISOR;
            cosi_short = cosi_short * table_short[1] / DIVISOR - sini_short * table_short[0] / DIVISOR;
            sini_short = tmp_short;
        }
        return sini_short;
    } else if (n <= 3 * N / 4) {
        if (n % M == 0) {
            return -sin_table_short[(n - N / 2) / M];
        }
        sini_short = -sin_table_short[(n - N / 2 - n % M) / M];
        cosi_short = -sin_table_short[(N - (n + N / 4) + n % M) / M];
        for (int i = 0; i < n % M; i++) {
            tmp_short = sini_short * table_short[1] / DIVISOR + cosi_short * table_short[0] / DIVISOR;
            cosi_short = cosi_short * table_short[1] / DIVISOR - sini_short * table_short[0] / DIVISOR;
            sini_short = tmp_short;
        }
        return sini_short;
    } else if (n <= N) {
        if (n % M == 0) {
            return -sin_table_short[(N - n) / M];
        }
        sini_short = -sin_table_short[(N - n + n % M) / M];
        cosi_short = sin_table_short[(((n + N / 4 - n % M) - N) / M)];
        for (int i = 0; i < n % M; i++) {
            tmp_short = sini_short * table_short[1] / DIVISOR + cosi_short * table_short[0] / DIVISOR;
            cosi_short = cosi_short * table_short[1] / DIVISOR - sini_short * table_short[0] / DIVISOR;
            sini_short = tmp_short;
        }
        return sini_short;
    } else {
        printf("Error!");
        return -1;
    }
}

short add_cos_short(int n) {
    i += N / 4;
    return add_sin_short(i);
}

// ビット反転並べ替え
void bit_reverse_short(short x_r[N], short x_i[N]) {
    for (int i = 0, j = 1; j < N; j++) {
        for (int k = N >> 1; k > (i ^= k); k >>= 1)
            ;
        if (i < j) {
            swap(x_r[i], x_r[j]); // 入れ替え
            swap(x_i[i], x_i[j]);
        }
    }
}

// 高速フーリエ変換
void fft_short(short x_r[N], short x_i[N]) {
    int m = N;
    while (m > 1) {
        for (int i = 0; i < N / m; i++) {
            for (int j = 0; j < m / 2; j++) {
                short a_r = x_r[i * m + j];
                short a_i = x_i[i * m + j];
                short b_r = x_r[i * m + j + m / 2];
                short b_i = x_i[i * m + j + m / 2];
                x_r[i * m + j] = a_r + b_r;
                x_i[i * m + j] = a_i + b_i;
                x_r[i * m + j + m / 2] =
                    (a_r - b_r) * add_cos_short(N / m * j) / DIVISOR + (a_i - b_i) * add_sin_short(N / m * j) / DIVISOR;
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-add_sin_short(N / m * j)) / DIVISOR +
                                         (a_i - b_i) * add_cos_short(N / m * j) / DIVISOR;
            }
        }
        m /= 2;
    }
    bit_reverse_short(x_r, x_i);
}

// ビット反転並べ替え
void bit_reverse_short_pointer(short *x_r, short *x_i) {
    for (int i = 0, j = 1; j < N; j++) {
        for (int k = N >> 1; k > (i ^= k); k >>= 1)
            ;
        if (i < j) {
            swap(x_r[i], x_r[j]); // 入れ替え
            swap(x_i[i], x_i[j]);
        }
    }
}

// 高速フーリエ変換
void fft_short_pointer(short *x_r, short *x_i) {
    int m = N;
    while (m > 1) {
        for (int i = 0; i < N / m; i++) {
            for (int j = 0; j < m / 2; j++) {
                short a_r = x_r[i * m + j];
                short a_i = x_i[i * m + j];
                short b_r = x_r[i * m + j + m / 2];
                short b_i = x_i[i * m + j + m / 2];
                x_r[i * m + j] = a_r + b_r;
                x_i[i * m + j] = a_i + b_i;
                x_r[i * m + j + m / 2] =
                    (a_r - b_r) * add_cos_short(N / m * j) / DIVISOR + (a_i - b_i) * add_sin_short(N / m * j) / DIVISOR;
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-add_sin_short(N / m * j)) / DIVISOR +
                                         (a_i - b_i) * add_cos_short(N / m * j) / DIVISOR;
            }
        }
        m /= 2;
    }
    bit_reverse_short_pointer(x_r, x_i);
}

// テーブル作成
void create_table_short() {
    table_short[0] = sin(2 * M_PI / N) * DIVISOR;
    table_short[1] = cos(2 * M_PI / N) * DIVISOR;

    for (int i = 1; i <= N / 4 / M; i++) {
        sin_table_short[i] = sin(2 * M_PI / N * M * i) * DIVISOR;
    }
}
