#include "CONST.hpp"
#include <iostream>

using namespace std;

int sin_table_int[N / 4 + 1];

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
                cout << (a_r - b_r) * add_cos_short(N / m * j) << endl;
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
                cout << (a_r - b_r) * add_cos_short(N / m * j) << endl;
            }
        }
        m /= 2;
    }
    bit_reverse_short_pointer(x_r, x_i);
}

short add_sin_short(int i) {
    int n = i % N;
    if (n <= N / 4) {
        return sin_table_int[n];
    } else if (n <= N / 2) {
        return sin_table_int[N / 2 - n];
    } else if (n <= 3 * N / 4) {
        return -sin_table_int[n - N / 2];
    } else if (n <= N) {
        return -sin_table_int[N - n];
    } else {
        printf("Error!");
        return -1;
    }
}

short add_cos_short(int i) {
    i += N / 4;
    return add_sin_short(i);
}

// テーブル作成
void create_table_short() {
    sin_table_int[0] = 0;

    for (int i = 1; i <= N / 4; i++) {
        sin_table_int[i] = sin(2 * M_PI / N * i) * DIVISOR;
    }
}
