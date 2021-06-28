#include "CONST.h"
#include <iostream>

using namespace std;

double sin_table[N / 4 + 1];
double sin_table[N / 4 / M + 1];
double table[2];
double sini, cosi, tmp;

// テーブルを用いたsin関数
double add_sin(int n) {
    n %= N;
    tmp = 0;

    if (n == 0) {
        return 0;
    } else if (n == 1 || n == N / 2 - 1) {
        return table[0];
    } else if (n == N / 4 - 1 || n == N / 2 - 1) {
        return table[1];
    } else if (n == N / 2 + 1 || n == N - 1) {
        return -table[0];
    } else if (n == 3 * N / 4 - 1) {
        return -table[1];
    } else if (n == N / 4) {
        return 1;
    } else if (n == 3 * N / 4) {
        return -1;
    } else if (n <= N / 4) {
        if (n % M == 0) {
            return sin_table[n / M];
        }
        sini = sin_table[(n - n % M) / M];
        cosi = sin_table[(N / 2 - (n + N / 4) + n % M) / M];
        for (int i = 0; i < n % M; i++) {
            tmp = sini * table[1] + cosi * table[0];
            cosi = cosi * table[1] - sini * table[0];
            sini = tmp;
        }
        return sini;
    } else if (n <= N / 2) {
        if (n % M == 0) {
            return sin_table[(N / 2 - n) / M];
        }
        sini = sin_table[(N / 2 - n + n % M) / M];
        cosi = -sin_table[((n + N / 4) - N / 2 - n % M) / M];
        for (int i = 0; i < n % M; i++) {
            tmp = sini * table[1] + cosi * table[0];
            cosi = cosi * table[1] - sini * table[0];
            sini = tmp;
        }
        return sini;
    } else if (n <= 3 * N / 4) {
        if (n % M == 0) {
            return -sin_table[(n - N / 2) / M];
        }
        sini = -sin_table[(n - N / 2 - n % M) / M];
        cosi = -sin_table[(N - (n + N / 4) + n % M) / M];
        for (int i = 0; i < n % M; i++) {
            tmp = sini * table[1] + cosi * table[0];
            cosi = cosi * table[1] - sini * table[0];
            sini = tmp;
        }
        return sini;
    } else if (n <= N) {
        if (n % M == 0) {
            return -sin_table[(N - n) / M];
        }
        sini = -sin_table[(N - n + n % M) / M];
        cosi = sin_table[(((n + N / 4 - n % M) - N) / M)];
        for (int i = 0; i < n % M; i++) {
            tmp = sini * table[1] + cosi * table[0];
            cosi = cosi * table[1] - sini * table[0];
            sini = tmp;
        }
        return sini;
    } else {
        printf("Error!");
        return -1;
    }
}

// テーブルを用いたcos関数
double add_cos(int n) {
    n += N / 4;
    return add_sin(n);
}

// テーブル作成
void create_table() {
    sin_table[0] = 0;

    for (int i = 1; i <= N / 4; i++) {
        sin_table[i] = sin(2 * M_PI / N * i);
    }
}

// ビット反転並べ替え
void bit_reverse_pointer(double *x_r, double *x_i) {
    for (int i = 0, j = 1; j < N; j++) {
        for (int k = N >> 1; k > (i ^= k); k >>= 1)
            ;
        if (i < j) {
            swap(x_r[i], x_r[j]); // 入れ替え
            swap(x_i[i], x_i[j]);
        }
    }
}

// ビット反転並べ替え
void bit_reverse(double x_r[N], double x_i[N]) {
    for (int i = 0, j = 1; j < N - 1; j++) {
        for (int k = N >> 1; k > (i ^= k); k >>= 1)
            ;
        if (i < j) {
            swap(x_r[i], x_r[j]); // 入れ替え
            swap(x_i[i], x_i[j]);
        }
    }
}

void fft(double x_r[N], double x_i[N]) {
    int m = N;
    while (m > 1) {
        for (int i = 0; i < N / m; i++) {
            for (int j = 0; j < m / 2; j++) {
                double a_r = x_r[i * m + j];
                double a_i = x_i[i * m + j];
                double b_r = x_r[i * m + j + m / 2];
                double b_i = x_i[i * m + j + m / 2];
                x_r[i * m + j] = a_r + b_r;
                x_i[i * m + j] = a_i + b_i;
                x_r[i * m + j + m / 2] = (a_r - b_r) * add_cos(N / m * j) + (a_i - b_i) * add_sin(N / m * j);
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-add_sin(N / m * j)) + (a_i - b_i) * add_cos(N / m * j);
            }
        }
        m /= 2;
    }
    bit_reverse(x_r, x_i);
}

void fft_pointer(double *x_r, double *x_i) {
    int m = N;
    while (m > 1) {
        for (int i = 0; i < N / m; i++) {
            for (int j = 0; j < m / 2; j++) {
                double a_r = x_r[i * m + j];
                double a_i = x_i[i * m + j];
                double b_r = x_r[i * m + j + m / 2];
                double b_i = x_i[i * m + j + m / 2];
                x_r[i * m + j] = a_r + b_r;
                x_i[i * m + j] = a_i + b_i;
                x_r[i * m + j + m / 2] = (a_r - b_r) * add_cos(N / m * j) + (a_i - b_i) * add_sin(N / m * j);
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-add_sin(N / m * j)) + (a_i - b_i) * add_cos(N / m * j);
            }
        }
        m /= 2;
    }
    bit_reverse_pointer(x_r, x_i);
}
