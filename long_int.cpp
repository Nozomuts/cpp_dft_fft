#include <chrono>
#include <fstream>
#include <iostream> // for cout

#define N 64    // 分割数
#define Fs 8000 // サンプリング周波数
#define A 1     // 振幅
#define F0 440  // 周波数
#define phi 0   // 初期位相
#define M 100

using namespace std;    // cout, endl, swap, ios, complex
using namespace chrono; // system_clock, duration_cast, microseconds, ofstream

short sin_table[N / 4 + 1];

short use_table_sin(short i) {
    short n = i % N;
    if (n <= N / 4) {
        return sin_table[n];
    } else if (n <= N / 2) {
        return sin_table[N / 2 - n];
    } else if (n <= 3 * N / 4) {
        return -sin_table[n - N / 2];
    } else if (n <= N) {
        return -sin_table[N - n];
    } else {
        printf("Error!");
        return -1;
    }
}

short use_table_cos(short i) {
    i += N / 4;
    return use_table_sin(i);
}

// テーブル作成
void create_table() {
    sin_table[0] = 0;

    for (short i = 1; i <= N / 4; i++) {
        sin_table[i] = sin(2 * M_PI / N * i) * M;
    }
}

// 高速フーリエ変換
void fft(short *x_r, short *x_i) {
    short m = N;
    while (m > 1) {
        for (short i = 0; i < N / m; i++) {
            for (short j = 0; j < m / 2; j++) {
                short a_r = x_r[i * m + j];
                short a_i = x_i[i * m + j];
                short b_r = x_r[i * m + j + m / 2];
                short b_i = x_i[i * m + j + m / 2];
                x_r[i * m + j] = a_r + b_r;
                x_i[i * m + j] = a_i + b_i;
                x_r[i * m + j + m / 2] =
                    (a_r - b_r) * use_table_cos(N / m * j) / M + (a_i - b_i) * use_table_sin(N / m * j) / M;
                x_i[i * m + j + m / 2] =
                    (a_r - b_r) * (-use_table_sin(N / m * j)) / M + (a_i - b_i) * use_table_cos(N / m * j) / M;
            }
        }
        m /= 2;
    }
}

// ビット反転並べ替え
void bit_reverse(short *x_r, short *x_i) {
    for (short i = 0, j = 1; j < N; j++) {
        for (short k = N >> 1; k > (i ^= k); k >>= 1)
            ;
        if (i < j) {
            swap(x_r[i], x_r[j]); // 入れ替え
            swap(x_i[i], x_i[j]);
        }
    }
}

int main() {
    short x_r[N], x_i[N];

    // 元データ作成
    for (short i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = (double)rand() / RAND_MAX * M; // t = i / Fs
        x_i[i] = (double)rand() / RAND_MAX * M;
    }

    create_table();

    fft(x_r, x_i);
    bit_reverse(x_r, x_i);

    // 結果をファイルに出力
    ofstream fft_ofs("long_int.csv");
    for (short i = 0; i < N; i++) {
        fft_ofs << x_r[i] / (double)M << "," << x_i[i] / (double)M << endl;
    }
    fft_ofs.close();

    return 0;
}
