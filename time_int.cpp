#include <chrono>
#include <fstream>
#include <iostream> // for cout

#define N 64   // 分割数
#define Fs 8000 // サンプリング周波数
#define A 1     // 振幅
#define F0 440  // 周波数
#define phi 0   // 初期位相
#define M 10000

using namespace std;    // cout, endl, swap, ios, complex
using namespace chrono; // system_clock, duration_cast, microseconds, ofstream

int sin_table[N / 4 + 1];

int use_table_sin(int i) {
    int n = i % N;
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

int use_table_cos(int i) {
    i += N / 4;
    return use_table_sin(i);
}

// テーブル作成
void create_table() {
    sin_table[0] = 0;

    for (int i = 1; i <= N / 4; i++) {
        sin_table[i] = sin(2 * M_PI / N * i) * M;
    }
}

// ビット反転並べ替え
void bit_reverse(int x_r[N], int x_i[N]) {
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
void fft(int x_r[N], int x_i[N]) {
    int m = N;
    while (m > 1) {
        for (int i = 0; i < N / m; i++) {
            for (int j = 0; j < m / 2; j++) {
                int a_r = x_r[i * m + j];
                int a_i = x_i[i * m + j];
                int b_r = x_r[i * m + j + m / 2];
                int b_i = x_i[i * m + j + m / 2];
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
    bit_reverse(x_r, x_i);
}

int main() {
    int x_r[N], x_i[N];
    system_clock::time_point start, end;

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        double rand1 = (double)rand() / RAND_MAX;
        double rand2 = (double)rand() / RAND_MAX;

        x_r[i] = rand1 * M; // t = i / Fs
        x_i[i] = rand2 * M;
    }

    create_table();

    start = system_clock::now();    // 計測スタート時刻を保存
    for (int i = 0; i < 1000; i++) { // N回繰り返して平均を求める
        fft(x_r, x_i);
    }
    end = system_clock::now(); // 計測終了時刻を保存
    // 要した時間を計算
    double fft_time = static_cast<double>(duration_cast<microseconds>(end - start).count() / 1000.0) / 1000;

    cout << fft_time << "ms" << endl;

    return 0;
}
