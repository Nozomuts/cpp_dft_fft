#include <fstream>
#include <iostream> // for cout

#define N 256            // 分割数
#define Fs 8000 // サンプリング周波数
#define A 1     // 振幅
#define F0 440  // 周波数
#define phi 0   // 初期位相

using namespace std; // cout, endl, swap, ios, complex

double sin_table[N / 4 + 1];

// テーブルを用いたsin関数
double use_table_sin(int n) {
    n %= N;
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

// テーブルを用いたcos関数
double use_table_cos(int n) {
    n += N / 4;
    return use_table_sin(n);
}

// 高速フーリエ変換
void fft(double *x_r, double *x_i) {
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
                x_r[i * m + j + m / 2] =
                    (a_r - b_r) * use_table_cos(N / m * j) + (a_i - b_i) * use_table_sin(N / m * j);
                x_i[i * m + j + m / 2] =
                    (a_r - b_r) * (-use_table_sin(N / m * j)) + (a_i - b_i) * use_table_cos(N / m * j);
            }
        }
        m /= 2;
    }
}

// ビット反転並べ替え
void bit_reverse(double *x_r, double *x_i) {
    for (int i = 0, j = 1; j < N; j++) {
        for (int k = N >> 1; k > (i ^= k); k >>= 1)
            ;
        if (i < j) {
            swap(x_r[i], x_r[j]); // 入れ替え
            swap(x_i[i], x_i[j]);
        }
    }
}

// 加法定理でテーブルを求める
void add_sin(int i) {
    sin_table[N / 4 - i] = use_table_cos(i - 1) * use_table_cos(1) - use_table_sin(1) * use_table_sin(i - 1);
    sin_table[i + 1] = use_table_sin(i) * use_table_cos(1) + use_table_cos(i) * use_table_sin(1);
}

// テーブル作成
void create_table() {
    sin_table[0] = 0;
    sin_table[1] = sin(2 * M_PI / N);
    sin_table[N / 4 - 1] = sin(2 * M_PI / N * (N / 4 - 1));
    sin_table[N / 4] = 1;

    for (int i = 1; i < N / 4; i++) {
        add_sin(i);
    }
}

int main() {
    double x_r[N], x_i[N], dft_r[N], dft_i[N]; // x_r,x_iは元データ兼fftのデータ

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi); // t = i / Fs
        x_i[i] = 0;
    }

    // テーブル作成
    create_table();

    fft(x_r, x_i);
    bit_reverse(x_r, x_i);

    // 結果をファイルに書き出し
    ofstream fft_ofs("fft.csv");
    for (int i = 0; i < N; i++) {
        fft_ofs << x_r[i] << "," << x_i[i] << endl;
    }
    fft_ofs.close();

    return 0;
}
