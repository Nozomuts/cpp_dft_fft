#include <chrono>
#include <fstream>
#include <iostream> // for cout

#define N 64    // 分割数
#define Fs 8000 // サンプリング周波数
#define A 1     // 振幅
#define F0 440  // 周波数
#define phi 0   // 初期位相

using namespace std; // cout, endl, swap, ios, complex

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

void dft(int x_r[N], int x_i[N], int *dft_r, int *dft_i) {
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            dft_r[k] += x_r[n] * use_table_cos(k * n) + x_i[n] * use_table_sin(k * n);
            dft_i[k] += x_r[n] * (-use_table_sin(k * n)) + x_i[n] * use_table_cos(k * n);
        }
    }
}

// 加法定理でテーブルを求める
// void add_sin(int i) {
//     sin_table[N / 4 - i] = use_table_cos(i - 1) * use_table_cos(1) - use_table_sin(1) * use_table_sin(i - 1);
//     sin_table[i + 1] = use_table_sin(i) * use_table_cos(1) + use_table_cos(i) * use_table_sin(1);
// }

// テーブル作成
void create_table() {
    sin_table[0] = 0;

    for (int i = 1; i <= N / 4; i++) {
        sin_table[i] = sin(2 * M_PI / N * i) * 100000;
    }
}

int main() {
    int x_r[N], x_i[N], dft_r[N], dft_i[N]; // x_r,x_iは元データ兼fftのデータ

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi) * 1000; // t = i / Fs
        x_i[i] = 0;
    }

    create_table();
    for (int i = 0; i < N; i++)
    {
       cout << use_table_sin(i) << endl;
    }
    

    dft(x_r, x_i, dft_r, dft_i);

    // 結果をファイルに出力
    ofstream dft_ofs("dft_int.csv");
    for (int i = 0; i < N; i++) {
        dft_ofs << dft_r[i] << "," << dft_i[i] << endl;
    }
    dft_ofs.close();

    return 0;
}
