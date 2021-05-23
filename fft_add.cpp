#include <fstream>
#include <iostream> // for cout

#define N 64            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数
#define phi (double)0   // 初期位相

using namespace std; // cout, endl, swap, ios, complex

double sin_table[N / 4 + 1];
double add_sin_table[4];
double sini, cosi, tmp;

double use_table_add_sin(int n) {
    n %= N;
    sini = 0;
    cosi = 1;
    tmp = 0;

    if (n == 0) {
        return add_sin_table[0];
    } else if (n == 1) {
        return add_sin_table[1];
    } else if (n == N / 4 - 1) {
        return add_sin_table[2];
    } else if (n == N / 4) {
        return add_sin_table[3];
    } else if (n <= N / 4) {
        for (int i = 0; i < n; i++) {
            tmp = sini * add_sin_table[2] + cosi * add_sin_table[1];
            cosi = cosi * add_sin_table[2] - sini * add_sin_table[1];
            sini = tmp;
        }
        return sini;
    } else if (n <= N / 2) {
        for (int i = 0; i < N / 2 - n; i++) {
            tmp = sini * add_sin_table[2] + cosi * add_sin_table[1];
            cosi = cosi * add_sin_table[2] - sini * add_sin_table[1];
            sini = tmp;
        }
        return sini;
    } else if (n <= 3 * N / 4) {
        for (int i = 0; i < n - N / 2; i++) {
            tmp = sini * add_sin_table[2] + cosi * add_sin_table[1];
            cosi = cosi * add_sin_table[2] - sini * add_sin_table[1];
            sini = tmp;
        }
        return -sini;
    } else if (n <= N) {
        for (int i = 0; i < N - n; i++) {
            tmp = sini * add_sin_table[2] + cosi * add_sin_table[1];
            cosi = cosi * add_sin_table[2] - sini * add_sin_table[1];
            sini = tmp;
        }
        return -sini;
    } else {
        printf("Error!");
        return -1;
    }
}

double use_table_add_cos(int n) {
    n += N / 4;
    return use_table_add_sin(n);
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

// 高速フーリエ変換(加法定理あり)
void add_sin_fft(double x_r[N], double x_i[N]) {
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
                    (a_r - b_r) * use_table_add_cos(N / m * j) + (a_i - b_i) * use_table_add_sin(N / m * j);
                x_i[i * m + j + m / 2] =
                    (a_r - b_r) * (-use_table_add_sin(N / m * j)) + (a_i - b_i) * use_table_add_cos(N / m * j);
            }
        }
        m /= 2;
    }
    bit_reverse(x_r, x_i);
}

void create_table() {
    add_sin_table[0] = 0;
    add_sin_table[1] = sin(2 * M_PI / N);
    add_sin_table[2] = sin(2 * M_PI / N * (N / 4 - 1));
    add_sin_table[3] = 1;
}

int main() {
    double x_r[N], x_i[N];

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi); // t = i / Fs
        x_i[i] = 0;
    }

    create_table();

    add_sin_fft(x_r, x_i);
    // 結果をファイルに書き出し
    ofstream fft_ofs("fft_add.csv");
    for (int i = 0; i < N; i++) {
        fft_ofs << x_r[i] << "," << x_i[i] << endl;
        cout << use_table_add_sin(i) << endl;
    }
    fft_ofs.close();

    return 0;
}