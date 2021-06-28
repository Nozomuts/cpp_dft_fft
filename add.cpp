#include <chrono>
#include <fstream>
#include <iostream> // for cout

#define N 64 // 分割数
#define M 2
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数
#define phi (double)0   // 初期位相

using namespace std;    // cout, endl, swap, ios, complex
using namespace chrono; // system_clock, duration_cast, microseconds, ofstream

double sin_table[N / 4 / M + 1];
double table[2];
double sini, cosi, tmp;

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

double add_cos(int n) {
    n += N / 4;
    return add_sin(n);
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
                x_r[i * m + j + m / 2] = (a_r - b_r) * add_cos(N / m * j) + (a_i - b_i) * add_sin(N / m * j);
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-add_sin(N / m * j)) + (a_i - b_i) * add_cos(N / m * j);
            }
        }
        m /= 2;
    }
    bit_reverse(x_r, x_i);
}

void create_table() {
    table[0] = sin(2 * M_PI / N);
    table[1] = cos(2 * M_PI / N);
    for (int i = 0; i <= N / 4 / M; i++) {
        sin_table[i] = sin(2 * M_PI / N * i * M);
    }
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

    fft(x_r, x_i);
    // 結果をファイルに書き出し
    ofstream fft_ofs("add.csv");
    for (int i = 0; i < N; i++) {
        fft_ofs << x_r[i] << "," << x_i[i] << endl;
        cout << add_sin(i) << endl;
    }
    fft_ofs.close();

    system_clock::time_point start, end;

    start = system_clock::now();    // 計測スタート時刻を保存
    for (int i = 0; i < 100; i++) { // N回繰り返して平均を求める
        fft(x_r, x_i);
    }
    end = system_clock::now(); // 計測終了時刻を保存
    // 要した時間を計算
    double fft_time = static_cast<double>(duration_cast<microseconds>(end - start).count() / 1000.0) / 100;

    cout << fft_time << "ms" << endl;

    return 0;
}
