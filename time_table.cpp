#include <chrono>
#include <iostream> // for cout

#define N 64            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数
#define phi (double)0   // 初期位相

using namespace std;    // cout, endl, swap, ios, complex
using namespace chrono; // system_clock, duration_cast, microseconds, ofstream

double sin_table[] = {
    0,       0.0980171, 0.19509,  0.290285, 0.382683, 0.471397, 0.55557,  0.634393, 0.707107,
    0.77301, 0.83147,   0.881921, 0.92388,  0.95694,  0.980785, 0.995185, 1,
};

double Sin(int i, int m) {
    int n = i % m;
    if (n <= m / 4) {
        return sin_table[n];
    } else if (n <= m / 2) {
        return sin_table[m / 2 - n];
    } else if (n <= 3 * m / 4) {
        return -sin_table[n - m / 2];
    } else if (n <= m) {
        return -sin_table[m - n];
    } else {
        printf("Error!");
        return -1;
    }
}

double Cos(int i, int m) {
    i += m / 4;
    return Sin(i, m);
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

// 高速フーリエ変換(テーブル無し)
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
                x_r[i * m + j + m / 2] = (a_r - b_r) * cos(2 * M_PI / m * j) + (a_i - b_i) * sin(2 * M_PI / m * j);
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-sin(2 * M_PI / m * j)) + (a_i - b_i) * cos(2 * M_PI / m * j);
            }
        }
        m /= 2;
    }
    bit_reverse(x_r, x_i);
}

// 高速フーリエ変換(テーブル有り)
void use_table_fft(double x_r[N], double x_i[N]) {
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
                x_r[i * m + j + m / 2] = (a_r - b_r) * Cos(j, m) + (a_i - b_i) * Sin(j, m);
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-Sin(j, m)) + (a_i - b_i) * Cos(j, m);
            }
        }
        m /= 2;
    }
    bit_reverse(x_r, x_i);
}

int main() {
    double x_r[N], x_i[N];
    double sum = 0;
    system_clock::time_point start, end;

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi); // t = i / Fs
        x_i[i] = 0;
    }

    start = system_clock::now();      // 計測スタート時刻を保存
    for (int i = 0; i < N * N; i++) { // N回繰り返して平均を求める
        fft(x_r, x_i);
    }
    end = system_clock::now(); // 計測終了時刻を保存
    // 要した時間を計算
    double fft_time = static_cast<double>(duration_cast<microseconds>(end - start).count() / 1000.0) / (N * N);

    start = system_clock::now();
    for (int i = 0; i < N * N; i++) {
        use_table_fft(x_r, x_i);
    }
    end = system_clock::now();
    double use_table_fft_time = static_cast<double>(duration_cast<microseconds>(end - start).count() / 1000.0) / (N * N);

    // 要した時間をミリ秒（1/1000秒）に変換して表示
    cout << "FFT(テーブル無し): " << fft_time << " ms" << endl;
    cout << "FFT(テーブル有り): " << use_table_fft_time << " ms" << endl;

    return 0;
}
