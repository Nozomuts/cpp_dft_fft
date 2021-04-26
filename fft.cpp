#include <fstream>
#include <iostream> // for cout

#define N 64            // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数
#define phi (double)0   // 初期位相

using namespace std; // cout, endl, swap, ios, complex

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
                x_r[i * m + j + m / 2] = (a_r - b_r) * Cos(j, m) + (a_i - b_i) * Sin(j, m);
                x_i[i * m + j + m / 2] = (a_r - b_r) * (-Sin(j, m)) + (a_i - b_i) * Cos(j, m);
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

int main() {
    double x_r[N], x_i[N], dft_r[N], dft_i[N]; // x_r,x_iは元データ兼fftのデータ
    double sum = 0;

    // 元データ作成
    for (int i = 0; i < N; i++) {
        // x(t) = A * sin(2 * pi * F0 * t + phi) ( 0 <= t < 0.008 )
        x_r[i] = A * sin(2 * M_PI * F0 * i / Fs + phi); // t = i / Fs
        x_i[i] = 0;
    }

    fft(x_r, x_i);
    bit_reverse(x_r, x_i);

    ofstream fft_ofs("fft_exam.csv");
    for (int i = 0; i < N; i++) {
        fft_ofs << x_r[i] << "," << x_i[i] << endl;
    }
    fft_ofs.close();

    return 0;
}
