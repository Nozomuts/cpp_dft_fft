#include <chrono>
#include <iostream> // for cout

#define N 9192          // 分割数
#define Fs (double)8000 // サンプリング周波数
#define A (double)1     // 振幅
#define F0 (double)440  // 周波数
#define phi (double)0   // 初期位相

using namespace std;    // cout, endl, swap, ios, complex
using namespace chrono; // system_clock, duration_cast, microseconds, ofstream

double sin_table[N / 4 + 1];

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

double use_table_cos(int n) {
    n += N / 4;
    return use_table_sin(n);
}

void add_sin(int i) {
    sin_table[N / 4 - i] = use_table_cos(i - 1) * use_table_cos(1) - use_table_sin(1) * use_table_sin(i - 1);
    sin_table[i + 1] = use_table_sin(i) * use_table_cos(1) + use_table_cos(i) * use_table_sin(1);
}

void create_table_add() {
    sin_table[0] = 0;
    sin_table[1] = sin(2 * M_PI / N);
    sin_table[N / 4 - 1] = sin(2 * M_PI / N * (N / 4 - 1));
    sin_table[N / 4] = 1;

    for (int i = 1; i < N / 4; i++) {
        add_sin(i);
    }
}

void create_table() {
    sin_table[0] = 0;

    for (int i = 1; i <= N / 4; i++) {
        sin_table[i] = sin(2 * M_PI / N * i);
    }
}

int main() {
    double x_r[N], x_i[N];
    double sum = 0;
    system_clock::time_point start, end;

    start = system_clock::now(); // 計測スタート時刻を保存
    for (int i = 0; i < 30; i++) {
        // テーブル作成
        create_table();
    }
    end = system_clock::now(); // 計測終了時刻を保存

    // 要した時間を計算
    double table_time = static_cast<double>(duration_cast<microseconds>(end - start).count() / 30) / N;

    start = system_clock::now(); // 計測スタート時刻を保存
    for (int i = 0; i < 30; i++) {
        // テーブル作成
    }
    create_table_add();
    end = system_clock::now(); // 計測終了時刻を保存

    // 要した時間を計算
    double add_table_time = static_cast<double>(duration_cast<microseconds>(end - start).count() / 30) / N;

    // 要した時間をミリ秒（1/1000秒）に変換して表示
    cout << "加法定理なし: " << table_time << " ms" << endl;
    cout << "加法定理あり: " << add_table_time << " ms" << endl;

    return 0;
}
