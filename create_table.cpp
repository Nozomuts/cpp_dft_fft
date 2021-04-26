#include <chrono>
#include <fstream>
#include <iostream> // for cout

#define N 64 // 分割数

using namespace std; // cout, endl, swap, ios, complex

double sin_table[N / 4 + 1];

double use_table_sin(int n) {
    if (n > N / 4) {
        n %= N;
    }
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
    cout << sin_table[i + 1] << endl;
}

int main() {
    sin_table[0] = 0;
    sin_table[1] = sin(2 * M_PI / N);
    sin_table[N / 4 - 1] = sin(2 * M_PI / N * (N / 4 - 1));
    sin_table[N / 4] = 1;

    for (int i = 1; i < N / 4; i++) {
        add_sin(i);
    }

    ofstream sin_ofs("sin.csv");
    for (int i = 0; i <= N / 4; i++) {
        sin_ofs << sin_table[i] << endl;
    }
    sin_ofs.close();

    return 0;
}
