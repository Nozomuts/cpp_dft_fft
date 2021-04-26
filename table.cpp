#include <chrono>
#include <fstream>
#include <iostream> // for cout

#define N 64 // 分割数

using namespace std; // cout, endl, swap, ios, complex

double sin_table[] = {
    0,       0.0980171, 0.19509,  0.290285, 0.382683, 0.471397, 0.55557,  0.634393, 0.707107,
    0.77301, 0.83147,   0.881921, 0.92388,  0.95694,  0.980785, 0.995185, 1,
};

double Sin(int i) {
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

double Cos(int i) {
    i += N / 4;
    return Sin(i);
}

int main() {
    ofstream sin_ofs("sin.csv");
    for (int i = 0; i < N; i++) {
        sin_ofs << Sin(i) << "," << endl;
    }
    sin_ofs.close();

    ofstream cos_ofs("cos_2.csv");
    for (int i = 0; i < N; i++) {
        cos_ofs << Cos(i) << "," << endl;
    }
    cos_ofs.close();

    return 0;
}
