#include <chrono>
#include <fstream>
#include <iostream> // for cout

#define N 64 // 分割数

using namespace std; // cout, endl, swap, ios, complex

int main() {

    ofstream sin_ofs("sin.csv");
    for (int i = 0; i <= N / 4; i++) {
        sin_ofs << sin(2 * M_PI / N * i) << "," << endl;
    }
    sin_ofs.close();

    return 0;
}
