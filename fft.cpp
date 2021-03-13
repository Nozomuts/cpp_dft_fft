#include <algorithm> // swap
#include <complex>
#include <iostream>

using namespace std;

typedef complex<double> Complex;
typedef unsigned int uint;

static Complex one(1.0, 0.0);
static Complex ione(0.0, 1.0);

class FFT {
  private:
    Complex *array;
    size_t size; // arrayのサイズ
    size_t use;  // コンストラクタで与えたサイズ

    bool verbose; // 配列ダンプのモード指定に使う
    size_t nextPow2(size_t s);
    void bitReverse();

  public:
    FFT(size_t s) : use(s), size(nextPow2(s)), verbose(false) { array = new Complex[size]; };
    ~FFT() { delete[] array; };

    Complex &operator[](int index) const {
        if (index < 0)
            return array[use - index];
        return array[index];
    };

    bool setVerbose(bool b) { return verbose = b; };
    void dump();

    void fft(bool isReverse = false);
};

size_t FFT::nextPow2(size_t s)
// s以上の最小の2のべき乗を返す
{
    size_t n = 1;
    while (n < s)
        n <<= 1;
    return n;
}

void FFT::bitReverse() {
    uint k, b, a;
    for (uint i = 0; i < size; i++) {
        k = 0;
        b = size >> 1;
        a = 1;
        while (b >= a) {
            if (b & i)
                k |= a;
            if (a & i)
                k |= b;
            b >>= 1;
            a <<= 1;
        }
        if (i < k)
            swap(array[i], array[k]);
    }
}

void FFT::dump() {
    uint end = verbose ? size : use;
    for (uint i = 0; i < end; i++)
        cout << array[i] << " ";
    cout << endl;
}

void FFT::fft(bool isReverse) {
    bitReverse();
    size_t m = 2;
    Complex w, ww, t;

    while (m <= size) {
        double arg = -2.0 * M_PI / m;
        w = Complex(cos(arg), sin(arg));
        if (isReverse)
            w = one / w; //-1乗 -(-2.0*PI/size) = 2.0*PI/size

        for (uint i = 0; i < size; i += m) {
            ww = 1.0;
            for (uint j = 0; j < m / 2; j++) {
                int a = i + j;
                int b = i + j + m / 2;

                t = ww * array[b];

                array[b] = array[a] - t;
                array[a] = array[a] + t;

                ww *= w;
            }
        }
        m *= 2;
    }
}

#define N (8)
#define RND ((double)rand() / (double)RAND_MAX)
int main() {
    srand(time(NULL));
    FFT box(N);

    for (int i = 0; i < N; i++)
        box[i] = Complex(RND, RND);
    box.setVerbose(true);

    Complex before[N], transform[N], itransform[N];
    // fft前
    for (int i = 0; i < N; i++)
        before[i] = box[i];
    // fft後
    box.fft();
    for (int i = 0; i < N; i++)
        transform[i] = box[i];
    cout << "before\t\t\ttransform" << endl;
    for (int i = 0; i < N; i++) {
        printf("(%9lf, %9lf) -> (%9lf, %9lf)\n", before[i].real(), before[i].imag(), transform[i].real(),
               transform[i].imag());
    }
    return 0;
}
