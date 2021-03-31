#include <fstream>
#include <iostream>
#include <string>

using namespace std; // string, ifstream, ios, cout, endl

void ifs_sum(string file_name) {
    ifstream ifs;
    string line;
    int n = 0;
    double sum = 0;
    ifs.open(file_name, ios::in);
    while (!ifs.eof()) {
        getline(ifs, line);
        sum += stod(line);
        n++;
    }
    ifs.close();
    cout << file_name << ": " << sum / n << endl;
}

int main() {
    ifs_sum("dft.txt");
    ifs_sum("fft.txt");

    return 0;
}
