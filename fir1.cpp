#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
/*
1403-11-07  
Fateme Ordikhani - FaOrdi
This code  created for N=fir1(n,[a,b])  in matlab that Calculate the filter coefficients.  n is  Filter order and a,b are Bandpass filter, normalyze  frequencies  
output is the same with MATLAB. (checked!)

*/
using namespace std;

// Function to generate a Hamming window
vector<double> hamming_window(int N) {
    vector<double> window(N);
    for (int n = 0; n < N; ++n) {
        window[n] = 0.54 - 0.46 * cos(2 * M_PI * n / (N - 1));
    }
    return window;
}

// Function to calculate the sinc function for filter design
double sinc(double x) {
    if (x == 0) return 1.0;
    return sin(M_PI * x) / (M_PI * x);
}

// Function to design the FIR filter using window method
vector<double> fir1(int N, const vector<double>& Wn) {
    assert(N > 0);
    assert(Wn.size() == 2 && Wn[0] < Wn[1]);  // Expecting a bandpass filter

    int L = N + 1;
    vector<double> filter(L, 0.0);
    
    // Frequency vector for bandpass filter design
    double low = Wn[0];
    double high = Wn[1];
    
    // Design the filter coefficients (bandpass FIR filter using sinc)
    for (int n = 0; n < L; ++n) {
        double k = n - N / 2;
        if (k == 0) {
            filter[n] = high - low; // The center frequency is at the middle
        } else {
            filter[n] = (sin(M_PI * high * k) - sin(M_PI * low * k)) / (M_PI * k);
        }
    }

    // Apply the window function (Hamming window)
    vector<double> window = hamming_window(L);
    for (int n = 0; n < L; ++n) {
        filter[n] *= window[n];
    }

    return filter;
}

int main() {
    int N = 48;  // Filter order
    vector<double> Wn = {0.001, 0.5};  // Bandpass filter, frequencies between 0.001 and 0.5
    
    // Design the filter
    vector<double> filter = fir1(N, Wn);
    
    // Output the filter coefficients
    cout << "Filter Coefficients: " << endl;
    for (double coeff : filter) {
        cout << coeff << " ";
    }
    cout << endl;

    return 0;
}


