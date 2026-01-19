#include <iostream>
#include <cmath>
#include <vector>
#include <numbers>
#include <iomanip>

double f(double x) { return std::sin(x);}

double romberg(double a, double b, double tol = 1e-8, int maxIter = 10) {
    std::vector<std::vector<double>> R(maxIter +1, std::vector<double>(maxIter + 1, 0.0));

    double h = b - a;
    R[0][0] = 0.5 * h * (f(a) + f(b));

    for (int k = 1; k <= maxIter; k++) {
        double h_k = (b - a) / (1 << k);

        double sum = 0.0;
        int newPoints = 1 << (k-1);
        for (int i = 1; i <= newPoints; i++) {
            double x = a + (2*i-1)*h_k;
            sum += f(x);
        }

        R[k][0] = 0.5 * R[k-1][0]+h_k*sum;

        for (int j = 1; j <= k; j++) {
            double pow4j = std::pow(4,j);
            R[k][j] = R[k][j-1]+(R[k][j-1]-R[k-1][j-1])/(pow4j -1);
        }

        if (std::fabs(R[k][k]-R[k-1][k-1]) < tol) {
            return R[k][k];
        }
    }

    return R[maxIter][maxIter];
}


int main() {
    double a = 0.0;
    double b = std::numbers::pi;
    double result = romberg(a, b);

    std::cout << std::fixed << std::setprecision(15);
    std::cout << "Rezultat Romberg-ove integracije je: " << result << std::endl;

    return 0;
}