#include <iostream>
#include <vector>
#include <cstdlib>  // For rand()
#include <ctime>    // For timing

// #include <Eigen/Dense>  // Include Eigen

std::vector<double> thomas_algorithm(const std::vector<double>& a,
                                     const std::vector<double>& b,
                                     const std::vector<double>& c,
                                     const std::vector<double>& d) {
    int n = b.size();
    std::vector<double> c_prime(n - 1);
    std::vector<double> d_prime(n);
    std::vector<double> x(n);

    // Forward sweep
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n - 1; ++i) {
        double denom = b[i] - a[i - 1] * c_prime[i - 1];
        c_prime[i] = c[i] / denom;
        d_prime[i] = (d[i] - a[i - 1] * d_prime[i - 1]) / denom;
    }
    d_prime[n - 1] = (d[n - 1] - a[n - 2] * d_prime[n - 2]) / (b[n - 1] - a[n - 2] * c_prime[n - 2]);

    // Back substitution
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}

int main() {
    const int n = 10;  // Size of system

    // Create tridiagonal system
    std::vector<double> a(n - 1, 4.0);  // sub-diagonal
    std::vector<double> b(n, 22.0);        // main diagonal
    std::vector<double> c(n - 1, 333.0);   // super-diagonal
    std::vector<double> d(n);

    // Random d vector
    for (int i = 0; i < n; ++i) {
        d[i] = i;
    }

    // Thomas Algorithm
    std::clock_t start = std::clock();
    std::vector<double> x = thomas_algorithm(a, b, c, d);
    double elapsed_thomas = double(std::clock() - start) / CLOCKS_PER_SEC;
    std::cout << "Thomas Algorithm finished in " << elapsed_thomas << " seconds.\n";

    /////////////////////////////////////
    //   This block was used to        //
    //   compare Thomas algo           //
    //   implementation with a         //                     
    //   standard result               //               
    /////////////////////////////////////  

    // // Build full Matrix A using Eigen
    // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    // for (int i = 0; i < n; ++i) {
    //     A(i, i) = b[i];  // Main diagonal
    //     if (i > 0) A(i, i - 1) = a[i - 1];  // Lower diagonal
    //     if (i < n - 1) A(i, i + 1) = c[i];  // Upper diagonal
    //     // std::cout << "done " << i << " row" << std::endl;
    // }

    // // Solve using Eigen Inverse
    // Eigen::VectorXd d_eigen = Eigen::VectorXd::Zero(n);
    // for (int i = 0; i < n; ++i) {
    //     d_eigen(i) = d[i];
    // }

    // // std::cout << "done d" << std::endl;

    // start = std::clock();
    // Eigen::VectorXd x2 = A.inverse() * d_eigen;  // Inverse method 
    // double elapsed_eigen = double(std::clock() - start) / CLOCKS_PER_SEC;
    // std::cout << "Eigen's Inverse method finished in " << elapsed_eigen << " seconds.\n";

    // // Print first few differences
    // std::cout << "First 5 differences (Thomas - Eigen): ";
    // for (int i = 0; i < 10; ++i) {
    //     std::cout << (x[i] - x2(i)) << " ";
    // }
    // std::cout << "\n";

    std::cout << "First 5 values (Thomas): ";
    for (int i = 0; i < 10; ++i) {
        std::cout << (x[i]) << " ";
    }
    std::cout << "\n";

    return 0;
}
