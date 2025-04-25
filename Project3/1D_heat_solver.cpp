#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

// Struct to store material properties
struct MaterialProperties {
    double thermalConductivity;  // W/m*K
    double specificHeatCapacity; // J/Kg*K
    double glassTransitionTemp;  // K
    double density;              // kg/m^3
};

// Define material properties as constants for different materials
const MaterialProperties THERMAL_PROTECTION {
    0.2,   // Thermal conductivity: 0.2 W/m*K
    1200,  // Specific heat capacity: 1200 J/Kg*K
    1200,  // Glass transition temperature: 1200K
    160    // Density: 160 kg/m^3
};

const MaterialProperties CARBON_FIBER {
    500,   // Thermal conductivity: 500 W/m*K
    700,   // Specific heat capacity: 700 J/Kg*K
    350,   // Glass transition temperature: 350K
    0      // Density not specified
};

const MaterialProperties GLUE {
    200,   // Thermal conductivity: 200 W/m*K
    900,   // Specific heat capacity: 900 J/Kg*K
    400,   // Glass transition temperature: 400K
    0      // Density not specified
};

const MaterialProperties STEEL {
    100,   // Thermal conductivity: 100 W/m*K
    500,   // Specific heat capacity: 500 J/Kg*K
    800,   // Maximum allowable temperature: 800K
    0      // Density not specified
};

// Normalize position from [0,2.5] meters to [0,1]
double normalizePosition(double z) {
    return z / 2.5;
}

// Carbon fiber thickness variation along z using sinusoidal profile
double carbonThickness(double z) {
    double normalizedZ = 1.0 - normalizePosition(z);
    const double f = 1.0;        // Frequency
    const double A = 1.5;        // Amplitude in cm
    return abs(A * sin(2.0 * M_PI * f * normalizedZ)) + 0.1;
}

// Glue thickness variation along z using logarithmic profile
double glueThickness(double z) {
    double normalizedZ = 1.0 - normalizePosition(z);
    const double A = 0.1;  // Curve steepness
    const double B = 20.0; // Slope modifier
    const double C = 0.01; // Offset
    return A * log(B * normalizedZ + 1.0) + C;
}

// Steel thickness variation along z using a sawtooth waveform
double steelThickness(double z) {
    double normalizedZ = 1.0 - normalizePosition(z);
    const double f = 5.0;        // Frequency
    const double A = 5.0;        // Amplitude in cm
    double phase = 2.0 * M_PI * f * normalizedZ;
    double sawtoothValue = 2.0 * (phase / (2.0 * M_PI) - floor(0.5 + phase / (2.0 * M_PI)));
    return (A / 2.0) * (sawtoothValue + 1.0) + 0.1;
}

// Initial external temperature profile as a function of z
double initialTemp(double z) {
    double normalizedZ = normalizePosition(z);
    const double A = -100.0;
    const double B = 8.0;
    const double C = 900.0;
    return A * log(B * normalizedZ + 1.0) + C;
}

// Thomas algorithm to efficiently solve tridiagonal systems Ax = d
std::vector<double> thomas_algorithm(const std::vector<double>& a,
                                     const std::vector<double>& b,
                                     const std::vector<double>& c,
                                     const std::vector<double>& d) {
    int n = b.size();
    std::vector<double> c_prime(n - 1);
    std::vector<double> d_prime(n);
    std::vector<double> x(n);

    // Forward elimination
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

// Solve 1D heat conduction through the thermal protection material
vector<double> calculateTempDistribution(double protectionThickness, double missionDuration, double dt, double dx, double initialExternalTemp) {
    double k = THERMAL_PROTECTION.thermalConductivity;
    double rho = THERMAL_PROTECTION.density;
    double cp = THERMAL_PROTECTION.specificHeatCapacity;
    double alpha = k / (rho * cp);  // Thermal diffusivity

    double thickness = protectionThickness/100; // Convert cm to meters
    int Nx = (int)(thickness / dx) + 1;          // Number of spatial nodes

    vector<double> T(Nx, 300.0); // Initialize temperatures to 300K (room temperature)

    // Apply boundary condition at external surface
    T[0] = initialExternalTemp;

    // BTCS method coefficients
    double lambda = alpha * dt / (dx * dx);
    vector<double> a(Nx - 1, -lambda);
    vector<double> b(Nx, 1 + (2 * lambda));
    vector<double> c(Nx - 1, -lambda);

    b[0] = 1.0;      // Fixed boundary at surface
    c[0] = 0.0;      // No forward connection at surface
    b.back() = 1+lambda;  // Slight modification at insulated end

    // Time stepping loop
    for (double t = 0; t < missionDuration; t += dt) {
        vector<double> d = T;
        T = thomas_algorithm(a, b, c, d);
    }

    return T;
}

// Calculate minimum required thickness of thermal protection at each z along the robot
vector<double> calculateRequiredProtectionThickness(double missionDuration, double dz) {
    bool DEBUG = false; 
    double totalSize = 2.50; // Robot height in meters
    vector<double> insulatorThickness((int)(totalSize/dz)+1, 0.0);

    int N = (int)(totalSize / dz + 1e-6);
    for (int i = 0; i <= N; i++) {
        double z = i*dz;
        double externalTemp = initialTemp(z);
        double maxAllowableTemp = CARBON_FIBER.glassTransitionTemp;
        double dx = 0.001;  // Spatial step size
        double dt = 1;      // Time step size
        double minThickness = 0.1;  // Minimum test thickness (cm)
        double maxThickness = 50.0; // Maximum allowable thickness (cm)
        double tolerance = 0.1;     // Search tolerance (cm)
        bool safeTemperature = false; 
        
        // Initialize maxThickness based on previous point for smoother variation
        if(i > 0) {
            maxThickness = insulatorThickness[i-1] + 2*tolerance;
        }
        
        double currentThickness = (minThickness + maxThickness) / 2.0;
        
        // Binary search to find minimum thickness satisfying temperature constraint
        while (maxThickness - minThickness > tolerance) {
            currentThickness = (minThickness + maxThickness) / 2.0;
            vector<double> T = calculateTempDistribution(currentThickness, missionDuration, dt, dx, externalTemp);
            double tempAtCarbonInterface = T.back();

            if (tempAtCarbonInterface <= CARBON_FIBER.glassTransitionTemp) {
                maxThickness = currentThickness;  // Try thinner protection
            } else {
                minThickness = currentThickness;  // Need thicker protection
            }
        }
        
        insulatorThickness[i] = currentThickness;
        cout << "Required thermal protection thickness at z = "<< z <<" m for " << (int)(missionDuration/3600) << " hour mission: " << insulatorThickness[i] << " cm" << endl;
    }

    return insulatorThickness;
}

// Main function
int main() {
    double dz = 0.1; // Spatial step size of 10 cm

    // Calculate required protection thickness for 1hr, 3hr, and 7hr missions
    vector<double> thickness1hr = calculateRequiredProtectionThickness(3600, dz);  
    vector<double> thickness3hr = calculateRequiredProtectionThickness(10800, dz); 
    vector<double> thickness7hr = calculateRequiredProtectionThickness(25200, dz); 
    
    cin.get();
    return 0;
}
