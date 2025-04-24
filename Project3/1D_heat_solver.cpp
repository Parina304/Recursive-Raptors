#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

// Material properties struct
struct MaterialProperties {
    double thermalConductivity;  // W/m*K
    double specificHeatCapacity; // J/Kg*K
    double glassTransitionTemp;  // K
    double density;              // kg/m^3
};

// Define material properties as constants
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
    0      // Density not specified in the file
};

const MaterialProperties GLUE {
    200,   // Thermal conductivity: 200 W/m*K
    900,   // Specific heat capacity: 900 J/Kg*K
    400,   // Glass transition temperature: 400K
    0      // Density not specified in the file
};

const MaterialProperties STEEL {
    100,   // Thermal conductivity: 100 W/m*K
    500,   // Specific heat capacity: 500 J/Kg*K
    800,   // Max temperature: 800K
    0      // Density not specified in the file
};

// Convert physical position to normalized position (0-1)
// Input: z in meters (0 = toe, 2.5 = head)
// Output: normalized position (0-1) where 0 is toe, 1 is head
double normalizePosition(double z) {
    return z / 2.5;
}

// Carbon fiber thickness function (sinusoidal profile)
double carbonThickness(double z) {
    // Convert from physical (0=toe, 2.5=head) to normalized and inverted (1=toe, 0=head)
    double normalizedZ = 1.0 - normalizePosition(z);
    
    // Parameters
    const double f = 1.0;        // Frequency in Hz
    const double A = 1.5;        // Amplitude in cm
    
    // Calculate thickness using sinusoidal function
    return abs(A * sin(2.0 * M_PI * f * normalizedZ)) + 0.1;
}

// Glue thickness function (logarithmic profile)
double glueThickness(double z) {
    // Convert from physical (0=toe, 2.5=head) to normalized and inverted (1=toe, 0=head)
    double normalizedZ = 1.0 - normalizePosition(z);
    
    // Parameters
    const double A = 0.1;  // Controls curve steepness
    const double B = 20.0; // Controls how fast it levels off
    const double C = 0.01; // Starting value
    
    // Calculate thickness using logarithmic function
    return A * log(B * normalizedZ + 1.0) + C;
}

// Steel thickness function (sawtooth profile)
double steelThickness(double z) {
    // Convert from physical (0=toe, 2.5=head) to normalized and inverted (1=toe, 0=head)
    double normalizedZ = 1.0 - normalizePosition(z);
    
    // Parameters
    const double f = 5.0;        // Frequency in Hz
    const double A = 5.0;        // Amplitude in cm
    
    // Calculate sawtooth wave and shift to positive range
    double phase = 2.0 * M_PI * f * normalizedZ;
    double sawtoothValue = 2.0 * (phase / (2.0 * M_PI) - floor(0.5 + phase / (2.0 * M_PI)));
    
    // Scale and shift to match the MATLAB implementation
    return (A / 2.0) * (sawtoothValue + 1.0) + 0.1; // Range: 0.1 to 5.1 cm
}

// Initial temperature profile
double initialTemp(double z) {
    // Convert from physical (0=toe, 2.5=head) to normalized (0=toe, 1=head)
    double normalizedZ = normalizePosition(z);
    
    // Parameters
    const double A = -100.0;
    const double B = 8.0;
    const double C = 900.0;
    
    // Using the same logic as in the MATLAB file
    return A * log(B * normalizedZ + 1.0) + C;
}

// Thomas algorithm for solving tridiagonal system (provided by you)
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

// Function to calculate temperature distribution through heat protection layer
vector<double> calculateTempDistribution(double protectionThickness, double missionDuration, double dt, double dx, double initialExternalTemp) {
    // Material properties
    double k = THERMAL_PROTECTION.thermalConductivity;
    double rho = THERMAL_PROTECTION.density;
    double cp = THERMAL_PROTECTION.specificHeatCapacity;
    double alpha = k / (rho * cp);  // thermal diffusivity

    double thickness = protectionThickness/100;
    int Nx = (int)(thickness / dx) + 1;

    // Initialize all temperatures to room temp (or some baseline)
    vector<double> T(Nx, 300.0);

    // Fix only the external boundary (first point)
    T[0] = initialExternalTemp;

    // Coefficients
    double lambda = alpha * dt / (dx * dx);
    cout << "Lambda:: "<< lambda << endl; 
    vector<double> a(Nx - 1, -lambda);
    vector<double> b(Nx, 1 + (2 * lambda));
    vector<double> c(Nx - 1, -lambda);

    // Apply Dirichlet condition at external surface
    b[0] = 1.0;
    c[0] = 0.0;

    b.back() = 1+lambda;

    // Time stepping
    for (int t = 0; t < missionDuration; t += dt) {
        vector<double> d = T;

        T = thomas_algorithm(a, b, c, d);

    }

    return T;
}


// Function to calculate minimum required protection layer thickness
vector<double> calculateRequiredProtectionThickness(double missionDuration, double dz) {

    bool DEBUG = true; 
    double totalSize = 2.50;
    vector<double> thickness(1, (int)(totalSize/dz)+1);
    // Initial conditions - worst case at toe where temp is 900K
    double externalTemp = 900.0;
    double maxAllowableTemp = CARBON_FIBER.glassTransitionTemp;
    double dx = 0.0001;  // 0.1 mm grid spacing
    double dt = 1;     // 1 second time step
    
    // Binary search to find minimum thickness
    double minThickness = 0.01;  // cm
    double maxThickness = 50.0; // cm
    double currentThickness = (minThickness + maxThickness) / 2.0;
    double tolerance = dx;    // Tolerance in cm
    bool safeTemperature = false; 
    
    while (maxThickness - minThickness > tolerance) {
        // Midpoint thickness to test
        currentThickness = (minThickness + maxThickness) / 2.0;
    
        // Simulate the temperature profile
        vector<double> T = calculateTempDistribution(currentThickness, missionDuration, dt, dx, externalTemp);
    
        // Get the temperature at the carbon interface
        double tempAtCarbonInterface = T.back();
    
        // Print debug info
        if (DEBUG) cout << "Testing thickness: " << currentThickness << " cm --> Temp at interface: " 
             << tempAtCarbonInterface << " K --> ";
    
        // Binary search update
        if (tempAtCarbonInterface <= CARBON_FIBER.glassTransitionTemp) {
            if (DEBUG) cout << "SAFE" << endl;
            maxThickness = currentThickness;  // Try thinner
        } else {
            if (DEBUG) cout << "UNSAFE" << endl;
            minThickness = currentThickness;  // Need thicker
        }
    }    
    thickness[0] = currentThickness;
    return thickness;
}

int main() {
    int dz = 0.1; // 1 cm step size
    // Calculate required protection thicknesses for different mission durations
    vector<double> thickness1hr = calculateRequiredProtectionThickness(3600, dz);  // 1 hour in seconds
    vector<double> thickness3hr = calculateRequiredProtectionThickness(10800, dz); // 3 hours in seconds
    vector<double> thickness7hr = calculateRequiredProtectionThickness(25200, dz); // 7 hours in seconds
    
    cout << "Required thermal protection thickness for 1 hour mission: " << thickness1hr[0] << " cm" << endl;
    cout << "Required thermal protection thickness for 3 hour mission: " << thickness3hr[0] << " cm" << endl;
    cout << "Required thermal protection thickness for 7 hour mission: " << thickness7hr[0] << " cm" << endl;
    cin.get();
    return 0;
}