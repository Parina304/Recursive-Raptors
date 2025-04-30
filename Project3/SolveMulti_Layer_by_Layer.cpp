#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

// Function to write two vectors into a CSV file with two columns
void writeVectorsToCSV(const string& filename, const vector<double>& vec1, const vector<double>& vec2) {
    // Check if vectors have the same size
    if (vec1.size() != vec2.size()) {
        cerr << "Error: Vectors must have the same size." << endl;
        return;
    }

    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    // Optional: Write a header
    file << "z value,thickness (cm)\n";

    // Write the data
    for (size_t i = 0; i < vec1.size(); ++i) {
        file << vec1[i] << "," << vec2[i] << "\n";
    }

    file.close();
}

// Struct to store material properties
struct MaterialProperties {
    double thermalConductivity;  // W/m*K
    double specificHeatCapacity; // J/kg*K
    double glassTransitionTemp;  // K
    double density;              // kg/m^3
    double alpha;                // thermal diffusivity (m^2/s)

    // Constructor
    MaterialProperties(double k, double cp, double Tg, double rho)
        : thermalConductivity(k),
          specificHeatCapacity(cp),
          glassTransitionTemp(Tg),
          density(rho),
          alpha(k / (rho * cp)) {} // calculate once at construction
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
    1600   // Density: 1600 kg/m^3
};

const MaterialProperties GLUE {
    200,   // Thermal conductivity: 200 W/m*K
    900,   // Specific heat capacity: 900 J/Kg*K
    400,   // Glass transition temperature: 400K
    1300   // Density: 1300 kg/m^3
};

const MaterialProperties STEEL {
    100,   // Thermal conductivity: 100 W/m*K
    500,   // Specific heat capacity: 500 J/Kg*K
    800,   // Maximum allowable temperature: 800K
    7850   // Density: 7850 kg/m^3
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
    double frequency = 5.0;
    double phase = std::fmod(frequency * normalizedZ, 1.0);
    if (phase < 0.0)
        phase += 1.0; // wrap negative times into [0,1)
    return 5.0 * phase + 0.1;
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

vector<double> calculateTempDistribution_layered(
    const vector<double>& T_prev,  // previous temp field
    double dt, double dx,
    const MaterialProperties& mat,
    double T_left_boundary)
{
    int Nx = T_prev.size();
    double lambda = mat.alpha * dt / (dx * dx);

    // Matrix
    vector<double> a(Nx - 1, -lambda);
    vector<double> b(Nx, 1 + 2 * lambda);
    vector<double> c(Nx - 1, -lambda);

    b[0] = 1.0;         // Dirichlet on outer side -- const temp boundary 
    c[0] = 0.0;
    b.back() = 1 + lambda;  // insulated inner side
    a[Nx - 2] = -lambda;

    // RHS
    vector<double> d = T_prev;
    d[0] = T_left_boundary;

    return thomas_algorithm(a, b, c, d);
}


vector<vector<double>> solveMultiLayerSequentialEvolving(
    const vector<double>& layerThicknesses_cm,
    const vector<MaterialProperties>& materials,
    double totalDuration, double dt,
    double dx, double T_ext)
{
    size_t numLayers = layerThicknesses_cm.size();
    vector<vector<double>> layerTemps(numLayers);

    // 🔹 Initialize all layers at t = 0
    for (size_t i = 0; i < numLayers; ++i) {
        int Nx = static_cast<int>(layerThicknesses_cm[i] / 100.0 / dx) + 1;
        layerTemps[i].assign(Nx, 300.0);  // Initial guess = 300K everywhere
    }

    int steps = static_cast<int>(totalDuration / dt);

    for (int step = 0; step < steps; ++step) {
        double T_left = T_ext;  // reapply external temp every step

        for (size_t i = 0; i < numLayers; ++i) {
            const auto& mat = materials[i];
            auto& T_current = layerTemps[i];

            // Solve 1 timestep starting from current temp field
            vector<double> T_new = calculateTempDistribution_layered(
                T_current, dt, dx, mat, T_left);

            T_current = T_new;
            T_left = T_new.back();  // pass interface temp forward
        }
    }

    return layerTemps;
}




// Main function
int main() {

    double dz = 0.1;       // Step in z-direction (meters)
    double dt = 1.0;       // Time step in seconds
    double dx = 0.0001;    // Space step in meters
    double totalDuration = 3600.0; // 1 hour

    vector<MaterialProperties> materials = {
        THERMAL_PROTECTION, GLUE, CARBON_FIBER, STEEL
    };

    cout << fixed << setprecision(2);
    cout << "z (m), T_CF (K), T_GLUE (K), T_STEEL (K)\n";

    for (double z = 0.0; z <= 2.5; z += dz) {
        // Dynamic layer thicknesses at current z (in cm)
        double ins = 18;  // cm — fixed for now
        double cf  = carbonThickness(z);  // cm
        double gl  = glueThickness(z);    // cm
        double st  = steelThickness(z);   // cm

        vector<double> layerThicknesses_cm = { ins, gl, cf, st };

        // External boundary condition
        double T_ext = initialTemp(z);  // K

        // Run full solver
        auto finalTemps = solveMultiLayerSequentialEvolving(
            layerThicknesses_cm, materials, totalDuration, dt, dx, T_ext);

        // Extract interface temps
        double T_cf   = finalTemps[1].back(); // after glue
        double T_glue = finalTemps[2].back(); // after carbon
        double T_steel = finalTemps[3].back(); // after steel

        cout << z << ", " << T_cf << ", " << T_glue << ", " << T_steel << "\n";
    }

    return 0;
}
    
