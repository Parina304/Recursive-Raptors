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

struct MaterialProperties {
    std::string name;            // Material name
    double thermalConductivity;  // W/m*K
    double specificHeatCapacity; // J/kg*K
    double glassTransitionTemp;  // K
    double density;              // kg/m^3
    double alpha;                // thermal diffusivity (m^2/s)

    // Constructor
    MaterialProperties(std::string n, double k, double cp, double Tg, double rho)
        : name(n),
          thermalConductivity(k),
          specificHeatCapacity(cp),
          glassTransitionTemp(Tg),
          density(rho),
          alpha(k / (rho * cp)) {} // precompute thermal diffusivity
};


// Define material properties as constants for different materials
const MaterialProperties THERMAL_PROTECTION {
    "THERMAL",  // name
    0.2,   // Thermal conductivity: 0.2 W/m*K
    1200,  // Specific heat capacity: 1200 J/Kg*K
    1200,  // Glass transition temperature: 1200K
    160    // Density: 160 kg/m^3
};

const MaterialProperties CARBON_FIBER {
    "CARBON",  // name
    500,   // Thermal conductivity: 500 W/m*K
    700,   // Specific heat capacity: 700 J/Kg*K
    350,   // Glass transition temperature: 350K
    1600      // Density not specified
};

const MaterialProperties GLUE {
    "GLUE",  // name
    200,   // Thermal conductivity: 200 W/m*K
    900,   // Specific heat capacity: 900 J/Kg*K
    400,   // Glass transition temperature: 400K
    1300      // Density not specified
};

const MaterialProperties STEEL {
    "STEEL",  // name
    100,   // Thermal conductivity: 100 W/m*K
    500,   // Specific heat capacity: 500 J/Kg*K
    800,   // Maximum allowable temperature: 800K
    7850      // Density not specified
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
    vector<double> zValues((int)(totalSize/dz)+1, 0.0);

    int N = (int)(totalSize / dz + 1e-6);
    for (int i = 0; i <= N; i++) {
        double z = i*dz;
        double externalTemp = initialTemp(z);
        double maxAllowableTemp = CARBON_FIBER.glassTransitionTemp;
        double dx = 0.001;  // Spatial step size
        double dt = 1;      // Time step size
        double minThickness = 0.1;  // Minimum test thickness (cm)
        double maxThickness = 50.0; // Maximum allowable thickness (cm)
        double tolerance = 0.001;     // Search tolerance (cm)
        bool safeTemperature = false; 
        
        // Initialize maxThickness based on previous point for smoother variation
        if(i > 0) {
            maxThickness = insulatorThickness[i-1] + 2;
            // maxThickness = 50;
        }
        
        double currentThickness = (minThickness + maxThickness) / 2.0;
        double tempAtCarbonInterface;
        // Binary search to find minimum thickness satisfying temperature constraint
        while (maxThickness - minThickness > tolerance) {
            currentThickness = (minThickness + maxThickness) / 2.0;
            vector<double> T = calculateTempDistribution(currentThickness, missionDuration, dt, dx, externalTemp);
            double tempAtCarbonInterface = T.back();

            if (tempAtCarbonInterface < CARBON_FIBER.glassTransitionTemp-3) 
            {
                maxThickness = currentThickness;  // Try thinner protection
            } 
            else if((tempAtCarbonInterface > CARBON_FIBER.glassTransitionTemp-3))
            {
                minThickness = currentThickness;  // Need thicker protection
            }
            else
            {
                maxThickness = currentThickness;
                minThickness = currentThickness;
            }
        }
        
        vector<double> T = calculateTempDistribution(currentThickness, missionDuration, dt, dx, externalTemp);
        tempAtCarbonInterface = T.back();

        insulatorThickness[i] = currentThickness;
        zValues[i] = z;

        cout << "Required thermal protection thickness at z = "<< z <<" m for " << (int)(missionDuration/3600) << " hour mission: " << insulatorThickness[i] << " cm -> "<< tempAtCarbonInterface << "K output" << endl;
    }
    
    // outputting thickness values in a csv
    string filename = "Thickness_" + to_string((int)(missionDuration/3600))+"_hr.csv";
    writeVectorsToCSV(filename, zValues, insulatorThickness);

    std::cout << "Output written in a csv!" << std::endl;
    return insulatorThickness;
}

vector<vector<double>> solveMultiLayer(
    const vector<double>& layerThicknesses_cm,
    const vector<int>& pointsPerLayer,
    const vector<MaterialProperties>& materials,
    double missionDuration,
    double dt,
    double T_ext,
    bool normalized
) {
    // Step 1: build global grid and property vectors
    vector<double> T, alpha, dx;
    vector<int> layerStartIdx;

    for (size_t i = 0; i < layerThicknesses_cm.size(); ++i) {
        int n = pointsPerLayer[i];
        double t_m = layerThicknesses_cm[i] / 100.0;
        double dx_i = t_m / (n - 1);

        if (i == 0) layerStartIdx.push_back(0);
        else layerStartIdx.push_back(T.size());

        for (int j = 0; j < n; ++j) {
            T.push_back(300.0);
            alpha.push_back(materials[i].alpha);
            dx.push_back(dx_i);
        }
        dx.pop_back(); // avoid double count between layers
        T.pop_back();
        alpha.pop_back();
    }
    dx.push_back(dx.back()); // restore final node
    T.push_back(300.0);
    alpha.push_back(alpha.back());

    int N = T.size();
    int steps = static_cast<int>(missionDuration / dt);
    vector<vector<double>> T_out(steps + 1, T);

    // Step 2: time marching
    for (int s = 1; s <= steps; ++s) {
        vector<double> a(N - 1), b(N), c(N - 1), d = T;

        for (int i = 0; i < N; ++i) {
            if (i == 0) {
                b[i] = 1.0; c[i] = 0.0; d[i] = T_ext;
            } else if (i == N - 1) {
                double dxL = dx[i - 1];
                double alpha_eff = alpha[i - 1];
                double lambda = alpha_eff * dt / (dxL * dxL);
                a[i - 1] = -lambda;
                b[i] = 1.0 + lambda;
            } else {
                double dxL = dx[i - 1];
                double dxR = dx[i];
                double alphaL = alpha[i - 1];
                double alphaR = alpha[i + 1];
                double alphaC = alpha[i];

                double a_eff = 2 * alphaL * alphaC / (alphaL + alphaC);
                double c_eff = 2 * alphaC * alphaR / (alphaC + alphaR);

                double lambdaL = a_eff * dt / (dxL * (dxL + dxR) / 2);
                double lambdaR = c_eff * dt / (dxR * (dxL + dxR) / 2);

                a[i - 1] = -lambdaL;
                b[i] = 1.0 + lambdaL + lambdaR;
                c[i] = -lambdaR;
            }
        }

        T = thomas_algorithm(a, b, c, d);
        T_out[s] = T;
    }

    return T_out;
}



// New function: Calculate required insulation thickness by solving full multilayer stack
vector<double> calculateRequiredInsulationThicknessMultiLayer(
    double missionDuration,
    double dz,
    bool normalized,
    const vector<int>& pointsPerLayer  // e.g. {15, 5, 5} for CF, glue, steel
) {
    double totalHeight = 2.50;
    int Nz = static_cast<int>(totalHeight / dz) + 1;
    vector<double> insThick(Nz), zVals(Nz);
    vector<MaterialProperties> mats = { THERMAL_PROTECTION, CARBON_FIBER, GLUE, STEEL };

    for (int i = 0; i < Nz; ++i) {
        double z = i * dz;
        zVals[i] = z;
        double dt = 1.0;

        // Thicknesses in cm
        double cf = carbonThickness(z);
        double gl = glueThickness(z);
        double st = steelThickness(z);

        // Binary search for required insulation
        double lo = 0.1, hi = 50.0, tol = 0.01;
        double mid = 0.0;

        while (hi - lo > tol) {
            mid = 0.5 * (lo + hi);
            vector<double> layers = { mid, cf, gl, st };
            vector<int> points = pointsPerLayer;
            for (int n : pointsPerLayer) points.push_back(n);

            auto Tdist = solveMultiLayer(layers, points, mats, missionDuration, dt, initialTemp(z), normalized);
            double T_if = Tdist.back()[9]; // interface with carbon (index 9)

            if (T_if <= CARBON_FIBER.glassTransitionTemp - 3.0) // why 3 K threshold? ... mention in report about safety 
                hi = mid;
            else
                lo = mid;
        }

        insThick[i] = 0.5 * (lo + hi);

        // Final printout
        if (normalized)
            cout << "z=" << normalizePosition(z) << " units, ins_thick=" << insThick[i] << " cm\n";
        else
            cout << "z=" << z << " m, ins_thick=" << insThick[i] << " cm\n";
    }

    writeVectorsToCSV("Thickness_" + to_string(int(missionDuration/3600)) + "_hr.csv", zVals, insThick);
    return insThick;
}


void MakeTimeSolution(double missionDuration,
                      const vector<double>& insulationThickness,
                      double dt, double dz, bool normalized,
                      const vector<int>& pointsPerLayer)
{
    int N = insulationThickness.size();
    vector<MaterialProperties> mats = {
        {"THERMAL_PROTECTION", 0.2, 1200, 1200, 160},
        {"CARBON_FIBER", 500, 700, 350, 1600},
        {"GLUE", 200, 900, 400, 1300},
        {"STEEL", 100, 500, 800, 7850}
    };

    ofstream fout("full_temp_profile.csv");

    int totalSteps = static_cast<int>(missionDuration / dt);

    for (int step = 0; step <= totalSteps; ++step) {
        double time = step * dt;

        for (int i = 0; i < N; ++i) {
            double z = i * dz;
            double ins = insulationThickness[i];
            double cf = carbonThickness(z);
            double gl = glueThickness(z);
            double st = steelThickness(z);

            vector<double> layers = {ins, cf, gl, st};
            vector<int> points = pointsPerLayer;
            double T_ext = initialTemp(z);

            vector<vector<double>> Tdist = solveMultiLayer(
                layers, points, mats,
                missionDuration, dt,
                T_ext, normalized
            );

            fout << "t=" << time << ",z=" << z << ",";

            int start = 0;
            for (size_t l = 0; l < mats.size(); ++l) {
                int n = points[l];
                double dx = layers[l] / (n - 1);
                fout << mats[l].name << "," << layers[l] << "," << dx;
                for (int j = 0; j < n; ++j) {
                    fout << "," << Tdist[step][start + j];
                }
                start += (n - 1); // avoid double-count at layer boundaries
            }

            fout << "\n";
        }
    }

    fout.close();
    cout << "full_temp_profile.csv generated.\n";
}


// Main function
int main() {
    double dz = 0.1;       // z‐step in meters
    double dt = 1.0;       // time step in seconds
    bool normalized = false;
    double missionDuration = 3600;  // 1 hour
    int Nz = static_cast<int>(2.5 / dz) + 1;

    // Resolution per layer: [insulator, carbon fiber, glue, steel]
    vector<int> pointsPerLayer = {10, 10, 10, 10};

    // Step 1: Calculate required insulation thickness at each z
    auto insulationThickness = calculateRequiredInsulationThicknessMultiLayer(missionDuration, dz, normalized, pointsPerLayer);

    // Step 2: Run full-time simulation using calculated insulation
    MakeTimeSolution(missionDuration, insulationThickness, dt, dz, normalized, pointsPerLayer);

    std::cin.get();
    return 0;
}
