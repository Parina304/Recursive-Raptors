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
    1600      // Density not specified
};

const MaterialProperties GLUE {
    200,   // Thermal conductivity: 200 W/m*K
    900,   // Specific heat capacity: 900 J/Kg*K
    400,   // Glass transition temperature: 400K
    1300      // Density not specified
};

const MaterialProperties STEEL {
    100,   // Thermal conductivity: 100 W/m*K
    500,   // Specific heat capacity: 500 J/Kg*K
    800,   // Maximum allowable temperature: 800K
    7850      // Density not specified
};

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

vector<double> solveMultiLayer(
    const vector<double>& layerThickness,
    const vector<MaterialProperties>& materials,
    double missionDuration,
    double dt,
    double dx,
    double T_ext,
    bool dummy
) 
{
    double totalThickness = 0;
    for (double t_cm : layerThickness) totalThickness += t_cm/100.0;
    int N = static_cast<int>(totalThickness / dx) + 1;

    vector<double> T(N, 300.0); // Initialize temperatures
    vector<double> alpha(N);    // thermal diffusivity at each node

    int idx = 0;
    for (size_t L = 0; L < materials.size(); ++L) {
        double t_m = layerThickness[L]/100.0;
        int nL = static_cast<int>(t_m / dx);
        for (int i = 0; i < nL && idx < N; ++i, ++idx) {
            alpha[idx] = materials[L].thermalConductivity / (materials[L].density * materials[L].specificHeatCapacity);
        }
    }
    while (idx < N) {
        alpha[idx++] = materials.back().thermalConductivity / (materials.back().density * materials.back().specificHeatCapacity);
    }

    vector<double> a(N-1, 0.0), b(N, 0.0), c(N-1, 0.0);

    for (int i = 0; i < N; ++i) {
        if (i == 0) {
            b[i] = 1.0;
            c[i] = 0.0;
        }
        else if (i == N-1) {
            double alphaL = alpha[i-1];
            double alphaR = alpha[i];
            double alpha_eff = 2.0 * alphaL * alphaR / (alphaL + alphaR);

            double lambda = alpha_eff * dt / (dx * dx);
            a[i-1] = -lambda;
            b[i] = 1.0 + lambda;
        }
        else {
            double alphaL = alpha[i-1];
            double alphaC = alpha[i];
            double alphaR = alpha[i+1];

            double alpha_eff_L = 2.0 * alphaL * alphaC / (alphaL + alphaC);
            double alpha_eff_R = 2.0 * alphaC * alphaR / (alphaC + alphaR);

            double lambda_L = alpha_eff_L * dt / (dx * dx);
            double lambda_R = alpha_eff_R * dt / (dx * dx);

            a[i-1] = -lambda_L;
            b[i]   = 1.0 + lambda_L + lambda_R;
            c[i]   = -lambda_R;
        }
    }

    T[0] = T_ext;
    for (double t = 0; t < missionDuration; t += dt) {
        vector<double> d = T;
        d[0] = T_ext;
        T = thomas_algorithm(a, b, c, d);
    }

    return T;
}


// New function: Calculate required insulation thickness by solving full multilayer ODE stack
vector<double> calculateRequiredInsulationThicknessMultiLayer(double missionDuration, double dz, bool normalized) 
{
    double totalHeight = 2.50;               // robot height in meters
    int Nz = static_cast<int>(totalHeight / dz) + 1;
    vector<double> insThick(Nz), zVals(Nz);
    vector<MaterialProperties> mats = { THERMAL_PROTECTION, CARBON_FIBER, GLUE, STEEL };
    vector<double> zValues((int)(totalHeight/dz)+1, 0.0);

    for (int i = 0; i < Nz; ++i) {
        double z = i * dz;
        zVals[i] = z;
        double dx = 0.0001;  // Spatial step size
        double dt = 1;      // Time step size

        // fixed layer thicknesses (cm)
        double cf = carbonThickness(z);
        double gl = glueThickness(z);
        double st = steelThickness(z);

        // binary search bounds (cm)
        double lo = 0.1, hi = 50.0;
        double tol = 0.01;
        double mid = 0;

        while (hi - lo > tol) {
            mid = 0.5 * (lo + hi);
            vector<double> layers = { mid, cf, gl, st };
            double T_ext = initialTemp(z);

            // solve full stack
            auto Tdist = solveMultiLayer(
                layers, mats,
                missionDuration,
                dt, dx,
                T_ext,
                true
            );

            // index at insulation/carbon interface
            int idx_if = static_cast<int>((mid/100.0) / dx + 0.5);
            double T_if = Tdist[idx_if];

            // check against carbon fiber limit
            if (T_if <= CARBON_FIBER.glassTransitionTemp-3) {
                hi = mid;
            } else {
                lo = mid;
            }
        }

        insThick[i] = 0.5 * (lo + hi);
        if(normalized)
        {
            cout << "z=" << normalizePosition(z) << " units, ins_thick=" << insThick[i]
                 << " cm -> T_if=";
        }
        else
        {
            cout << "z=" << z << " m, ins_thick=" << insThick[i]
                 << " cm -> T_if=";
        }
        // final solve for reporting
        auto Tfin = solveMultiLayer({insThick[i], cf, gl, st}, mats,
                                     missionDuration, dt, dx, initialTemp(z), normalized);
        int idx = static_cast<int>((insThick[i]/100.0) / dx + 0.5);
        cout << Tfin[idx] << " K\n";
        zValues[i] = z;
    }

    // outputting thickness values in a csv
    string filename = "Thickness_" + to_string((int)(missionDuration/3600))+"_hr.csv";
    writeVectorsToCSV(filename, zValues, insThick);

    std::cout << "Output written in a csv!" << std::endl;
    return insThick;
}


int main() {
    double dz = 0.1;       // z-step in meters
    double dt = 1.0;       // time step in seconds
    double dx = 0.0001;    // x-step in meters
    double totalTime = 3600.0;
    bool normalized = true;

    // 1) Compute insulation thickness using your existing method
    auto thick1 = calculateRequiredProtectionThickness(totalTime, dz);
    int expectedLines = thick1.size();

    // 2) Open CSV with precomputed thickness profile
    ifstream file("thickness_profile.csv");
    if (!file.is_open()) {
        cerr << "Failed to open thickness_profile.csv\n";
        return 1;
    }

    string line;
    getline(file, line); // skip header

    int i = 0;
    while (getline(file, line)) {
        if (i >= expectedLines) {
            cerr << "CSV has more rows than expected (" << expectedLines << "). Aborting.\n";
            break;
        }

        stringstream ss(line);
        string token;

        double z_val, T_ext, cf, gl, st;
        getline(ss, token, ','); z_val = stod(token);
        getline(ss, token, ','); T_ext = stod(token);
        getline(ss, token, ','); cf = stod(token);
        getline(ss, token, ','); gl = stod(token);
        getline(ss, token, ','); st = stod(token);

        double expected_z = i * dz;
        if (abs(z_val - expected_z) > 1e-4) {
            cerr << "Mismatch: expected z=" << expected_z << ", got z=" << z_val << "\n";
            break;
        }

        double ins = thick1[i]; // cm

        vector<double> layers = { ins, cf, gl, st };
        vector<MaterialProperties> mats = {
            THERMAL_PROTECTION,
            CARBON_FIBER,
            GLUE,
            STEEL
        };

        // Solve through layers
        auto Tdist = solveMultiLayer(
            layers, mats,
            totalTime, dt, dx,
            T_ext,
            true // normalized
        );

        // Compute interface indices
        int i_ins_end  = int(((ins + cf)       / 100.0) / dx);
        int i_cf_end   = int(((ins + cf + gl)  / 100.0) / dx);
        int i_glue_end = int(((ins + cf + gl + st) / 100.0) / dx);

        double T_cf_if   = Tdist[i_ins_end];
        double T_glue_if = Tdist[i_cf_end];
        double T_stl_if  = Tdist[i_glue_end];

        double z_norm = normalizePosition(z_val);

        if (normalized) {
            cout << fixed << setprecision(2)
                 << "z=" << z_norm << " units -> "
                 << "T_CF="   << T_cf_if   << "K, "
                 << "T_GLUE=" << T_glue_if << "K, "
                 << "T_STL="  << T_stl_if  << "K\n";
        } else {
            cout << fixed << setprecision(2)
                 << "z=" << z_val << " m -> "
                 << "T_CF="   << T_cf_if   << "K, "
                 << "T_GLUE=" << T_glue_if << "K, "
                 << "T_STL="  << T_stl_if  << "K\n";
        }

        i++;
    }

    if (i != expectedLines) {
        cerr << "Warning: Processed " << i << " lines, but expected " << expectedLines << "\n";
    }

    file.close();
    cin.get();
    return 0;
}
