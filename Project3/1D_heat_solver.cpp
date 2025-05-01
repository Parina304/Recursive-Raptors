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
    double thermalConductivity;
    double specificHeatCapacity;
    double glassTransitionTemp;
    double density;
    double alpha;

    // Now only four parameters
    MaterialProperties(double k, double cp, double Tg, double rho)
      : thermalConductivity(k),
        specificHeatCapacity(cp),
        glassTransitionTemp(Tg),
        density(rho),
        alpha(k/(rho*cp)) {}
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
    623.15,   // Glass transition temperature: 350K --> 623.15K
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
    double alpha = THERMAL_PROTECTION.alpha;  // Thermal diffusivity

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
vector<double> calculateRequiredProtectionThickness(double missionDuration, double dz, double dt) {
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
    int timePoints = static_cast<int>(missionDuration / dt) + 1;

    vector<double> T(N, 300.0); // Initialize temperatures
    vector<vector<double>> T_out(timePoints, vector<double>(N, 300.0));
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
    for (int i = 0; i < timePoints; i++) {
        vector<double> d = T;
        d[0] = T_ext;
        T = thomas_algorithm(a, b, c, d);
        T_out[i] = T;  // Safe: using integer loop index
    }

    return T_out;
}


// New function: Calculate required insulation thickness by solving full multilayer ODE stack
vector<double> calculateRequiredInsulationThicknessMultiLayer(double missionDuration, double dz, double dt, bool normalized) 
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
            int idx_c = static_cast<int>(((mid+cf)/100.0) / dx + 0.5);
            int idx_g = static_cast<int>(((mid+cf+gl)/100.0) / dx + 0.5);
            int idx_s = static_cast<int>(((mid+cf+gl+st)/100.0) / dx + 0.5);
            double T_if = Tdist.back()[idx_if];
            double T_c = Tdist.back()[idx_c];
            double T_g = Tdist.back()[idx_g];
            double T_s = Tdist.back()[idx_s];

            // check against carbon fiber limit
            if (T_if <= CARBON_FIBER.glassTransitionTemp-3 && T_c <=GLUE.glassTransitionTemp-3 && T_g <=STEEL.glassTransitionTemp-3) {
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
        cout << Tfin.back()[idx] << " K\n";
        zValues[i] = z;
    }

    // outputting thickness values in a csv
    string filename = "Thickness_" + to_string((int)(missionDuration/3600))+"_hr.csv";
    writeVectorsToCSV(filename, zValues, insThick);

    std::cout << "Output written in a csv!" << std::endl;
    return insThick;
}

void MakeTimeSolution(double missionDuration, vector<double> insulationThickness, double dt, double dx, double dz, bool normalized)
{
    ofstream fout("full_temp_profile_1hr_shorter.csv");
    // 2) Now simulate the full 4‐layer stack at each z for, say, the 1 hr mission:
    int N = insulationThickness.size();  // == totalSize/dz + 1
    for (int i = 0; i < N; ++i) {
        double z   = i * dz;
        double ins = insulationThickness[i];        // in cm
        double cf  = carbonThickness(z);            // in cm
        double gl  = glueThickness(z);              // in cm
        double st  = steelThickness(z);             // in cm

        vector<double>   layers = { ins, cf, gl, st };
        vector<MaterialProperties> mats = {
            THERMAL_PROTECTION,
            CARBON_FIBER,
            GLUE,
            STEEL
        };
        double T_ext = initialTemp(z);
        // call the multi‐layer solver:
        auto Tdist = solveMultiLayer(
            layers, mats,
            missionDuration,   // missionDuration (1 hr)
            dt, dx,
            T_ext,
            true
        );

        // cout << Tdist.size() << endl;
        // cout << Tdist[1].size() << endl;


        int totalXPoints = static_cast<int>((ins+cf+gl+st)/(100*dx))+1;
        int totalTPoints = static_cast<int>((missionDuration)/dt)+1;
        // cout << totalXPoints << endl;
        // cout << totalTPoints << endl;
        
        for (int i = 0; i < totalTPoints; ++i) 
        {
            fout << z << ", ";
            fout << i << ", ";
            for (int j = 0; j < totalXPoints; ++j) 
            {
                fout << Tdist[i][j] << ", ";
            }
            fout << "\n";
        }
        fout << "\n";

        // extract interface indices
        int i_ins_end  = int(((ins+cf)/100.0)/dx);
        int i_cf_end   = int(((ins+cf+gl)/100.0)/dx);
        int i_glue_end = int(((ins+cf+gl+st)/100.0)/dx);

        
        double T_cf_if   = Tdist.back()[i_ins_end];
        double T_glue_if = Tdist.back()[i_cf_end];
        double T_stl_if  = Tdist.back()[i_glue_end];
        
        double z_norm = normalizePosition(z);

        if(normalized)
        {
            cout << fixed << setprecision(2)
                 << "z=" << z_norm << " units -> "
                 << "T_CF="   << T_cf_if   << "K, "
                 << "T_GLUE=" << T_glue_if << "K, "
                 << "T_STL="  << T_stl_if  << "K\n";
        }
        else
        {
            cout << fixed << setprecision(2)
                 << "z=" << z << " m -> "
                 << "T_CF="   << T_cf_if   << "K, "
                 << "T_GLUE=" << T_glue_if << "K, "
                 << "T_STL="  << T_stl_if  << "K\n";
        }
    }
    std::cout << "Done" << std::endl; 
}
// Main function
int main() {
    double dz = 0.1;    // z‐step in meters
    double dt = 5.0;    // time step in seconds
    double dx = 0.0001;  // x‐step in meters 
    bool normalized = false;
    vector<double> testThickness(static_cast<int>(2.50/dz)+1,18.0);
    std::cout << "3hr Prof change of glass temp to 350K to 623.15K run" << std::endl; 
    // 1) Compute the minimum insulation thickness at each z for 1 hr, 3 hr, 7 hr
    // auto thick1 = calculateRequiredInsulationThicknessMultiLayer( 3600, dz, dt, normalized);
    auto thick3 = calculateRequiredInsulationThicknessMultiLayer(10800, dz, dt, normalized);
    // auto thick7 = calculateRequiredInsulationThicknessMultiLayer(25200, dz, dt, normalized);
    // auto thick1 = calculateRequiredProtectionThickness(3600, dz);
    // auto thick3 = calculateRequiredProtectionThickness(10800, dz);
    // auto thick7 = calculateRequiredProtectionThickness(25200, dz);

    MakeTimeSolution(10800, thick3, dt, dx, dz, normalized);

    std::cin.get();
    return 0;
}