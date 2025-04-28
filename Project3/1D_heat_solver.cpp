#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
using namespace std;

// z is 0 to 2.5 m from toe to head and thicknesses are defined in cm 

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
    7850   // Density: 7850 kg/m^3
};

const MaterialProperties STEEL {
    100,   // Thermal conductivity: 100 W/m*K
    500,   // Specific heat capacity: 500 J/Kg*K
    800,   // Maximum allowable temperature: 800K
    1300   // Density: 1300 kg/m^3
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
    
    double totalSize = 2.50; // Robot height in meters
    vector<double> insulatorThickness((int)(totalSize/dz)+1, 0.0);
    vector<double> zValues((int)(totalSize/dz)+1, 0.0);

    int N = (int)(totalSize / dz + 1e-6);
    for (int i = 0; i <= N; i++) {
        double z = i*dz;
        double externalTemp = 900; // ---- initialTemp(z) needs to work with the CSV format! --- sorry I changed it such that it is wrong rn 
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
        zValues[i] = z;

        cout << "Required thermal protection thickness at z = "<< z <<" m for " << (int)(missionDuration/3600) << " hour mission: " << insulatorThickness[i] << " cm" << endl;
    }
    
    // outputting thickness values in a csv
    string filename = "Thickness_" + to_string((int)(missionDuration/3600))+"_hr.csv";
    writeVectorsToCSV(filename, zValues, insulatorThickness);

    std::cout << "Output written in a csv!" << std::endl;
    return insulatorThickness;
}

// Solve 1D heat conduction through all layers of material
vector<double> calculateMultiLayerTempDistribution(double missionDuration, double dt,
    const vector<int>& nodesPerLayer,
    const vector<double>& layerDx_m,
    const vector<MaterialProperties>& layerMaterials) {
    
    // Step 1: Build the mesh
    vector<double> x_positions;
    vector<double> alpha;  // thermal diffusivity at each node

    for (size_t i = 0; i < nodesPerLayer.size(); ++i) {
        for (int j = 0; j < nodesPerLayer[i]; ++j) {
            double xpos = (j * layerDx_m[i]);
            if (!x_positions.empty()) xpos += x_positions.back();
            x_positions.push_back(xpos);

            alpha.push_back(layerMaterials[i].alpha);
        }
        if (!x_positions.empty()) x_positions.pop_back(); // remove duplicate interface node
        if (!alpha.empty()) alpha.pop_back();
    }

    int N = x_positions.size();
    vector<double> T(N, 300.0); // Initial temperature
    T[0] = 900.0;               // Boundary condition: fixed 900K at outer surface

    // Step 2: Setup tridiagonal matrix
    vector<double> a(N-1, 0.0); // sub-diagonal
    vector<double> b(N, 0.0);   // main diagonal
    vector<double> c(N-1, 0.0); // super-diagonal

    double dx = (x_positions[1] - x_positions[0]); // Assume uniform dx inside each layer
    double lambda = 0.0; 

    for (int i = 0; i < N; ++i) {
        lambda = alpha[i] * dt / (dx * dx); 

        if (i > 0) a[i-1] = -lambda;        // sub-diagonal
        b[i] = 1.0 + 2.0 * lambda;           // main diagonal
        if (i < N-1) c[i] = -lambda;         // super-diagonal
    }

    // Apply boundary conditions
    b[0] = 1.0;
    c[0] = 0.0;
    a[N-2] = -lambda;
    b[N-1] = 1.0 + lambda;

    // Step 3: Time stepping loop
    for (double time = 0; time < missionDuration; time += dt) {
        vector<double> d = T;
        d[0] = 900.0; 
        T = thomas_algorithm(a, b, c, d);
    }

    return T;
}


//Generates dx values per layer based on number of nodes 
vector<double> generateDxList(const vector<double>& layerThicknesses_cm, const vector<int>& nodes_per_layer) {
    vector<double> dx_per_layer_m;

    for (size_t i = 0; i < layerThicknesses_cm.size(); ++i) {
        double dx = layerThicknesses_cm[i] / 100.0 / (nodes_per_layer[i] - 1); // -1 because we include start and end node & div by 100 for cm to m  
        dx_per_layer_m.push_back(dx);
    }

    return dx_per_layer_m;
}

//Resolves dx values per layer such that the nodes are equally spaced 
pair<vector<double>, vector<int>> resolveDxList(const vector<double>& layerThicknesses_cm, const vector<double>& dx_guess_per_layer_cm) {
    
    vector<double> resolved_dx_per_layer_m;
    vector<int> nodes_per_layer;

    for (size_t i = 0; i < layerThicknesses_cm.size(); ++i) {
        double thickness = layerThicknesses_cm[i] / 100.0; // cm → meters
        double dx_guess = dx_guess_per_layer_cm[i] / 100.0; // cm → meters

        // Calculate approximate number of intervals
        int intervals = (int) round(thickness / dx_guess);
        if (intervals < 1) intervals = 1; // Avoid zero division
            
        double dx_exact = thickness / intervals;

        resolved_dx_per_layer_m.push_back(dx_exact);
        nodes_per_layer.push_back(intervals + 1); // Nodes = intervals + 1
    }

    return {resolved_dx_per_layer_m, nodes_per_layer};
}


void printTemperatureProfile(const vector<double>& temperatures,
    const vector<double>& dx_per_layer_m,
    const vector<MaterialProperties>& materials,
    const vector<int>& nodes_per_layer) {

    cout << "\n=== Material Glass Transition Temperatures ===" << endl;
    for (size_t i = 0; i < materials.size(); ++i) {
        cout << "Material " << i + 1 << ": Glass Transition Temp = "
        << materials[i].glassTransitionTemp << " K" << endl;
    }
    cout << "==============================================" << endl;
    cout << "\nDistance (m)     Temperature (K)    Material" << endl;
    cout << "-----------------------------------------------" << endl;

    double position = 0.0;
    int nodeIndex = 0;

    for (size_t layer = 0; layer < materials.size(); ++layer) {
        double dx = dx_per_layer_m[layer];
        int nodes = nodes_per_layer[layer];

        for (int n = 0; n < nodes; ++n) {
        if (nodeIndex >= temperatures.size()) break;

        cout << fixed << position << "           "
        << temperatures[nodeIndex] << "           "
        << "Layer " << layer + 1 << endl;

        position += dx;
        nodeIndex++;
        }
    cout << "----- End of Material " << layer + 1 << " -----" << endl;
    }
}




// Main function
int main() {

    // cout << "At z=0 [toes]: " << endl; // use csv now 
    // cout << "Initial temp " << initialTemp(0) << endl;
    // cout << "Carbon thickness " << carbonThickness(0) << endl;
    // cout << "Glue thickness " << glueThickness(0) << endl;
    // cout << "Steel thickness " << steelThickness(0) << endl; cout << endl;

    // Define layers (thicknesses in cm)
    vector<double> layerThicknesses_cm = {17, 0.314, 0.1, 0.315, 2.6}; // [thermal protection, glue, carbon fiber, glue, steel] --- usually from the functions, this is just the first layer at toes 
    //props of the first layer

    // Define material list
    vector<MaterialProperties> layerMaterials = {THERMAL_PROTECTION, GLUE, CARBON_FIBER, GLUE, STEEL};

    // Define mission settings
    double missionDuration = 3600.0; // 1 hour
    double dt = 1.0;                 // 1 second

    // ------------ Define dx HEREEEEEEEE via dx guess or num of nodes! ------------ // // LOOK HERE TO CHANGE STUFFFF // // 
    vector<double> guessed_dx_cm = {0.1, 0.005, 0.001, 0.005, 0.01}; // cm
    vector<int> nodesPerLayer = {10, 10, 10, 10, 10}; // Example: 10 nodes per layer
    // Decide: define based on dx guesses or numNodes?
    bool use_dx_guesses = true; // <<< switch this true/false depending on case

    vector<double> layerDx_m; // in meters
    if (use_dx_guesses) {
        // If starting with rough dx guesses, call resolveDxList()
        auto [resolved_dx, resolved_nodes] = resolveDxList(layerThicknesses_cm, guessed_dx_cm);
        layerDx_m = resolved_dx;
        nodesPerLayer = resolved_nodes;
        
    } else {
        // If starting with fixed number of nodes
        layerDx_m = generateDxList(layerThicknesses_cm, nodesPerLayer);
    }

    // Call the solver
    vector<double> finalTemps = calculateMultiLayerTempDistribution(
        missionDuration, dt, nodesPerLayer, layerDx_m, layerMaterials);
    
        // Print results
    printTemperatureProfile(finalTemps, layerDx_m, layerMaterials, nodesPerLayer);


    /*
    // Calculate required protection thickness for 1hr, 3hr, and 7hr missions
    vector<double> thickness1hr = calculateRequiredProtectionThickness(3600, dz);  
    vector<double> thickness3hr = calculateRequiredProtectionThickness(10800, dz); 
    vector<double> thickness7hr = calculateRequiredProtectionThickness(25200, dz); 
    
    cin.get(); */

    return 0;
}
