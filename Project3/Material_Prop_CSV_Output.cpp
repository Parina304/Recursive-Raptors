#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

// Constants
const double TOTAL_HEIGHT = 2.5; // meters
const double FIXED_THERMAL_PROTECTION_THICKNESS_CM = 15.0; // uniform 15 cm

// Normalize position from [0,2.5] meters to [0,1]
double normalizePosition(double z) {
    return z / TOTAL_HEIGHT;
}

// Carbon fiber thickness variation along z using sinusoidal profile
double carbonThickness(double z) {
    double normalizedZ = 1.0 - normalizePosition(z);
    const double f = 1.0;        
    const double A = 1.5;        
    return abs(A * sin(2.0 * M_PI * f * normalizedZ)) + 0.1; // in cm
}

// Glue thickness variation along z using logarithmic profile
double glueThickness(double z) {
    double normalizedZ = 1.0 - normalizePosition(z);
    const double A = 0.1;  
    const double B = 20.0; 
    const double C = 0.01;
    return A * log(B * normalizedZ + 1.0) + C; // in cm
}

// Steel thickness variation along z using a sawtooth waveform
double steelThickness(double z) {
    double normalizedZ = 1.0 - normalizePosition(z);
    const double f = 5.0;
    const double A = 5.0;
    double phase = 2.0 * M_PI * f * normalizedZ;
    double sawtoothValue = 2.0 * (phase / (2.0 * M_PI) - floor(0.5 + phase / (2.0 * M_PI)));
    return (A / 2.0) * (sawtoothValue + 1.0) + 0.1; // in cm
}

// Initial external temperature profile as a function of z
double initialTemp(double z) {
    double normalizedZ = normalizePosition(z);
    const double A = -100.0;
    const double B = 8.0;
    const double C = 900.0;
    return A * log(B * normalizedZ + 1.0) + C; // in Kelvin
}

// Write all data to a CSV
void writeThicknessProfileToCSV(const string& filename, double dz) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    file << "z (m),Initial Temp (K),Thermal Protection Thickness (cm),Glue Thickness (cm),Carbon Fiber Thickness (cm),Glue Thickness (cm),Steel Thickness (cm)\n";

    int numSteps = static_cast<int>(TOTAL_HEIGHT / dz) + 1;

    for (int i = 0; i < numSteps; ++i) {
        double z = i * dz;
        if (z > TOTAL_HEIGHT) z = TOTAL_HEIGHT; // ensure it doesn't slightly overshoot

        double temp = initialTemp(z);
        double thermalProtection = FIXED_THERMAL_PROTECTION_THICKNESS_CM;
        double glue1 = glueThickness(z);
        double carbon = carbonThickness(z);
        double glue2 = glueThickness(z); // second glue layer
        double steel = steelThickness(z);

        file << fixed << setprecision(6)
             << z << ","
             << setprecision(3) << temp << ","
             << thermalProtection << ","
             << glue1 << ","
             << carbon << ","
             << glue2 << ","
             << steel << "\n";
    }

    file.close();
    cout << "CSV file written successfully: " << filename << endl;
}

int main() {
    double dz;
    cout << "Enter dz value (meters) for dividing 2.5m height (recommended: 0.01 to 0.05 m): ";
    cin >> dz;

    if (dz <= 0.0 || dz > TOTAL_HEIGHT) {
        cerr << "Invalid dz value. Must be > 0 and <= 2.5 meters." << endl;
        return 1;
    }

    writeThicknessProfileToCSV("thickness_profile.csv", dz);

    return 0;
}
