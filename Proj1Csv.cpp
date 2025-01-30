#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream> // For string stream to parse CSV
#include <cstring>

using namespace std;

/**
* Reads a CSV file and stores its contents in a 2D vector.
* @param filename The name of the CSV file to read.
* @return A 2D vector containing the contents of the CSV file.
*/
vector<vector<string>> readCSV(const string &filename) {
     vector<vector<string>> data; // 2D vector to store CSV content
    ifstream file(filename);    // Input file stream to read the file

    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return data; // Return empty vector if file cannot be opened
    }

    string line;
    // Read the file line by line
    while (getline(file, line)) {
        vector<string> row; // Vector to store one row of the CSV
        stringstream ss(line); // String stream to split the line
        string cell;

         // Split the line by commas and store in the row vector
        while (getline(ss, cell, ',')) {
            row.push_back(cell);
        }

        data.push_back(row); // Add the row to the 2D vector
    }

    file.close(); // Close the file
    return data;
}

/**
* Prints a 2D vector to the console.
* @param data The 2D vector to print.
*/
void print2DVector(const vector<vector<string>> &data) {
    for (const auto &row : data) {
        for (const auto &cell : row) {
            cout << cell << " "; // Print each cell separated by a space
        }
        cout << endl; // Move to the next line after each row
    }
}

int main() {
    // Prompt user for the CSV file name
    string filename;
    cout << "Enter the CSV file name: ";
    cin >> filename;

// Read the CSV file into a 2D vector
    vector<vector<string>> csvData = readCSV(filename);

// Check if data was successfully read
    if (csvData.empty()) {
        cout << "No data found or failed to read the file." << endl;
    } else {
        cout << "CSV file contents:" << endl;
    // Print the contents of the CSV file
        print2DVector(csvData);
    }

    return 0;
}
