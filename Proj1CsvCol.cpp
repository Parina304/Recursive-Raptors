
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
     stringstream ss("");
    ifstream file(filename);    // Input file stream to read the file
    int ncol = 0, flag = 0;
    int ridx;

    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return data; // Return empty vector if file cannot be opened
    }

    string line;

    while (getline(file, line)) {
        // cout << line;
        // ss << line;
        vector<string> col; // Vector to store one column of the CSV
        stringstream ss(line); // String stream to split the line
        string cell;

        ridx = 0;
        while (getline(ss, cell, ',')) {
            if(!flag){
                // if first line, init data
                data.push_back(vector<string>());
                data[ncol++].push_back(cell);
            }
            else{
                // else push back to columns
                data[ridx++].push_back(cell);
            }
        }
        // turn off first line flag
        flag = 1;
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

void printDataFormat(const vector<vector<string>> &data, const int &cidx){
    printf("%s data (units in %s): ", data[cidx][0].c_str(), data[cidx][1].c_str());
    for (int i = 2; i < data[cidx].size(); i++)printf("%s ", data[cidx][i].c_str());
    printf("\n");
}

int main() {
    // Prompt user for the CSV file name
    string filename;
    int cidx;
    cout << "Enter the CSV file name: ";
    cin >> filename;

// Read the CSV file into a 2D vector
    vector<vector<string>> csvData = readCSV(filename);

// Check if data was successfully read
    if (csvData.empty()) {
        cout << "No data found or failed to read the file." << endl;
    } else {
        cout << "CSV file loaded." << endl;
    // Print the contents of the CSV file
        // print2DVector(csvData);
    }
    while(1){
        printf("Enter column number: ");
        cin >> cidx;
        printDataFormat(csvData, cidx);
    }
    return 0;
}