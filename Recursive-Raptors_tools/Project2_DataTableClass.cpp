#include <vector> 
#include <iostream> 
using namespace std; 

class MyDataTable {

    private: 
        vector<vector <string>> table; // stores the data  

    public: 
        // Constructor 
        MyDataTable (char*** dataInput, int numCols, int* colSizes) {

            table.resize(numCols); 
            for (int col = 0; col < numCols; ++col) {
                for (int row = 0; row < colSizes[col]; ++row) {
                    table[col].push_back(dataInput[col][row]); 
                }
            }
        }

        ostream& printTable(ostream& os = cout) const {
            int maxRows = 0; // so I know how many empty spaces I need to have everything work 
            for (const auto& col : table) {
                maxRows = max(maxRows, (int)col.size()); // casted as an int bc it is currently a sizet datatype 
            }

            // fill out the data in the correct places 
            // cycle thru cols within the same row to print nicely 
            for (int row = 0; row < maxRows; ++row) {
                for (const auto& col : table) { // auto means it figures out the datatype during compile-time 
                    if (row < col.size()) {
                        os << col[row] << "\t\t";
                    } else {
                        os << " \t\t"; // Empty cell
                    }
                }
                os << endl;
            }
        return os; 
    }

    friend ostream& operator<<(ostream& os, const MyDataTable& table) {
        os << "In the overloading function" << endl;
        return table.printTable();
    }

};



int main () {

    char* col1[] = { (char*)"Temperature", (char*)"30", (char*)"15", (char*)"60", (char*)"25" };
    char* col2[] = { (char*)"Humidity", (char*)"80", (char*)"50", (char*)"70" };
    char* col3[] = { (char*)"Pressure", (char*)"1013", (char*)"1005" };

    char** columns[] = { col1, col2, col3 };
    int sizes[] = { 5, 4, 3 };

    MyDataTable test_table(columns, 3, sizes);
    cout << "Table Output:\n";
    test_table.printTable();
    
    cout << test_table; 

    return 0; 
} 