#include<iostream>
#include<typeinfo>
using namespace std;

int main(void){
    int var1 = 0, var2 = 1;
    int* p1, p2;
    cout << "addr1: " << &var1 << endl;
    p1 = &var1;
    cout << "p1: " << p1 << endl;
    cout << "*p1: " << *p1 << endl;
    cout << typeid(p2).name();
    return 0;
}