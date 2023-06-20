#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <string.h>
#include <fstream>
#include <time.h>
#include "Polynomial.h"

using namespace std;

// Extended Euclidean Algorithm
void EEA(Polynomial P, Polynomial Q, Polynomial& s0, Polynomial& s1, Polynomial& t0, Polynomial& t1, Polynomial& W, int e0);

int main(void){
    int x;
    int e0, e1;
    int count;
    bool failure = false;
    ifstream inFile;
    ofstream outFile;
    Polynomial generator(generator_coeff);          // generator polynomial g(x)
    Polynomial information;                         // Information bits I(x)
    Polynomial codeword;                            // Codeword C(x) (with error and erasure)
    Polynomial Answer;                              // correct codewrod A(x)          
    Polynomial received_vector;                     // received vector R(x)
    Polynomial error;                               // Error E(x)
    Polynomial syndrome;                            // Syndrome S(x)
    Polynomial s0, s1, t0, t1;                      // Polynomial for EEA
    Polynomial tmp;                                 // temporary polynomial
    Polynomial I0, I1;                              // sigma0(x) and its formal derivative
    Polynomial W;                                   // error evaluator polynomial omega(x)
    Polynomial xr;                                  // x^r

// Encoding
    // Initilization
    information.degree = k - 1;
    Answer.degree = n - 1;
    codeword.degree = n - 1;
    received_vector.degree = n - 1;
    error.degree = n - 1;
    syndrome.degree = r - 1;
    I0.degree = 1;
    tmp.degree = 1;
    srand(time(NULL));

    // Generate testcase
    outFile.open("testcase.txt");
    for(int testcase = 0; testcase < 100; testcase++){
        for(int i = 0; i <= k - 1; i++) information.data[i]  = rand() % (n + 1);
        codeword = information * generator;                                         // C(x) = I(x) * g(x)
        for(int i = 0; i <= n - 1; i++) outFile << codeword.data[i] << " ";
        outFile << endl;
        for(int i = 0; i <= n - 1; i++) {                                           // Add error and erasure
            int randi = rand() % 10;
            if(randi == 0) outFile << rand() % 64 << " ";
            else if(randi == 1) outFile << "*" << " ";
            else outFile << codeword.data[i] << " ";
        }
        outFile << endl;
    }
    outFile.close();

// Decoding
    inFile.open("testcase.txt");
    string s;
    for(int testcase = 0; testcase < 10; testcase++){
        // Initialization
        e0 = 0; e1 = 0;                                                             
        Answer.degree = n - 1;
        codeword.degree = n - 1;
        received_vector.degree = n - 1;
        error.degree = n - 1;
        syndrome.degree = r - 1;
        xr.degree = r;
        tmp.degree = 1;
        I0 = 1;
        s0 = 1; s1 = 0; 
        t0 = 0; t1 = 1;
        //s0.Print(); s1.Print(); t0.Print(); t1.Print(); I0.Print();
        for(int i = 0; i <= n - 1; i++){
            error.data[i] = 0;
            syndrome.data[i] = 0;
            xr.data[i] = 0;
        }
        xr.data[r] = 1;
        for(int i = 0; i <= n - 1; i++) inFile >> Answer.data[i]; Answer.Print();
        for(int i = 0; i <= n - 1; i++){                                            // Compute R'(x) and sigma0(x)
            inFile >> s;
            cout << s << " ";
            if(s == "*") {
                tmp.data[0] = 1; tmp.data[1] = pow_table[i];                        
                I0 = I0 * tmp;
                e0++;                                                               // e0 = # of erasure
                codeword.data[i] = 0;                                               
            }
            else {
                codeword.data[i] = stoi(s);
                if(codeword.data[i] != Answer.data[i]) e1++;
            }
        }
        cout << endl;
        cout << "e0 = " << e0 << "; e1 = " << e1 << endl;
        if(e0 > r) {
            cout << "failure" << endl;
            continue;
        }
        received_vector = codeword;                                                 
        received_vector = received_vector % generator;                              // modulo g(x) for faster calculation of R(alpha^i)
        
        // Compute Syndrome
        for(int i = 0; i <= r - 1; i++) {                                           // Sj = R(alpha^j)
            syndrome.data[i] = received_vector.get_value(pow_table[i + 1]);     
        }
        syndrome.degree = r - 1;
        while(syndrome.data[syndrome.degree] == 0) syndrome.degree--;
        syndrome = (syndrome * I0) % xr;                                            // S0(x) = sigma0(x) * S(x) (mod x^r)
        EEA(xr, syndrome, s0, s1, t0, t1, W, e0);                                   // Perform EEA to find sigma1(x) and omega(x)
                                                                                    // t1 = sigma1(x) and W = omega(x)
        I0 = I0 * t1;                                                               // sigma(x) = sigma0(x) * sigma1(x)
        I1 = I0.formal_derivative();                                                // Compute the formal derivative

        // Time Domain Approach
        failure = false;                                                            // boolean variable for decode failure or not
        if(I0.data[0] == 0 || W.degree >= e0 + t1.degree) failure = true;
        else{
            count = 0;
            for(int i = 0; i <= n - 1; i++){
                x = GF64_div(1, pow_table[i]);
                if(I0.get_value(x) == 0 && I1.get_value(x) != 0){
                    count++;
                    error.data[i] = GF64_div(W.get_value(x), I1.get_value(x));      // Ei = -omega(alpha ^ -i) / sigma'(alpha ^ -i)
                }
                else error.data[i] = 0;
            }
            if(count != I0.degree) failure = true;
        }
        if(failure) cout << "failure!" << endl;
        else{                                                                       
            codeword = codeword + error;                                            // C(x) = R(x) - E(x)
            for(int i = 0; i <= n - 1; i++){                                        // Compare C(x) with A(x)
                if(codeword.data[i] != Answer.data[i]){
                    codeword.Print();
                    Answer.Print();
                    system("pause");
                }
            }
            cout << "Testcase " << testcase << " pass!" << endl;
        }
        system("pause");
    }
    inFile.close();
    return 0;
}

void EEA(Polynomial P, Polynomial Q, Polynomial& s0, Polynomial& s1, Polynomial& t0, Polynomial& t1, Polynomial& W, int e0){
    int u = (r - e0) / 2;                       // u = ceil(r - e0 / 2)
    int v = r - 1 - u;                          // v = r - 1 - u
    Polynomial q, tmp, tmps, tmpt;
    while(Q.degree > v || t1.degree > u){       // terminate condition : deg(rj(x) <= v) and deg(vj(x)) <= u
        q = P / Q;                              // q = P / Q
        tmp = P % Q;                            // r_j+1 = r_j-1 - q * r_j
        P = Q;                              
        Q = tmp;
        tmps = s0 + q * s1;                     // u_j+1 = u_j-1 - q * u_j
        s0 = s1; s1 = tmps;
        tmpt = t0 + q * t1;                     // v_j+1 = v_j-1 - q * v_j
        t0 = t1; t1 = tmpt;
    }
    W = Q;
}