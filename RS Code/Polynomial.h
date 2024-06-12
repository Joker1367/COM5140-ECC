#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <math.h>
#include <string.h>

using namespace std;

#define MAX_Bit 1000
int n;
int k;
int r;

// GF(64) with a is a primitive element satisfying a^6 + a + 1

// pow_table[i] = a ^ i
/*vector<int> pow_table = {1, 2, 4, 8, 16, 32, 3, 6, 
                        12, 24, 48, 35, 5, 10, 20, 40, 
                        19, 38, 15, 30, 60, 59, 53, 41, 
                        17, 34, 7, 14, 28, 56, 51, 37, 
                        9, 18, 36, 11, 22, 44, 27, 54, 
                        47, 29, 58, 55, 45, 25, 50, 39, 
                        13, 26, 52, 43, 21, 42, 23, 46, 
                        31, 62, 63, 61, 57, 49, 33};

// log_table[i] = log_a i with log_a 0 = -1
vector<int> log_table = {-1, 0, 1, 6, 2, 12, 7, 26, 
                        3, 32, 13, 35, 8, 48, 27, 
                        18, 4, 24, 33, 16, 14, 52, 
                        36, 54, 9, 45, 49, 38, 28, 
                        41, 19, 56, 5, 62, 25, 11, 
                        34, 31, 17, 47, 15, 23, 53, 
                        51, 37, 44, 55, 40, 10, 61, 
                        46, 30, 50, 22, 39, 43, 29, 
                        60, 42, 21, 20, 59, 57, 58};*/

// coefficient of generator polynomial
/*vector<int> generator_coeff = {58, 62, 59, 7, 35, 58, 63, 47, 51, 6, 33, 
                               43, 44, 27, 7, 53, 39, 62, 52, 41, 44, 1};*/

vector<int> pow_table;
vector<int> log_table;
vector<int> generator_coeff;

void Initialize(int n, int k, int r, int primitive_poly) {

    pow_table.assign(n, 0);
    log_table.assign(n + 1, 0);
    generator_coeff.assign(r + 1, 0);

    log_table[0] = -1;

    pow_table[0] = 1;
    log_table[1] = 0;
    for(int i = 1; i < n; i++) {
        pow_table[i] = pow_table[i - 1] << 1;
        if(pow_table[i] & (n + 1)) {
            pow_table[i] ^= primitive_poly;
        }
        log_table[pow_table[i]] = i;
    }

    generator_coeff[0] = 1;
    for(int i = 1; i <= r; i++) {
        for(int j = i; j >= 1; j--) {
            if(generator_coeff[j] == 0) {
                generator_coeff[j] = generator_coeff[j - 1];
            }
            else {
                generator_coeff[j] = generator_coeff[j - 1] ^ pow_table[(log_table[generator_coeff[j]] + i) % n];
            }
        }
        generator_coeff[0] = pow_table[(log_table[generator_coeff[0]] + i) % n];
    }
}

// Addition in GF(n+1)
int GF_add(int a, int b){
    return a ^ b;
}

// multiplication in GF(n+1)
int GF_mul(int a, int b){
    if(a == 0 || b == 0) return 0;
    else{
        return pow_table[(log_table[a] + log_table[b]) % n];
    }
}

// Division in GF(64)
int GF_div(int a, int b){
    if(a == 0) return 0;
    else if(b == 0) {
        cout << "Divide by zero!!" << endl;
        return -1;
    }
    else return pow_table[(log_table[a] - log_table[b] + n) % n];
}

// Polynomial in GF64 : representing P(x) = P0 + P1 * x + P2 * x^2 ..... Pn * x^n
class Polynomial{
public:
    int degree;                         // Degree of P(x)
    vector<int> data;                   // data[i] = Pi
//constructors
    Polynomial(); 
    Polynomial(int);
    Polynomial(vector<int>);
    
//overloaded arithmetic operators as member functions
    Polynomial operator+(Polynomial);
    Polynomial operator*(Polynomial);
    Polynomial operator/(Polynomial);
    Polynomial operator%(Polynomial);
    Polynomial formal_derivative();
    int get_value(int);                 // compute P(a) if a is the input
    void left_shift();                  // P(x) -> P(x) * x
    void right_shift();                 // P(x) -> P(x) / x
    void Print();                       // Print P(x) (only coefficient)
};

Polynomial::Polynomial(){
    degree = 0;
    data.assign(MAX_Bit, 0);
    for(int i = 0; i < MAX_Bit; i++) data[i] = 0;
}

Polynomial::Polynomial(int x){
    degree = 0;
    data.assign(MAX_Bit, 0);
    data[0] = x;
}

Polynomial::Polynomial(vector<int> d){
    degree = d.size() - 1;
    data.assign(MAX_Bit, 0);
    for(int i = 0; i <= degree; i++) data[i] = d[i];
}

Polynomial Polynomial::operator+(Polynomial y){                             // res(x) = A(x) + B(x)
    Polynomial res;         
    int degree;
    int x_len = this->degree;
    int y_len = y.degree;
    for(degree = 0; degree <= x_len || degree <= y_len; degree++){
        res.data[degree] = GF_add(this->data[degree], y.data[degree]);    // res[i] = A[i] + B[i]    
    }
    while(degree >= 1 && res.data[degree] == 0) {                           // check prefix zero and update degree
        degree--;
    }
    res.degree = degree;
    return res;
}

Polynomial Polynomial::operator*(Polynomial y){                             // res(x) = A(x) * B(x)
    Polynomial res;
    int x_len = this->degree;
    int y_len = y.degree;
    int degree = x_len + y_len;
    for(int i = 0; i <= y_len; i++){                                        // res[i] = sum(A[j] * B[i - j])
        for(int j = 0; j <= x_len; j++) {
            res.data[i + j] = GF_add(res.data[i + j], GF_mul(this->data[j], y.data[i]));
        }
    }
    while(degree >= 1 && res.data[degree] == 0) {                           // check prefix zero and update degree
        degree--;
    }
    res.degree = degree;
    return res;
}

Polynomial Polynomial::operator/(Polynomial y){                             // A(x) = B(x) * res(x) + t(x)
    Polynomial t, tmp, res;
    if(this->degree < y.degree) return res;                                 // if deg(A(x)) < deg(B(x)) then res(x) = 0

    int i;
    int r_len = 0;
    t.degree = y.degree;
    for(i = 0; i <= y.degree; i++){                                         // 長除法 (long division)
        t.data[y.degree - i] = this->data[this->degree - i];
    }
    while(true){
        if(t.degree == y.degree){
            res.data[0] = GF_div(t.data[t.degree], y.data[y.degree]);
            t = t + y * res.data[0];
        }
        if(i <= this->degree){
            t.left_shift();
            t.data[0] = this->data[this->degree - i];
            res.left_shift();
            i++;
        }
        else break;
    }
    return res;
}

Polynomial Polynomial::operator%(Polynomial y){                             // A(x) = B(x) * q(x) + res(x)
    Polynomial t, tmp, res;
    if(this->degree < y.degree) return *this;                               // if deg(A(x)) < deg(B(x)) then res(x) = A(x) 

    int i;
    int q;
    int r_len = 0;
    t.degree = y.degree;
    for(i = 0; i <= y.degree; i++){                                         // 長除法 (long division)
        t.data[y.degree - i] = this->data[this->degree - i];
    }
    while(true){
        if(t.degree == y.degree){
            res.data[0] = GF_div(t.data[t.degree], y.data[y.degree]);
            t = t + y * res.data[0];
        }
        if(i <= this->degree){
            t.left_shift();
            t.data[0] = this->data[this->degree - i];
            res.left_shift();
            i++;
        }
        else break;
    }
    return t;
}

Polynomial Polynomial::formal_derivative(){                                 // A(x) -> A'(x)
    Polynomial res;
    if(this->degree == 0) return res;
    else{
        res.degree = this->degree - 1;
        for(int i = 0; i <= res.degree; i++){                               // An * x^n -> n * An * x^(n - 1)
            if(i % 2 == 0) res.data[i] = this->data[i + 1];                 // GF(64) has characteristic 2    
            else res.data[i] = 0;
        }
        while(res.data[res.degree] == 0) res.degree--;
        return res;
    }
}

int Polynomial::get_value(int alpha){                                       // Compute A(alpha)
    int pow = alpha;
    int res = data[degree];                                                 // Horner's rule
    for(int i = degree - 1; i >= 0; i--){                                   
        res = GF_add(data[i], GF_mul(res, pow));
    }
    return res;
}

void Polynomial::Print(){
    for(int i = 0; i <= degree; i++) cout << setw(2) << setfill(' ') << data[i] << " ";
    cout << endl;
}

void Polynomial::left_shift(){
    if(this->degree == 0 && this->data[0] == 0) return;
    for(int i = this->degree; i >= 0; i--) this->data[i + 1] = this->data[i];
    this->data[0] = 0;
    this->degree++;
}

void Polynomial::right_shift(){
    for(int i = 1; i <= this->degree; i++) this->data[i - 1] = this->data[i];
    this->data[this->degree] = 0;
    if(this->degree > 0) this->degree--;
}