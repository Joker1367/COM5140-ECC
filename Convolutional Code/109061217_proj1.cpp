#include <iostream>
#include <iomanip>
#include <bitset>
#include <cmath>
#include <limits>
#include <float.h>
#include <stdlib.h>
#include <vector>
#include <fstream>

#define TR_LEN 32
#define INFO_PERIOD 63

using namespace std;

const long long int para_1 = 4101842887655102017LL;
const long long int para_2 = 2685821657736338717LL;
const double para_3 = 5.42101086242752217E-20;
unsigned long long int SEED = 14;
unsigned long long int RANV;
int RANI = 0;

double Ranq1(){
    if(RANI == 0){
        RANV = SEED ^ para_1;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV *= para_2;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;
    return RANV * para_2 * para_3;
}
void Normal(double& n1, double& n2, double std_dev){
    double x1, x2, s;
    do{
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while (s >= 1.0);
    n1 = std_dev * x1 * sqrt(-2 * log(s) / s);
    n2 = std_dev * x2 * sqrt(-2 * log(s) / s);
}

// each node represent a state, storing the current information and previous information
typedef struct Node{
    double prev_metric = DBL_MAX;               // metric in previous state
    double cur_metric = DBL_MAX;                // metric in current state
    int prev_node = 100000;                     // index of previous node
    vector<vector<int>> next;                   // next[0][0] = index of next state when u = 0
                                                // next[1][0] = index of next state when u = 1
    bitset<TR_LEN> pre_path {0};                    // each bit represent u that was sent previously, in previous state 
    bitset<TR_LEN> path {0};                        // each bit represent u that was sent previously, in current state
} node;

double SNR[20] = {1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8};

int main(void){
    int ERROR = 0;                              // # of bit error
    int u, tmp;                                 // temporary storage
    int m1, m2;                                 // [m1, m2] = uG
    int num_of_0, num_of_1;                     // # of 0 and 1 while we using the Majority Vote
    int truncated_len = 0;
    double x1, x2;                              // Signal point after BPSK and adding noise
    double STD_DEV;                             // standard deviation correspond to SNR to generate Normal r.v
    double bias;                                // the average of the metric to prevent overflow
    int index = 0;                              
    long long int information_index = 0;        // # of transmitted information bit
    double n1, n2;                              // Noise
    double M1, M2, M;
    ofstream outFile;
    node* trellis = new node[64];
    vector<int> information;                    // information bits
    vector<int> state;                          // current state

    // [m1, m2] = uG
    // m1 = u + s[2] + s[3] + s[5] + s[6]
    // m2 = u + s[1] + s[2] + s[3] + s[6]
    for(int i = 0; i < 64; i++){
        m1 = 0;
        m2 = 0;
        if(i & 1) m2 = m2 + 1;
        if(i & 2) {
            m1 = m1 + 1;
            m2 = m2 + 1;
        }
        if(i & 4) {
            m1 = m1 + 1;
            m2 = m2 + 1;
        }
        if(i & 16) m1 += 1;
        if(i & 32) {
            m1 = m1 + 1;
            m2 = m2 + 1;
        }
        index = i * 2;
        if(index >= 64) index -= 64;
        trellis[i].next.push_back({index, m1 % 2, m2 % 2});
        trellis[i].next.push_back({index + 1, (m1 + 1) % 2, (m2 + 1) % 2});
    }
    state.assign(7, 0);
    information.assign(63, 0);

    // generate information
    // u[l + 6] = u[l + 1] + u[l] for l >= 0. The period of the sequence is 63
    information[0] = 1;
    for(int i = 6; i < INFO_PERIOD; i++) information[i] = information[i - 6] ^ information[i - 5];

    // Start Testing for each SNR
    for(int testcase = 0; testcase < 15; testcase++){
        STD_DEV = sqrt(pow(10, -SNR[testcase] / 10));            // Calculate standard deviation correspond to SNR
        ERROR = 0;                                               // Initialize # of error
        information_index = 0;                                   // Initialize # of transmitted information bit
        for(int i = 0; i < 64; i++){                             // Initialize each state node
            trellis[i].prev_metric = DBL_MAX;
            trellis[i].cur_metric = DBL_MAX;
            trellis[i].path = 0;
            trellis[i].pre_path = 0;
        }
        trellis[0].prev_metric = 0;                              // Start from all zero state
        cout << "SNR = " << SNR[testcase] << "dB\n";                
        while(ERROR < 1000){
            u = information[information_index % INFO_PERIOD];               // Transmit an information bit
            m1 = (u + state[2] + state[3] + state[5] + state[6]) % 2;       // Encode 
            m2 = (u + state[1] + state[2] + state[3] + state[6]) % 2;
            state[6] = state[5];                                            // Update State
            state[5] = state[4];
            state[4] = state[3];
            state[3] = state[2];
            state[2] = state[1];
            state[1] = u;
            Normal(n1, n2, STD_DEV);                              // Perform BPSK and add noise 
            x1 = -2 * m1 + 1 + n1;
            x2 = -2 * m2 + 1 + n2;
            /*if(x1 >= 0) x1 = 1;                                 // Receiving the Signal with Demodulation (Soft decision / Hard Decision)
            else x1 = -1;
            if(x2 >= 0) x2 = 1;
            else x2 = -1;*/

            // Updating the information in trellis diagram
            for(int i = 0; i < 64; i++){
                if(trellis[i].prev_metric == DBL_MAX) continue;   // Check whether a state node is already reached

                // Calculate the metric send from the previous state node while u = 0
                index = trellis[i].next[0][0];    
                if(trellis[i].next[0][1] == 0) M1 = (x1 - 1) * (x1 - 1);
                else M1 = (x1 + 1) * (x1 + 1);
                if(trellis[i].next[0][2] == 0) M2 = (x2 - 1) * (x2 - 1);
                else M2 = (x2 + 1) * (x2 + 1);
                M = trellis[i].prev_metric + M1 + M2;
                if(trellis[index].cur_metric > M){                // Updating information in the state node
                    trellis[index].cur_metric = M;
                    trellis[index].prev_node = i;
                }

                // Calculate the metric send from the previous state node while u = 0
                index = trellis[i].next[1][0];                     
                if(trellis[i].next[1][1] == 0) M1 = (x1 - 1) * (x1 - 1);
                else M1 = (x1 + 1) * (x1 + 1);
                if(trellis[i].next[1][2] == 0) M2 = (x2 - 1) * (x2 - 1);
                else M2 = (x2 + 1) * (x2 + 1);
                M = trellis[i].prev_metric + M1 + M2;                       
                if(trellis[index].cur_metric > M){                // Updating information in the state node 
                    trellis[index].cur_metric = M;
                    trellis[index].prev_node = i;
                }
            }
            bias = 0;                                             // Initialize bias
            tmp = 0;                                              // Initialize # of reached state node
            for(int i = 0; i < 64; i++){
                if(trellis[i].cur_metric == DBL_MAX) continue;    // check whether the state node is reached
                tmp++;
                trellis[i].prev_metric = trellis[i].cur_metric;   // Update Metric
                bias += trellis[i].cur_metric;                              
                trellis[i].cur_metric = DBL_MAX;                            
                index = trellis[i].prev_node;                     // Update the path store in the state node
                if(i == trellis[index].next[0][0]) trellis[i].path = trellis[index].pre_path << 1;
                else if(i == trellis[index].next[1][0]) {
                    trellis[i].path = trellis[index].pre_path << 1;
                    trellis[i].path[0] = 1;
                }
                else cout << "ERROR!\n";
            }
            for(int i = 0; i < 64; i++) trellis[i].pre_path = trellis[i].path;
            bias = bias / tmp;
            for(int i = 0; i < 64; i++) {
                if(trellis[i].prev_metric == DBL_MAX) continue;   // remove the dc term(average) of metric in each state node
                trellis[i].prev_metric -= bias;
            }

            if(information_index >= TR_LEN - 1){
                // Best state
                index = 0;
                double best_metric = trellis[0].prev_metric;       // Output the best state result
                for(int i = 1; i < 64; i++){
                    if(trellis[i].prev_metric < best_metric){
                        index = i;
                        best_metric = trellis[i].prev_metric;
                    }
                }
                // Check if there's error
                if(trellis[index].path[TR_LEN - 1] != information[(information_index - TR_LEN + 1) % INFO_PERIOD]) ERROR++;

                // Fixed State
                //if(trellis[0].path[TR_LEN - 1] != information[(information_index - TR_LEN + 1) % INFO_PERIOD]) ERROR++;

                // Majority vote
                /*num_of_0 = 0;
                num_of_1 = 1;
                for(int i = 0; i < 64; i++){
                    if(trellis[i].path[31] == 0) num_of_0++;
                    else num_of_1++;
                }
                if(num_of_0 >= num_of_1 && information[(information_index - 31) % INFO_PERIOD] == 1) ERROR++;
                else if(num_of_0 < num_of_1 && information[(information_index - 31) % INFO_PERIOD] == 0) ERROR++;*/
            }
            information_index = (information_index + 1);
        }
        cout << "ERROR = " << ERROR << " Number of bits = " << information_index - 30 << endl;
        cout << "BER = " << ERROR * 1.0 / (information_index - 30) << endl;
    }
    return 0;
}
