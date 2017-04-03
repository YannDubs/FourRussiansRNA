//
//  FourRussiansRNA.cpp
//  FourRussiansRNA
//
//  Created by Yann Dubois on 01.04.17.
//  Copyright © 2017 Yann Dubois. All rights reserved.
//

#include "FourRussiansRNA.hpp"
#include <fstream>
#include <iostream>
#include <math.h>

//CONSTANTS
const string folder("/Users/yanndubois/Desktop/GitHub/FourRussiansRNA/Data/");

//CODE

// scoring function
int scoreB(char a, char b){
    // returns 1 if match 0 else
    if ((a == 'a' && b == 'u') || (a == 'u' && b == 'a') || (a == 'c' && b == 'g') || (b == 'c' && a == 'g')){
        //cout << " return " << 1 << endl;
        return 1;
    }
    return 0;
}

// preprocessing helper
int maxVal(ull x, ull y, const size_t q) {
    int max = 0, sum1 = 0, sum2 = 0;
    for (int k = int(q-1); k >= 0; k--) {
        if ((x & (1 << k)) != 0) sum1 = sum1 + 1;
        if ((y & (1 << k)) != 0) sum2 = sum2 - 1;
        if (sum1 + sum2 > max) max = sum1 + sum2;
    }
    return max;
}

void debugPrint(vector<vector<double> > &m) {
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[0].size(); j++) {
            cout << m[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

//Four russians running
double nussinovFourRussians(const string& x){
    // INITIALIZATION
    size_t n (x.size());
    vector<vector<double>> D (n,vector<double> (n, -1));
    const size_t q (round(log(n)));
    //Index
    //vector<vector<size_t>> Index (m,vector<size_t> (n));
    
    // preprocessing step for table R
    size_t qsq = pow(2, q);
    vvi R(qsq, vi(qsq));
    for (ull x = 0; x < (1 << q); ++x) {
        for (ull y = 0; y < (1 << q); ++y) {
            // x and y are horizontal and vertical difference bit vectors
            // represented by unsigned long longs
            R[x][y] = maxVal(x, y, q);
        }
    }
    
    int max_q = int(ceil(n/q))+1;
    vvu hvs(n, vu(max_q)); // horizontal diff vector store
    vvu vvs(n, vu(max_q)); // vertical diff vector store
    
    // ITERATION
    for (size_t j(0); j < n; ++j){
        for (long int i(j); i >= 0; --i){
            if (j - i <= 0) {
                D[i][j] = 0;
                continue;
            }
            D[i][j] = max(max(
                              D[i+1][j],
                              D[i][j-1]),
                          scoreB(x[i],x[j]) + D[i+1][j-1]);
            size_t groupI ((i)/q);
            size_t groupJ ((j)/q);
            // the matrix is diagonal so groupI doesn't start always from the same position but from
            // the diagonal
            int nGroupsBetween (int(groupJ - groupI - 1));
            // + q to get right most. -1 to put back in 0 index
            size_t iI = min(q*groupI + q - 1, n-2);
            size_t jJ(q*groupJ);
            // for all cells in the first group
            for (size_t k(i); k <= iI; ++k){
                D[i][j] = max (D[i][j], D[i][k] + D[k+1][j]);
            }
            //for all cells in last group
            for (size_t k(jJ+1); k <= j; ++k){
                D[i][j] = max (D[i][j], D[i][k-1] + D[k][j]);
            }
            for (int K(0); K < nGroupsBetween; ++K){
                // take right most element of group I, adds 1 to get left most of next
                // then adds K*q to shift depedning on the block
                size_t l(iI + 1 + K*q);
                size_t t(l+1);
                ull hdiff = hvs[i][(iI + q*K + 1)/q];
                ull vdiff = vvs[j][t/q];
                D[i][j] = max(D[i][j], D[i][l]+D[t][j]+R[hdiff][vdiff]);
            }
            
            // compute the vertical difference vector
            if (i % q == 1 && i + q <= n) {
                // compute and store the v¯ vector i/qth group for column j
                ull vdiff = 0; size_t c = q-2;
                for (size_t k(i); k < i+q-1; ++k) {
                    if (D[k][j] - D[k+1][j] == 1) {
                        vdiff = (vdiff | (1 << c));
                    } c--;
                }
                vvs[j][groupI] = vdiff; // i/qth group for column j
            }
            
            // compute the horizontal difference vector
            if (j % q == q - 1) {
                // compute and store the v vector (j − 1)/qth group for row i
                ull hdiff = 0; size_t c = q-2;
                for (size_t k(j+1-q); k <= j; ++k) {
                    if (D[i][k+1] - D[i][k] == 1) {
                        hdiff = (hdiff | (1 << c));
                    } c--;
                }
                hvs[i][groupJ] = hdiff; // (j − 1)/qth group for row i
            }                   
        }
    }
    cout << D[0][n-1] << endl;
    return D[0][n-1];
    //TRACEBACK
}

string LoadSeq(string file){
    
    ifstream in(folder + file);
    
    if (!in) {
        cout << "Cannot open file.\n";
        return "1";
    }
    string str { istreambuf_iterator<char>(in), istreambuf_iterator<char>() };
    
    in.close();
    
    return str;
}
