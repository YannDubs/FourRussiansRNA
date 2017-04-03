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
<<<<<<< HEAD:FourRussiansRNA/old yann/FourRussiansRNA.cpp
    int max = -1e5, sum1 = 0, sum2 = 0;
    for (int k = 0; k < q; k++) {
=======
    int max = 0, sum1 = 0, sum2 = 0;
    for (int k = int(q-1); k >= 0; k--) {
>>>>>>> origin/master:FourRussiansRNA/FourRussiansRNA.cpp
        if ((x & (1 << k)) != 0) sum1 = sum1 + 1;
        cout << "SUM1  :" << sum1 << endl;
        if ((y & (1 << k)) != 0) sum2 = sum2 - 1;
        cout << "SUM2  :" << sum2 << endl;
        cout << "SUM12  :" << (sum1+sum2) << endl;
        if (sum1 + sum2 > max) max = sum1 + sum2;
    }
    cout << "MAX  :" << max << endl << endl;
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
<<<<<<< HEAD:FourRussiansRNA/old yann/FourRussiansRNA.cpp
    int n (x.size());
    vector<vector<int>> D (n,vector<int> (n));
    const size_t q (2); //(round(log(n)));
    //cout << " q " << q << endl;
    //cout << " n " << n << endl;
=======
    size_t n (x.size());
    vector<vector<double>> D (n,vector<double> (n, -1));
    const size_t q (round(log(n)));
>>>>>>> origin/master:FourRussiansRNA/FourRussiansRNA.cpp
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
<<<<<<< HEAD:FourRussiansRNA/old yann/FourRussiansRNA.cpp
    for (int j(0); j < n; ++j){
        D[j][j] = 0;
        if (j < n-1){
            D[j+1][j] =-1;
        }
        //cout << " n " << n << endl;
        //cout << " j " << j << endl;
        for (int i(j-1); i >= 0; --i){
            //cout << " i " << i << endl;
            //cout << " j " << j << endl;
            D[i][j] = scoreB(x[i],x[j]) + D[i+1][j-1];
            int groupI ((i)/q);
            int groupJ ((j)/q);
            
            cout << endl <<   " i " << i << " j " << j << endl << " groupI " << groupI << " groupJ " << groupJ << endl ;
            
=======
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
>>>>>>> origin/master:FourRussiansRNA/FourRussiansRNA.cpp
            // the matrix is diagonal so groupI doesn't start always from the same position but from
            // the diagonal
            int nGroupsBetween (int(groupJ - groupI - 1));
            // + q to get right most. -1 to put back in 0 index
<<<<<<< HEAD:FourRussiansRNA/old yann/FourRussiansRNA.cpp
            int iI(q*groupI + q - 1);
            int jJ(q*groupJ );
            cout << " Ii " << iI << " jJ " << jJ << endl;
            cout << " nGroup " << nGroupsBetween << endl;
            
            //cout << "II" << endl;
            // for all cells in the first group
            // needs to take min iI and n because last group could have less columns than q
            for (size_t k(i+1); k <= min(iI,n-1); ++k){
                
                //cout << " iI " << iI << " i " << i << " j " << j << " k " << k << endl;
                D[i][j] = max (D[i][j], D[i][k-1] + D[k][j]);
=======
            size_t iI = min(q*groupI + q - 1, n-2);
            size_t jJ(q*groupJ);
            // for all cells in the first group
            for (size_t k(i); k <= iI; ++k){
                D[i][j] = max (D[i][j], D[i][k] + D[k+1][j]);
>>>>>>> origin/master:FourRussiansRNA/FourRussiansRNA.cpp
            }
            //for all cells in last group
<<<<<<< HEAD:FourRussiansRNA/old yann/FourRussiansRNA.cpp
            cout << "D[i][j] before block c: " <<  D[i][j] <<endl;
            for (int k(jJ); k <= j; ++k){
                //cout << "D[i][k-1] : " << D[i][k-1] <<endl;
                //cout << "D[k][j] : " << D[k][j] <<endl;
                //cout << " i " << i << " j " << j << " k " << k << endl;
                D[i][j] = max(D[i][j], D[i][k-1] + D[k][j]);
            }
            cout << "D[i][j] after block c: " <<  D[i][j] <<endl;
            
            /*
            cout << "I " << groupI << endl;
            cout << "J " << groupJ << endl;
            cout << "nGroup " << nGroupsBetween << endl;*/
            //for all groups in between
=======
            for (size_t k(jJ+1); k <= j; ++k){
                D[i][j] = max (D[i][j], D[i][k-1] + D[k][j]);
            }
>>>>>>> origin/master:FourRussiansRNA/FourRussiansRNA.cpp
            for (int K(0); K < nGroupsBetween; ++K){
                // take right most element of group I, adds 1 to get left most of next
                // then adds K*q to shift depedning on the block
                size_t l(iI + 1 + K*q);
<<<<<<< HEAD:FourRussiansRNA/old yann/FourRussiansRNA.cpp
                //not sure
                size_t t(l);
                cout << " l " << l << " t " << t << endl;
                cout << "D[i][l] : " << D[i][l]  <<endl;
                cout << "D[t][j] : " << D[t][j] <<endl;
                cout << "R[hvs[i][K]][vvs[j][K]] : " << R[hvs[i][K]][vvs[j][K]] <<endl;
                D[i][j] = max(D[i][j],D[i][l] + D[t][j] + R[hvs[i][K]][vvs[j][K]]);
=======
                size_t t(l+1);
                ull hdiff = hvs[i][(iI + q*K + 1)/q];
                ull vdiff = vvs[j][t/q];
                D[i][j] = max(D[i][j], D[i][l]+D[t][j]+R[hdiff][vdiff]);
>>>>>>> origin/master:FourRussiansRNA/FourRussiansRNA.cpp
            }
            cout <<"D[i][j] after block b : " <<D[i][j] <<endl;;
            
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
<<<<<<< HEAD:FourRussiansRNA/old yann/FourRussiansRNA.cpp
    
    cout << endl << "D : " << endl;
    for (int i = 0; i < n; i++){
        //cout << "Row number " << i << endl;
        for (int j = 0; j < n; j++){
            cout << D[i][j] << " ";
        }
        cout << endl;
    }
    
=======
>>>>>>> origin/master:FourRussiansRNA/FourRussiansRNA.cpp
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
