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
    int max = -1e5, sum1 = 0, sum2 = 0;
    for (int k = 0; k < q; k++) {
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

//Four russians running
void nussinovFourRussians(const string& x){
    // INITIALIZATION
    int n (x.size());
    vector<vector<int>> D (n,vector<int> (n));
    const size_t q (2); //(round(log(n)));
    //cout << " q " << q << endl;
    //cout << " n " << n << endl;
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
            // cout << R[x][y] << " " << x << " " << y << endl;
        }
    }
    
    int max_q = int(ceil(n/q))+1;
    vvu hvs(n, vu(max_q)); // horizontal diff vector store
    vvu vvs(n, vu(max_q)); // vertical diff vector store
    
    
    // ITERATION
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
            
            // the matrix is diagonal so groupI doesn't start always from the same position but from
            // the diagonal
            int nGroupsBetween (int(groupJ - groupI - 1));
            // + q to get right most. -1 to put back in 0 index
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
            }
            
            //cout << "JJ" << endl;
            //for all cells in last group
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
            for (int K(0); K < nGroupsBetween; ++K){
                /*cout << "K " << K << endl;
                cout << "nGroup " << nGroupsBetween << endl;
                cout << "K<nGroup " << (K < nGroupsBetween) << endl;*/
                // take right most element of group I, adds 1 to get left most of next
                // then adds K*q to shift depedning on the block
                size_t l(iI + 1 + K*q);
                //not sure
                size_t t(l);
                cout << " l " << l << " t " << t << endl;
                cout << "D[i][l] : " << D[i][l]  <<endl;
                cout << "D[t][j] : " << D[t][j] <<endl;
                cout << "R[hvs[i][K]][vvs[j][K]] : " << R[hvs[i][K]][vvs[j][K]] <<endl;
                D[i][j] = max(D[i][j],D[i][l] + D[t][j] + R[hvs[i][K]][vvs[j][K]]);
            }
            cout <<"D[i][j] after block b : " <<D[i][j] <<endl;;
            
            // compute the vertical difference vector
            if (i % q == 1) {
                // compute and store the v¯ vector i/qth group for column j
                ull vdiff = 0; size_t c = q-2;
                for (size_t k(i); k < i+q-1; ++k) {
                    if (D[k-1][j] - D[k][j] == 1) {
                        vdiff = (vdiff | (1 << c));
                    }
                    c--;
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
                    }
                    c--;
                }
                hvs[i][groupJ] = hdiff; // (j − 1)/qth group for row i
            }
        }
    }
    
    cout << endl << "D : " << endl;
    for (int i = 0; i < n; i++){
        //cout << "Row number " << i << endl;
        for (int j = 0; j < n; j++){
            cout << D[i][j] << " ";
        }
        cout << endl;
    }
    
    cout << D[0][n-1] << endl;
    
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


/*
 // these are the groups in which (i,j) would be if the matrix wasn't diagonal
 size_t groupRow ((j+1)/q);
 // the -1 shifts j so that j = 0 is in group 1 (we use groups in base 1)
 size_t groupColumn ((i+1)/q + 1);



 // the matrix is diagonal so groupI doesn't start always from the same position but from
 // the diagonal
 size_t groupHorizontal (groupRow - groupColumn + 1);
 size_t groupVertical (groupColumn);

 // minus 1 is to put back in 0 index. + q at the end i to get right most entry
 size_t iI (q*groupRow -1 + q);
 size_t jJ (q*groupRow -1 );
 */

/*
 void bitsetTimingExperiment() {
   int q = 12;
   int qsq = pow(2, q);
   vvi R(qsq, vi(qsq));
   clock_t start;
   start = clock();
   for (int x = 0; x < (1 << q); ++x) {
     for (int y = 0; y < (1 << q); ++y) {
       bitset<12> xb(x); bitset<12> yb(y);
       int max = 0, sum1 = 0, sum2 = 0;
       for (int k = 0; k < q; k++) {
         sum1 = sum1 + xb[k];
         sum2 = sum2 - yb[k];
         if (sum1 + sum2 > max) max = sum1 + sum2;
       }
       R[x][y] = max;
     }
   }
   cout << "Time: " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;
   //  Time: 4624.03 ms
   start = clock();
   for (int x = 0; x < (1 << q); ++x) {
     for (int y = 0; y < (1 << q); ++y) {
       R[x][y] = maxval(x, y, q);
       // cout << R[x][y] << " " << x << " " << y << endl;
     }
   }
   cout << "Time: " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;
   // Time: 1427.21 ms
 }
 */

/*
void exampleBitComputationTest() {
  int q = 5;
  int n = 1000;
  int max_q = int(ceil(n/q))+1;
  vvu hvs(n, vu(max_q)); // horizontal diff vector store
  vvu vvs(n, vu(max_q)); // vertical diff vector store

  vector<vector<double>> D (n,vector<double> (n));

  int i; int j;

  j = 10;
  // vertical difference
  D[4][j] = 8;
  D[5][j] = 8;
  D[6][j] = 7;
  D[7][j] = 6;
  // binary output should be 0 1 1

  i = 10;
  // horizontal difference
  D[i][8] = 5;
  D[i][9] = 6;
  D[i][10] = 7;
  D[i][11] = 7;
  // binary output should be 1 1 0

  int groupI = 20;
  int groupJ = 20;

  i = 5; j = 10;
  q = 4;
  // compute the vertical difference vector
  // i = 5. q = 4. 5 % 4 == 1
  if (i % q == 1) {
    // compute and store the v¯ vector i/qth group for column j
    ull vdiff = 0; size_t c = q-2;
    for (size_t k(i); k < i+q-1; ++k) {
      if (D[k-1][j] - D[k][j] == 1) {
        vdiff = (vdiff | (1 << c));
      }
      c--;
    }
    vvs[j][groupI] = vdiff; // i/qth group for column j
    cout << vvs[j][groupI] << endl;
    // 3 == 011
    if (vvs[j][groupI] != 3) {
      cout << "incorrect vertical bit vector computed" << endl;
    }
  }

  i = 10; j = 11;

  // compute the horizontal difference vector
  if (j % q == q - 1) {
    // compute and store the v vector (j − 1)/qth group for row i
    ull hdiff = 0; size_t c = q-2;
    for (size_t k(j+1-q); k <= j; ++k) {
      if (D[i][k+1] - D[i][k] == 1) {
        hdiff = (hdiff | (1 << c));
      }
      c--;
    }
    hvs[i][groupJ] = hdiff; // (j − 1)/qth group for row i
    cout << hvs[i][groupJ] << endl;
    // 6 == 110
    if (hvs[i][groupJ] != 6) {
      cout << "incorrect horizontal bit vector computed" << endl;
    }
  }
}
*/
