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
    //cout << " q " << q << endl;
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
    for (size_t j(0); j < n; ++j){
        //        cout << " n " << n << endl;
        //cout << " j " << j << endl;
        for (long int i(j); i >= 0; --i){
            //            cout << " i " << i << endl;
            //cout << " j " << j << endl;
            if (j - i <= 0) {
                D[i][j] = 0;
                continue;
            }
            D[i][j] = max(max(
                              D[i+1][j],
                              D[i][j-1]),
                          scoreB(x[i],x[j]) + D[i+1][j-1]);
            //            if (j == n-4 && i == n-4){
            //            cout << " score " <<D[i][j] << endl;
            //            }
            size_t groupI ((i)/q);
            size_t groupJ ((j)/q);
            // the matrix is diagonal so groupI doesn't start always from the same position but from
            // the diagonal
            int nGroupsBetween (int(groupJ - groupI - 1));
            // + q to get right most. -1 to put back in 0 index
            size_t iI = min(q*groupI + q - 1, n-1);
            size_t jJ(q*groupJ);
            //cout << " Ii " << iI << endl;
            
            //cout << "II" << endl;
            // for all cells in the first group
            for (size_t k(i+1); k < iI; ++k){
                
                //                cout << " i " << i << " j " << j << " k " << k << endl;
                //                cout << D[i][k-1] + D[k][j] << "a" <<endl;
                D[i][j] = max (D[i][j], D[i][k-1] + D[k][j]);
            }
            
            //cout << "JJ" << endl;
            //for all cells in last group
            for (size_t k(jJ+1); k < j; ++k){
                //                cout << " i " << i << " j " << j << " k " << k << endl;
                //                if (i == 0 && j == 10) {
                //                cout << D[i][k-1] + D[k][j] << "c" << endl;
                //                }
                D[i][j] = max (D[i][j], D[i][k-1] + D[k][j]);
            }
            
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
                size_t t(l+1);
//                if (i == 0 && j == 9) {
//                cout << i << " " << (iI + q*K + 1)/q << endl;
//                cout << j << " " << t/q << endl;
//                }
                ull hdiff = hvs[i][(iI + q*K + 1)/q];
                ull vdiff = vvs[j][t/q];
                double maxdiff = R[hdiff][vdiff];
                double baseH = D[i][l];
                double baseV = D[t][j];
                double base = baseH + baseV;
//                if (i == 0 && j == 9) {
//                    cout << "i: " << i << " l: " << l << " t: " << t << " j: " << j << " q: " << q << " K: " << K << " nGroupsBetween: " << nGroupsBetween << endl;
//                    cout << "base H: " << baseH << endl;
//                    cout << "base V: " << baseV << endl;
//                    cout << "hdiff: " << bitset<8>(hdiff) << endl;
//                    cout << "vdiff: " << bitset<8>(vdiff) << endl;
//                    cout << "maxdiff: " << maxdiff << endl;
//                }
                D[i][j] = max(D[i][j], base + maxdiff);
                
//                 D[i][j] = max(D[i][j],D[i][l] + D[t][j] + R[hvs[i][K]][vvs[j][K]]);
            }
            
            // compute the vertical difference vector
            if (i % q == 1 && i + q <= n) {
                // compute and store the v¯ vector i/qth group for column j
                ull vdiff = 0; size_t c = q-2;
                for (size_t k(i); k < i+q-1; ++k) {
                    if (D[k][j] - D[k+1][j] == 1) {
                        vdiff = (vdiff | (1 << c));
                    }
                    c--;
                }
//                cout << "calculating vdiff for " << i << " to " << i+q-1 << ": " << bitset<8>(vdiff) << endl;
//                cout << "group: " << groupI << "for j: " << j << endl;
//                cout << "vvs: " << j << " " << groupI << " " << bitset<8>(vdiff) << endl;
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
//                cout << "calculating hdiff for " << j+1-q << " to " << j << ": " << bitset<8>(hdiff) << endl;
//                cout << "group: " << groupJ << " for i: "<< i  << endl;
//                cout << "hvs: " << i << " " << groupJ << " " << bitset<8>(hdiff) << endl;
                hvs[i][groupJ] = hdiff; // (j − 1)/qth group for row i
            }
//            cout << "i: " <<  i << " j: " << j << endl;
//            debugPrint(D);                   
        }
    }
    
//    debugPrint(D);
//    cout << D[n][n] << endl;
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
