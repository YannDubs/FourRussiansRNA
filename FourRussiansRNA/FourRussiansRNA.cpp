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
int b(char a, char b){
    // returns 1 if match 0 else
    if ((a == 'a' && b == 'u') || (a == 'u' && b == 'a') || (a == 'c' && b == 'g') || (b == 'c' && a == 'g')){
        return 1;
    }
    return 0;
}

//Four russians running
void nussimovFourRussioans(const string& x){
    // INITIALIZATION
    size_t n (x.size());
    vector<vector<double>> D (n,vector<double> (n));
    const size_t q (round(log(n)));
    //Index
    //vector<vector<size_t>> Index (m,vector<size_t> (n));
    
    // ITERATION
    for (size_t j(0); j < n; ++j){
        for (size_t i(j-1); i <= 0; --i){
            D[i][j] = b(x[i],x[j]) + D[i+1][j-1];
            size_t groupI ((i+1)/q);
            size_t groupJ ((j+1)/q);
            // the matrix is diagonal so groupI doesn't start always from the same position but from
            // the diagonal
            size_t nGroupsBetween (groupJ - groupI - 1);
            // + q to get right most. -1 to put back in 0 index
            size_t iI(q*groupI + q - 1);
            size_t jJ(q*groupJ - 1);
            
            // for all cells in the first group
            for (size_t k(i+1); i <= iI; ++k){
                D[i][j] = max (D[i][j], D[i][k-1] + D[k][j]);
            }
            
            //for all cells in last group
            for (size_t k(jJ); k <= j; ++k){
                D[i][j] = max (D[i][j], D[i][k-1] + D[k][j]);
            }
            
            //for all groups in between
            for (size_t K(0); K < nGroupsBetween; ++K){
                // take right most element of group I, adds 1 to get left most of next
                // then adds K*q to shift depedning on the block
                size_t l(iI + 1 + K*q);
                //not sure
                size_t t(l);
                D[i][j] = max(D[i][j],D[i][l] + D[t][j] + getR(hvstore(i,K),vvstore(j,K)));
            }
        }
    }
    
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
