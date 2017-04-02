//
//  nussinov.cpp
//  FourRussiansRNA
//
//  Created by Leo Cho on 2017-04-02.
//  Copyright © 2017 Yann Dubois. All rights reserved.
//

#include "nussinov.hpp"

#include <vector>
#include <iostream>

typedef vector<int> vi;
typedef vector<vi> vvi;

// scoring function
int b(char a, char b){
    // returns 1 if match 0 else
    if ((a == 'a' && b == 'u') || (a == 'u' && b == 'a') || (a == 'c' && b == 'g') || (b == 'c' && a == 'g')){
        return 1;
    }
    return 0;
}

int foldScore(int i, int j, const string &seq, vvi &matrix) {
    if (j-i < 2) {
        return 0;
    }
    if (matrix[i][j] == -1) {
        matrix[i][j] = max(max(
                           foldScore(i+1, j, seq, matrix),
                           foldScore(i, j-1, seq, matrix)),
                           foldScore(i+1, j-1, seq, matrix)+b(seq[i], seq[j]));
        for (int k = i + 1; k < j; k++) {
            matrix[i][j] = max(
                               matrix[i][j],
                               foldScore(i, k, seq, matrix)+foldScore(k+1, j, seq, matrix));
        }
    }
    return matrix[i][j];
}

void nussinovScore(const string& x) {
    size_t n (x.size());
    vvi matrix(n, vi(n));
    cout << foldScore(0, int(n)-1, x, matrix) << endl;
}
