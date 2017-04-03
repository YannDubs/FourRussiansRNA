//
//  nussinov.cpp
//  FourRussiansRNA
//
//  Created by Leo Cho on 2017-04-02.
//  Copyright © 2017 Yann Dubois. All rights reserved.
//

#include "nussinov.hpp"

#include <iostream>


// scoring function
int B(char a, char b){
    // returns 1 if match 0 else
    if ((a == 'a' && b == 'u') || (a == 'u' && b == 'a') || (a == 'c' && b == 'g') || (b == 'c' && a == 'g')){
        return 1;
    }
    return 0;
}

void debugPrint(vvi &m) {
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[0].size(); j++) {
            cout << m[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}


int foldScore(int i, int j, const string &seq, vvi &matrix) {
    if (j-i < 2) {
        return 0;
    }
    if (matrix[i][j] == -1) {
        matrix[i][j] = max(max(
                           foldScore(i+1, j, seq, matrix),
                           foldScore(i, j-1, seq, matrix)),
                           foldScore(i+1, j-1, seq, matrix)+B(seq[i], seq[j]));
        for (int k = i + 1; k < j; k++) {
            matrix[i][j] = max(
                               matrix[i][j],
                               foldScore(i, k, seq, matrix)+foldScore(k+1, j, seq, matrix));
        }
    }
    return matrix[i][j];
}

int foldScoreIterative(const string &seq, vvi &matrix) {
    int n = int(seq.size());
    for (int j = 0; j < n; j++){
         for (int i = j; i >= 0; i--) {
             if (j - i <= 0) {
                 matrix[i][j] = 0;
                 continue;
             }
             matrix[i][j] = max(max(
                                    matrix[i+1][j],
                                    matrix[i][j-1]),
                                matrix[i+1][j-1]+B(seq[i], seq[j]));
             for (int k = i + 1; k < j; k++) {
                 matrix[i][j] = max(
                                    matrix[i][j],
                                    matrix[i][k]+matrix[k+1][j]);
             }
         }
    }
    return matrix[0][n-1];
}

int nussinovScore(const string& x) {
    size_t n (x.size());
    vvi matrix(n, vi(n, -1));
    int score = foldScoreIterative(x, matrix);
    cout << score << endl;
    return score;
}
