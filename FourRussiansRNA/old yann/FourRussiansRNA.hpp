//
//  FourRussiansRNA.hpp
//  FourRussiansRNA
//
//  Created by Yann Dubois on 01.04.17.
//  Copyright © 2017 Yann Dubois. All rights reserved.
//

#include <string>
#include <vector>

using namespace std;

typedef unsigned long long ull;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<ull> vu;
typedef vector<vu> vvu;

#ifndef FourRussiansRNA_h
#define FourRussiansRNA_h

// scoring function
int scoreB(char a, char b);

void nussinovFourRussians(const string& x);
string LoadSeq(string file);
int maxVal(ull x, ull y, const size_t q);
#endif /* FourRussiansRNA_h */
