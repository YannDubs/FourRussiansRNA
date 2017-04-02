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

#ifndef FourRussiansRNA_h
#define FourRussiansRNA_h

// scoring function
int b(char a, char b);

void nussinovFourRussians(const string& x);
string LoadSeq(string file);
#endif /* FourRussiansRNA_h */
