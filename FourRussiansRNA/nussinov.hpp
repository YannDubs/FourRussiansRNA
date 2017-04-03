//
//  nussinov.hpp
//  FourRussiansRNA
//
//  Created by Leo Cho on 2017-04-02.
//  Copyright © 2017 Yann Dubois. All rights reserved.
//

#include <string>
#include <vector>

using namespace std;

#ifndef nussinov_hpp
#define nussinov_hpp

typedef vector<int> vi;
typedef vector<vi> vvi;

int nussinovScore(const string& x);
int foldScore(int i, int j, const string &seq, vvi &matrix);
int B(char a, char b);

#endif /* nussinov_hpp */
