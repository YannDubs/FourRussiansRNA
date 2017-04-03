//
//  main.cpp
//  FourRussiansRNA
//
//  Created by Yann Dubois on 01.04.17.
//  Copyright © 2017 Yann Dubois. All rights reserved.
//

#include <iostream>
#include <string>

#include "FourRussiansRNA.hpp"
#include "nussinov.hpp"

using namespace std;

clock_t timerStart() {
    clock_t start;
    start = clock();
    return start;
}

void timerEnd(clock_t start, string context) {
    cout << context << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;
}


int main(int argc, const char * argv[]) {
// string x(LoadSeq("sequence.txt"));
    // 994 bp
//    string seq = "gcgaggcuagcgcuacccgugcgccugcguggaacgauucuguggcgagugccggccgaaagcuagguccggauugcacguggagggccgcccgaagggcacucucggacauuaacccgcauucuguaccauggggcgcaaguuggacccuacgaaggagaagcgggggccaggccgaaaggcccggaagcagaagggugccgagacagaacucgucagauucuugccugcaguaagugacgaaaauuccaagaggcugucuagucgugcucgaaagagggcagccaagaggagauugggcucuguugaagccccuaagacaaauaagucuccugaggccaaaccauugccuggaaagcuaccaaaaggagcuguccagacagcugguaagaagggaccccagucccuauuuaaugcuccucgaggcaagaagcgcccagcaccuggcagugaugaggaagaggaggaggaagacucugaagaagaugguauggugaaccacggggaccucuggggcuccgaggacgaugcugauacgguagaugacuauggagcugacuccaacucugaggaugaggaggaaggugaagcguugcugcccauugaaagagcugcucggaagcagaaggcccgggaagcugcugcugggauccaguggagugaagaggagaccgaggacgaggaggaagagaaagaagugaccccugagucaggccccccaaagguggaagaggcagaugggggccugcagaucaauguggaugaggaaccauuugugcugcccccugcuggggagauggagcaggaugcccaggcuccagaccugcaacgaguucacaagcggauccaggauauugugggaauucugcgugauuuuggggcucagcgggaggaagggcggucucguucugaauaccugaaccggcucaagaaggaucuggccauuuacuacuccuauggagacuuccugcuuggcaagcucauggaccucuuc";
    string seq = "gcgagg";

    clock_t clock = timerStart();
    nussinovFourRussians(seq);
    timerEnd(clock, "Nussinov Four Russians: ");
    
    clock = timerStart();
    nussinovScore(seq);
    timerEnd(clock, "Nussinov Time: ");
    
    return 0;
}

