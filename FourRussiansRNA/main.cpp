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

void compareNussinov(string seq) {
    cout << "Seq: " << seq << endl;
    
    clock_t clock = timerStart();
    double frscore = nussinovFourRussians(seq);
    timerEnd(clock, "Four Russians Time: ");
    
    clock = timerStart();
    double nscore = double(nussinovScore(seq));
    timerEnd(clock, "Nussinov Time: ");
    
    if (frscore != nscore) {
        cout << "******* Score Mismatch!" << endl;
    }
}


int main(int argc, const char * argv[]) {
    // string x(LoadSeq("sequence.txt"));
    string seq = "gcgaggcuagcgcuacccgugcgccugcguggaacgauucuguggcgagugccggccgaaagcuagguccggauugcacguggagggccgcccgaagggcacucucggacauuaacccgcauucuguaccauggggcgcaaguuggacccuacgaaggagaagcgggggccaggccgaaaggcccggaagcagaagggugccgagacagaacucgucagauucuugccugcaguaagugacgaaaauuccaagaggcugucuagucgugcucgaaagagggcagccaagaggagauugggcucuguugaagccccuaagacaaauaagucuccugaggccaaaccauugccuggaaagcuaccaaaaggagcuguccagacagcugguaagaagggaccccagucccuauuuaaugcuccucgaggcaagaagcgcccagcaccuggcagugaugaggaagaggaggaggaagacucugaagaagaugguauggugaaccacggggaccucuggggcuccgaggacgaugcugauacgguagaugacuauggagcugacuccaacucugaggaugaggaggaaggugaagcguugcugcccauugaaagagcugcucggaagcagaaggcccgggaagcugcugcugggauccaguggagugaagaggagaccgaggacgaggaggaagagaaagaagugaccccugagucaggccccccaaagguggaagaggcagaugggggccugcagaucaauguggaugaggaaccauuugugcugcccccugcuggggagauggagcaggaugcccaggcuccagaccugcaacgaguucacaagcggauccaggauauugugggaauucugcgugauuuuggggcucagcgggaggaagggcggucucguucugaauaccugaaccggcucaagaaggaucuggccauuuacuacuccuauggagacuuccugcuuggcaagcucauggaccucuuc";
    string bcyrn1 = "ggccgggcgcgguggcucacgccuguaaucccagcucucagggaggcuaagaggcgggaggauagcuugagcccaggaguucgagaccugccugggcaauauagcgagaccccguucuccagaaaaaggaaaaaaaaaaacaaaagacaaaaaaaaaauaagcguaacuucccucaaagcaacaaccccccccccccuuu";
    string seq1 = "gcgagg";
    string seq2 = "cagac";
    string seq3 = "caagaacaag";
    string seq4 = "cagcag";
    string seq5 = "gggcauuaacccg";
    string seq6 = "cgccugcguggaacgauucuguccggccgaaagcuaaga";
    vector<string> groups{seq, bcyrn1, seq1, seq2, seq3, seq4, seq5, seq6}; 
    for (auto s : groups) {
        compareNussinov(s);
        cout << endl;
    }
    
    return 0;
}

