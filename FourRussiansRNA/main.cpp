//
//  main.cpp
//  FourRussiansRNA
//
//  Created by Yann Dubois on 01.04.17.
//  Copyright © 2017 Yann Dubois. All rights reserved.
//

#include <iostream>
#include <string>
#include <math.h>

#include "FourRussiansRNA.hpp"
#include "nussinov.hpp"

using namespace std;

clock_t timerStart() {
    clock_t start;
    start = clock();
    return start;
}

double timerEnd(clock_t start) {
    return (clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
}

void compareNussinov(string seq) {
    
    clock_t clock = timerStart();
    double frscore = nussinovFourRussians(seq);
    cout << "Four Russians Time: " << timerEnd(clock) << " ms" << endl;
    
    clock = timerStart();
    double nscore = double(nussinovScore(seq));
    cout << "Nussinov: " << timerEnd(clock) << " ms" << endl;
    
    if (frscore != nscore) {
        cout << "******* Score Mismatch!" << endl;
    }
}

void crossValidateBase(string seq) {
    double currentTime(1e6);
    double bestTime(1e6);
    double bestBase(0);
    vector<double> bases({2,3,5,7,10});
    size_t n (seq.size());
    
    for (auto baseLog : bases){
        int q (round(log(n)/log(baseLog)));
        cout << "baseLog " << baseLog << " Q " << q << endl;
        
        clock_t clock = timerStart();
        double frscore = nussinovFourRussians(seq, q);
        //best seems to be base 2
        currentTime = timerEnd(clock);
        cout << "Four Russians Time: " << currentTime  << " ms" << endl << endl;
        
        if (currentTime < bestTime){
            bestTime = currentTime;
            bestBase = baseLog ;
        }
        
    }
    
    cout << " bestBase " << bestBase << endl << endl;
}

void crossValidateQ(string seq) {
    double currentTime(1e6);
    double bestTime(1e6);
    double bestQ(0);
    vector<double> qs({3,4,5,6,7,8,9});
    
    for (auto q : qs){
        
        cout << "Q: " << q << endl;
        
        clock_t clock = timerStart();
        double frscore = nussinovFourRussians(seq, q);
        //best seems to be base 2
        currentTime = timerEnd(clock);
        cout << "Four Russians Time: " << currentTime  << " ms" << endl << endl;
        
        if (currentTime < bestTime){
            bestTime = currentTime;
            bestQ = q ;
        }
        
    }
    
    cout << " bestQ " << bestQ << endl << endl;
}


int main(int argc, const char * argv[]) {
    
    //const vector<string> fileSequences({"sequence2kbp.txt","sequence3kbp.txt","sequence4kbp.txt","sequence5kbp.txt","sequence6kbp.txt"});
    
    const vector<string> fileSequences({"sequence200bp.txt"});//,"sequence1kbp.txt","sequence3kbp.txt","sequence6kbp.txt"});
    //const vector<string> fileSequences({"sequence3kbp.txt"});
    //const vector<string> fileSequences{"sequence200bp.txt","sequence1kbp.txt"};
    
    for (auto file : fileSequences) {
        string seq(LoadSeq(file));
        cout << "Sequence from file : " << file << endl;
        crossValidateBase(seq);
        //crossValidateQ(seq);
        //compareNussinov(seq);
        cout << endl;
    }
    
    return 0;
}

/*string seq = "gcgaggcuagcgcuacccgugcgccugcguggaacgauucuguggcgagugccggccgaaagcuagguccggauugcacguggagggccgcccgaagggcacucucggacauuaacccgcauucuguaccauggggcgcaaguuggacccuacgaaggagaagcgggggccaggccgaaaggcccggaagcagaagggugccgagacagaacucgucagauucuugccugcaguaagugacgaaaauuccaagaggcugucuagucgugcucgaaagagggcagccaagaggagauugggcucuguugaagccccuaagacaaauaagucuccugaggccaaaccauugccuggaaagcuaccaaaaggagcuguccagacagcugguaagaagggaccccagucccuauuuaaugcuccucgaggcaagaagcgcccagcaccuggcagugaugaggaagaggaggaggaagacucugaagaagaugguauggugaaccacggggaccucuggggcuccgaggacgaugcugauacgguagaugacuauggagcugacuccaacucugaggaugaggaggaaggugaagcguugcugcccauugaaagagcugcucggaagcagaaggcccgggaagcugcugcugggauccaguggagugaagaggagaccgaggacgaggaggaagagaaagaagugaccccugagucaggccccccaaagguggaagaggcagaugggggccugcagaucaauguggaugaggaaccauuugugcugcccccugcuggggagauggagcaggaugcccaggcuccagaccugcaacgaguucacaagcggauccaggauauugugggaauucugcgugauuuuggggcucagcgggaggaagggcggucucguucugaauaccugaaccggcucaagaaggaucuggccauuuacuacuccuauggagacuuccugcuuggcaagcucauggaccucuuc";
 string bcyrn1 = "ggccgggcgcgguggcucacgccuguaaucccagcucucagggaggcuaagaggcgggaggauagcuugagcccaggaguucgagaccugccugggcaauauagcgagaccccguucuccagaaaaaggaaaaaaaaaaacaaaagacaaaaaaaaaauaagcguaacuucccucaaagcaacaaccccccccccccuuu";
 string seq1 = "gcgagg";
 string seq2 = "cagac";
 string seq3 = "caagaacaag";
 string seq4 = "cagcag";
 string seq5 = "gggcauuaacccg";
 string seq6 = "cgccugcguggaacgauucuguccggccgaaagcuaaga";*/
