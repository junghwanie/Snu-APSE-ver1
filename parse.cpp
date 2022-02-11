#include "parse.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstdio>
#include <vector>

#define ONESEQ 30000

vector<string> Parse::parsing_c(string c_path) {
    ifstream tmp_query;
    tmp_query.open(c_path);
    
    // fasta file format
    if (!tmp_query.is_open()) {
        cout << "Failed to open file: " << c_path << endl;
    }
    
    int seq_size = 0;
    int acc_size = 0;
    int per_seq = 0;
    int carriage = 0;
    
    while (getline(tmp_query, buff)) {
        if (buff[0] != '>') {
            tmp_seq += buff;
            seq_size += buff.size();
        }
        else {
            tmp_acc += buff;
            acc_size += 1;
        }
    }
    
    for (int i = 0; i < seq_size; i++) {
        if (tmp_seq[i] == '\r') {
            carriage += 1;
        }
    }
    
    per_seq = (seq_size + 1) / acc_size;
    
    size_t seq_len = tmp_seq.length();
    int q = 0;
    while (seq_len > q) { // each seq -> vector
        sub_seq = tmp_seq.substr(q, per_seq);
        or_seq.push_back(sub_seq);
        
        q += sub_seq.length();
    }
    
    return or_seq;
}

vector<string> Parse::parsing_a(string a_path) {
    ifstream tmp_acc;
    tmp_acc.open(a_path);
    
    if (!tmp_acc.is_open()) {
        cout << "Failed to open file: " << a_path << endl;
    }
    
    int seq_size = 0;
    string acc_buff;
    
    while (getline(tmp_acc, buff)) {
        if (buff[0] != '>') {
            tmp_seq += buff;
            seq_size += buff.size();
        }
        else {
            buff.erase(0, 1); // erase ">"
            or_acc.push_back(buff); 
        }
    }
    // store accession ID

    return or_acc;
}

// position-based trim for multiple sequences
vector<string> Parse::ParseGAPTrimFunc(string m_TrimPath, int m_TrimPos) {
    vector<string> m_TrimMaxGapChar;
    vector<string> m_TrimMaxGapTaxon;
    vector<string> m_TrimGapMultiple;
    m_TrimMaxGapChar = parsing_c(m_TrimPath);
    m_TrimMaxGapTaxon = parsing_a(m_TrimPath);
    int m_TrimSize = m_TrimMaxGapTaxon.size();
   
    m_TrimMaxGapChar.erase(m_TrimMaxGapChar.begin() + m_TrimPos);
    m_TrimMaxGapTaxon.erase(m_TrimMaxGapTaxon.begin() + m_TrimPos);
    
    for (int i = 0; i < m_TrimSize; i++) {
        m_TrimGapMultiple.push_back(">");
        m_TrimGapMultiple.push_back(m_TrimMaxGapTaxon[i]);
    }
    
    return m_TrimGapMultiple;
}

vector<string> Parse::ParseDefiChar(string m_tmpPath) {
    
    ifstream m_tmpAll;
    m_tmpAll.open(m_tmpPath);
    
    if (!m_tmpAll.is_open()) {
        cout << "Failed to open file: " << m_tmpPath << endl;
    }
    
    string m_tmpDefChar;
    vector<string> m_TrimDefChar;
    string m_buff;
    
    while (getline(m_tmpAll, m_buff)) {
        m_tmpDefChar += m_buff;
    }
    
    size_t prev = m_tmpDefChar.find_first_not_of('>', 0);
    size_t cont = m_tmpDefChar.find_first_of('>', prev);
    while (string::npos != prev || string::npos != cont) {
        m_TrimDefChar.push_back(m_tmpDefChar.substr(prev, cont - prev));
        prev = m_tmpDefChar.find_first_not_of('>', cont);
        cont = m_tmpDefChar.find_first_of('>', prev);
    }
    
    vector<string> m_RecDefChar;
    size_t m_TrimAllSize = m_TrimDefChar.size();
    int d = 0;
    size_t m_AccPos = 0;
    size_t prevPos = 0;
    string m_TrimBuff;
    vector<string> m_RecAcc;
    m_RecAcc = parsing_a(m_tmpPath);
    
    while (d < m_TrimAllSize) {
        m_TrimBuff = m_TrimDefChar[d];
        m_TrimBuff.insert(0, ">");
        m_AccPos = m_RecAcc[d].find_last_of(m_TrimDefChar[d]) + 2; // m_RecAcc delimit removed +1
        m_TrimBuff.insert(m_AccPos, "\n");
        prevPos = m_RecAcc[d].find_first_of(m_TrimDefChar[d]);
        m_TrimBuff.insert(prevPos, "\n");
        
        m_RecDefChar.push_back(m_TrimBuff);
        d++;
        m_AccPos = 0;
    }
    m_RecDefChar[0].erase(0, 1);
    
    return m_RecDefChar;
}

