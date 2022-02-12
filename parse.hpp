#ifndef parse_hpp
#define parse_hpp

#include <iostream>
#include <vector>
#include <string>
using namespace std;

class Parse {
    string buff;
    string tmp_seq;
    string tmp_acc;
    string sub_seq;
    string sub_acc;
    
    vector<string> or_seq;
    vector<string> or_acc;
    vector<double> m_COutputFreq;
    vector<string> m_AOutputFreq;
    
public:
    vector<string> parsing_c(string c_path);
    vector<string> parsing_a(string a_path);
    vector<string> ParseGAPTrimFunc(string m_TrimPath, int m_TrimPos);
    vector<string> ParseDefiChar(string m_tmpPath);
};

#endif /* parse_hpp */
