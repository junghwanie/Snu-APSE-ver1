//
//  parse.hpp
//  master0.0
//
//  Created by 이정환 on 2020/03/20.
//  Copyright © 2020 이정환. All rights reserved.
//

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
