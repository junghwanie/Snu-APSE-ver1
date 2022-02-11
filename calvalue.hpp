#ifndef cvalue_hpp
#define cvalue_hpp

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "console.hpp"
#include "Spread.hpp"

using namespace std;

vector<vector<int>> CalNumberBase(vector<string> m_RefSeq, int m_RefSize, int m_OneSize);

class Calvalue {

    string e_file;
    string u_skew;
    string path;
    string c_file;
    string v_path;
    
    int a_len;
    int a_sum = 0, g_sum = 0, c_sum = 0, t_sum = 0;
    double r_sum, y_sum;
    double a_avg = 0, g_avg = 0, c_avg = 0, t_avg = 0;
    double r_avg, y_avg;
    
    vector<double> a_abs, g_abs, c_abs, t_abs;
    vector<string> pre_cval;
    vector<string> pre_skew;
    
    vector<string> pre_freq;
    vector<string> m_CFreqDescription;
    
    vector<string> pre_rcfv;
    vector<string> pre_mdata;
    vector<double> AT_skew, GC_skew, AG_skew, CT_skew;
    vector<double> r_skew;
    vector<double> meanA, meanG, meanC, meanT;
    
    vector<double> dist;
    vector<string> c_acc;
   
protected:
    
public:
    vector<double> c_value(string c_path, vector<string> m_CvalChar);
    vector<vector<double>> v_skew(string s_path, vector<string> m_SkewChar);
    vector<vector<double>> nuc_freq(string f_path, vector<string> m_FracChar);
    vector<double> rcfv_value(string r_path, vector<string> m_RcfvChar);
    vector<vector<double>> shared_mdata(string m_path, vector<string> m_MissChar);
    
    void GapTrimFunc(vector<vector<double>> m_CFreqFile, vector<string> m_CFreqAcc, int m_CFracASize, int m_GapMaxVal);
};

#endif /* cvalue_hpp */
