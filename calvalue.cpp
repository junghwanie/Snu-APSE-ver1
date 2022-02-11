#include "calvalue.hpp"
#include "console.hpp"
#include "parse.hpp"
#include "return.hpp"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cmath>

#define NUC_TYPE 5
#define FREQ_TYPE 9
#define SKEW_TYPE 4

using namespace std;

vector<vector<int>> CalNumberBase(vector<string> m_RefSeq, int m_RefSize, int m_OneSize) {
    string tmp_freq;
    vector<vector<int>> m_RefBase(m_RefSize, vector<int>(NUC_TYPE, 0));
    
    for (int i = 0; i < m_RefSize; i++) {
        tmp_freq = m_RefSeq[i];
        for (int j = 0; j < m_OneSize; j++) {
            if (tmp_freq[j] == 'A')
                m_RefBase[i][0]++;
            else if (tmp_freq[j] == 'G')
                m_RefBase[i][1]++;
            else if (tmp_freq[j] == 'C')
                m_RefBase[i][2]++;
            else if (tmp_freq[j] == 'T')
                m_RefBase[i][3]++;
            else
                m_RefBase[i][4]++;
        }
    }
    return m_RefBase;
}

vector<double> Calvalue::c_value(string c_path, vector<string> m_CvalChar) {
    
    pre_cval = m_CvalChar;
    
    map<string, string> m_ti;
    m_ti.insert(make_pair("A", "G"));
    m_ti.insert(make_pair("G", "A"));
    m_ti.insert(make_pair("C", "T"));
    m_ti.insert(make_pair("T", "C"));
    
    map<string, string> m_tv1;
    m_tv1.insert(make_pair("A", "C"));
    m_tv1.insert(make_pair("C", "A"));
    m_tv1.insert(make_pair("T", "A"));
    m_tv1.insert(make_pair("G", "C"));
    map<string, string> m_tv2;
    m_tv2.insert(make_pair("A", "T"));
    m_tv2.insert(make_pair("C", "G"));
    m_tv2.insert(make_pair("T", "G"));
    m_tv2.insert(make_pair("G", "T"));
    
    map<string, string> m_match;
    m_match.insert(make_pair("A", "A"));
    m_match.insert(make_pair("G", "G"));
    m_match.insert(make_pair("C", "C"));
    m_match.insert(make_pair("T", "T"));
    
    string ext_acc;
    string ext_comp;
    
    int m = 0;
    int n = 0;
    int ti_cnt; // transition
    int tv_cnt; // transversion
    int mat_cnt; // match (p-distance)
    
    vector<double> tri_cnt;
    vector<double> trv_cnt;
    vector<double> vi_vt;
    
    string n_carr;
    int c_cnt = 0;
    n_carr = pre_cval[0];
    
    size_t c_size = pre_cval.size();
    size_t ext_len = 0;
    
    size_t ext_acc_len = ext_acc.length(); // # of characeter 764
    for (int i = 0; i < ext_acc_len; i++) {
        if (n_carr[i] == '\r' || '\0')
            c_cnt++;
    }
    vector<vector<int>> ti_nuc(c_size, vector<int>(c_size, 0));
    vector<vector<int>> tv_nuc(c_size, vector<int>(c_size, 0));
    vector<vector<double>> p_dist(c_size, vector<double>(c_size, 0));
    vector<vector<double>> new_dist(c_size, vector<double>(c_size, 0));
    vector<vector<double>> titv(c_size, vector<double>(c_size, 0));
    
    while (n < c_size) {
        m = 0;
        while (m < c_size) {

            ti_cnt = 0;
            tv_cnt = 0;
            mat_cnt = 0;
            ext_len = pre_cval[n].length();
       
            for (int i = 0; i < ext_len; i++) {
                if (m_ti[pre_cval[n].substr(i, 1)] == pre_cval[m].substr(i, 1)) {
                    ti_cnt++;
                }
                if (m_tv1[pre_cval[n].substr(i, 1)] == pre_cval[m].substr(i, 1)) {
                    tv_cnt++;
                }
                if (m_tv2[pre_cval[n].substr(i, 1)] == pre_cval[m].substr(i, 1)) {
                    tv_cnt++;
                }
                if (m_match[pre_cval[n].substr(i, 1)] == pre_cval[m].substr(i, 1)) {
                    mat_cnt++;
                }
            }
            ti_nuc[n][m] = ti_cnt;
            tv_nuc[n][m] = tv_cnt;
            p_dist[n][m] = mat_cnt;
            m++;
        }
        n++;
    }
    
    for (int i = 0; i < c_size; i++) {
        for (int j = 0; j < c_size; j++) {
            if (i == j) {
                titv[i][j] = 0;
            }
            else
                titv[i][j] = static_cast<double>(ti_nuc[i][j])
                / tv_nuc[i][j];
        }
    }
    
    // sum of ti_tv vector
    double t_sum[c_size];
    double v_sum[c_size];
    double p_sum[c_size];
    double pv_sum[c_size];
    double t_avg[c_size];
    double v_avg[c_size];
    double p_avg[c_size];
    double pv_avg[c_size];
    
    for (int k = 0; k < c_size; k++) {
        for (int i = 0; i < c_size; i++) {
            t_sum[k] += titv[k][i];
        }
    }
    for (int i = 0; i < c_size; i++) {
        t_avg[i] = t_sum[i] / c_size;
    }
    for (int i = 0; i < c_size; i++) {
        for (int j = 0; j < c_size; j++) {
            titv[i][j] = pow(titv[i][j] - t_avg[i], 2);
        }
    }
    for (int i = 0; i < c_size; i++) {
        for (int j = 0; j < c_size; j++) {
            v_sum[i] += titv[i][j];
        }
    }
    for (int i = 0; i < c_size; i++) {
        v_avg[i] = sqrt(v_sum[i] / c_size);
    }
    
    size_t w_cnt = pre_cval[0].length();
    for (int i = 0; i < c_size; i++) {
        for (int j = 0; j < c_size; j++) {
            p_dist[i][j] = static_cast<double>(abs(p_dist[i][j] - w_cnt)) / w_cnt;
        }
    }
    for (int k = 0; k < c_size; k++) {
        for (int i = 0; i < c_size; i++) {
            p_sum[k] += p_dist[k][i];
        }
    }
    for (int i = 0; i < c_size; i++) {
        p_avg[i] = p_sum[i] / c_size;
    }
    
    for (int i = 0; i < c_size; i++) {
        for (int j = 0; j < c_size; j++) {
            p_dist[i][j] = pow(p_dist[i][j] - p_avg[i], 2);
        }
    }
    for (int i = 0; i < c_size; i++) {
        for (int j = 0; j < c_size; j++) {
            pv_sum[i] += p_dist[i][j];
        }
    }
    for (int i = 0; i < c_size; i++) {
        pv_avg[i] = sqrt(pv_sum[i] / c_size);
    }
    
    for (int i = 0; i < c_size; i++)
        dist.push_back(static_cast<double>(v_avg[i]) / pv_avg[i]);
    
    return dist;
}

vector<vector<double>> Calvalue::v_skew(string s_path, vector<string> m_SkewChar) {
    
    pre_skew = m_SkewChar;
   
    string tmp_cnt;
    size_t s_size = pre_skew.size();
    size_t s_len = pre_skew[0].length();
    
    vector<vector<int>> s_nuc(s_size, vector<int>(NUC_TYPE, 0));
    s_nuc = CalNumberBase(pre_skew, (int)s_size, (int)s_len);
    
    vector<vector<double>> skew_nuc(s_size, vector<double>(SKEW_TYPE, 0));
    for (int i = 0; i < s_size; i++) {
        skew_nuc[i][0] = static_cast<double>(s_nuc[i][0] - s_nuc[i][3])
        / (s_nuc[i][0] + s_nuc[i][3]);
        skew_nuc[i][1] = static_cast<double>(s_nuc[i][1] - s_nuc[i][2])
        / (s_nuc[i][1] + s_nuc[i][2]);
        skew_nuc[i][2] = static_cast<double>(s_nuc[i][0] - s_nuc[i][1])
        / (s_nuc[i][0] + s_nuc[i][1]);
        skew_nuc[i][3] = static_cast<double>(s_nuc[i][2] - s_nuc[i][3])
        / (s_nuc[i][2] + s_nuc[i][3]);
    }
    
    return skew_nuc;
    
}

vector<vector<double>> Calvalue::nuc_freq(string f_path, vector<string> m_FracChar) {

    pre_freq = m_FracChar;
    
    size_t f_size = pre_freq.size();
    size_t f_len = pre_freq[0].length();

    vector<vector<int>> f_nuc(f_size, vector<int>(NUC_TYPE, 0));
    f_nuc = CalNumberBase(pre_freq, (int)f_size, (int)f_len);

    vector<int> m_nuc;
    vector<int> m_aNuc; // all
    for (int i = 0; i < f_size; i++) {
        m_nuc.push_back(f_nuc[i][0] + f_nuc[i][1] + f_nuc[i][2] + f_nuc[i][3]);
        m_aNuc.push_back(f_nuc[i][0] + f_nuc[i][1] + f_nuc[i][2] + f_nuc[i][3] + f_nuc[i][4]);
    }
    
    vector<vector<double>> freq_nuc(f_size, vector<double>(FREQ_TYPE, 0));
    for (int i = 0; i < f_size; i++) {
        freq_nuc[i][0] = fabs(static_cast<double>(f_nuc[i][0]) / m_nuc[i]);
        freq_nuc[i][1] = fabs(static_cast<double>(f_nuc[i][1]) / m_nuc[i]);
        freq_nuc[i][2] = fabs(static_cast<double>(f_nuc[i][2]) / m_nuc[i]);
        freq_nuc[i][3] = fabs(static_cast<double>(f_nuc[i][3]) / m_nuc[i]);
        
        freq_nuc[i][4] = static_cast<double>(f_nuc[i][1] + f_nuc[i][2]) /
        (f_nuc[i][0] + f_nuc[i][1] + f_nuc[i][2] + f_nuc[i][3]);
        
        freq_nuc[i][5] = static_cast<double>(f_nuc[i][0] + f_nuc[i][3]) /
        (f_nuc[i][0] + f_nuc[i][1] + f_nuc[i][2] + f_nuc[i][3]);
        
        freq_nuc[i][6] = static_cast<double>(f_nuc[i][0] + f_nuc[i][1]) /
        (f_nuc[i][0] + f_nuc[i][1] + f_nuc[i][2] + f_nuc[i][3]);
        
        freq_nuc[i][7] = static_cast<double>(f_nuc[i][3] + f_nuc[i][2]) /
        (f_nuc[i][0] + f_nuc[i][1] + f_nuc[i][2] + f_nuc[i][3]);
        
        freq_nuc[i][8] = static_cast<double>(1 - ((double)m_nuc[i] / m_aNuc[i]));
    }
    
    return freq_nuc;
}

void Calvalue::GapTrimFunc(vector<vector<double>> m_CFreqFile, vector<string> m_CFreqAcc, int m_CFracASize, int m_GapMaxVal) {
    vector<vector<double>> m_CFreqElement(m_CFracASize, vector<double>(FREQ_TYPE, 0));
    m_CFreqElement = m_CFreqFile;
    
    unsigned int m_RowDelete = m_GapMaxVal;
    m_CFreqElement.erase(m_CFreqElement.begin() + m_RowDelete);
}

vector<double> Calvalue::rcfv_value(string r_path, vector<string> m_RcfvChar) {
    
    pre_rcfv = m_RcfvChar;
    string tmp_rcfv;
    size_t r_size = pre_rcfv.size();
    size_t r_len = pre_rcfv[0].length();
    
    vector<vector<int>> r_nuc(r_size, vector<int>(NUC_TYPE, 0));
    vector<vector<double>> m_Rnuc(r_size, vector<double>(NUC_TYPE, 0));
    r_nuc = CalNumberBase(pre_rcfv, (int)r_size, (int)r_len);

    for (int i = 0; i < r_size; i++) {
        for (int j = 0; j < NUC_TYPE; j++) {
            m_Rnuc[i][j] =  r_nuc[i][j] / (double)r_len;
        }
    }
    
    for (int i = 0; i < r_size; i++) {
        a_sum += m_Rnuc[i][0];
        g_sum += m_Rnuc[i][1];
        c_sum += m_Rnuc[i][2];
        t_sum += m_Rnuc[i][3];
    }
    
    a_avg = (double)a_sum / r_size;
    g_avg = (double)g_sum / r_size;
    c_avg = (double)c_sum / r_size;
    t_avg = (double)t_sum / r_size;
  
    for (int i = 0; i < r_size; i++) {
        a_abs.push_back(fabs(m_Rnuc[i][0] - (double)a_avg));
        g_abs.push_back(fabs(m_Rnuc[i][1] - (double)g_avg));
        c_abs.push_back(fabs(m_Rnuc[i][2] - (double)c_avg));
        t_abs.push_back(fabs(m_Rnuc[i][2] - (double)t_avg));
    }
    
    vector<double> r_standard;
    vector<double> m_RCFV;
    for (int i = 0; i < r_size; i++) {
        r_standard.push_back((double)(a_abs[i] + g_abs[i] + c_abs[i]
                                      + t_abs[i]) / (double)r_size);
        
    }

    return r_standard;
}

vector<vector<double>> Calvalue::shared_mdata(string m_path, vector<string> m_MissChar) {
    
    pre_mdata = m_MissChar;
    
    map<string, string> m_amb;
    m_amb.insert(make_pair("N", "N"));
    m_amb.insert(make_pair("-", "-"));
    
    int m = 0, n = 0;
    int m_cnt;
    string m_data;
    m_data = pre_mdata[0];
    
    size_t m_size = pre_mdata.size();
   
    vector<vector<double>> shared_mdata(m_size, vector<double>(m_size, 0));
    size_t amb_len;
    while (n < m_size) {
        m = 0;
        while (m < m_size) {
            m_cnt = 0;
            amb_len = pre_mdata[n].length();
            for (int i = 0; i < amb_len; i++) {
                if (m_amb[pre_mdata[n].substr(i, 1)] == pre_mdata[m].substr(i, 1)) {
                    m_cnt++;
                }
            }
            shared_mdata[n][m] = static_cast<double>(m_cnt) / amb_len;
            m++;
        }
        n++;
    }

    return shared_mdata;
}




