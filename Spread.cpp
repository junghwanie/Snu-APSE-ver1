#include "Spread.hpp"
#define FREQ_TYPE 9
#define SKEW_TYPE 4


void Spread::OSpreadDestroy() {
    m_OSpath.clear();
}

void Spread::OSpreadCval() {
    
    Console con;
    Parse par;
    Calvalue cal;
    Ret r;
    
    vector<string> m_Pcval;
    m_OSpath = con.CLinkInputPath();
    m_Pcval = par.parsing_c(m_OSpath);
    int m_CvalSize = m_Pcval.size();
    
    vector<double> m_IOuputCval;
    m_IOuputCval = cal.c_value(m_OSpath, m_Pcval);
    
    vector<string> m_Cacc;
    m_Cacc = par.parsing_a(m_OSpath);
    
    std::cout << std::endl;
    std::cout << "***** Organism-to Organism C Factor *****" << std::endl << std::endl;
    std::cout << "File: " << m_OSpath << std::endl;
    std::cout << "Ex) 1~30 taxon -> C Factor between taxon 1 and other taxon" << std::endl << std::endl;
   
    static vector<pair<string, double>> c; // utility
    for (auto i = 0; i < m_CvalSize; i++) {
        c.push_back(make_pair(m_Cacc[i], m_IOuputCval[i]));
    }
    
    for (auto i = 0; i < m_CvalSize; i++) {
        cout << c[i].first << ' ';
        cout.setf(ios::fixed);
        cout.setf(ios::showpoint);
        cout.precision(4);
        cout << c[i].second << endl;
    }
    cout << endl;
    
    // p-distance vector with jukes-cantor distance model
    /*
    vector<double> p_dist;
    for (int i = 0; i < a_len; i++) {
         p_dist.push_back((double)(-3.0f / 4.0f) * (double)log(1.0f - ((4.0f / 3.0f) * ((double)ext_acc_len - (double)match_cnt[i]) / (double)ext_acc_len)));
    }
    */
    
    char outsel;
    char retrn;
    cout << "Output file ([y]/[n])? ";
    cin >> outsel;

    if (outsel == 'y') {
        string ABSOLUTE_FILE_PATH; // class x
        cout << "Enter the storage path >> ";
        cin >> ABSOLUTE_FILE_PATH;
        
        ofstream outfile;
        ABSOLUTE_FILE_PATH += "CValue.txt";
        outfile.open(ABSOLUTE_FILE_PATH, ios::app);
        if (outfile.is_open()) {
            for (int i = 0; i < m_CvalSize; i++) {
                outfile << c[i].first << " ";
                outfile.setf(ios::fixed);
                outfile.setf(ios::showpoint);
                outfile.precision(4);
                outfile << c[i].second << endl;
            }
        }
        outfile.close();
    
        cout << "Return ([y]/[n])? ";
        cin >> retrn;
        if (retrn != 'y')
            exit(0);
        else {
            r.exitNow();
        }
    }
    else {
        cout << "Return ([y]/[n])? ";
        cin >> retrn;
        if (retrn != 'y')
            exit(0);
        else {
            r.exitNow();
        }
    }
    OSpreadDestroy();
}

void Spread::OSpreadSkew() {
    
    Console con;
    Parse par;
    Calvalue cal;
    Ret r;
    
    vector<string> m_Pskew;
    m_OSpath = con.CLinkInputPath();
    m_Pskew = par.parsing_c(m_OSpath);
    
    vector<string> v_acc;
    v_acc = par.parsing_a(m_OSpath);
    
    int m_CSkewSize;
    m_CSkewSize = m_Pskew.size();
    
    vector<vector<double>> m_IOutputSkew(m_CSkewSize, vector<double>(SKEW_TYPE, 0));
    m_IOutputSkew = cal.v_skew(m_OSpath, m_Pskew);
    
    double s_sum[SKEW_TYPE];
    double s_avg[SKEW_TYPE];
    for (int i = 0; i < m_CSkewSize; i++) {
        s_sum[0] += m_IOutputSkew[i][0];
        s_sum[1] += m_IOutputSkew[i][1];
        s_sum[2] += m_IOutputSkew[i][2];
        s_sum[3] += m_IOutputSkew[i][3];
    }
    
    s_avg[0] = static_cast<double>(s_sum[0]) / m_CSkewSize;
    s_avg[1] = static_cast<double>(s_sum[1]) / m_CSkewSize;
    s_avg[2] = static_cast<double>(s_sum[2]) / m_CSkewSize;
    s_avg[3] = static_cast<double>(s_sum[3]) / m_CSkewSize;
    
    std::cout << std::endl;
    std::cout << ">> Skew Value: AT / GC / AG / CT" << std::endl;
    static std::map<string, tuple<double,double,double,double>> skew_map;
    for (auto i = 0; i < m_CSkewSize; i++) {
        skew_map[v_acc[i]] = make_tuple(m_IOutputSkew[i][0],m_IOutputSkew[i][1],m_IOutputSkew[i][2],m_IOutputSkew[i][3]);
    }
    for (auto i = 0; i < m_CSkewSize; i++) {
        std::cout << v_acc[i] << ' ';
        std::cout.setf(ios::fixed);
        std::cout.setf(ios::showpoint);
        std::cout.precision(4);
        std::cout << std::get<0>(skew_map[v_acc[i]]) << ' ' << std::get<1>(skew_map[v_acc[i]]) << ' ' << std::get<2>(skew_map[v_acc[i]]) << ' ' << std::get<3>(skew_map[v_acc[i]]) << std::endl;;
    }
    std::cout << std::endl;
    std::vector<double> skew_av;

    for (auto i = 0; i < SKEW_TYPE; i++) {
        skew_av.push_back(s_avg[i]);
    }
    std::cout << ">> Mean: AT / GC / AG / CT" << std::endl;
    for (int i = 0; i < SKEW_TYPE; i++) {
        std::cout << skew_av[i] << ' ';
    }
    std::cout << std::endl;

    char outsel;
    char retrn;
    std::cout << std::endl;
    cout << "Output file ([y]/[n])? ";
    cin >> outsel;

    if (outsel == 'y') {
        string ABSOLUTE_FILE_PATH; // class x
        cout << "Enter the storage path >> ";
        cin >> ABSOLUTE_FILE_PATH;
        
        ofstream outfile;
        ABSOLUTE_FILE_PATH += "SKEW.txt";
        outfile.open(ABSOLUTE_FILE_PATH, ios::app);
        if (outfile.is_open()) {
            for (auto i = 0; i < m_CSkewSize; i++) {
                outfile << v_acc[i] << ' ';
                outfile.setf(ios::fixed);
                outfile.setf(ios::showpoint);
                outfile.precision(4);
                outfile << std::get<0>(skew_map[v_acc[i]]) << ' ' << std::get<1>(skew_map[v_acc[i]]) << ' ' << std::get<2>(skew_map[v_acc[i]]) << ' ' << std::get<3>(skew_map[v_acc[i]]) << std::endl;
            }
            outfile << std::endl;
        }
        outfile.close();
    
        cout << "Return ([y]/[n])? ";
        cin >> retrn;
        if (retrn != 'y')
            exit(0);
        else {
            r.exitNow();
        }
    }
    else {
        cout << "Return ([y]/[n])? ";
        cin >> retrn;
        if (retrn != 'y')
            exit(0);
        else {
            r.exitNow();
        }
    }
    OSpreadDestroy();
}

void Spread::OSpreadFreq() {
    
    Console con;
    Parse par;
    Calvalue cal;
    Ret r;
    
    vector<string> m_Pfreq;
    m_OSpath = con.CLinkInputPath();
    
    m_Pfreq = par.parsing_c(m_OSpath);
    
    int m_CFracSize;
    m_CFracSize = m_Pfreq.size();
    
    vector<vector<double>> m_IOutputFreq(m_CFracSize, vector<double>(FREQ_TYPE, 0));
    m_IOutputFreq = cal.nuc_freq(m_OSpath, m_Pfreq);
    
    vector<string> m_Facc;
    m_Facc = par.parsing_a(m_OSpath);
    
    std::cout << std::endl;
    std::cout << ">> Nucleotide Frequecies: A / G / C / T / GC / AT / R / Y / GAP" << std::endl;
    static std::map<string, tuple<double,double,double,double,double,double,double,double,double>> freq_map;
    for (auto i = 0; i < m_CFracSize; i++) {
        freq_map[m_Facc[i]] = make_tuple(m_IOutputFreq[i][0],m_IOutputFreq[i][1],m_IOutputFreq[i][2],m_IOutputFreq[i][3],m_IOutputFreq[i][4],m_IOutputFreq[i][5],m_IOutputFreq[i][6],m_IOutputFreq[i][7],m_IOutputFreq[i][8]);
    }
    for (auto i = 0; i < m_CFracSize; i++) {
        std::cout << m_Facc[i] << ' ';
        std::cout.setf(ios::fixed);
        std::cout.setf(ios::showpoint);
        std::cout.precision(4);
        std::cout << std::get<0>(freq_map[m_Facc[i]]) << ' ' << std::get<1>(freq_map[m_Facc[i]]) << ' ' << std::get<2>(freq_map[m_Facc[i]]) << ' ' << std::get<3>(freq_map[m_Facc[i]]) << ' ' << std::get<4>(freq_map[m_Facc[i]]) << ' ' << std::get<5>(freq_map[m_Facc[i]]) << ' ' << std::get<6>(freq_map[m_Facc[i]]) << ' ' << std::get<7>(freq_map[m_Facc[i]]) << ' ' <<
        std::get<8>(freq_map[m_Facc[i]]) << std::endl;
    }
    std::cout << std::endl;
    
    for (int i = 0; i < m_CFracSize; i++) {
        m_Asum += m_IOutputFreq[i][0];
        m_Gsum += m_IOutputFreq[i][1];
        m_Csum += m_IOutputFreq[i][2];
        m_Tsum += m_IOutputFreq[i][3];
        m_GCsum += m_IOutputFreq[i][4];
        m_ATsum += m_IOutputFreq[i][5];
        m_Rsum += m_IOutputFreq[i][6];
        m_Ysum += m_IOutputFreq[i][7];
        m_GAPsum += m_IOutputFreq[i][8];
    }
    vector<double> m_CFracAve;
    m_CFracAve.push_back(m_Asum / m_CFracSize);
    m_CFracAve.push_back(m_Gsum / m_CFracSize);
    m_CFracAve.push_back(m_Csum / m_CFracSize);
    m_CFracAve.push_back(m_Tsum / m_CFracSize);
    m_CFracAve.push_back(m_GCsum / m_CFracSize);
    m_CFracAve.push_back(m_ATsum / m_CFracSize);
    m_CFracAve.push_back(m_Rsum / m_CFracSize);
    m_CFracAve.push_back(m_Ysum / m_CFracSize);
    m_CFracAve.push_back(m_GAPsum / m_CFracSize);
    
    std::cout << ">> Mean: A / G / C / T / GC / AT / R / Y / GAP" << std::endl;
    for (int i = 0; i < FREQ_TYPE; i++) {
        std::cout << m_CFracAve[i] << ' ';
    }
    std::cout << std::endl;
    /* Graph -> summary로 제공, 다른 파라미터도 추가
    cout << "#Mean: R(Purine) Bias " << r_avg << endl;
    cout << "#Mean: Y(Pyrimidine) Bias " << y_avg << endl << endl;
    */
    
    // Max Frequency of the GAP code (cutoff = 0.2)
    vector<double> m_CGapFrequencies;
    for (int i = 0; i < m_CFracSize; i++) {
        m_CGapFrequencies.push_back(m_IOutputFreq[i][8]);
    }
    double m_CGapCutOff = 0.2;
    int m_CGapCutOffPos = 0;
    vector<int> m_CGapCutOffIndex;
    vector<string> m_CGapCutOffTaxa;
    string m_CGapMaxTaxon;
    
    for (int i = 0; i < m_CFracSize; i++) {
        if (m_CGapFrequencies[i] >= m_CGapCutOff) {
            m_CGapCutOffPos++;
            m_CGapCutOffTaxa.push_back(m_Facc[i]);
            m_CGapCutOffIndex.push_back(i); // ith Gap acc
        }
    }
    
    char m_bTrimGap;
    char m_bTrimCutoff;
    char m_bUserCutoff;
    char m_bUserCont = true;
    char m_bOutsel;
    char m_bRetrn;
    vector<string> m_TrimDefCharMultiple;
    int m_CGapCutOffSize = m_CGapCutOffIndex.size(); // # of cutoff
    double m_CUserCutOff;
    
    std::cout << std::endl;
    std::cout << "Trim GAP frequencies? ([y]/[n])? ";
    std::cin >> m_bTrimGap;
    
    if (m_bTrimGap == 'y') {
        std::cout << std::endl;
        std::cout << "********** Flitering System **********" << std::endl << std::endl;
        std::cout << "                For GAP" << std::endl << std::endl;
        
        std::cout << "--------------- SUMMARY ---------------" << std::endl;
        std::cout << "The Cutoff Value of GAP Frequencies: 0.2" << std::endl;
        std::cout << "Taxon: ";
        for (int i = 0; i < m_CGapCutOffPos; i++) {
            std::cout << m_CGapCutOffTaxa[i] << ' ';
        }
        std::cout << "(is/are) trimmed AND" << std::endl;
        std::cout << m_CGapCutOffPos << " of " << m_CFracSize << " taxa (is/are) trimmed, Right ([y]/[n])? ";
        std::cin >> m_bTrimCutoff;
        if (m_bTrimCutoff == 'y') {
            std::cout << std::endl;
            std::cout << "******* Modified Multiple Aligned File *******" << std::endl;
            m_TrimDefCharMultiple = par.ParseDefiChar(m_OSpath);
            for (int i = 0; i < m_CGapCutOffSize; i++) {
                m_TrimDefCharMultiple.erase(m_TrimDefCharMultiple.begin() + m_CGapCutOffIndex[i]);
            }
            for (int i = 0; i < m_CFracSize - m_CGapCutOffSize; i++) {
                std::cout << m_TrimDefCharMultiple[i];
            }
            std::cout << std::endl;
            // multiple alignment output function
            std::cout << std::endl;
            std::cout << "Output file ([y]/[n])? ";
            cin >> m_bOutsel;
            if (m_bOutsel == 'y') {
                string ABSOLUTE_FILE_PATH; // class x
                cout << "Enter the storage path >> ";
                cin >> ABSOLUTE_FILE_PATH;
                
                ofstream outfile;
                ABSOLUTE_FILE_PATH += "Multiple.fas";
                outfile.open(ABSOLUTE_FILE_PATH, ios::app);
                if (outfile.is_open()) {
                    for (auto i = 0; i < m_CFracSize - m_CGapCutOffSize; i++) {
                        outfile << m_TrimDefCharMultiple[i];
                    }
                }
                outfile.close();
            
                cout << "Return ([y]/[n])? ";
                cin >> m_bRetrn;
                if (m_bRetrn != 'y')
                    exit(0);
                else {
                    r.exitNow();
                }
            }
            else {
                cout << "Return ([y]/[n])? ";
                cin >> m_bRetrn;
                if (m_bRetrn != 'y')
                    exit(0);
                else {
                    r.exitNow();
                }
            }
        }
        else {
            while (m_bUserCont) {
                std::cout << std::endl;
                std::cout << ">> Enter the Cutoff Value: ";
                std::cin >> m_CUserCutOff;
            
                int m_CGapCutOffUserPos = 0;
                vector<int> m_CGapCutOffUserIndex;
                vector<string> m_CGapCutOffUserTaxa;
                vector<string> m_TrimDefCharRevised;
                string m_CGapMaxTaxon;
            
                for (int i = 0; i < m_CFracSize; i++) {
                    if (m_CGapFrequencies[i] >= m_CUserCutOff) {
                        m_CGapCutOffUserPos++;
                        m_CGapCutOffUserTaxa.push_back(m_Facc[i]);
                        m_CGapCutOffUserIndex.push_back(i); // ith Gap acc
                    }
                }
                std::cout << std::endl;
                std::cout << "--------------- SUMMARY ---------------" << std::endl;
                std::cout << "The Cutoff Value of GAP Frequencies: " << m_CUserCutOff << std::endl;
                std::cout << "Taxon: ";
                for (int i = 0; i < m_CGapCutOffUserPos; i++) {
                    std::cout << m_CGapCutOffUserTaxa[i] << ' ';
                }
                std::cout << "(is/are) trimmed AND" << std::endl;
                std::cout << m_CGapCutOffUserPos << " of " << m_CFracSize << " taxa (is/are) trimmed, Right ([y]/[n])? ";
                std::cin >> m_bUserCutoff;
                if (m_bUserCutoff == 'y') {
                    m_bUserCont = false;
                    std::cout << "******* Modified Multiple Aligned File *******" << std::endl;
                    m_TrimDefCharRevised = par.ParseDefiChar(m_OSpath);
                    int m_CGapCutOffUserSize = m_CGapCutOffUserIndex.size();
                    for (int i = 0; i < m_CGapCutOffUserSize; i++) {
                        m_TrimDefCharRevised.erase(m_TrimDefCharRevised.begin() + m_CGapCutOffUserIndex[i] - i);
                    }
                    for (int i = 0; i < m_CFracSize - m_CGapCutOffUserSize; i++) {
                        std::cout << m_TrimDefCharRevised[i];
                    }
                    std::cout << std::endl;
                    // multiple alignment output function
                    std::cout << std::endl;
                    std::cout << "Output file ([y]/[n])? ";
                    cin >> m_bOutsel;
                    if (m_bOutsel == 'y') {
                        string ABSOLUTE_FILE_PATH; // class x
                        cout << "Enter the storage path >> ";
                        cin >> ABSOLUTE_FILE_PATH;
                
                        ofstream outfile;
                        ABSOLUTE_FILE_PATH += "Revised_Multiple.fas";
                        outfile.open(ABSOLUTE_FILE_PATH, ios::app);
                        if (outfile.is_open()) {
                            for (auto i = 0; i < m_CFracSize - m_CGapCutOffUserSize; i++) {
                                outfile << m_TrimDefCharRevised[i];
                            }
                        }
                        outfile.close();
            
                        cout << "Return ([y]/[n])? ";
                        cin >> m_bRetrn;
                        if (m_bRetrn != 'y')
                            exit(0);
                        else {
                            r.exitNow();
                        }
                    }
                    else {
                        cout << "Return ([y]/[n])? ";
                        cin >> m_bRetrn;
                        if (m_bRetrn != 'y')
                            exit(0);
                        else {
                            r.exitNow();
                        }
                    }
                }
            }
        } // m_bUserCont loop
        
        cout << "Output file ([y]/[n])? ";
        cin >> m_bOutsel;

        if (m_bOutsel == 'y') {
            string ABSOLUTE_FILE_PATH; // class x
            cout << "Enter the storage path >> ";
            cin >> ABSOLUTE_FILE_PATH;
            
            ofstream outfile;
            ABSOLUTE_FILE_PATH += "Base_Frequencies.txt";
            outfile.open(ABSOLUTE_FILE_PATH, ios::app);
            if (outfile.is_open()) {
                for (auto i = 0; i < m_CFracSize-1; i++) {
                    
                }
            }
            outfile.close();
        
            cout << "Return ([y]/[n])? ";
            cin >> m_bRetrn;
            if (m_bRetrn != 'y')
                exit(0);
            else {
                r.exitNow();
            }
        }
        else {
            cout << "Return ([y]/[n])? ";
            cin >> m_bRetrn;
            if (m_bRetrn != 'y')
                exit(0);
            else {
                r.exitNow();
            }
        }
    }
    else {
        cout << "Output file ([y]/[n])? ";
        cin >> m_bOutsel;

        if (m_bOutsel == 'y') {
            string ABSOLUTE_FILE_PATH; // class x
            cout << "Enter the storage path >> ";
            cin >> ABSOLUTE_FILE_PATH;
        
            ofstream outfile;
            ABSOLUTE_FILE_PATH += "Base_Frequencies.txt";
            outfile.open(ABSOLUTE_FILE_PATH, ios::app);
            if (outfile.is_open()) {
                for (auto i = 0; i < m_CFracSize; i++) {
                    outfile << m_Facc[i] << ' ';
                    outfile.setf(ios::fixed);
                    outfile.setf(ios::showpoint);
                    outfile.precision(4);
                    outfile << std::get<0>(freq_map[m_Facc[i]]) << ' ' << std::get<1>(freq_map[m_Facc[i]]) << ' ' << std::get<2>(freq_map[m_Facc[i]]) << ' ' << std::get<3>(freq_map[m_Facc[i]]) << ' ' << std::get<4>(freq_map[m_Facc[i]]) << ' ' << std::get<5>(freq_map[m_Facc[i]]) << ' ' << std::get<6>(freq_map[m_Facc[i]]) << ' ' << std::get<7>(freq_map[m_Facc[i]]) << ' ' <<
                    std::get<8>(freq_map[m_Facc[i]]) << std::endl;
                }
            }
            outfile.close();
    
            cout << "Return ([y]/[n])? ";
            cin >> m_bRetrn;
            if (m_bRetrn != 'y')
                exit(0);
            else {
                r.exitNow();
            }
        }
        else {
            cout << "Return ([y]/[n])? ";
            cin >> m_bRetrn;
            if (m_bRetrn != 'y')
                exit(0);
            else {
                r.exitNow();
            }
        }
    }
    OSpreadDestroy();
}

void Spread::OSpreadRcfv() {
    
    Console con;
    Parse par;
    Calvalue cal;
    Ret r;

    vector<string> m_Prcfv;
    m_OSpath = con.CLinkInputPath();
    
    m_Prcfv = par.parsing_c(m_OSpath);
    
    int m_CRcfvSize;
    m_CRcfvSize = m_Prcfv.size();

    vector<double> m_IOutputRcfv;
    m_IOutputRcfv = cal.rcfv_value(m_OSpath, m_Prcfv);
    
    vector<string> m_Racc;
    m_Racc = par.parsing_a(m_OSpath);
    
    double m_Rsum = 0;
    double m_Ravr;
    for (int i = 0; i < m_CRcfvSize; i++)
        m_Rsum += m_IOutputRcfv[i];
    m_Ravr = static_cast<double>(m_Rsum) / m_CRcfvSize;

    static vector<pair<string, double>> rcfv; // utility
    std::cout << std::endl;
    std::cout << "***** Taxon Specific RCFV Value *****" << std::endl << std::endl;
    std::cout << "File: " << m_OSpath << std::endl << std::endl;
    for (auto i = 0; i < m_CRcfvSize; i++) {
        rcfv.push_back(make_pair(m_Racc[i], m_IOutputRcfv[i]));
    }
    
    for (auto i = 0; i < m_CRcfvSize; i++) {
        std::cout << rcfv[i].first << ' ';
        cout.setf(ios::fixed);
        cout.setf(ios::showpoint);
        cout.precision(4);
        std::cout << rcfv[i].second << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Mean: RCFV" << ' ' << m_Ravr << std::endl << std::endl;
    
    //return r_standard;
    char outsel;
    char retrn;
    std::cout << "Output file ([y]/[n])? ";
    std::cin >> outsel;

    if (outsel == 'y') {
        string ABSOLUTE_FILE_PATH; // class x
        cout << "Enter the storage path >> ";
        cin >> ABSOLUTE_FILE_PATH;
        
        ofstream outfile;
        ABSOLUTE_FILE_PATH += "RCFV.txt";
        outfile.open(ABSOLUTE_FILE_PATH, ios::app);
        if (outfile.is_open()) {
            for (int i = 0; i < m_CRcfvSize; i++) {
                outfile << rcfv[i].first << ' ';
                outfile.setf(ios::fixed);
                outfile.setf(ios::showpoint);
                outfile.precision(4);
                outfile << rcfv[i].second << endl;
            }
        }
        outfile.close();
    
        cout << "Return ([y]/[n])? ";
        cin >> retrn;
        if (retrn != 'y')
            exit(0);
        else {
            r.exitNow();
        }
    }
    else {
        cout << "Return ([y]/[n])? ";
        cin >> retrn;
        if (retrn != 'y')
            exit(0);
        else {
            r.exitNow();
        }
    }
    OSpreadDestroy();
}

void Spread::OSpreadMiss() {
    Console con;
    Parse par;
    Calvalue cal;
    Ret r;

    vector<string> m_Pmiss;
    m_OSpath = con.CLinkInputPath();
    m_Pmiss = par.parsing_c(m_OSpath);
    
    int m_CMissSize;
    m_CMissSize = m_Pmiss.size();
    vector<vector<double>> m_IOutputMiss(m_CMissSize, vector<double>(m_CMissSize, 0));
    m_IOutputMiss = cal.shared_mdata(m_OSpath, m_Pmiss);

    double m_sum[m_CMissSize];
    double m_avg[m_CMissSize];
    double m_ssum = 0;
    long double m_aavg;
    
    for (int i = 0; i < m_CMissSize; i++)
        for (int j = 0; j < m_CMissSize; j++) {
            m_sum[i] += m_IOutputMiss[i][j];
        }
    
    for (int i = 0; i < m_CMissSize; i++)
        m_avg[i] = static_cast<double>(m_sum[i]) / m_CMissSize;
    
    for (int i = 0; i < m_CMissSize; i++)
        m_ssum += m_avg[i];
    
    m_aavg = static_cast<double>(m_ssum) / m_CMissSize;
    
    vector<string> m_acc;
    m_acc = par.parsing_a(m_OSpath);
    std::cout << std::endl;
    std::cout << "***** Shared Missing Data *****" << std::endl << std::endl;
    std::cout << "File: " << m_OSpath << std::endl;
    std::cout << "Ex) Organism-to-Organism Negative Overlap Missing Data" << std::endl;
    std::cout << "Ex) Column and Row of 2D Vector for Taxa" << std::endl << std::endl;
    
    // negative overlapping calculation of ambiguity characters
    for (int i = 0; i < m_CMissSize; i++) {
        for (int j = 0; j < m_CMissSize; j++) {
            cout.setf(ios::fixed);
            cout.setf(ios::showpoint);
            cout.precision(4);
            cout << m_IOutputMiss[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
    static vector<pair<string, double>> md;
    for (auto i = 0; i < m_CMissSize; i++) {
        md.push_back(make_pair(m_acc[i], m_avg[i]));
    }
    cout << ">> Mean: Taxon-Specific Missing Data" << endl;
    for (auto i = 0; i < m_CMissSize; i++) {
        cout << md[i].first << ' ';
        cout.setf(ios::fixed);
        cout.setf(ios::showpoint);
        cout.precision(4);
        cout << md[i].second << endl;
    }
    cout << endl;
    cout << ">> Mean (" << m_OSpath << "): " << m_aavg << endl << endl;
    
    char outsel;
    char retrn;
    cout << "Output file ([y]/[n])? ";
    cin >> outsel;

    if (outsel == 'y') {
        string ABSOLUTE_FILE_PATH; // class x
        cout << "Enter the storage path >> ";
        cin >> ABSOLUTE_FILE_PATH;
        
        ofstream outfile;
        ABSOLUTE_FILE_PATH += "Shared_Missing.txt";
        outfile.open(ABSOLUTE_FILE_PATH, ios::app);
        if (outfile.is_open()) {
            for (int i = 0; i < m_CMissSize; i++) {
                for (int j = 0; j < m_CMissSize; j++) {
                    outfile.setf(ios::fixed);
                    outfile.setf(ios::showpoint);
                    outfile.precision(4);
                    outfile << m_IOutputMiss[i][j] << ' ';
                }
                outfile << endl;
            }
        }
        outfile.close();
    
        cout << "Return ([y]/[n])? ";
        cin >> retrn;
        if (retrn != 'y')
            exit(0);
        else {
            r.exitNow();
        }
    }
    else {
        cout << "Return ([y]/[n])? ";
        cin >> retrn;
        if (retrn != 'y')
            exit(0);
        else {
            r.exitNow();
        }
    }
    OSpreadDestroy();
}

