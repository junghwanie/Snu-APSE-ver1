#include "console.hpp"
#include "calvalue.hpp"
#include "parse.hpp"
#include <iostream>
#define FREQ_TYPE 8
using namespace std;

int menuDraw() {
    std::cout << std::endl;
    std::cout << "************** Connection System  **************" << std::endl << std::endl;
    std::cout << "                    For APSE                    " << std::endl << std::endl;
    std::cout << ">> Choose The Menu To Calculate Systematic Error" << std::endl << std::endl;
    std::cout << "1. C Value" << std::endl;
    std::cout << "2. Skew Value" << std::endl;
    std::cout << "3. Base Frequencies" << std::endl;
    std::cout << "4. RCFV value" << std::endl;
    std::cout << "5. Shared Missing Data" << std::endl;
    std::cout << "9. Exit" << std::endl;
    std::cout << "Enter the Number and Press Return >>";
    
    int s_value;
    cin >> s_value;
    
    return s_value;
}

void Console::linkInit() {
    Calvalue cv;
    int user_s;
    
    user_s = menuDraw();
    
    if (user_s == 1) {
        CValueOutput();
    }
    else if (user_s == 2){
        CSkewOutput();
    }
    else if (user_s == 3) {
        CFreqOutput();
    }
    else if (user_s == 4) {
        CRcfvOutput();
    }
    else if (user_s == 5) {
        CMissOutput();
    }
    else if (user_s == 9) {
        exit(0);
    }
    else {
        std::cerr << "ERROR::PROGRAM::LINKING::FAILED\n";
    }
}

string Console::CLinkInputPath() {
    
    Calvalue cv;
    Parse pr;
    
    string f_path, m_linkfile;
    m_linkfile = "/0";
       
    cout << "Enter your multiple alignment file(.fasta, .aln, .txt) >>";
    cin >> m_linkfile;
    
    return m_linkfile;
}

void Console::CValueOutput() {
    Spread sp;
    sp.OSpreadCval();
}

void Console::CSkewOutput() {
    Spread sp;
    sp.OSpreadSkew();
}

void Console::CFreqOutput() {
    Spread sp;
    sp.OSpreadFreq();
}

void Console::CRcfvOutput() {
    Spread sp;
    sp.OSpreadRcfv();
}

void Console::CMissOutput() {
    Spread sp;
    sp.OSpreadMiss();
}
