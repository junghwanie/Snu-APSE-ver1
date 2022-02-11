#ifndef console_hpp
#define console_hpp

#include <iostream>
using namespace std;

int menuDraw();

class Console {
    
public:
    void linkInit();
    string CLinkInputPath(); // path 하나만 있으면 됨

    void CValueOutput();
    void CSkewOutput();
    void CRcfvOutput();
    void CFreqOutput();
    void CMissOutput();
};

#endif /* console_hpp */
