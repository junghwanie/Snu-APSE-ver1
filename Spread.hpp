#ifndef Spread_hpp
#define Spread_hpp

#include <stdio.h>
#include <fstream>
#include <string>
#include "parse.hpp"
#include "calvalue.hpp"
#include "return.hpp"

class Spread {
    double m_Asum = 0, m_Gsum = 0, m_Csum = 0, m_Tsum = 0,
    m_GCsum = 0, m_ATsum = 0, m_Rsum = 0, m_Ysum = 0,
    m_GAPsum = 0;
    
    string m_OSpath;
    
public:
    void OSpreadDestroy();
    void OSpreadFreq();
    void OSpreadRcfv();
    void OSpreadMiss();
    void OSpreadCval();
    void OSpreadSkew();
};

#endif /* Spread_hpp */
