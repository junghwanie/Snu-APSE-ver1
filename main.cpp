/*
    Version 1.0 2020/01 Lab of Computational Biology and Bioinformatics
    Seoul National University All right reserved.

    Created by Junghwan Lee on 2020/03/15.
    Copyright © 2020 Lab of Computational Biology and Bioinformatics. All rights reserved.
*/

#include "console.hpp"
#include "parse.hpp"
#include <iostream>
#include <stdlib.h>
using namespace std;

extern void exitNow() {
    Console con;
    con.linkInit();
}

int main() {
    
    std::cout << "************** Assessment Program for Systematic Error **************" << std::endl;
    std::cout << "                            For APSE v1.0                            " << std::endl;
    std::cout << "                                                                     " << std::endl;
    std::cout << " Systematic Bias Analyzer for Phylogenetic Estimation.               " << std::endl;
    std::cout << " Version 1.0 2020/01 Lab of Computational Biology and Bioinformatics " << std::endl;
    std::cout << " Seoul National University All right reserved.                       " << std::endl;
    std::cout << " Created by Junghwan Lee on 2020/03/15.                              " << std::endl;
    std::cout << " Copyright 2020 Lab of Computational Biology and Bioinformatics      " << std::endl;
    std::cout << " All rights reserved.                                                " << std::endl;
    
    Console con;
    con.linkInit();
    
    
    return 0;
    
}


