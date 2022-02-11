//
//  return.cpp
//  master0.0
//
//  Created by 이정환 on 2020/09/28.
//  Copyright © 2020 이정환. All rights reserved.
//

#include "return.hpp"
#include "console.hpp"
#include <iostream>

void Ret::exitNow() {
    system("clear"); // run terminal
    Console con;
    con.linkInit();
}

