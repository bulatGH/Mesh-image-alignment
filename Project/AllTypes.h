#ifndef ALLTYPES_H
#define ALLTYPES_H

#include <ppl.h>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include "stdafx.h"
#define NOMINMAX
#include <Windows.h>
#include <algorithm>
#include <string>
#include "linalg.h"
#include "stdafx.h"
#include <ctime>
#include <cstdlib>
#include <random>
#include <time.h>
#include <thread>
#include <mutex> 


using namespace concurrency;
using namespace std;

typedef unsigned char byte;

const short minValueS = -32000;
const short maxValueS =  32000;

const int minValueI = -2000000000;
const int maxValueI =  2000000000;

const float minValueF = -10000000;
const float maxValueF =  10000000;

const double PI = 3.14159265359;


#endif