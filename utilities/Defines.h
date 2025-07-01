#pragma once
#ifndef DEFINES_H_
#define DEFINES_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <random>
#include <assert.h> //stardand header file for c.
#include <omp.h>


#include <sys/types.h>
#include <sys/mman.h>  //memory management
#include <sys/stat.h>
#include <sys/time.h>  //gettimeofday function
#include <fcntl.h>  //IO in linux or unix : open and close function
#include <unistd.h>  //IO in linux or unix : read and write function

#include <dirent.h>
#include <limits>

#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <chrono>
#include <iomanip>
#include <cstdarg>

#include <fcntl.h> 
#include <sys/socket.h>  
#include <netinet/in.h>  
#include <arpa/inet.h>  
#include <functional>
#include <tuple>

#include <cstdarg>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <functional>
#include <future>
#include <stdexcept>
#include <type_traits>  // 包含 std::invoke_result

#include "CoreIndex.h"


#define NDEBUG // must precede cassert to disable assert.

using ui = unsigned int;
typedef long long unsigned ll_ui;
using std::string;
using std::ifstream;
using std::ofstream;

using std::cerr;
using std::cout;

using std::stoi;
using std::sort;
using std::ceil;
using std::max;
using std::abs;
using std::endl;
using std::to_string;

using std::pair;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::queue;

using std::hash;
using std::mt19937;
using std::random_device;
using std::uniform_int_distribution;


#define pb push_back
#define mp make_pair

const unsigned INF = INT32_MAX;
const uint DEFAULT_EDGE_BUF_SIZE = 5000;

enum GraphStore { uncompressed, byte_compressed, nibble_compressed };
enum GraphOrientation { original_graph, degree_oriented };

#endif /* DEFINES_H_ */