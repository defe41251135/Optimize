#ifndef INCLUDE_3000WORLD_H
#define INCLUDE_3000WORLD_H
/* linux */
#include <dirent.h>                 // linux系统文件夹操作
#include <sys/time.h>               // linux系统时间操作      gettimeofday()
#include <time.h>                   // clock_t
#include <unistd.h>                 // 提供对 POSIX 操作系统 API 的访问功能的头文件

/* SSE */
#include <xmmintrin.h>              // SSE
#include <emmintrin.h>              // SSE2

/* C */
#include <stdlib.h>                 // C标准库函数       rand()随机数发生器
#include <stdio.h>                  // C标准输入输出使用
#include <locale.h>                 // C本地化函数。 这些函数用于在处理多种自然语言的软件编程设计时，把程序调整到特定的区域设置.。这些区域设置影响到C语言标准库的输入/输出函数。
#include <cstring>                  // std::memset():数值数组的初始化   std::memcpy(): 源指针指向的内容复制到目的指针指向的内容

/* STL */
#include <random>
#include <list>
#include <vector>                   // STL 动态数组容器
#include <queue>                    // STL 队列容器
#include <deque>                    // STL 双端队列容器
#include <set>
#include <unordered_set>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <memory>                   // 内存操作相关, 智能指针
#include <iterator>                 // 迭代器      std::advance
#include <numeric>                  // 数值计算
#include <cmath>                    // 数学计算
#include <chrono>                   // 计时
#include <assert.h>                 // 断言
#include <thread>                   // 多线程 std::thread
#include <atomic>                   // std::atomic_flag可用于多线程之间的同步操作，类似于linux中的信号量。使用atomic_flag可实现mutex.
#include <mutex>                    // 锁    std::mutex, std::unique_lock
#include <condition_variable>       // 条件变量 std::condition_variable
#include <functional>               // C++标准库bind()函数   bind函数的最根本的作用就是可以把一个参数较多的函数给封装成参数较少的函数
#include <iomanip>                  // setprecision() 设置浮点数的精度
#include <utility>                  //STL 通用模板类
#include <exception>                // 异常类
using namespace std::placeholders;  // 占位符命名空间



/* Eigen */
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>



/* SSE2 */
#include <emmintrin.h>              // SSE2头文件，此头文件里包含SSE头文件

///* g2o */
//#include <g2o/core/block_solver.h>
//#include <g2o/core/optimization_algorithm_levenberg.h>
//#include <g2o/core/base_vertex.h>
//#include <g2o/core/base_unary_edge.h>
//#include <g2o/core/base_binary_edge.h>
//#include <g2o/core/base_multi_edge.h>
//#include <g2o/core/robust_kernel_impl.h>

//#include <g2o/types/sba/types_six_dof_expmap.h>
//#include <g2o/types/sba/types_sba.h>
////#include <g2o/types/slam3d/vertex_se3.h>    // VertexSE3
//#include <g2o/types/slam3d/types_slam3d.h>
//#include <g2o/types/slam3d/isometry3d_mappings.h>

//#include <g2o/solvers/eigen/linear_solver_eigen.h>



/* openCV */
//#include <opencv/cv.h>
#include <opencv2/core/core.hpp>
//#include <opencv2/core/types_c.h>
//#include <opencv2/imgproc/imgproc.hpp>          //cv::cvtColor
//#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>



/* Pangolin */
#include <pangolin/pangolin.h>



/* sophus */
#include "sophus/sim3.hpp"
#include "sophus/se3.hpp"



/* ziplib */
#include <zip.h>            // 读取.zip压缩文件

#endif // INCLUDE_3000WORLD_H
