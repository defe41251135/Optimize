#ifndef SETTINGS_H
#define SETTINGS_H

#include "util/all_util_include.h"
using namespace world3000;

// 用全局变量保存系统的参数设置
namespace SLAMSystem {
/* GN迭代优化中求解增量方程H　* deltaX = -b 的模式 */
#define SOLVER_SVD (int)1
#define SOLVER_ORTHOGONALIZE_SYSTEM (int)2
#define SOLVER_ORTHOGONALIZE_POINTMARG (int)4
#define SOLVER_ORTHOGONALIZE_FULL (int)8
#define SOLVER_SVD_CUT7 (int)16
#define SOLVER_REMOVE_POSEPRIOR (int)32
#define SOLVER_USE_GN (int)64
#define SOLVER_FIX_LAMBDA (int)128
#define SOLVER_ORTHOGONALIZE_X (int)256
#define SOLVER_MOMENTUM (int)512
#define SOLVER_STEPMOMENTUM (int)1024
#define SOLVER_ORTHOGONALIZE_X_LATER (int)2048


/* pyramid 图片金字塔 */
#define PyramidLevels 8
/* 每层的宽和高 */
extern int widthPyramid[PyramidLevels];
extern int heightPyramid[PyramidLevels];
/* 每层的相机内参 */
extern float fxPyramid[PyramidLevels];
extern float fyPyramid[PyramidLevels];
extern float cxPyramid[PyramidLevels];
extern float cyPyramid[PyramidLevels];
extern Mat33f KPyramid[PyramidLevels];
// 阈值,用于判断重投影后像素坐标是否出界
extern float width_thres;
extern float height_thres;

void setPicturePyramid(int width, int height, const Mat33f& K);


#define MULTITHREAD

extern int settings_ThreadNum;          //
#define ThreadNum 8
extern bool settings_UseMultiThread;

extern int settings_OptimizeIterNum;

extern float settings_huberTH;
extern float settings_pixelCoordTH;

extern float settings_outlierTHSumComponent;


// frameEnergyTH
extern float settings_frameEnergyTHConstWeight;
extern float settings_frameEnergyTH;
extern float settings_frameEnergyTHFacMedian;
extern float settings_overallEnergyTHWeight;
extern float settings_coarseCutoffTH;

extern float settings_affineOptModeA;
extern float settings_affineOptModeB;

extern bool settings_forceAceptStep;
extern int settings_solverMode;
extern double settings_solverModeDelta;

extern float settings_thresGNOptimize;

// #define scale 逆深度idepth，仿射亮度参数(a,b)，相机内参，Ｒ，ｔ都需要单位换算成统一的单位，来进行计算，这里定义相应的换算比例
#define Idepth_scale 1.0f
#define aff_a_scale 1.0f
#define aff_b_scale 1.0f
#define cam_f_scale 1.0f
#define cam_c_scale 1.0f
#define R_scale 1.0f
#define t_scale 1.0f



#define res_pattern_pointNum 8      // 直接法residual pattern 中像素点的个数
extern int staticPattern[10][40][2];
#define patternP staticPattern[8]
}
#endif // SETTINGS_H
