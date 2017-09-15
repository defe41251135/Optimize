#include "FullSystem/Settings.h"


namespace SLAMSystem {

int widthPyramid[PyramidLevels];
int heightPyramid[PyramidLevels];
float fxPyramid[PyramidLevels];
float fyPyramid[PyramidLevels];
float cxPyramid[PyramidLevels];
float cyPyramid[PyramidLevels];
Mat33f KPyramid[PyramidLevels];

float width_thres;
float height_thres;

void setPicturePyramid(int width, int height, const Mat33f &K)      // used in FullSystem
{
    widthPyramid[0] = width;
    heightPyramid[0] = height;
    KPyramid[0] = K;
    fxPyramid[0] = K(0,0);
    fyPyramid[0] = K(1,1);
    cxPyramid[0] = K(0,2);
    cyPyramid[0] = K(1,2);

    assert(pow(2,(PyramidLevels-1)) <= width && pow(2,(PyramidLevels-1)) <= height);
    for(int level = 1; level < PyramidLevels; level++)
    {
        widthPyramid[level] = width / pow(2,level);
        heightPyramid[level] = height / pow(2,level);

        fxPyramid[level] = fxPyramid[level-1] / 2;
        fyPyramid[level] = fyPyramid[level-1] / 2;
        cxPyramid[level] = (cxPyramid[0] + 0.5f) / pow(2, level) - 0.5f;
        cyPyramid[level] = (cyPyramid[0] + 0.5f) / pow(2, level) - 0.5f;
        KPyramid[level] << fxPyramid[level], 0.0f, cxPyramid[level], 0.0f, fyPyramid[level], cyPyramid[level], 0.0f, 0.0f, 1;

    }

    width_thres = width - 1.0f;
    height_thres = height - 1.0f;
}

int settings_ThreadNum = 8;
bool settings_UseMultiThread = true;

//int settings_PyramidLevels = 8;

int settings_OptimizeIterNum = 5;

float settings_huberTH = 9; //  Huber Threshold   huber核函数的阈值
float settings_pixelCoordTH = 3; // 特征点法中像素点的像素坐标观测值与预测值之间的差值的阈值
/* Outlier Threshold on photometric energy */   //光度误差的外点阈值     直接法使用，特征点法不需要
float settings_outlierTHSumComponent = 50*50; 		// higher -> less strong gradient-based reweighting .   基于梯度的权重调整

// parameters controlling adaptive energy threshold computation.
float settings_frameEnergyTHConstWeight = 0.5;
float settings_frameEnergyTH = 0.7f;
float settings_frameEnergyTHFacMedian = 1.5;
float settings_overallEnergyTHWeight = 1;
float settings_coarseCutoffTH = 20;

float settings_affineOptModeA = 1e12; //-1: fix. >=0: optimize (with prior, if > 0).     光度参数a是否优化, 小于0不优化，大于等于0优化
float settings_affineOptModeB = 1e8; //-1: fix. >=0: optimize (with prior, if > 0).

// GN迭代优化时强制使用给定步长,先只考虑为ture的情况
bool settings_forceAceptStep = true;

// 设置求解线性系统Hx=b的模式
/* some modes for solving the resulting linear system (e.g. orthogonalize wrt. unobservable dimensions) */
int settings_solverMode = SOLVER_FIX_LAMBDA | SOLVER_ORTHOGONALIZE_X_LATER;
double settings_solverModeDelta = 0.00001;  // SVD分解求解deltaX用到

float settings_thresGNOptimize = 1.2;    // GN迭代优化终止阈值,当某次优化前后状态差值小于阈值时停止优化

int staticPattern[10][40][2] = {
        {{0,0}, 	  {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},	// .    只取该点
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

        {{0,-1},	  {-1,0},	   {0,0},	    {1,0},	     {0,1}, 	  {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},	// +    取该点及该点上下左右共5点
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

        {{-1,-1},	  {1,1},	   {0,0},	    {-1,1},	     {1,-1}, 	  {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},	// x    取该点及该点左下右下左上右上共5点
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

        {{-1,-1},	  {-1,0},	   {-1,1},		{-1,0},		 {0,0},		  {0,1},	   {1,-1},		{1,0},		 {1,1},       {-100,-100},	// full-tight   取该点以及紧邻该点的共9个点
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

        {{0,-2},	  {-1,-1},	   {1,-1},		{-2,0},		 {0,0},		  {2,0},	   {-1,1},		{1,1},		 {0,2},       {-100,-100},	// full-spread-9
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

        {{0,-2},	  {-1,-1},	   {1,-1},		{-2,0},		 {0,0},		  {2,0},	   {-1,1},		{1,1},		 {0,2},       {-2,-2},   // full-spread-13
         {-2,2},      {2,-2},      {2,2},       {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

        {{-2,-2},     {-2,-1}, {-2,-0}, {-2,1}, {-2,2}, {-1,-2}, {-1,-1}, {-1,-0}, {-1,1}, {-1,2}, 										// full-25
         {-0,-2},     {-0,-1}, {-0,-0}, {-0,1}, {-0,2}, {+1,-2}, {+1,-1}, {+1,-0}, {+1,1}, {+1,2},
         {+2,-2}, 	  {+2,-1}, {+2,-0}, {+2,1}, {+2,2}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

        {{0,-2},	  {-1,-1},	   {1,-1},		{-2,0},		 {0,0},		  {2,0},	   {-1,1},		{1,1},		 {0,2},       {-2,-2},   // full-spread-21
         {-2,2},      {2,-2},      {2,2},       {-3,-1},     {-3,1},      {3,-1}, 	   {3,1},       {1,-3},      {-1,-3},     {1,3},
         {-1,3},      {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

        {{0,-2},	  {-1,-1},	   {1,-1},		{-2,0},		 {0,0},		  {2,0},	   {-1,1},		{0,2},		 {-100,-100}, {-100,-100},	// 8 for SSE efficiency     论文中的figure4
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
         {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

        {{-4,-4},     {-4,-2}, {-4,-0}, {-4,2}, {-4,4}, {-2,-4}, {-2,-2}, {-2,-0}, {-2,2}, {-2,4}, 										// full-45-SPREAD
         {-0,-4},     {-0,-2}, {-0,-0}, {-0,2}, {-0,4}, {+2,-4}, {+2,-2}, {+2,-0}, {+2,2}, {+2,4},
         {+4,-4}, 	  {+4,-2}, {+4,-0}, {+4,2}, {+4,4}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200},
         {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}},
};
}
