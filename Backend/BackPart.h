#ifndef BACKPART_H
#define BACKPART_H

#include "Frontend/Part.h"
#include "Backend/BackEnd.h"
using namespace std;
using namespace world3000;

/*
 * 系统的后端优化的组件
 *
 */

namespace SLAMSystem {

class Frame;
class BackPixelPoint;
class AffLight;
class BackFrame
{
public:
    BackFrame(Frame* frame);    // 由前端建立相应的后端
private:
    inline void init();
public:
    void setTcw0(const SE3& Tcw0);
    void setX0(const SE3& Tcw0, const Vec8& x0);
    void setX0(const SE3 &Tcw0, const AffLight& aff_g2l);

    void set_state_x0(const Vec8& state_x0);

    void setX(const Vec8& x);

    AffLight getAffineParams();     // return (a,b)
    AffLight getAffineParams0();    //return (a0,b0)
public:
    Frame* frontData;
    vector<BackPixelPoint*> backPixelPoints;

//*******************************位姿
    SE3 Tcw0;   // evaluation point
    SE3 Tcw;
    SE3 Twc;

//*******************************1)关键帧的状态量Rtab
    Vec8 state_x0;  // state at evaluation point x0     [0-5: worldToCam-leftEps. 6-7: a,b]
    Vec8 state_x;   // state Now
    Vec8 state_x_backup;//R,t,(a,b), 直接法有光度参数(a,b)，加上位姿，先前状态的待优化变量可以表示为8维向量．(state_backup)
    Vec8 step;//
    Vec8 step_backup;//
    Vec8 delta_0ToNow;//状态量的累积增量    从当前状态减去最初状态
    Vec8 state_prior;   // 其值与settings相关
//******************************2)对于特征点法，此处关键帧的状态量只有R,t
    Vec6 state_x0_Feature;
    Vec6 state_x_Feature;
    Vec6 state_x_backup_Feature;//R,t,(a,b), 直接法有光度参数(a,b)，加上位姿，先前状态的待优化变量可以表示为8维向量．(state_backup)
    Vec6 step_Feature;//
    Vec6 step_backup_Feature;//
    Vec6 delta_0ToNow_Feature;//状态量的累积增量
    Vec6 state_prior_Feature;


//*****************************Nullspace
    Mat66 nullspace_pose;
    Mat42 nullspace_affine;
    Vec6  nullspace_scale;

// ID
    int keyframeId;    //       滑动窗口中该帧的序号

};

class PixelPoint;
class BackResidual;
class BackPixelPoint
{
public:
    BackPixelPoint(PixelPoint* point);

    float idepth0;
    float idepth;
    float idepth_backup;
    float step;
    float step_backup;
    float delta_idepth;     // deltaF    = idepth - idepth0

    float idepth_hessian;   // hessian value (inverse variance) of inverse depth.   ???
    float maxRelBaseline;   //                                                      ???

    vector<BackResidual*> backResiduals;

    float b_Sum;
    float H_depth;  // HdiF

    float prior;    // priorF


/* 线性增量方程H * deltaX = -b 中与点深度相关的H和b */
    // 线性化的相关量
/* 分为是否优化相机内参，不优化时H_C_idepth_accLF，H_C_idepth_accAF应设置为0 */
    float H_idepth_idepth_accLF;
    Vec4f H_C_idepth_accLF;
    float b_idepth_accLF;
    // 激活的相关量
/* 分为是否优化相机内参，不优化时H_C_idepth_accLF，H_C_idepth_accAF应设置为0 */
    float H_idepth_idepth_accAF;
    Vec4f H_C_idepth_accAF;
    float b_idepth_accAF;


//frontend
    PixelPoint* frontData;

};

class Residual;
class ResidualJacobian;

class BackResidual              //残差的雅可比J放到前端好呢还是后端好呢???
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    BackResidual(Residual* res, BackPixelPoint* backPnt, BackFrame* backHost, BackFrame* backTarget)
        : frontData(res), backPoint(backPnt), backHostFrame(backHost), backTargetFrame(backTarget),
          J(NULL/* todo*/), isLinearized(false), backResidualIdxInBackPoint(0)
    {

    }
    ~BackResidual() {}

public:
    void getJacobianPart1AndCalculateJacobianPart2();

public:
    BackPixelPoint* backPoint;
    BackFrame* backHostFrame;
    BackFrame* backTargetFrame;
    size_t backHostIndex, backTargetIndex;

    ResidualJacobian* J;

/* 直接法 */
    Vec8f res_toZeroF_direct;
/* 特征点法 */
    Vec2f res_toZeroF_feature;  //  <2,1>

    // Jacobian part2
/* 直接法 */
    Vec8f JpJdF;        // Jacobian矩阵中与逆深度相关的部分   (1.205)
/* 特征点法 */
    Vec6f JpJdF_Feature;    // Jacobian矩阵中与逆深度相关的部分   (1.205)


    bool isLinearized;     //

    // if residual is not OOB & not OUTLIER & should be used during accumulations
    bool isActiveAndIsGoodNEW;  // 标志位，残差的状态不是OOB,不是OUTLIER且在加速计算过程中被使用

    size_t backResidualIdxInBackPoint;
// frontend
    Residual* frontData;

};

class Camera;
class BackCamera
{
public:
    BackCamera(Camera* camera)
        : frontData(camera), intriParas(Vec4f::Constant(0)), step(Vec4f::Constant(0))
    {

    }
    ~BackCamera()
    {

    }

public:
    Vec4f intriParas_0;
    Vec4f intriParas;       // 高斯牛顿迭代优化在每一次迭代结束后相机内参的量
    float fx,fy,cx,cy;

/* 不优化相机内参时该项固定为0 */
    Vec4f step;             // 高斯牛顿迭代优化在每一次迭代结束后deltaX中相机内参的增量

    Vec4f intriParas_backup;    // 备份
    Vec4f step_backup;          // 备份

    Vec4f delta_0ToNow;      // 从初始的内参到当前内参的变化量

// frontend
    Camera* frontData;
};
}
#endif // BACKPART_H
