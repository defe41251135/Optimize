#ifndef PART_H
#define PART_H
#include "util/all_util_include.h"
#include "FullSystem/Settings.h"
#include "Backend/BackPart.h"
#include "Backend/BackEnd.h"
#include "util/projection.h"

using namespace std;
using namespace world3000;

/*
 * 构成SLAM系统的基本元素，包括帧，像素点，相机, residual
 */

namespace SLAMSystem {

class BackCamera;
class Camera
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    Camera() {}
    ~Camera() {}
public:
//    float getfx(){return fx;}
//    float getfy(){return fy;}
//    float getcx(){return cx;}
//    float getcy(){return cy;}

    Vec4f intriParas0;      // fx,fy,cx,cy
//    float fx,fy,cx,cy;


    // backend
    BackCamera* backData;
};



class Frame;
class RawPixelPoint
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    RawPixelPoint(int u, int v, Frame* phostFrame, float type, Camera* pCamera);
    ~RawPixelPoint();

    RawPixelPointState trackRawPixelPoints();// todo

public:
    float u, v;
    float idepth;
    Frame* hostFrame;
};
class Frame;
class Residual;
class BackPixelPoint;
class PixelPoint
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    PixelPoint(const RawPixelPoint* const rawPixelPoint);
    ~PixelPoint();

public:
    Frame* hostFrame;   // 生成该点的Frame   host-frame of point
    float u,v;      // 该点在hostframe中的像素坐标
    float idepth;   // 该点在hostframe中的逆深度

    vector<Frame*> observeFrames;   // 可观测到该点的关键帧
    vector<Vec3f> observePixelPoints;    // 该像素点在可观测到该点的关键帧中的像素坐标及逆深度    (这两个量在跟踪过程中不断更新，像素坐标及逆深度为观测值，用于计算特征点法的重投影误差)

    vector<Residual*> pointResiduals;   // 包含了该点所有的点帧残差 only contains good residuals (not OOB and not OUTLIER). Arbitrary order.
/* 直接法参数 */
    // static values    直接法中每个residual pattern中的灰度和权重   (对于特征点法，就是一个个单独的点，没有pattern一说)
    float gray[res_pattern_pointNum];
    float weight[res_pattern_pointNum];

    // backend
    BackPixelPoint* backData;
};

class KeyFrame
{
    // todo
};

class AffLight
{
public:
    AffLight(double a_, double b_)
        : a(a_), b(b_)
    {}
    AffLight()
        : a(0), b(0)
    {}
    static Vec2f getHostToTargetAffineParams(float hostExposureTime, float targetExposureTime, AffLight g2host, AffLight g2target)
    {
        double a = exp(g2host.a - g2target.a) * targetExposureTime/hostExposureTime;
        double b = g2target.b - a*g2host.b;
        return Vec2(a,b).cast<float>();
    }

public:
    double a;
    double b;
};

class Image
{
public:
    Image(int w_, int h_, double timestamp_=0)
        : w(w_), h(h_), timestamp(timestamp_), imageGray(NULL)
    {
        imageGray = new float [w*h];
    }
    Image(const Image& other) = delete;
    Image& operator = (const Image& other) = delete;
    ~Image()
    {
        delete [] imageGray;
    }

public:
    float* imageGray;
    int w;
    int h;
    double timestamp;
    float exposureTime;
};

class FrameShell
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    FrameShell()
        : id(0), incoming_id(0), timestamp(0.0), Trefc(SE3()), trackingRef(NULL), Twc(SE3()), marginalizedAt(-1)
    {}
public:
    size_t id;
    size_t incoming_id;
    double timestamp;

    // 只更新一次
    SE3 Trefc;                  // 相对于参考帧(最新的KF)的位姿的逆
    FrameShell* trackingRef;    // 该帧的参考帧(最新的KF)

    // 持续更新
    SE3 Twc;                    // 该帧相对于世界坐标系的位姿的逆
    AffLight aff_g2l;           // 从世界到该帧的仿射参数(a,b)
    int marginalizedAt;                 // 记录该帧由哪一帧进行的边缘化
};
class FrametoFrame;
class BackFrame;
class Frame                                     // 先把帧和关键帧写在一起，后来完善系统的时候再把它们分开
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Frame()
        : keyframeId(-1), keyframeShellId(-1), exposureTime(-1), image(NULL), Tcw0(SE3()), toBeMarginalization(false), shell(NULL), backData(NULL)
    {
        initPyramid();
    }
private:
    inline void initPyramid()
    {
        for(int i=0; i<PyramidLevels; i++)
        {
            imagePyramid[i] = NULL;
            pixelGradientNorm22[i] = NULL;
        }
    }
public:
    Frame(const Frame& other)=delete;
    Frame& operator=(const Frame& other)=delete;
    ~Frame()
    {

    }

    /**
     * @brief 构建图像金字塔，确定每层的像素灰度值以及像素梯度值
     * @param inputImage
     */
    void makeFrameFromImage(Image* inputImage)
    {
        for(unsigned i = 0; i < PyramidLevels; i++)
        {
            imagePyramid[i] = new Vec3f [widthPyramid[i]*heightPyramid[i]];
            pixelGradientNorm22[i] = new float [widthPyramid[i]*heightPyramid[i]];
        }
        // 确定最底层，即原始图像的灰度
        image = imagePyramid[0];
        imageWidth = widthPyramid[0];
        imageHeight = heightPyramid[0];
        for(unsigned i = 0; i < imageWidth*imageHeight; i++)
            image[i](0) = inputImage->imageGray[i];
        // 确定最底层以外的其他层的灰度值
        for(int i = 1; i < PyramidLevels; i++)
        {
            Vec3f* imagei = imagePyramid[i];
            Vec3f* imageireduce1 = imagePyramid[i-1];
            int imageireduce1width = widthPyramid[i-1];                               // imageireduce1: 2x2   ---->  imagei : 1x1
            for(int j = 0; j < heightPyramid[i]; j++)
                for(int k = 0; k < widthPyramid[i]; k++)
                {
                    imagei[k + j*widthPyramid[i]](0) = 0.25f * (imageireduce1[2*k + 2*j*imageireduce1width](0) +       //(二倍降采样的四个点取均值)
                            imageireduce1[2*k+1 + 2*j*imageireduce1width](0) +
                            imageireduce1[2*k + 2*j*imageireduce1width+imageireduce1width](0) +
                            imageireduce1[2*k+1 + 2*j*imageireduce1width+imageireduce1width](0))  ;
                }
        }
        // 确定每一层的每个像素点的像素梯度值
        for(unsigned i = i; i < PyramidLevels; i++)
        {
            Vec3f* imagei = imagePyramid[i];

            for(int idx = widthPyramid[i]; idx < widthPyramid[i]*(heightPyramid[i]-1); idx++)
            {
                float dpdy = 0.5f * (imagei[idx+widthPyramid[i]](0) - imagei[idx-widthPyramid[i]](0));
                float dpdx = 0.5f * (imagei[idx+1](0) - imagei[idx-1](0));                                  // 此处有待完善，图片边缘的像素点的像素梯度不应该这么算
                imagei[idx](1) = dpdx;
                imagei[idx](2) = dpdy;

                pixelGradientNorm22[i][idx] = (dpdx*dpdx + dpdy*dpdy);
            }
        }
    }


public:
// ID
    size_t keyframeId;
    size_t keyframeShellId;

//**********************************图片信息
    // 光度校准的相关变量
    float frameEnergyTH;	// 阈值，用于决定某一点帧残差的状态, 根据系统当前所有的点帧残差值设定     set dynamically depending on tracking residual
    double exposureTime;
/*    // 特征点法
    float* imagePyramid[PyramidLevels];
    float* image;       // == imagePyramid[0]   */
    // 直接法
    Vec3f* imagePyramid[PyramidLevels];
    Vec3f* image;       // == imagePyramid[0]
    int imageWidth, imageHeight;
    float* pixelGradientNorm22[PyramidLevels];      // 像素梯度的二范数的平方，用于稀疏直接法中的像素点选择

//************************** 帧中的不同状态的稀疏的像素点
    vector<PixelPoint*> pixelPointsActive;
    vector<PixelPoint*> pixelPointsMarginalized;
    vector<PixelPoint*> PixelPointsDiscard;
    vector<RawPixelPoint*> rawPixelPoints;

//***************************相机位姿
    SE3 Tcw0;

    SE3 Tcw;
    SE3 Twc;

//***********帧与帧之间的预计算量 R,t, a,b...tobe extend
    vector<FrametoFrame*> thisToTargets;

// 边缘化标志位，该帧是否被边缘化
    bool toBeMarginalization;

// 帧的外壳
    FrameShell* shell;


//backend
    BackFrame* backData;

};

class Camera;

class FrametoFrame
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    FrametoFrame(Frame* host, Frame* target, Camera* cam, TrackMode T);

    FrametoFrame() { hostframe = targetframe = NULL;}

    void PreCalculate(Frame* host, Frame* target, Camera* cam);

    ~FrametoFrame() {}

public:
    Frame* hostframe;
    Frame* targetframe;
    TrackMode trackMode;

    Mat33f R0_h2t;          //R hostframe to targetframe at evaluation point
    Vec3f  t0_h2t;

    Mat33f R_h2t;          //R hostframe to targetframe
    Vec3f  t_h2t;           //t

    Mat33f KRKinv;          // K * R * K.inverse()
    Mat33f RKinv;           // R * K.inverse()
    Vec3f  Kt;              // K * t

    // 直接法参数(a,b),在特征点法中不需要
    Vec2f afflight;          //  从hostframe到targetframe的仿射亮度参数(a,b)
    float b0;              // hostframe在评估点处的仿射亮度参数(a0,b0)中的b0       at evaluation point
};



class ResidualJacobian
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    ResidualJacobian()
    {
        init();
    }
    ~ResidualJacobian()
    {}
private:
    void init()
    {
        res_direct = Vec8f::Constant(0);
        J_p_d_xi[0] = J_p_d_xi[1] = Vec6f::Constant(0);
    }

public:
// 残差
/* 对于直接法，residual pattern中有8个点 */
    Vec8f res_direct;
/* 对于特征点法，residual就是该点的重投影误差 */
    float res_feature_Norm;
    Vec2f res_feature_uv;

/* 对于直接法，有灰度值对像素坐标的偏导，即像素梯度 */
    Vec8f J_I_d_p[2];
/* 对于特征点法，没有像素梯度这一项 */

/* 对于直接法，有残差关于仿射亮度参数(a,b)的偏导 */
    Vec8f J_res_d_aff[2];   // Jphoto
/* 对于特征点法，没有仿射亮度参数这一项 */

    Vec6f J_p_d_xi[2];      // 投影点p'(u',v')对位姿(R,t)的偏导
    Vec4f J_p_d_cam[2];     // 投影点p'(u',v')对相机内参(焦距，主点)的偏导
    Vec2f J_p_d_idepth;     // 投影点p'(u',v')对idepth的偏导

/* 直接法残差的雅可比中的一些中间计算量 */
    Mat22f J_I_d_p_2;           // J_I_d_p^T * J_I_d_p
    Mat22f J_res_d_aff_J_I_d_p; // J_res_d_aff^T * J_I_d_p
    Mat22f J_res_d_aff_2;       // J_res_d_aff^T * J_res_d_aff
};

class BackResidual;

class RawPixelPointTempResidual     // RawPoint对应的临时的残差项  17.8.14
{
public:
    ResState tempRes_State;
    ResState tempRes_NewState;

    double tempRes_Energy;
    double tempRes_NewEnergy;

    Frame* hostFrame;
    Frame* targetFrame;
};

/*
 * 残差项，hostframe与targetframe之间在某一点(该点属于hostframe)的残差
 */
class Residual
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Residual() {assert(false);}
    Residual(PixelPoint* pnt, Frame* host, Frame* target, TrackMode T1)
        : point(pnt), hostFrame(host), targetFrame(target), resState(Res_INIT),
          residualIdxInPoint(0)
    {
        Jaco == new ResidualJacobian();
//        if(T == Direct)
//            Jaco = new ResidualJacobian_Direct();
//        else if(T == Feature)
//            Jaco = new ResidualJacobian_Feature();
    }

    ~Residual()
    {
        if (Jaco != NULL) {delete Jaco; Jaco=NULL;}
    }

public:
    double calcuResidualAndJacobianPart1(Camera* cam, TrackMode T1, CameraMode T2);

    void applyRes();
    void resetOOB();


public:
    TrackMode trackMode;

    PixelPoint* point;
    Frame* hostFrame;
    Frame* targetFrame;
    ResidualJacobian* Jaco;

    ResState resState;
    double energy;

    ResState newResState;
    double newEnergy;

/* 直接法，该变量用于targetframe能量阈值的确定 */
    double newEnergyForTH;


    size_t residualIdxInPoint;


//backend
    BackResidual* backData;
};



}
#endif // PART_H
