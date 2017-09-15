#ifndef BACKAND_H
#define BACKAND_H

#include "FullSystem/FullSystem.h"
#include "Frontend/Tracking.h"
#include "Backend/BackPart.h"
#include "util/all_util_include.h"
#include "Backend/CalculateHessian/AccumulatedTopHessian.h"
#include "Backend/CalculateHessian/AccumulatedSCHessian.h"
using namespace world3000;


namespace SLAMSystem {



class Tracking;
class BackFrame;
class BackPixelPoint;
class Residual;
class Camera;
//class Frame;
//class PixelPoint;
//class Camera;

//class EFFrame;
//class EFPixelPoint;
//class EFResidual;
class AccumulatedSCHessian;

//template<TrackMode T1, CameraMode T2> class AccumulatedTopHessian;
class AccumulatedTopHessian;

/*
 * 系统的能量函数类，Tracking过程中后端滑动优化
 */
class EnergyFunction
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    EnergyFunction();
    EnergyFunction(TrackMode T1, CameraMode T2);
    void setTracking(Tracking* pTracking);      // 传递指针的函数，在FullSystem中完成指针的传递操作    (仿照orbSLAM)
    ~EnergyFunction();

    void makeIndex();


    float backGNOptimize(int IterNum);
private:
/* called by backGNOptimize() */
    Vec3 accumulateAllRes(bool fixLinearization);  // Vec3 linearizeAll(bool fixLinearization);
    void accumulateAllRes_Reductor(bool fixLinearization, vector<Residual*>* res_toRemove, int min, int max, float* basicUnit, int threadId);
        /* 用于直接法 */
    void setNewFrameEnergyTH();     // 对于直接法，更新帧的能量阈值

    double calculateLinearizeEnergy();
    double calculateMarginalizeEnergy();

    void applyResidual();
    void applyResidual_Reductor(int min, int max, Vec10* basicUnit, int threadId);

    void backupState(bool backupLastStep);  // 备份每次优化前的量

    void solveSystem(int iteration, double lambda); //
        void getNullspace(vector<VecX>& nullspace_pose, vector<VecX>& nullspace_scale,
                      vector<VecX>& nullspace_affA, vector<VecX>& nullspace_affB);      /* scale 考虑去掉，特征点法不更新光度参数 */
        void solve(int iteration, double lambda, Camera* cam);
            void accumulateActivatePart(MatXX& H, VecX& b, bool MT);
            void accumulateLinearPart(MatXX& H, VecX& b, bool MT);
            void accumulateSchurPart(MatXX& H, VecX& b, bool MT);
            VecX getStitchDelta();
            void calculateDeltaXBySVD(MatXX& H, VecX& b, VecX& deltaX);
            void feedBackDeltaX(VecX delta_X, Camera* cam, bool MT);
                void feedBackDeltaIdepthDirect(const Vec4f& deltaX_cam, Mat18f* deltaX_xi_adjoint, int min, int max, Vec10* basicUnit, int thd_idx);
                void feedBackDeltaIdepthFeature(const Vec4f &deltaX_cam, Mat16f *deltaX_xi_adjoint, int min, int max, Vec10 *basicUnit, int thd_idx);

    bool addStepFromBackup();
public:
    Tracking* m_pTracking;      // frontend
    MultiThread<float>* threadReduce1;  // 多线程求解器,基本单元是float浮点型, 用于累加所有残差值（对于重投影误差的光度误差都是float值）
    MultiThread<Vec10>* threadReduce2;

    TrackMode trackMode;
    CameraMode cameraMode;      // 需要赋初值
public:


    vector<BackFrame*> allBackKeyFrames;
    vector<BackPixelPoint*> allBackPixelPoints;
    vector<BackPixelPoint*> allBackPixelPointsToMarg;


    MatXX H_M;  // 系统累积的边缘化后的像素点对应的H的拼接HM
    VecX b_M;

    int allBackKeyFramesNum;
    int allBackPixelPointsNum;
    int allBackResidualsNum;

    int resNum_Active;
    int resNum_Linearized;
    int resNum_Marginalize;


    // delta
/* 直接法 */
    Mat18f* adHTdeltaF; // 1x8的数组, 表示从host到target的状态变化量
/* 特征点法 */
    Mat16f* adHTdeltaF_Feature; //

/* 是否优化内参, 不优化时delta为0 */
    Vec4f deltaF_C; // delta_C

    Vec4 cPrior;

/* 直接法 */
    Mat88* adHost;// hostframe 的伴随矩阵  double类型
    Mat88* adTarget;// targetframe 的伴随矩阵  double类型

/* 特征点法 */
    Mat66* adHost_Feature;// hostframe 的伴随矩阵  double类型
    Mat66* adTarget_Feature;// targetframe 的伴随矩阵  double类型

    // Nullspace
    vector<VecX> lastNullspace_pose;
    vector<VecX> lastNullspace_scale;
    vector<VecX> lastNullspace_affA;
    vector<VecX> lastNullspace_affB;

    // 累加的全局的H,b
    AccumulatedTopHessian* acc_top_Activate;
    AccumulatedTopHessian* acc_top_Linearized;
    AccumulatedSCHessian* acc_bot_Marginalize;

    // 求解整体的线性方程Ｈ * deltaｘ ＝ ｂ的 H,deltax,b的最新值
    MatXX last_H;
    VecX last_b;
    VecX last_deltaX;

};
}

#endif // BACKAND_H
