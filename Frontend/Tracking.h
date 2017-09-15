#ifndef TRACKING_H
#define TRACKING_H

#include "util/all_util_include.h"
#include "FullSystem/Settings.h"
#include "FullSystem/FullSystem.h"
#include "Frontend/Part.h"
#include "Backend/BackPart.h"
#include "Backend/BackEnd.h"

using namespace world3000;

using namespace std;

namespace SLAMSystem{



class EnergyFunction;
class Image;
class PixelPoint;
class RawPixelPoint;
class RawPixelPointTempResidual;
class FrameShell;
class Frame;
class Camera;
class Residual;

//class Residual<Feature>;

class Tracking
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Tracking()=default;
    Tracking(Image* image, int id, TrackMode T1, CameraMode T2); // tracking所有的流程
    ~Tracking();
    Tracking(const Tracking&)=delete;
    Tracking& operator =(const Tracking&)=delete;

public:
    void setEnergyFunction(EnergyFunction* pEF);
public:
// tracking流程(包括了前后端的操作)
    // 1. adds a new frame, and creates point & residual structs.
    Frame* addFrame(Image* image, int id);
    // 2. initialize
    void doInitialize();
    // 3. 确定永远跟踪新帧
    void trackNewestKeyframe();
    // 4. 更新shell及评估点
    void updateShellAndEvalPT(Frame* frame);    //
    // 5. 跟踪先前关键帧中的rawPixelPoint在当前帧中的状态
    void traceRawPixelPointInNowFrame(Frame* frame);
    // 6. 判断是否需要添加新的关键帧
    bool decideMakeNewKeyframe();
    // 7. 根据判断结果执行deliver操作
    void makeNonKeyframe(Frame* frame);
    void makeKeyframe(Frame* frame);
    /* 8~end called by makeKeyFrame() */
    // 8. 将某一帧设置为可被边缘化 use by makeKeyFrame
    void flagKeyframesForMarg(Frame* frame);
    // 9. 将该关键帧加入到KFtracking记录中(包括前后端)
    void addNewKeyframe(Frame* frame);
    // 10. 进行一些关键帧之间的预计算量的计算(包括前后端)
    void setPreCalcValues();
    // 11. add new residuals for old points 为系统添加新的残差
    void addNewResidualsForPixelPoints();
    // 12. 根据当前关键帧与先前关键帧的关系激活先前关键帧中的rawPoint, 激活后的Point用于之后的后端优化
    void activateRawPixelPoints();
    // 12.1 多线程激活点，点由rawPoints变成Points
    void activate(vector<PixelPoint*>* const activatedPoints, const vector<RawPixelPoint*>* const toActivatePoints);
    // 12.2 激活一个像素点
    PixelPoint* activateARawPixelPoint(RawPixelPoint* aRawPoint, RawPixelPointTempResidual* tempresiduals);
    // 13. 在后端中对激活后的所有PixelPoint以及关键帧进行统一编号
    void makeAllIdxForBackend();
    // 14. 在后端中设定优化次数，迭代终止阈值，开始G-N优化
    float runBackendOptimize(int IterNum);
    // 15. 优化后会有包含的残差项为0的pixelpoint,移除这些点(前后端)
    void removeOutliers();
    // 16. 对剩下的包含残差项不为0的点进行标记，判断这些点的状态(要被边缘化还是被直接丢弃)，然后执行相应的边缘化及丢弃操作．(前后端)
    void flagPixelPointsAndMargOrDrop();
    // 17. 为当前帧新建新的rawPixelPoint
    void makeRawPixelPointsForNowKeyframe();
    // 18. 边缘化边缘化标志被置位的关键帧
    void margFlagedKeyframe();
//*******************************************************************************************************************************


// backend
    Vec3 linearizeAll(bool b);
    void accumulateAllResidual(bool b, vector<Residual*>* todo, int first, int last , Vec10* stats, int tidx);

public:
    // 2. initialize
    bool hasInitialized;            //  dso系统标志位, 系统初始化标志位，已经初始化了就置位

    MultiThread<float> threadReduce1;


    vector<Frame*> allKeyframes;            // 系统所有的关键帧 (现在存在的)
    vector<FrameShell*> allKeyframeShells;  // 系统所有关键帧的Shell
    vector<Residual*> allResiduals;       // 系统所有激活的残差项

    vector<FrameShell*> allFrameShells;// track过的所有帧

    Camera* camera; // 相机内参类

/* 直接法，用于设置当前帧的残差的阈值(特征点法不需要) */
    vector<float> allResVec;    // targetframe为当前帧的所有点帧残差的能量值(state_NewEnergyWithOutlier)

// backend
    EnergyFunction* EF;  // 里程计实时优化的目标函数

    TrackMode trackMode;
    CameraMode cameraMode;
};
}
#endif // TRACKING_H
