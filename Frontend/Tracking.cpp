#include "Frontend/Tracking.h"

namespace SLAMSystem {


Tracking::Tracking(Image *image, int id, TrackMode T1, CameraMode T2)
    : trackMode(T1), cameraMode(T2)
{
    /* tracking class object初始化 todo */


//************************Run********************************
    // 1. adds a new frame
    Frame* frame = addFrame(image, id);
    // 2. initialize
    if(!hasInitialized)     doInitialize();
    // 3. 确定永远跟踪新帧
    trackNewestKeyframe();
    // 4. 更新shell及评估点
    updateShellAndEvalPT(frame);
    // 5. 跟踪先前关键帧中的rawPixelPoint在当前帧中的状态
    traceRawPixelPointInNowFrame(frame);
    // 6. 判断是否需要添加新的关键帧
    bool make = decideMakeNewKeyframe();
    // 7. 根据判断结果执行deliver操作
    if(!make) makeNonKeyframe(frame);
    else makeKeyframe(frame);
    // 8.
    //todo
//**********************End**********************************
}

Tracking::~Tracking()
{
    //todo
}

void Tracking::setEnergyFunction(EnergyFunction *pEF)
{
    EF = pEF;
}

Frame *Tracking::addFrame(Image *image, int id)     // should be called by Tracking conductor
{
    // front
    Frame* frame = new Frame();
    FrameShell* shell = new FrameShell();

    shell->id = (size_t)allFrameShells.size();
    shell->incoming_id = id;
    shell->timestamp = image->timestamp;//
    allFrameShells.push_back(shell);
    frame->shell = shell;
    frame->exposureTime = image->exposureTime;
    frame->makeFrameFromImage(image);
    // back
    BackFrame* backFrame = new BackFrame(frame);
    frame->backData = backFrame;

    return frame;
}

void Tracking::doInitialize()
{
    // todo
}

void Tracking::trackNewestKeyframe()
{
    // todo
}
void Tracking::updateShellAndEvalPT(Frame *frame)
{
    frame->shell->Twc = frame->shell->trackingRef->Twc * frame->shell->Trefc;
//    frame->backData->setTcw0(frame->shell->Twc.inverse());
    frame->backData->setX0(frame->shell->Twc.inverse(), frame->shell->aff_g2l);
}

void Tracking::traceRawPixelPointInNowFrame(Frame *frame)
{
    Mat33f K = Mat33f::Identity();
    K(0,0) = camera->backData->fx; K(1,1) = camera->backData->fy; K(0,2) = camera->backData->cx; K(1,2) = camera->backData->cy;
    // 计算从先前关键帧到当前帧的位姿变换Tji，以及一些中间量K*R*Kinv，K*t，从先前帧到当前帧的仿射亮度参数(a,b)
    for(Frame* hostFrame : allKeyframes)
    {
        SE3 T_hostToNow = frame->Tcw * hostFrame->Twc;
        Vec2f affineParameters = AffLight::getHostToTargetAffineParams(hostFrame->exposureTime,frame->exposureTime,hostFrame->backData->getAffineParams(),frame->backData->getAffineParams());
        // 判断先前关键帧中的每个raw像素点在当前帧中的状态
        for(RawPixelPoint* rawPoint : hostFrame->rawPixelPoints)    // 对于稀疏的特征点法，对应未激活的特征点
        {
            rawPoint->trackRawPixelPoints();
        }
    }
}
bool Tracking::decideMakeNewKeyframe()
{
    // todo
}

//void Tracking::deliverTrackedFrame(Frame *frame, bool make)
//{
//    if(make) makeKeyframe(frame);
//    else makeNonKeyframe(frame);
//    // todo
//}

void Tracking::makeKeyframe(Frame *frame)
{
    addNewKeyframe(frame);// todo 7.11

    setPreCalcValues();

    addNewResidualsForPixelPoints();

    activateRawPixelPoints();

    makeAllIdxForBackend();

    //Begin Optimize
    runBackendOptimize(settings_OptimizeIterNum);

    removeOutliers();

    flagPixelPointsAndMargOrDrop();

    makeRawPixelPointsForNowKeyframe();

    margFlagedKeyframe();
}

void Tracking::makeNonKeyframe(Frame *frame)
{
    delete frame;
}

void Tracking::flagKeyframesForMarg(Frame *frame)
{
    //1. 如果setting_minFrameAge超过了设定的滑动窗口中帧的最大个数N，那么对所有帧中除了最后N帧的可被边缘化标志位置位

    //2. flag marginalize all frames that have not enough points.
    /* 进行某些帧的边缘化标志位置位操作，这些帧的选取满足以下条件：
    *      1. 该帧观测到的3D点中，in/(in+out)<0.05 或者 从最后一帧到该帧的仿射亮度参数小于阈值
    *      2. 剩余的未被边缘化的帧的个数必须大于滑动窗口中帧最少个数阈值
    */

    //3. marginalize one.
    // 如果上一步的边缘化置位操作以后，剩余的帧数大于滑动窗口中帧的最大个数N，那么继续边缘化置位操作
    /* 要被边缘化的帧的选取依据：论文公式(20)
     * 找到要边缘化的帧toMarginalize,该帧使得(20)得分最大
     */
}

void Tracking::addNewKeyframe(Frame *frame)
{
    // frontend
    frame->keyframeId = allKeyframes.size();
    allKeyframes.push_back(frame);
    frame->keyframeShellId = allKeyframeShells.size();
    allKeyframeShells.push_back(frame->shell);
    // 前端加进去以后，后端也要加进去
// =========================== 在系统的后端中加入该KF ======================================
//    ef->insertFrame(fh, &Hcalib);               // energy function插入当前帧(当前帧fh已经作为关键帧)
    EF->allBackKeyFrames.push_back(frame->backData);
    EF->allBackKeyFramesNum++;
    // H_M, b_M  矩阵扩容
    EF->H_M.conservativeResize(EF->allBackKeyFramesNum*8+4, EF->allBackKeyFramesNum*8+4);
    EF->b_M.conservativeResize(EF->allBackKeyFramesNum*8+4);
    EF->H_M.rightCols<8>().setZero();
    EF->H_M.bottomRows<8>().setZero();
    EF->b_M.tail<8>().setZero();
    // 设置相应的伴随矩阵 todo
}

void Tracking::setPreCalcValues()
{
    // 1. front: FrametoFrame
    for(size_t i = 0; i < allKeyframes.size(); i++)
    {
        Frame* frame = allKeyframes[i];
        frame->thisToTargets.resize(allKeyframes.size());
    }
    for(size_t i = 0; i < allKeyframes.size(); i++)
    {
        for(size_t j = 0; j < allKeyframes.size(); j++)
        {
            Frame* frame = allKeyframes[i];
            if(frame->thisToTargets[j]!=NULL) continue;
            frame->thisToTargets[j] = new FrametoFrame(frame,allKeyframes[j],camera,trackMode);
        }
    }

// 2. back: delta
    if(trackMode == Direct)
    {
        // 2.1 adHTdeltaF(从host到target的状态变化量)
        if(EF->adHTdeltaF!=NULL) delete [] EF->adHTdeltaF;
        EF->adHTdeltaF = new Mat18f [EF->allBackKeyFramesNum * EF->allBackKeyFramesNum];
        for(size_t i = 0; i < EF->allBackKeyFramesNum; i++)
            for(size_t j = 0; j < EF->allBackKeyFramesNum; j++)
            {
                int idx = i+j*EF->allBackKeyFramesNum;
                EF->adHTdeltaF[ idx ] = EF->allBackKeyFrames[i]->delta_0ToNow.cast<float>().transpose() * EF->adHost[idx].cast<float>() +
                        EF->allBackKeyFrames[j]->delta_0ToNow.cast<float>().transpose() * EF->adTarget[idx].cast<float>();
            }
        // 2.2 deltaF_C
        if(cameraMode == Optimize)
            EF->deltaF_C = camera->backData->delta_0ToNow.cast<float>();
        else if(cameraMode == Fix)
            EF->deltaF_C = Vec4f::Zero();
        // 2.3 deltaF_di
        for(BackFrame* frame : EF->allBackKeyFrames)
            for(BackPixelPoint* point : frame->backPixelPoints)
            {
                point->delta_idepth = point->idepth - point->idepth0;
            }
    } else if(trackMode == Feature)
    {
        if(EF->adHTdeltaF_Feature!=NULL) delete [] EF->adHTdeltaF_Feature;
        EF->adHTdeltaF_Feature = new Mat16f [EF->allBackKeyFramesNum * EF->allBackKeyFramesNum];
        for(size_t i = 0; i < EF->allBackKeyFramesNum; i++)
            for(size_t j = 0; j < EF->allBackKeyFramesNum; j++)
            {
                int idx = i+j*EF->allBackKeyFramesNum;
                EF->adHTdeltaF_Feature[idx] = EF->allBackKeyFrames[i]->delta_0ToNow_Feature.cast<float>().transpose() * EF->adHost_Feature[idx].cast<float>() +
                        EF->allBackKeyFrames[j]->delta_0ToNow_Feature.cast<float>().transpose() * EF->adTarget_Feature[idx].cast<float>();
            }
        if(cameraMode == Optimize)
            EF->deltaF_C = camera->backData->delta_0ToNow.cast<float>();
        else if(cameraMode == Fix)
            EF->deltaF_C = Vec4f::Zero();
        for(BackFrame* frame : EF->allBackKeyFrames)
            for(BackPixelPoint* point : frame->backPixelPoints)
            {
                point->delta_idepth = point->idepth - point->idepth0;
            }
    }
}
/**
 * @brief 新建以nowFrame为targetFrame的残差项, 放入对应的点中
 */
void Tracking::addNewResidualsForPixelPoints()
{

    Frame* nowFrame = allKeyframes.back();
    for(Frame* frame : allKeyframes)
    {
        if(frame == nowFrame) continue;
        for(PixelPoint* point : frame->pixelPointsActive)
        {
            // 1. 新建以nowFrame为targetFrame的残差项, 放入对应的点中
            Residual* res = new Residual(point, frame, nowFrame, trackMode);
            res->resState = Res_IN;
            res->residualIdxInPoint = point->pointResiduals.size();
            point->pointResiduals.push_back(res);

            // 2. 后端相应地加入该残差项, 放入对应的点中 (后端由前端建立)
            BackResidual* backRes = new BackResidual(res, point->backData, frame->backData, nowFrame->backData);
            backRes->backResidualIdxInBackPoint = point->backData->backResiduals.size();
            point->backData->backResiduals.push_back(backRes);
            EF->allBackResidualsNum++;

            res->backData = backRes;

        }
    }
}

void Tracking::activateRawPixelPoints()
{
    Frame* nowFrame = allKeyframes.back();

    vector<RawPixelPoint*> toActivate; toActivate.reserve(20000);
    for(Frame* frame : allKeyframes)
    {
        if(frame == nowFrame) continue;
        //经过一些条件筛选以及中间量的计算 todo 8.14
        for(size_t i=0; i < frame->rawPixelPoints.size(); i++)
        {
            //经过一些条件筛选以及中间量的计算 todo 8.14
            // 选出需要激活的点，放到toActivate中 todo      toActivate.push_back() 8.14
        }
    }

    vector<PixelPoint*> activated; activated.resize(toActivate.size()); // 激活后的点放到activated中

    //调用Tracking::activate多线程处理,激活点
    auto newfunction = std::bind(&Tracking::activate, this, &activated, &toActivate);
    threadReduce1.reduce(newfunction, 0, toActivate.size(), 50);
}

void Tracking::activate(vector<PixelPoint *> * const activatedPoints, const vector<RawPixelPoint *> * const toActivatePoints)
{
    RawPixelPointTempResidual* residuals = new RawPixelPointTempResidual [toActivatePoints->size()];
    for(size_t i = 0; i < toActivatePoints->size(); i++)
    // 调用FullSystem::activateARawPixelPoint,激活一个点
        (*activatedPoints)[i] = activateARawPixelPoint((*toActivatePoints)[i], residuals);
}

// 激活一个像素点　RawPoint-->Point
PixelPoint * Tracking::activateARawPixelPoint(RawPixelPoint * aRawPoint, RawPixelPointTempResidual* residuals)
{

    size_t num = 0;
    for(Frame* frame : allKeyframes)
    {
        if(frame = aRawPoint->hostFrame) continue;

        residuals[num].tempRes_Energy=0;
        residuals[num].tempRes_NewEnergy=0;
        residuals[num].tempRes_State=ResState::Res_IN;//???
        residuals[num].tempRes_NewState=ResState::Res_OUT;//???
        residuals[num].targetFrame = frame;
        num++;
    }
    assert(num == allKeyframes.size()-1);

    PixelPoint* point = new PixelPoint(aRawPoint);
    // todo 8.14
}

void Tracking::makeAllIdxForBackend()
{
    EF->makeIndex();
}

float Tracking::runBackendOptimize(int IterNum)
{
    return EF->backGNOptimize(settings_OptimizeIterNum);
}

void Tracking::removeOutliers()
{
    // todo 8.14
}

void Tracking::flagPixelPointsAndMargOrDrop()
{
    // todo 8.14
}

void Tracking::makeRawPixelPointsForNowKeyframe()
{
    // todo 8.14
}


void Tracking::margFlagedKeyframe()
{
    // 后端边缘化掉一帧，目标函数中的总帧数减少
    EF->allBackKeyFramesNum--;
    // todo 8.14
}





//******************************************************************************************88
Vec3 Tracking::linearizeAll(bool b)
{
    double lastEnergy = 0;

    vector<Residual*> presidualsToRemove[settings_ThreadNum];
    for(int i = 0; i < settings_ThreadNum; i++)
    {
        for(Residual* presidual:presidualsToRemove[i])
            presidual = NULL;
        presidualsToRemove[i].clear();
    }

#ifdef MULTITHREAD
//    threadReduce.reduce();
#else

#endif


}

//void Tracking::accumulateAllResidual(bool b, vector<Residual *> *todo, int first, int last, Vec10 *stats, int tidx)
//{
//    for(int i = first; i < last; i++)
//    {
//        Residual* pres = allResiduals[i];
////        (*stats)[0] += pres
//    }
//}



}
