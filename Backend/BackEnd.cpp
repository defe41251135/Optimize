#include "Backend/BackEnd.h"

namespace SLAMSystem {

EnergyFunction::EnergyFunction()
{
    // todo 8.14
}

EnergyFunction::EnergyFunction(TrackMode T1, CameraMode T2)
    : trackMode(T1), cameraMode(T2)
{
    acc_top_Activate = new AccumulatedTopHessian(T1, T2);
    acc_top_Linearized = new AccumulatedTopHessian(T1, T2);
    acc_bot_Marginalize = new AccumulatedSCHessian(T1, T2);
}

void EnergyFunction::setTracking(Tracking *pTracking)
{
    m_pTracking = pTracking;
}

EnergyFunction::~EnergyFunction()
{
    // todo 8.14
}
/**
 * @brief 为后端的帧，残差进行统一编号
 */
void EnergyFunction::makeIndex()
{
    // 为后端的每一个关键帧标号
    for(size_t index = 0; index < allBackKeyFrames.size(); index++)
    {
        allBackKeyFrames[index]->keyframeId = index;
    }

    // 为后端的每个残差项标号
    allBackPixelPoints.clear();
    for(BackFrame* backFrame : allBackKeyFrames)
        for(BackPixelPoint* backPoint : backFrame->backPixelPoints)
        {
            allBackPixelPoints.push_back(backPoint);
            for(BackResidual* backResidual : backPoint->backResiduals)
            {
                backResidual->backHostIndex = backResidual->backHostFrame->keyframeId;
                backResidual->backTargetIndex = backResidual->backTargetFrame->keyframeId;
            }
        }
}


























}
