#include "Frontend/Part.h"

namespace SLAMSystem {

//template<TrackMode T> double Residual::calcuResidualAndJacobianPart1(Camera *cam)      // 8.21 17:00 finished.


//void Residual::applyRes()       // finished 8.23


//void Residual::resetOOB()










//void FrametoFrame::PreCalculate(Frame *host, Frame *target, Camera *cam)    // done On 17.8.16




RawPixelPointState RawPixelPoint::trackRawPixelPoints()
{
    //todo, track，修改该raw点的状态
}

PixelPoint::PixelPoint(const RawPixelPoint * const rawPixelPoint)
{
    // todo, 17.8.14
}

double Residual::calcuResidualAndJacobianPart1(Camera *cam, TrackMode T1, CameraMode T2)      // 求解雅可比矩阵的核心函数
{
    // 1. 如果点帧残差的状态为OOB，那么新的状态也为OOB,返回
    if(resState == Res_OOB)
    {
        newResState = Res_OOB;
        return energy;
    }
    // 2. 如果投影点的像素坐标从hostframe变到targetframe后不在targetframe的区域内，那么新的残差状态设置为OOB
    //        if(T == Direct)
    //            FrametoFrame_Direct* h2t = hostFrame->thisToTargets[targetFrame->keyframeId];
    //        else if(T == Feature)
    //            FrametoFrame_Feature* h2t = hostFrame->thisToTargets[targetFrame->keyframeId];
    //        else { assert(false && "No such mode!"); }
    FrametoFrame* h2t = hostFrame->thisToTargets[targetFrame->keyframeId];

    const Vec3f* image = targetFrame->image;
    const Mat33f& KRKinv = h2t->KRKinv;
    const Vec3f& Kt = h2t->Kt;
    const Mat33f& R0_h2t = h2t->R0_h2t;
    const Vec3f& t0_h2t = h2t->t0_h2t;

    float normalizedFactor,X_new, Y_new, u_new, v_new, idepth_new;
    Vec3f normalized_3Dcoord_old;// 3D点的相机归一化坐标(在生成该3D点的hostframe坐标系中) [X/Z, Y/Z, 1].transpose()

    if(projectPoint(point->u,point->v,point->idepth,0,0,cam->backData->fx,cam->backData->fy,cam->backData->cx,cam->backData->cy,R0_h2t,t0_h2t,normalizedFactor,X_new,Y_new,u_new,v_new,normalized_3Dcoord_old,idepth_new)==false)
    {
        newResState = Res_OOB;
        return energy;
    }

    if(T1 == Direct)
    {
        const Vec2f& afflight = h2t->afflight;
        const float& b0 = h2t->b0;
        const float* const gray = point->gray;
        const float* const weight = point->weight;

        //3. 残差的雅可比求导
        {
            float d_u_new_d_idepth, d_v_new_d_idepth; // 像素坐标(u',v')对像素点逆深度的偏导
            Vec6f d_u_new_d_xi, d_v_new_d_xi;   // 像素坐标(u',v')对位姿R,t的偏导
            /* 分情况，是否优化相机内参，不优化就不用求偏导 */
            Vec4f d_u_new_d_cam, d_v_new_d_cam; // 像素坐标(u',v')对相机内参的偏导



            d_u_new_d_idepth = normalizedFactor * (t0_h2t[0] - t0_h2t[2]*X_new*normalizedFactor) * cam->backData->fx;
            d_v_new_d_idepth = normalizedFactor * (t0_h2t[1] - t0_h2t[2]*Y_new*normalizedFactor) * cam->backData->fy;

            d_u_new_d_xi[0] = idepth_new * cam->backData->fx;
            d_u_new_d_xi[1] = 0;
            d_u_new_d_xi[2] = -idepth_new * X_new * normalizedFactor * cam->backData->fx;
            d_u_new_d_xi[3] = -X_new * Y_new * normalizedFactor * normalizedFactor * cam->backData->fx;
            d_u_new_d_xi[4] = (1 + X_new * X_new * normalizedFactor * normalizedFactor) * cam->backData->fx;
            d_u_new_d_xi[5] = -Y_new * normalizedFactor * cam->backData->fx;
            d_v_new_d_xi[0] = 0;
            d_v_new_d_xi[1] = idepth_new * cam->backData->fy;
            d_v_new_d_xi[2] = -idepth_new * Y_new * normalizedFactor * cam->backData->fy;
            d_v_new_d_xi[3] = -(1 + X_new * X_new * normalizedFactor * normalizedFactor) * cam->backData->fy;
            d_v_new_d_xi[4] = X_new * Y_new * normalizedFactor * normalizedFactor * cam->backData->fy;
            d_v_new_d_xi[5] = X_new * normalizedFactor * cam->backData->fy;

            /* 是否优化相机内参，不优化就不用求偏导 */
            if(T2 == Optimize)
            {
                d_u_new_d_cam[0] = normalized_3Dcoord_old[0] * normalizedFactor * (R0_h2t(2,0) * X_new * normalizedFactor - R0_h2t(0,0)) + X_new * normalizedFactor;
                d_u_new_d_cam[1] = normalized_3Dcoord_old[1] * cam->backData->fx * normalizedFactor * (R0_h2t(2,0) * X_new * normalizedFactor -R0_h2t(0,1)) / cam->backData->fy;
                d_u_new_d_cam[2] = normalizedFactor * (R0_h2t(2,0) * X_new * normalizedFactor - R0_h2t(0,0));
                d_u_new_d_cam[3] = cam->backData->fx * normalizedFactor * (R0_h2t(2,1) * X_new * normalizedFactor - R0_h2t(0,1)) / cam->backData->fy;

                d_v_new_d_cam[0] = normalized_3Dcoord_old[0] * cam->backData->fy * normalizedFactor * (R0_h2t(2,0) * Y_new * normalizedFactor - R0_h2t(1,0)) / cam->backData->fx;
                d_v_new_d_cam[1] = normalized_3Dcoord_old[1] * normalizedFactor * (R0_h2t(2,0) * Y_new * normalizedFactor - R0_h2t(1,1)) + Y_new * normalizedFactor;
                d_v_new_d_cam[2] = cam->backData->fy * normalizedFactor * (R0_h2t(2,0) * Y_new * normalizedFactor - R0_h2t(1,0)) / cam->backData->fx;
                d_v_new_d_cam[3] = normalizedFactor * (R0_h2t(2,0) * Y_new * normalizedFactor - R0_h2t(1,1)) + 1;

                /* 是否优化相机内参，不优化就不用求偏导 */
                Jaco->J_p_d_cam[0] = d_u_new_d_cam;
                Jaco->J_p_d_cam[1] = d_v_new_d_cam;
            }
            else if(T2 == Fix)
            {
                /* 是否优化相机内参，不优化就不用求偏导 */
                Jaco->J_p_d_cam[0] = Vec4f::Zero();
                Jaco->J_p_d_cam[1] = Vec4f::Zero();
            }

            Jaco->J_p_d_xi[0] = d_u_new_d_xi;
            Jaco->J_p_d_xi[1] = d_v_new_d_xi;

            Jaco->J_p_d_idepth = Vec2f(d_u_new_d_idepth, d_v_new_d_idepth);


        }
        // 4. 对于直接法和特征点法分别计算光度误差和重投影误差
        /* 直接法中的光度误差 */
        {

            // 将点的residual　pattern 中的每个点从hostframe重建再投影到targetframe，得到在targetframe对应的像素坐标
            //  4.1 如果像素坐标不满足阈值，直接设置点帧残差状态为OOB，返回state_energy
            //  4.2 将像素坐标在targetframe的dI中插值，得到targetframe对应的像素坐标的灰度值以及像素梯度构成的3维向量grayAndGrayGradient
            //  4.3 targetframe中的对应的像素坐标的灰度值减去hostframe中的像素坐标的灰度值，得到了灰度残差residual,即直接法中的光度误差

            // 直接法残差的雅可比中的一些中间计算量
            float J_I_d_p_2_00 = 0, J_I_d_p_2_01 = 0, J_I_d_p_2_10 = 0, J_I_d_p_2_11 = 0;
            float J_res_d_aff_J_I_d_p_00 = 0, J_res_d_aff_J_I_d_p_01 = 0, J_res_d_aff_J_I_d_p_10 = 0, J_res_d_aff_J_I_d_p_11 = 0;
            float J_res_d_aff_2_00 = 0, J_res_d_aff_2_01 = 0, J_res_d_aff_2_10 = 0, J_res_d_aff_2_11 = 0;

            float w_J_I_d_p_2_sum = 0;

            Vec2f projectedPoint[res_pattern_pointNum];
            float energyLeft = 0;
            for(int idx = 0; idx < res_pattern_pointNum; idx++)
            {
                // 4.1
                float u_direct_new, v_direct_new;
                if(projectPoint(point->u+patternP[idx][0], point->v+patternP[idx][1], point->idepth, KRKinv, Kt, u_direct_new, v_direct_new) == false)
                { newResState = ResState::Res_OOB; return energy; }
                projectedPoint[idx][0] = u_direct_new; projectedPoint[idx][1] = v_direct_new;

                // 4.2 插值
                Vec3f grayAndGrayGradient = (world3000::getInterpolatedElement33(image, u_direct_new, v_direct_new, widthPyramid[0]));
                if(!std::isfinite((float)grayAndGrayGradient[0]))
                { newResState = ResState::Res_OOB; return energy; }

                // 4.3 计算光度误差，计算残差的雅可比的中间量．
                float residual_direct = grayAndGrayGradient[0] - (float)(afflight[0] * gray[idx] + afflight[1]);
                float drdA = (gray[idx]-b0);
                //settings_outlierTHSumComponent: 基于梯度的权重调整
                float w = sqrt(settings_outlierTHSumComponent / (settings_outlierTHSumComponent + grayAndGrayGradient.tail<2>().squaredNorm()));
                w = 0.5f*(w + weight[idx]);

                float huberKernelWeight = fabs(residual_direct) < settings_huberTH ? 1 : settings_huberTH / fabs (residual_direct);
                energyLeft += w * w * huberKernelWeight * residual_direct * residual_direct * (2-huberKernelWeight);
                // 4.3.2 计算残差的雅可比的中间量．
                {
                    if(huberKernelWeight < 1) huberKernelWeight = sqrt(huberKernelWeight);
                    huberKernelWeight = huberKernelWeight * w;

                    grayAndGrayGradient[1] *= huberKernelWeight;
                    grayAndGrayGradient[2] *= huberKernelWeight;

                    Jaco->res_direct[idx] = residual_direct * huberKernelWeight;

                    Jaco->J_I_d_p[0][idx] = grayAndGrayGradient[1];
                    Jaco->J_I_d_p[1][idx] = grayAndGrayGradient[2];     // 图片对像素点的偏导，即像素梯度

                    Jaco->J_res_d_aff[0][idx] = drdA * huberKernelWeight;
                    Jaco->J_res_d_aff[1][idx] = huberKernelWeight;                     // 残差对光度参数(a,b)的偏导 Jphoto

                    J_I_d_p_2_00 += pow(grayAndGrayGradient[1],2);
                    J_I_d_p_2_01 += grayAndGrayGradient[1] * grayAndGrayGradient[2];
                    J_I_d_p_2_10 = J_I_d_p_2_01;
                    J_I_d_p_2_11 += pow(grayAndGrayGradient[2],2);

                    J_res_d_aff_J_I_d_p_00 += drdA * huberKernelWeight * grayAndGrayGradient[1];
                    J_res_d_aff_J_I_d_p_01 += drdA * huberKernelWeight * grayAndGrayGradient[2];
                    J_res_d_aff_J_I_d_p_10 += huberKernelWeight * grayAndGrayGradient[1];
                    J_res_d_aff_J_I_d_p_11 += huberKernelWeight * grayAndGrayGradient[2];

                    J_res_d_aff_2_00 += drdA * drdA * huberKernelWeight * huberKernelWeight;
                    J_res_d_aff_2_01 += drdA * huberKernelWeight *huberKernelWeight;
                    J_res_d_aff_2_10 = J_res_d_aff_2_01;
                    J_res_d_aff_2_11 += huberKernelWeight * huberKernelWeight;

                    w_J_I_d_p_2_sum += huberKernelWeight*huberKernelWeight*(grayAndGrayGradient[1]*grayAndGrayGradient[1]+grayAndGrayGradient[2]*grayAndGrayGradient[2]);
                    if(settings_affineOptModeA < 0) Jaco->J_res_d_aff[0][idx] = 0;
                    if(settings_affineOptModeB < 0) Jaco->J_res_d_aff[1][idx] = 0;
                }


            }
            // 4.3.3 更新残差的雅可比的相应量
            Jaco->J_I_d_p_2(0,0) = J_I_d_p_2_00;
            Jaco->J_I_d_p_2(0,1) = J_I_d_p_2_01;
            Jaco->J_I_d_p_2(1,0) = J_I_d_p_2_10;
            Jaco->J_I_d_p_2(1,1) = J_I_d_p_2_11;

            Jaco->J_res_d_aff_J_I_d_p(0,0) = J_res_d_aff_J_I_d_p_00;
            Jaco->J_res_d_aff_J_I_d_p(0,1) = J_res_d_aff_J_I_d_p_01;
            Jaco->J_res_d_aff_J_I_d_p(1,0) = J_res_d_aff_J_I_d_p_10;
            Jaco->J_res_d_aff_J_I_d_p(1,1) = J_res_d_aff_J_I_d_p_11;

            Jaco->J_res_d_aff_2(0,0) = J_res_d_aff_2_00;
            Jaco->J_res_d_aff_2(0,1) = J_res_d_aff_2_01;
            Jaco->J_res_d_aff_2(1,0) = J_res_d_aff_2_10;
            Jaco->J_res_d_aff_2(1,1) = J_res_d_aff_2_11;

            // energyleft大于阈值，则newstate为outlier,frameEnergyTH为光度相关参数
            if(energyLeft > std::max<float>(hostFrame->frameEnergyTH, targetFrame->frameEnergyTH) || w_J_I_d_p_2_sum < 2)
            {
                //        energyLeft = std::max<float>(hostFrame->frameEnergyTH, m_target->frameEnergyTH);
                newResState = ResState::Res_OUT;
            }
            // 否则newstate为IN
            else
            {
                newResState = ResState::Res_IN;
            }

            newEnergy = energyLeft;
            newEnergyForTH = energyLeft;

            return energyLeft;
        }

    }
    else if(T1 == Feature)
    {
//        Vec2f afflight = Vec2f::Zero();
//        float b0 = 0;

        //3. 残差的雅可比求导
        {
            float d_u_new_d_idepth, d_v_new_d_idepth; // 像素坐标(u',v')对像素点逆深度的偏导
            Vec6f d_u_new_d_xi, d_v_new_d_xi;   // 像素坐标(u',v')对位姿R,t的偏导
            /* 分情况，是否优化相机内参，不优化就不用求偏导 */
            Vec4f d_u_new_d_cam, d_v_new_d_cam; // 像素坐标(u',v')对相机内参的偏导



            d_u_new_d_idepth = normalizedFactor * (t0_h2t[0] - t0_h2t[2]*X_new*normalizedFactor) * cam->backData->fx;
            d_v_new_d_idepth = normalizedFactor * (t0_h2t[1] - t0_h2t[2]*Y_new*normalizedFactor) * cam->backData->fy;

            d_u_new_d_xi[0] = idepth_new * cam->backData->fx;
            d_u_new_d_xi[1] = 0;
            d_u_new_d_xi[2] = -idepth_new * X_new * normalizedFactor * cam->backData->fx;
            d_u_new_d_xi[3] = -X_new * Y_new * normalizedFactor * normalizedFactor * cam->backData->fx;
            d_u_new_d_xi[4] = (1 + X_new * X_new * normalizedFactor * normalizedFactor) * cam->backData->fx;
            d_u_new_d_xi[5] = -Y_new * normalizedFactor * cam->backData->fx;
            d_v_new_d_xi[0] = 0;
            d_v_new_d_xi[1] = idepth_new * cam->backData->fy;
            d_v_new_d_xi[2] = -idepth_new * Y_new * normalizedFactor * cam->backData->fy;
            d_v_new_d_xi[3] = -(1 + X_new * X_new * normalizedFactor * normalizedFactor) * cam->backData->fy;
            d_v_new_d_xi[4] = X_new * Y_new * normalizedFactor * normalizedFactor * cam->backData->fy;
            d_v_new_d_xi[5] = X_new * normalizedFactor * cam->backData->fy;

            /* 是否优化相机内参，不优化就不用求偏导 */
            if(T2 == Optimize)
            {
                d_u_new_d_cam[0] = normalized_3Dcoord_old[0] * normalizedFactor * (R0_h2t(2,0) * X_new * normalizedFactor - R0_h2t(0,0)) + X_new * normalizedFactor;
                d_u_new_d_cam[1] = normalized_3Dcoord_old[1] * cam->backData->fx * normalizedFactor * (R0_h2t(2,0) * X_new * normalizedFactor -R0_h2t(0,1)) / cam->backData->fy;
                d_u_new_d_cam[2] = normalizedFactor * (R0_h2t(2,0) * X_new * normalizedFactor - R0_h2t(0,0));
                d_u_new_d_cam[3] = cam->backData->fx * normalizedFactor * (R0_h2t(2,1) * X_new * normalizedFactor - R0_h2t(0,1)) / cam->backData->fy;

                d_v_new_d_cam[0] = normalized_3Dcoord_old[0] * cam->backData->fy * normalizedFactor * (R0_h2t(2,0) * Y_new * normalizedFactor - R0_h2t(1,0)) / cam->backData->fx;
                d_v_new_d_cam[1] = normalized_3Dcoord_old[1] * normalizedFactor * (R0_h2t(2,0) * Y_new * normalizedFactor - R0_h2t(1,1)) + Y_new * normalizedFactor;
                d_v_new_d_cam[2] = cam->backData->fy * normalizedFactor * (R0_h2t(2,0) * Y_new * normalizedFactor - R0_h2t(1,0)) / cam->backData->fx;
                d_v_new_d_cam[3] = normalizedFactor * (R0_h2t(2,0) * Y_new * normalizedFactor - R0_h2t(1,1)) + 1;

                /* 是否优化相机内参，不优化就不用求偏导 */
                Jaco->J_p_d_cam[0] = d_u_new_d_cam;
                Jaco->J_p_d_cam[1] = d_v_new_d_cam;
            }
            else if(T2 == Fix)
            {
                /* 是否优化相机内参，不优化就不用求偏导 */
                Jaco->J_p_d_cam[0] = Vec4f::Zero();
                Jaco->J_p_d_cam[1] = Vec4f::Zero();
            }

            Jaco->J_p_d_xi[0] = d_u_new_d_xi;
            Jaco->J_p_d_xi[1] = d_v_new_d_xi;

            Jaco->J_p_d_idepth = Vec2f(d_u_new_d_idepth, d_v_new_d_idepth);

        }
        // 4. 对于直接法和特征点法分别计算光度误差和重投影误差
        /* 特征点法的重投影误差 */
        {

            //得到预测值predictPoint，即点从hostframe变换到targetframe应有的像素坐标及逆深度.
            Vec3f predictPoint = Vec3f(u_new, v_new, idepth_new);
            //观测值observePoint, point在targerframe中的像素坐标观测值应由tracking前端给出
            auto idx = find(point->observeFrames.begin(), point->observeFrames.end(), targetFrame);
            assert(idx != point->observeFrames.end());
            int index = distance(point->observeFrames.begin(), idx);
            Vec3f observePoint = point->observePixelPoints[index];
            //观测值减去预测值，得到重投影误差，即特征点法中的残差
            float residual_feature_Norm = sqrt(pow((observePoint[0] - predictPoint[0]),2) + pow((observePoint[1] - predictPoint[1]),2));
            float residual_feature_u = observePoint[0] - predictPoint[0];
            float residual_feature_v = observePoint[0] - predictPoint[1];

            float residual_feature_uvmax = max(residual_feature_u, residual_feature_v);
            float huberKernelWeight = fabs(residual_feature_uvmax) < settings_huberTH ? 1 : settings_huberTH / fabs (residual_feature_uvmax);

            Jaco->res_feature_uv = Vec2f(residual_feature_u * huberKernelWeight, residual_feature_v * huberKernelWeight);
            Jaco->res_feature_Norm = Jaco->res_feature_uv.norm();

            float energyLeft = huberKernelWeight * residual_feature_Norm * residual_feature_Norm * (2-huberKernelWeight);


            if(abs(observePoint[0] - predictPoint[0]) > settings_pixelCoordTH || abs(observePoint[1] - predictPoint[1]) > settings_pixelCoordTH)
            {
                newResState = ResState::Res_OUT;
            }
            else
            {
                newResState = ResState::Res_IN;
            }
            newEnergy = energyLeft;
            return energyLeft;
        }

    }
    else assert(false && "No such mode!");
    //        const Vec2f& afflight = h2t->afflight;
}

void Residual::applyRes()
{

    if(newResState == ResState::Res_IN && resState != ResState::Res_OOB)
    {
        backData->isActiveAndIsGoodNEW = true;
        backData->getJacobianPart1AndCalculateJacobianPart2();

        resState = newResState;
        energy = newEnergy;
    }
    else
        backData->isActiveAndIsGoodNEW = false;

}

void Residual::resetOOB()
{
    energy = newEnergy = 0;
    resState = Res_IN;
    newResState = Res_OUT;
}



FrametoFrame::FrametoFrame(Frame *host, Frame *target, Camera *cam, TrackMode T)
    : hostframe(host), targetframe(target), trackMode(T)
{
    SE3 T0_h2t = target->Tcw0 * host->Tcw0.inverse();  // 公式(6)
    R0_h2t = T0_h2t.rotationMatrix().cast<float>();
    t0_h2t = T0_h2t.translation().cast<float>();

    SE3 T_h2t = target->Tcw * host->Tcw.inverse();
    R_h2t = T_h2t.rotationMatrix().cast<float>();
    t_h2t = T_h2t.translation().cast<float>();

    Mat33f K = Mat33f::Identity();
    K(0,0) = cam->backData->fx;
    K(1,1) = cam->backData->fy;
    K(0,2) = cam->backData->cx;
    K(1,2) = cam->backData->cy;


    //todo
}

void FrametoFrame::PreCalculate(Frame *host, Frame *target, Camera *cam)
{
    hostframe = host; targetframe = target;

    SE3 T0_h2t = target->Tcw0 * host->Tcw0.inverse();
    R0_h2t = T0_h2t.rotationMatrix().cast<float>();
    t0_h2t = T0_h2t.translation().cast<float>();

    SE3 T_h2t = target->Tcw * host->Tcw.inverse();
    R_h2t = T_h2t.rotationMatrix().cast<float>();
    t_h2t = T_h2t.translation().cast<float>();

    Mat33f K = Mat33f::Identity();
    K(0,0) = cam->backData->fx;
    K(1,1) = cam->backData->fy;
    K(0,2) = cam->backData->cx;
    K(1,2) = cam->backData->cy;

    KRKinv = K * R_h2t * K.inverse();
    RKinv = R_h2t * K.inverse();
    Kt = K * t_h2t;

    if(trackMode == Direct)
    {
        // 直接法参数(a,b),在特征点法中不需要
        afflight = AffLight::getHostToTargetAffineParams(hostframe->exposureTime, targetframe->exposureTime,
                                                         hostframe->backData->getAffineParams(), targetframe->backData->getAffineParams());
        b0 = hostframe->backData->getAffineParams0().b;
    } else if(trackMode == Feature)
    {
        afflight = Vec2f::Zero();
        b0 = 0;
    } else { assert(false && "No such mode!");}

}





}
