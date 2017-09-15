#include "Backend/BackPart.h"

namespace SLAMSystem {

BackFrame::BackFrame(Frame *frame)
    : frontData(frame), Tcw0(frame->Tcw0), Tcw(frame->Tcw), Twc(frame->Twc),
      state_x0(Vec8::Constant(NAN)), state_x(Vec8::Constant(NAN)), state_x_backup(Vec8::Constant(NAN)),
      step(Vec8::Constant(0)), step_backup(Vec8::Constant(0)), delta_0ToNow(Vec8::Constant(0)), keyframeId(frame->keyframeId)
{

}

void BackFrame::init()
{

}

void BackFrame::setTcw0(const SE3 &Tcw0)
{
    this->Tcw0 = Tcw0;
}

void BackFrame::setX0(const SE3 &Tcw0, const Vec8 &x0)
{
    state_x0 = x0;
    //    this->Tcw0 = Tcw0;
}

void BackFrame::setX0(const SE3 &Tcw0, const AffLight &aff_g2l)
{
    setTcw0(Tcw0);
    state_x0 = Vec8::Zero();
    state_x0[6] = aff_g2l.a;
    state_x0[7] = aff_g2l.b;

}

void BackFrame::set_state_x0(const Vec8 &state_x0)
{
    assert(state_x0.head<6>().squaredNorm() < 1e-20);
    this->state_x0 = state_x0;

    for(int i = 0; i <6; i++)
    {
        Vec6 epsilon = Vec6::Constant(0);
        epsilon[i] = 1e-3;
        SE3 EepsP = Sophus::SE3::exp(epsilon);
        SE3 EepsM = Sophus::SE3::exp(-epsilon);
        SE3 w2c_leftEps_P_x0 = Tcw0 * EepsP * Tcw0.inverse();
        SE3 w2c_leftEps_M_x0 = Tcw0 * EepsM * Tcw0.inverse();

        nullspace_pose.col(i) = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log()) / (2e-3);
    }

    // scale change
    SE3 w2c_leftEps_P_x0 = Tcw0;
    w2c_leftEps_P_x0.translation() *= 1.00001;
    w2c_leftEps_P_x0 = w2c_leftEps_P_x0 * Tcw0.inverse();
    SE3 w2c_leftEps_M_x0 = (Tcw0);
    w2c_leftEps_M_x0.translation() /= 1.00001;
    w2c_leftEps_M_x0 = w2c_leftEps_M_x0 * Tcw0.inverse();
    nullspace_scale = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log())/(2e-3);


    nullspace_affine.setZero();
    nullspace_affine.topLeftCorner<2,1>()  = Vec2(1,0);
    assert(frontData->exposureTime > 0);
    nullspace_affine.topRightCorner<2,1>() = Vec2(0, expf(state_x0[6])*frontData->exposureTime);
}





void BackFrame::setX(const Vec8 &x)
{
    state_x = x;
    //    Tcw =  SE3::exp(state_xNow.segment<6>(0)) * Tcw0;
}
// return (a,b)
AffLight BackFrame::getAffineParams()
{
    return AffLight(state_x[6],state_x[7]);       // 此处的(a,b)都是单位统一后用来程序计算的
}

AffLight BackFrame::getAffineParams0()
{
    return AffLight(state_x0[6],state_x0[7]);
}

void BackResidual::getJacobianPart1AndCalculateJacobianPart2()      // to be modified 8.23
{
    std::swap<ResidualJacobian*> (J, frontData->Jaco);  // getJacobianPart1

/* 直接法，有像素梯度J_I_d_p及图片对仿射亮度参数的偏导J_p_d_idepth, 需要计算这些量到增量方程的一些中间量*/
    {
        Vec2f J_I_d_p_2_J_p_d_idepth = J->J_I_d_p_2 * J->J_p_d_idepth;
        for(int i = 0; i < 6; i++)
        {
            JpJdF[i] = J->J_p_d_xi[0][i] * J_I_d_p_2_J_p_d_idepth[0] + J->J_p_d_xi[1][i] * J_I_d_p_2_J_p_d_idepth[1];
        }
        JpJdF.segment<2>(6) = J->J_res_d_aff_J_I_d_p * J->J_p_d_idepth;
    }
/* 特征点法，没有像素梯度以及仿射亮度参数的偏导，中间量JpJdF按如下方式计算 */
    /* 不确定，有待后续的修改 (09.11认为是正确的 )*/
    {
        for(int i = 0; i < 6; i++)
        {
            JpJdF[i] = J->J_p_d_xi[0][i] * J->J_p_d_idepth[0] + J->J_p_d_xi[1][i] * J->J_p_d_idepth[1];
        }
        JpJdF.segment<2>(6) = J->J_p_d_idepth;
    }

}

BackPixelPoint::BackPixelPoint(PixelPoint *point)
{

}


//BackFrame::BackFrame(Frame *frame)
//    : frontData(frame), state_x0(setPrior())
//{
//    prior = 0;
//}

//Vec8 BackFrame::setPrior()
//{
//    Vec8 piror = Vec8::Zero();
//    return piror;
//}



}
