#ifndef PROJECTION_H
#define PROJECTION_H

#include "util/all_util_include.h"
#include "FullSystem/Settings.h"
#include "Frontend/Part.h"

using namespace std;
using namespace SLAMSystem;

namespace world3000 {

inline bool projectPoint(
        const float& u_old, const float& v_old,
        const float& idepth,
        const Mat33f& KRKinv, const Vec3f& Kt,
        float& u_new, float& v_new)
{
    Vec3f point_3D = KRKinv * Vec3f(u_old, v_old, 1) + Kt * idepth;   // 重建出像素点对应的三维点的空间坐标
    u_new = point_3D[0] / point_3D[2];
    v_new = point_3D[1] / point_3D[2];

    return u_new > 1.0f && v_new > 1.0f && u_new < width_thres && v_new < height_thres;

}
/**
 * @brief projectPoint
 * @param u_old                 3D点对应的像素坐标(u,v)(在hostframe中)
 * @param v_old
 * @param idepth                3D点的逆深度(相对于host帧)
 * @param u_old_drift
 * @param v_old_drift           像素坐标的偏移量，一般为0
 * @param cam
 * @param R                     R   从hostframe到targetframe的位姿变换的R
 * @param t                     t   从hostframe到targetframe的位姿变换的t
 * 输出量
 * @param normalizedFactor      坐标归一化因子 1/Z'
 * @param u_new
 * @param v_new                 该3D点在target帧中的像素坐标(u',v')
 * @param normalized_3Dcoord_old    3D点的相机归一化坐标(在生成该3D点的hostframe坐标系中) [X/Z, Y/Z, 1].transpose()
 * @param idepth_new            3D点的逆深度(相对于target帧)
 * @return 投影点的像素坐标满足阈值要求，返回true,否则返回false
 */

inline bool projectPoint(
        const float& u_old, const float& v_old,
        const float& idepth,
        const int& u_old_drift, const int& v_old_drift,
        const float& fx, const float& fy, const float& cx, const float& cy,
        const Mat33f& R, const Vec3f& t,
        float& normalizedFactor, float& X_new, float& Y_new,
        float& u_new, float& v_new, Vec3f& normalized_3Dcoord_old,
        float& idepth_new)
{
    normalized_3Dcoord_old = Vec3f((u_old + u_old_drift - cx) / fx,
                                   (v_old + v_old_drift - cy) / fy,
                                   1);

    Vec3f point_3Dcoord_new = R * normalized_3Dcoord_old + t * idepth;

    normalizedFactor = 1.0f / point_3Dcoord_new[2];
    X_new = point_3Dcoord_new[0];
    Y_new = point_3Dcoord_new[1];
    assert(normalizedFactor>0);
    idepth_new = normalizedFactor * idepth;

    u_new = normalizedFactor * point_3Dcoord_new[0] * fx + cx;
    v_new = normalizedFactor * point_3Dcoord_new[1] * fy + cy;

    return u_new > 1.0f && v_new > 1.0f && u_new < width_thres && v_new < height_thres;
}

}

#endif // PROJECTION_H
