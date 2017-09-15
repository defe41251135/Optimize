#ifndef ACCUMULATEDHESSIAN_H
#define ACCUMULATEDHESSIAN_H

#include "util/all_util_include.h"
#include "FullSystem/Settings.h"
#include "Backend/BackEnd.h"
#include "Backend/BackPart.h"

using namespace std;
using namespace world3000;

namespace SLAMSystem
{
/* 这部分求解的是刚激活的点及线性化了的点的残差及它们对应的host,target帧的H和b */


// C, xi, a,b res   4, 6, 2, 1                  // 考虑多种情形，此处应建立继承类，根据具体情况进行子类构建
/* 是否优化相机内参, C4 */
/* 直接法OR特征点法，(a,b) */
/* 得到局部的H,b */
class AccumulateMatrixXXf
{
public:
    AccumulateMatrixXXf()
    {

    }
    virtual ~AccumulateMatrixXXf()
    {

    }
    virtual void init() = 0;
    virtual void upDate(const float* const C4_u, const float* const xi6_u,
                        const float* const C4_v, const float* const xi6_v,
                        const Mat22f& J_I_d_p_2,
                        const Mat22f& J_res_d_aff_J_I_d_p, const Vec2f& J_I_d_p_J_r,
                        const Mat22f& J_res_d_aff_2, const Vec2f& J_res_d_aff_J_r, float res2) {}
    virtual void upDate(const float* const xi6_u, const float* const xi6_v,
                        const Mat22f& J_I_d_p_2,                                                // topLeft
                        const Mat22f& J_res_d_aff_J_I_d_p, const Vec2f& J_I_d_p_J_r,            // topRight
                        const Mat22f& J_res_d_aff_2, const Vec2f& J_res_d_aff_J_r, float res2) {}   // botRight
    virtual void upDate(const float* const C4_u, const float* const xi6_u,
                        const float* const C4_v, const float* const xi6_v,
                        const Vec2f& feature_res) {}
    virtual void upDate(const float* const xi6_u, const float* const xi6_v,
                        const Vec2f& feature_res) {}
    virtual void done() = 0;
public:
    size_t hNum;
    MatXXf H;
protected:
    int hNumL0, hNumL1, hNumL2;

};

/* 直接法且优化相机内参   4 + 6 + 2 + 1 = 13 */
class AccumulateMatrix1313f : public AccumulateMatrixXXf
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    void init()
    {
        memset(topLeftDataL0, 0, sizeof(topLeftDataL0));
        memset(topLeftDataL1, 0, sizeof(topLeftDataL1));
        memset(topLeftDataL2, 0, sizeof(topLeftDataL2));

        memset(topRightDataL0, 0, sizeof(topRightDataL0));
        memset(topRightDataL1, 0, sizeof(topRightDataL1));
        memset(topRightDataL2, 0, sizeof(topRightDataL2));

        memset(botRightDataL0, 0, sizeof(botRightDataL0));
        memset(botRightDataL1, 0, sizeof(botRightDataL1));
        memset(botRightDataL2, 0, sizeof(botRightDataL2));

        hNum = hNumL0 = hNumL1 = hNumL2 = 0;
    }
    // 分层级更新H右上半侧每一个元素
    void upDate(const float* const C4_u, const float* const xi6_u,
                const float* const C4_v, const float* const xi6_v,
                const Mat22f& J_I_d_p_2,
                const Mat22f& J_res_d_aff_J_I_d_p, const Vec2f& J_I_d_p_J_r,
                const Mat22f& J_res_d_aff_2, const Vec2f& J_res_d_aff_J_r, float res2)
    {
        float a = J_I_d_p_2(0, 0);
        float b = J_I_d_p_2(0, 1);
        float c = J_I_d_p_2(1, 1);

        topLeftDataL0[0] += a * C4_u[0] * C4_u[0] + b * (C4_u[0] * C4_v[0] + C4_v[0] * C4_u[0]) + c * C4_v[0] * C4_v[0];
        topLeftDataL0[1] += a * C4_u[1] * C4_u[0] + b * (C4_u[1] * C4_v[0] + C4_v[1] * C4_u[0]) + c * C4_v[1] * C4_v[0];
        topLeftDataL0[2] += a * C4_u[2] * C4_u[0] + b * (C4_u[2] * C4_v[0] + C4_v[2] * C4_u[0]) + c * C4_v[2] * C4_v[0];
        topLeftDataL0[3] += a * C4_u[3] * C4_u[0] + b * (C4_u[3] * C4_v[0] + C4_v[3] * C4_u[0]) + c * C4_v[3] * C4_v[0];
        topLeftDataL0[4] += a * xi6_u[0] * C4_u[0] + b * (xi6_u[0] * C4_v[0] + xi6_v[0] * C4_u[0]) + c * xi6_v[0] * C4_v[0];
        topLeftDataL0[5] += a * xi6_u[1] * C4_u[0] + b * (xi6_u[1] * C4_v[0] + xi6_v[1] * C4_u[0]) + c * xi6_v[1] * C4_v[0];
        topLeftDataL0[6] += a * xi6_u[2] * C4_u[0] + b * (xi6_u[2] * C4_v[0] + xi6_v[2] * C4_u[0]) + c * xi6_v[2] * C4_v[0];
        topLeftDataL0[7] += a * xi6_u[3] * C4_u[0] + c * xi6_v[3] * C4_v[0] +  b * (xi6_u[3] * C4_v[0] + xi6_v[3] * C4_u[0]);
        topLeftDataL0[8] += a * xi6_u[4] * C4_u[0] + c * xi6_v[4] * C4_v[0] +  b * (xi6_u[4] * C4_v[0] + xi6_v[4] * C4_u[0]);
        topLeftDataL0[9] += a * xi6_u[5] * C4_u[0] + c * xi6_v[5] * C4_v[0] +  b * (xi6_u[5] * C4_v[0] + xi6_v[5] * C4_u[0]);

        topLeftDataL0[10] += a * C4_u[1] * C4_u[1] + c * C4_v[1] * C4_v[1] +  b * (C4_u[1] * C4_v[1] + C4_v[1] * C4_u[1]);
        topLeftDataL0[11] += a * C4_u[2] * C4_u[1] + c * C4_v[2] * C4_v[1] +  b * (C4_u[2] * C4_v[1] + C4_v[2] * C4_u[1]);
        topLeftDataL0[12] += a * C4_u[3] * C4_u[1] + c * C4_v[3] * C4_v[1] +  b * (C4_u[3] * C4_v[1] + C4_v[3] * C4_u[1]);
        topLeftDataL0[13] += a * xi6_u[0] * C4_u[1] + c * xi6_v[0] * C4_v[1] +  b * (xi6_u[0] * C4_v[1] + xi6_v[0] * C4_u[1]);
        topLeftDataL0[14] += a * xi6_u[1] * C4_u[1] + c * xi6_v[1] * C4_v[1] +  b * (xi6_u[1] * C4_v[1] + xi6_v[1] * C4_u[1]);
        topLeftDataL0[15] += a * xi6_u[2] * C4_u[1] + c * xi6_v[2] * C4_v[1] +  b * (xi6_u[2] * C4_v[1] + xi6_v[2] * C4_u[1]);
        topLeftDataL0[16] += a * xi6_u[3] * C4_u[1] + c * xi6_v[3] * C4_v[1] +  b * (xi6_u[3] * C4_v[1] + xi6_v[3] * C4_u[1]);
        topLeftDataL0[17] += a * xi6_u[4] * C4_u[1] + c * xi6_v[4] * C4_v[1] +  b * (xi6_u[4] * C4_v[1] + xi6_v[4] * C4_u[1]);
        topLeftDataL0[18] += a * xi6_u[5] * C4_u[1] + c * xi6_v[5] * C4_v[1] +  b * (xi6_u[5] * C4_v[1] + xi6_v[5] * C4_u[1]);

        topLeftDataL0[19] += a * C4_u[2] * C4_u[2] + c * C4_v[2] * C4_v[2] +  b * (C4_u[2] * C4_v[2] + C4_v[2] * C4_u[2]);
        topLeftDataL0[20] += a * C4_u[3] * C4_u[2] + c * C4_v[3] * C4_v[2] +  b * (C4_u[3] * C4_v[2] + C4_v[3] * C4_u[2]);
        topLeftDataL0[21] += a * xi6_u[0] * C4_u[2] + c * xi6_v[0] * C4_v[2] +  b * (xi6_u[0] * C4_v[2] + xi6_v[0] * C4_u[2]);
        topLeftDataL0[22] += a * xi6_u[1] * C4_u[2] + c * xi6_v[1] * C4_v[2] +  b * (xi6_u[1] * C4_v[2] + xi6_v[1] * C4_u[2]);
        topLeftDataL0[23] += a * xi6_u[2] * C4_u[2] + c * xi6_v[2] * C4_v[2] +  b * (xi6_u[2] * C4_v[2] + xi6_v[2] * C4_u[2]);
        topLeftDataL0[24] += a * xi6_u[3] * C4_u[2] + c * xi6_v[3] * C4_v[2] +  b * (xi6_u[3] * C4_v[2] + xi6_v[3] * C4_u[2]);
        topLeftDataL0[25] += a * xi6_u[4] * C4_u[2] + c * xi6_v[4] * C4_v[2] +  b * (xi6_u[4] * C4_v[2] + xi6_v[4] * C4_u[2]);
        topLeftDataL0[26] += a * xi6_u[5] * C4_u[2] + c * xi6_v[5] * C4_v[2] +  b * (xi6_u[5] * C4_v[2] + xi6_v[5] * C4_u[2]);

        topLeftDataL0[27] += a * C4_u[3] * C4_u[3] + c * C4_v[3] * C4_v[3] +  b * (C4_u[3] * C4_v[3] + C4_v[3] * C4_u[3]);
        topLeftDataL0[28] += a * xi6_u[0] * C4_u[3] + c * xi6_v[0] * C4_v[3] +  b * (xi6_u[0] * C4_v[3] + xi6_v[0] * C4_u[3]);
        topLeftDataL0[29] += a * xi6_u[1] * C4_u[3] + c * xi6_v[1] * C4_v[3] +  b * (xi6_u[1] * C4_v[3] + xi6_v[1] * C4_u[3]);
        topLeftDataL0[30] += a * xi6_u[2] * C4_u[3] + c * xi6_v[2] * C4_v[3] +  b * (xi6_u[2] * C4_v[3] + xi6_v[2] * C4_u[3]);
        topLeftDataL0[31] += a * xi6_u[3] * C4_u[3] + c * xi6_v[3] * C4_v[3] +  b * (xi6_u[3] * C4_v[3] + xi6_v[3] * C4_u[3]);
        topLeftDataL0[32] += a * xi6_u[4] * C4_u[3] + c * xi6_v[4] * C4_v[3] +  b * (xi6_u[4] * C4_v[3] + xi6_v[4] * C4_u[3]);
        topLeftDataL0[33] += a * xi6_u[5] * C4_u[3] + c * xi6_v[5] * C4_v[3] +  b * (xi6_u[5] * C4_v[3] + xi6_v[5] * C4_u[3]);

        topLeftDataL0[34] += a * xi6_u[0] * xi6_u[0] + c * xi6_v[0] * xi6_v[0] +  b * (xi6_u[0] * xi6_v[0] + xi6_v[0] * xi6_u[0]);
        topLeftDataL0[35] += a * xi6_u[1] * xi6_u[0] + c * xi6_v[1] * xi6_v[0] +  b * (xi6_u[1] * xi6_v[0] + xi6_v[1] * xi6_u[0]);
        topLeftDataL0[36] += a * xi6_u[2] * xi6_u[0] + c * xi6_v[2] * xi6_v[0] +  b * (xi6_u[2] * xi6_v[0] + xi6_v[2] * xi6_u[0]);
        topLeftDataL0[37] += a * xi6_u[3] * xi6_u[0] + c * xi6_v[3] * xi6_v[0] +  b * (xi6_u[3] * xi6_v[0] + xi6_v[3] * xi6_u[0]);
        topLeftDataL0[38] += a * xi6_u[4] * xi6_u[0] + c * xi6_v[4] * xi6_v[0] +  b * (xi6_u[4] * xi6_v[0] + xi6_v[4] * xi6_u[0]);
        topLeftDataL0[39] += a * xi6_u[5] * xi6_u[0] + c * xi6_v[5] * xi6_v[0] +  b * (xi6_u[5] * xi6_v[0] + xi6_v[5] * xi6_u[0]);



        topLeftDataL0[40] += a * xi6_u[1] * xi6_u[1] + c * xi6_v[1] * xi6_v[1] +  b * (xi6_u[1] * xi6_v[1] + xi6_v[1] * xi6_u[1]);
        topLeftDataL0[41] += a * xi6_u[2] * xi6_u[1] + c * xi6_v[2] * xi6_v[1] +  b * (xi6_u[2] * xi6_v[1] + xi6_v[2] * xi6_u[1]);
        topLeftDataL0[42] += a * xi6_u[3] * xi6_u[1] + c * xi6_v[3] * xi6_v[1] +  b * (xi6_u[3] * xi6_v[1] + xi6_v[3] * xi6_u[1]);
        topLeftDataL0[43] += a * xi6_u[4] * xi6_u[1] + c * xi6_v[4] * xi6_v[1] +  b * (xi6_u[4] * xi6_v[1] + xi6_v[4] * xi6_u[1]);
        topLeftDataL0[44] += a * xi6_u[5] * xi6_u[1] + c * xi6_v[5] * xi6_v[1] +  b * (xi6_u[5] * xi6_v[1] + xi6_v[5] * xi6_u[1]);


        topLeftDataL0[45] += a * xi6_u[2] * xi6_u[2] + c * xi6_v[2] * xi6_v[2] +  b * (xi6_u[2] * xi6_v[2] + xi6_v[2] * xi6_u[2]);
        topLeftDataL0[46] += a * xi6_u[3] * xi6_u[2] + c * xi6_v[3] * xi6_v[2] +  b * (xi6_u[3] * xi6_v[2] + xi6_v[3] * xi6_u[2]);
        topLeftDataL0[47] += a * xi6_u[4] * xi6_u[2] + c * xi6_v[4] * xi6_v[2] +  b * (xi6_u[4] * xi6_v[2] + xi6_v[4] * xi6_u[2]);
        topLeftDataL0[48] += a * xi6_u[5] * xi6_u[2] + c * xi6_v[5] * xi6_v[2] +  b * (xi6_u[5] * xi6_v[2] + xi6_v[5] * xi6_u[2]);


        topLeftDataL0[49] += a * xi6_u[3] * xi6_u[3] + c * xi6_v[3] * xi6_v[3] +  b * (xi6_u[3] * xi6_v[3] + xi6_v[3] * xi6_u[3]);
        topLeftDataL0[50] += a * xi6_u[4] * xi6_u[3] + c * xi6_v[4] * xi6_v[3] +  b * (xi6_u[4] * xi6_v[3] + xi6_v[4] * xi6_u[3]);
        topLeftDataL0[51] += a * xi6_u[5] * xi6_u[3] + c * xi6_v[5] * xi6_v[3] +  b * (xi6_u[5] * xi6_v[3] + xi6_v[5] * xi6_u[3]);


        topLeftDataL0[52] += a * xi6_u[4] * xi6_u[4] + c * xi6_v[4] * xi6_v[4] +  b * (xi6_u[4] * xi6_v[4] + xi6_v[4] * xi6_u[4]);
        topLeftDataL0[53] += a * xi6_u[5] * xi6_u[4] + c * xi6_v[5] * xi6_v[4] +  b * (xi6_u[5] * xi6_v[4] + xi6_v[5] * xi6_u[4]);

        topLeftDataL0[54] += a * xi6_u[5] * xi6_u[5] + c * xi6_v[5] * xi6_v[5] +  b * (xi6_u[5] * xi6_v[5] + xi6_v[5] * xi6_u[5]);

        float d = J_res_d_aff_J_I_d_p(0,0);
        float e = J_res_d_aff_J_I_d_p(0,1);
        float f = J_res_d_aff_J_I_d_p(1,0);
        float g = J_res_d_aff_J_I_d_p(1,1);

        float h = J_I_d_p_J_r[0];
        float i = J_I_d_p_J_r[1];

        topRightDataL0[0] += d * C4_u[0] + e * C4_v[0];
        topRightDataL0[1] += f * C4_u[0] + g * C4_v[0];
        topRightDataL0[2] += h * C4_u[0] + i * C4_v[0];

        topRightDataL0[3] += d * C4_u[1] + e * C4_v[1];
        topRightDataL0[4] += f * C4_u[1] + g * C4_v[1];
        topRightDataL0[5] += h * C4_u[1] + i * C4_v[1];

        topRightDataL0[6] += d * C4_u[2] + e * C4_v[2];
        topRightDataL0[7] += f * C4_u[2] + g * C4_v[2];
        topRightDataL0[8] += h * C4_u[2] + i * C4_v[2];

        topRightDataL0[9]  += d * C4_u[3] + e * C4_v[3];
        topRightDataL0[10] += f * C4_u[3] + g * C4_v[3];
        topRightDataL0[11] += h * C4_u[3] + i * C4_v[3];

        topRightDataL0[12] += d * xi6_u[0] + e * xi6_v[0];
        topRightDataL0[13] += f * xi6_u[0] + g * xi6_v[0];
        topRightDataL0[14] += h * xi6_u[0] + i * xi6_v[0];

        topRightDataL0[15] += d * xi6_u[1] + e * xi6_v[1];
        topRightDataL0[16] += f * xi6_u[1] + g * xi6_v[1];
        topRightDataL0[17] += h * xi6_u[1] + i * xi6_v[1];

        topRightDataL0[18] += d * xi6_u[2] + e * xi6_v[2];
        topRightDataL0[19] += f * xi6_u[2] + g * xi6_v[2];
        topRightDataL0[20] += h * xi6_u[2] + i * xi6_v[2];

        topRightDataL0[21] += d * xi6_u[3] + e * xi6_v[3];
        topRightDataL0[22] += f * xi6_u[3] + g * xi6_v[3];
        topRightDataL0[23] += h * xi6_u[3] + i * xi6_v[3];

        topRightDataL0[24] += d * xi6_u[4] + e * xi6_v[4];
        topRightDataL0[25] += f * xi6_u[4] + g * xi6_v[4];
        topRightDataL0[26] += h * xi6_u[4] + i * xi6_v[4];

        topRightDataL0[27] += d * xi6_u[5] + e * xi6_v[5];
        topRightDataL0[28] += f * xi6_u[5] + g * xi6_v[5];
        topRightDataL0[29] += h * xi6_u[5] + i * xi6_v[5];

        float j = J_res_d_aff_2(0,0);
        float k = J_res_d_aff_2(0,1);
        float l = J_res_d_aff_J_r[0];
        float m = J_res_d_aff_2(1,1);
        float n = J_res_d_aff_J_r[1];
        float o = res2;

        botRightDataL0[0] += j;
        botRightDataL0[1] += k;
        botRightDataL0[2] += l;
        botRightDataL0[3] += m;
        botRightDataL0[4] += n;
        botRightDataL0[5] += o;



        hNum++;
        hNumL0++;
        shiftUp(false);
    }

    // 得到最终的H(13,13)
    void done()
    {
        H.setZero(13,13);
        shiftUp(true);
        assert(hNumL0 == hNumL1 == 0);

        int topLeftDataIdx = 0;
        for(int i = 0; i < 10; i++)
            for(int j = 0; j < 10; j++)
            {
                H(i,j) = H(j,i) = topLeftDataL2[topLeftDataIdx++];
            }

        int topRightDataIdx = 0;
        for(int i = 0; i < 10; i++)
            for(int j = 10; j < 13; j++)
            {
                H(i,j) = H(j,i) = topRightDataL2[topRightDataIdx++];
            }

        int botRightDataIdx = 0;
        for(int i = 10; i <13; i++)
            for(int j = 10; j < 13; j++)
            {
                H(i,j) = H(j,i) = botRightDataL2[botRightDataIdx++];
            }
    }

public:
//    Mat1313f H;

private:
    float topLeftDataL0[56];    // 13*13的H矩阵的左上部分,10*10大小的对称方阵，需要保存的矩阵元素个数为(10*10-10)/2+10 =55
    float topLeftDataL1[56];        // SSE2指令集使用8个128位(16字节)寄存器，对于float数据(32位，4字节)，一个寄存器可以储存4个float数据
    float topLeftDataL2[56];

    float topRightDataL0[32];   // H右上，3*10 = 30
    float topRightDataL1[32];
    float topRightDataL2[32];

    float botRightDataL0[8];    // H右下，3*3-3/2+3 = 6
    float botRightDataL1[8];
    float botRightDataL2[8];

    void shiftUp(bool force)    // L0 --> L1 --> L2  多层级向上移位
    {
        if(hNumL0 > 1000 || force)
        {
             // L1 = L1 + L0, L0 = 0
            for(int i = 0; i < 56; i=i+4)   //  i=i+4 一次对4组float数据进行处理
                _mm_store_ps(topLeftDataL1 + i, _mm_add_ps(_mm_load_ps(topLeftDataL1 + i), _mm_load_ps(topLeftDataL0 + i)));
            for(int i = 0; i < 32; i = i + 4)
                _mm_store_ps(topRightDataL1 + i, _mm_add_ps(_mm_load_ps(topRightDataL1 + i), _mm_load_ps(topRightDataL0 + i)));
            for(int i = 0; i < 8; i = i + 4)
                _mm_store_ps(botRightDataL1 + i, _mm_add_ps(_mm_load_ps(botRightDataL1 + i), _mm_load_ps(botRightDataL0 + i)));

            hNumL1 += hNumL0;
            hNumL0 = 0;

            memset(topLeftDataL0, 0, sizeof(topLeftDataL0));
            memset(topRightDataL0, 0, sizeof(topRightDataL0));
            memset(botRightDataL0, 0, sizeof(botRightDataL0));
        }

        if(hNumL1 > 1000 || force)
        {
            // L2 = L2 + L1, L1 = 0
           for(int i = 0; i < 56; i=i+4)   //  i=i+4 一次对4组float数据进行处理
               _mm_store_ps(topLeftDataL2 + i, _mm_add_ps(_mm_load_ps(topLeftDataL2 + i), _mm_load_ps(topLeftDataL1 + i)));
           for(int i = 0; i < 32; i = i + 4)
               _mm_store_ps(topRightDataL2 + i, _mm_add_ps(_mm_load_ps(topRightDataL2 + i), _mm_load_ps(topRightDataL1 + i)));
           for(int i = 0; i < 8; i = i + 4)
               _mm_store_ps(botRightDataL2 + i, _mm_add_ps(_mm_load_ps(botRightDataL2 + i), _mm_load_ps(botRightDataL1 + i)));

           hNumL2 += hNumL1;
           hNumL1 = 0;

           memset(topLeftDataL1, 0, sizeof(topLeftDataL1));
           memset(topRightDataL1, 0, sizeof(topRightDataL1));
           memset(botRightDataL1, 0, sizeof(botRightDataL1));
        }
    }
};

/* 直接法不优化相机内参  6 + 2 + 1 = 9 */
class AccumulateMatrix99f : public AccumulateMatrixXXf
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    AccumulateMatrix99f() = default;
    ~AccumulateMatrix99f() = default;
    AccumulateMatrix99f(const AccumulateMatrix99f& other) = delete;
    AccumulateMatrix99f& operator =(const AccumulateMatrix99f& other) = delete;
public:
    void init()
    {
        memset(topLeftDataL0, 0, sizeof(topLeftDataL0));
        memset(topLeftDataL1, 0, sizeof(topLeftDataL1));
        memset(topLeftDataL2, 0, sizeof(topLeftDataL2));

        memset(topRightDataL0, 0, sizeof(topRightDataL0));
        memset(topRightDataL1, 0, sizeof(topRightDataL1));
        memset(topRightDataL2, 0, sizeof(topRightDataL2));

        memset(botRightDataL0, 0, sizeof(botRightDataL0));
        memset(botRightDataL1, 0, sizeof(botRightDataL1));
        memset(botRightDataL2, 0, sizeof(botRightDataL2));

        hNum = hNumL0 = hNumL1 = hNumL2 = 0;
    }
    void upDate(const float* const xi6_u, const float* const xi6_v,
                const Mat22f& J_I_d_p_2,                                                // topLeft
                const Mat22f& J_res_d_aff_J_I_d_p, const Vec2f& J_I_d_p_J_r,            // topRight
                const Mat22f& J_res_d_aff_2, const Vec2f& J_res_d_aff_J_r, float res2   // botRight
                )
    {
        float a = J_I_d_p_2(0, 0);
        float b = J_I_d_p_2(0, 1);
        float c = J_I_d_p_2(1, 1);

        topLeftDataL0[0] += a * xi6_u[0] * xi6_u[0] + b * xi6_u[0] * xi6_v[0] + b * xi6_v[0] * xi6_u[0] + c * xi6_v[0] * xi6_v[0];
        topLeftDataL0[1] += a * xi6_u[0] * xi6_u[1] + b * xi6_u[0] * xi6_v[1] + b * xi6_v[0] * xi6_u[1] + c * xi6_v[0] * xi6_v[1];
        topLeftDataL0[2] += a * xi6_u[0] * xi6_u[2] + b * xi6_u[0] * xi6_v[2] + b * xi6_v[0] * xi6_u[2] + c * xi6_v[0] * xi6_v[2];
        topLeftDataL0[3] += a * xi6_u[0] * xi6_u[3] + b * xi6_u[0] * xi6_v[3] + b * xi6_v[0] * xi6_u[3] + c * xi6_v[0] * xi6_v[3];
        topLeftDataL0[4] += a * xi6_u[0] * xi6_u[4] + b * xi6_u[0] * xi6_v[4] + b * xi6_v[0] * xi6_u[4] + c * xi6_v[0] * xi6_v[4];
        topLeftDataL0[5] += a * xi6_u[0] * xi6_u[5] + b * xi6_u[0] * xi6_v[5] + b * xi6_v[0] * xi6_u[5] + c * xi6_v[0] * xi6_v[5];

        topLeftDataL0[6]  += a * xi6_u[1] * xi6_u[1] + b * xi6_u[1] * xi6_v[1] + b * xi6_v[1] * xi6_u[1] + c * xi6_v[1] * xi6_v[1];
        topLeftDataL0[7]  += a * xi6_u[1] * xi6_u[2] + b * xi6_u[1] * xi6_v[2] + b * xi6_v[1] * xi6_u[2] + c * xi6_v[1] * xi6_v[2];
        topLeftDataL0[8]  += a * xi6_u[1] * xi6_u[3] + b * xi6_u[1] * xi6_v[3] + b * xi6_v[1] * xi6_u[3] + c * xi6_v[1] * xi6_v[3];
        topLeftDataL0[9]  += a * xi6_u[1] * xi6_u[4] + b * xi6_u[1] * xi6_v[4] + b * xi6_v[1] * xi6_u[4] + c * xi6_v[1] * xi6_v[4];
        topLeftDataL0[10] += a * xi6_u[1] * xi6_u[5] + b * xi6_u[1] * xi6_v[5] + b * xi6_v[1] * xi6_u[5] + c * xi6_v[1] * xi6_v[5];

        topLeftDataL0[11] += a * xi6_u[2] * xi6_u[2] + b * xi6_u[2] * xi6_v[2] + b * xi6_v[2] * xi6_u[2] + c * xi6_v[2] * xi6_v[2];
        topLeftDataL0[12] += a * xi6_u[2] * xi6_u[3] + b * xi6_u[2] * xi6_v[3] + b * xi6_v[2] * xi6_u[3] + c * xi6_v[2] * xi6_v[3];
        topLeftDataL0[13] += a * xi6_u[2] * xi6_u[4] + b * xi6_u[2] * xi6_v[4] + b * xi6_v[2] * xi6_u[4] + c * xi6_v[2] * xi6_v[4];
        topLeftDataL0[14] += a * xi6_u[2] * xi6_u[5] + b * xi6_u[2] * xi6_v[5] + b * xi6_v[2] * xi6_u[5] + c * xi6_v[2] * xi6_v[5];

        topLeftDataL0[15] += a * xi6_u[3] * xi6_u[3] + b * xi6_u[3] * xi6_v[3] + b * xi6_v[3] * xi6_u[3] + c * xi6_v[3] * xi6_v[3];
        topLeftDataL0[16] += a * xi6_u[3] * xi6_u[4] + b * xi6_u[3] * xi6_v[4] + b * xi6_v[3] * xi6_u[4] + c * xi6_v[3] * xi6_v[4];
        topLeftDataL0[17] += a * xi6_u[3] * xi6_u[5] + b * xi6_u[3] * xi6_v[5] + b * xi6_v[3] * xi6_u[5] + c * xi6_v[3] * xi6_v[5];

        topLeftDataL0[18] += a * xi6_u[4] * xi6_u[4] + b * xi6_u[4] * xi6_v[4] + b * xi6_v[4] * xi6_u[4] + c * xi6_v[4] * xi6_v[4];
        topLeftDataL0[19] += a * xi6_u[4] * xi6_u[5] + b * xi6_u[4] * xi6_v[5] + b * xi6_v[4] * xi6_u[5] + c * xi6_v[4] * xi6_v[5];

        topLeftDataL0[20] += a * xi6_u[5] * xi6_u[5] + b * xi6_u[5] * xi6_v[5] + b * xi6_v[5] * xi6_u[5] + c * xi6_v[5] * xi6_v[5];

        float d = J_res_d_aff_J_I_d_p(0,0);
        float e = J_res_d_aff_J_I_d_p(0,1);
        float f = J_res_d_aff_J_I_d_p(1,0);
        float g = J_res_d_aff_J_I_d_p(1,1);

        float h = J_I_d_p_J_r[0];
        float i = J_I_d_p_J_r[1];

        topRightDataL0[0] += d * xi6_u[0] + e * xi6_v[0];
        topRightDataL0[1] += f * xi6_u[0] + g * xi6_v[0];
        topRightDataL0[2] += h * xi6_u[0] + i * xi6_v[0];

        topRightDataL0[3] += d * xi6_u[1] + e * xi6_v[1];
        topRightDataL0[4] += f * xi6_u[1] + g * xi6_v[1];
        topRightDataL0[5] += h * xi6_u[1] + i * xi6_v[1];

        topRightDataL0[6] += d * xi6_u[2] + e * xi6_v[2];
        topRightDataL0[7] += f * xi6_u[2] + g * xi6_v[2];
        topRightDataL0[8] += h * xi6_u[2] + i * xi6_v[2];

        topRightDataL0[9]  += d * xi6_u[3] + e * xi6_v[3];
        topRightDataL0[10] += f * xi6_u[3] + g * xi6_v[3];
        topRightDataL0[11] += h * xi6_u[3] + i * xi6_v[3];

        topRightDataL0[12] += d * xi6_u[4] + e * xi6_v[4];
        topRightDataL0[13] += f * xi6_u[4] + g * xi6_v[4];
        topRightDataL0[14] += h * xi6_u[4] + i * xi6_v[4];

        topRightDataL0[15] += d * xi6_u[5] + e * xi6_v[5];
        topRightDataL0[16] += f * xi6_u[5] + g * xi6_v[5];
        topRightDataL0[17] += h * xi6_u[5] + i * xi6_v[5];

        float j = J_res_d_aff_2(0,0);
        float k = J_res_d_aff_2(0,1);
        float l = J_res_d_aff_J_r[0];
        float m = J_res_d_aff_2(1,1);
        float n = J_res_d_aff_J_r[1];
        float o = res2;

        botRightDataL0[0] += j;
        botRightDataL0[1] += k;
        botRightDataL0[2] += l;
        botRightDataL0[3] += m;
        botRightDataL0[4] += n;
        botRightDataL0[5] += o;


        hNum++;
        hNumL0++;
        shiftUp(false);
    }
    // 得到最终的H(9,9)
    void done()
    {
        H.setZero(9,9);
        shiftUp(true);
        assert(hNumL0 == hNumL1 == 0);

        int topLeftDataIdx = 0;
        for(int i = 0; i < 6; i++)
            for(int j = 0; j < 6; j++)
            {
                H(i,j) = H(j,i) = topLeftDataL2[topLeftDataIdx++];
            }

        int topRightDataIdx = 0;
        for(int i = 0; i < 6; i++)
            for(int j = 6; j < 9; j++)
            {
                H(i,j) = H(j,i) = topRightDataL2[topRightDataIdx++];
            }

        int botRightDataIdx = 0;
        for(int i = 6; i < 9; i++)
            for(int j = 6; j < 9; j++)
            {
                H(i,j) = H(j,i) = botRightDataL2[botRightDataIdx++];
            }
    }

public:
//    Mat99f H;
private:
    float topLeftDataL0[24];    // 9*9的H矩阵的左上部分,6*6大小的对称方阵，需要保存的矩阵元素个数为(6*6-6)/2+6 = 21
    float topLeftDataL1[24];        // SSE2指令集使用8个128位(16字节)寄存器，对于float数据(32位，4字节)，一个寄存器可以储存4个float数据
    float topLeftDataL2[24];

    float topRightDataL0[24];   // H右上，3*6 = 18
    float topRightDataL1[24];
    float topRightDataL2[24];

    float botRightDataL0[8];    // H右下，3*3-3/2+3 = 6
    float botRightDataL1[8];
    float botRightDataL2[8];

    void shiftUp(bool force)
    {
        if(hNumL0 > 1000 || force)
        {
             // L1 = L1 + L0, L0 = 0
            for(int i = 0; i < 56; i=i+4)   //  i=i+4 一次对4组float数据进行处理
                _mm_store_ps(topLeftDataL1 + i, _mm_add_ps(_mm_load_ps(topLeftDataL1 + i), _mm_load_ps(topLeftDataL0 + i)));
            for(int i = 0; i < 32; i = i + 4)
                _mm_store_ps(topRightDataL1 + i, _mm_add_ps(_mm_load_ps(topRightDataL1 + i), _mm_load_ps(topRightDataL0 + i)));
            for(int i = 0; i < 8; i = i + 4)
                _mm_store_ps(botRightDataL1 + i, _mm_add_ps(_mm_load_ps(botRightDataL1 + i), _mm_load_ps(botRightDataL0 + i)));

            hNumL1 += hNumL0;
            hNumL0 = 0;

            memset(topLeftDataL0, 0, sizeof(topLeftDataL0));
            memset(topRightDataL0, 0, sizeof(topRightDataL0));
            memset(botRightDataL0, 0, sizeof(botRightDataL0));
        }

        if(hNumL1 > 1000 || force)
        {
            // L2 = L2 + L1, L1 = 0
           for(int i = 0; i < 56; i=i+4)   //  i=i+4 一次对4组float数据进行处理
               _mm_store_ps(topLeftDataL2 + i, _mm_add_ps(_mm_load_ps(topLeftDataL2 + i), _mm_load_ps(topLeftDataL1 + i)));
           for(int i = 0; i < 32; i = i + 4)
               _mm_store_ps(topRightDataL2 + i, _mm_add_ps(_mm_load_ps(topRightDataL2 + i), _mm_load_ps(topRightDataL1 + i)));
           for(int i = 0; i < 8; i = i + 4)
               _mm_store_ps(botRightDataL2 + i, _mm_add_ps(_mm_load_ps(botRightDataL2 + i), _mm_load_ps(botRightDataL1 + i)));

           hNumL2 += hNumL1;
           hNumL1 = 0;

           memset(topLeftDataL1, 0, sizeof(topLeftDataL1));
           memset(topRightDataL1, 0, sizeof(topRightDataL1));
           memset(botRightDataL1, 0, sizeof(botRightDataL1));
        }
    }
};

/* 特征点法优化相机内参 4 + 6 + 1 = 11 */
class AccumulateMatrix1111f : public AccumulateMatrixXXf
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    AccumulateMatrix1111f() = default;
    ~AccumulateMatrix1111f() = default;
    AccumulateMatrix1111f(const AccumulateMatrix1111f& other) = delete;
    AccumulateMatrix1111f& operator =(const AccumulateMatrix1111f& other) = delete;
public:
    void init()
    {
        memset(topLeftDataL0, 0, sizeof(topLeftDataL0));
        memset(topLeftDataL1, 0, sizeof(topLeftDataL1));
        memset(topLeftDataL2, 0, sizeof(topLeftDataL2));

        memset(topRightDataL0, 0, sizeof(topRightDataL0));
        memset(topRightDataL1, 0, sizeof(topRightDataL1));
        memset(topRightDataL2, 0, sizeof(topRightDataL2));

        memset(botRightDataL0, 0, sizeof(botRightDataL0));
        memset(botRightDataL1, 0, sizeof(botRightDataL1));
        memset(botRightDataL2, 0, sizeof(botRightDataL2));

        hNum = hNumL0 = hNumL1 = hNumL2 = 0;
    }

    void upDate(const float* const C4_u, const float* const xi6_u,
                const float* const C4_v, const float* const xi6_v,
                const Vec2f& feature_res)
    {

        topLeftDataL0[0] += C4_u[0] * C4_u[0] + C4_v[0] * C4_v[0];
        topLeftDataL0[1] += C4_u[0] * C4_u[1] + C4_v[0] * C4_v[1];
        topLeftDataL0[2] += C4_u[0] * C4_u[2] + C4_v[0] * C4_v[2];
        topLeftDataL0[3] += C4_u[0] * C4_u[3] + C4_v[0] * C4_v[3];
        topLeftDataL0[4] += C4_u[0] * xi6_u[0] + C4_v[0] * xi6_v[0];
        topLeftDataL0[5] += C4_u[0] * xi6_u[1] + C4_v[0] * xi6_v[1];
        topLeftDataL0[6] += C4_u[0] * xi6_u[2] + C4_v[0] * xi6_v[2];
        topLeftDataL0[7] += C4_u[0] * xi6_u[3] + C4_v[0] * xi6_v[3];
        topLeftDataL0[8] += C4_u[0] * xi6_u[4] + C4_v[0] * xi6_v[4];
        topLeftDataL0[9] += C4_u[0] * xi6_u[5] + C4_v[0] * xi6_v[5];

        topLeftDataL0[10] += C4_u[1] * C4_u[1] + C4_v[1] * C4_v[1];
        topLeftDataL0[11] += C4_u[1] * C4_u[2] + C4_v[1] * C4_v[2];
        topLeftDataL0[12] += C4_u[1] * C4_u[3] + C4_v[1] * C4_v[3];
        topLeftDataL0[13] += C4_u[1] * xi6_u[0] + C4_v[1] * xi6_v[0];
        topLeftDataL0[14] += C4_u[1] * xi6_u[1] + C4_v[1] * xi6_v[1];
        topLeftDataL0[15] += C4_u[1] * xi6_u[2] + C4_v[1] * xi6_v[2];
        topLeftDataL0[16] += C4_u[1] * xi6_u[3] + C4_v[1] * xi6_v[3];
        topLeftDataL0[17] += C4_u[1] * xi6_u[4] + C4_v[1] * xi6_v[4];
        topLeftDataL0[18] += C4_u[1] * xi6_u[5] + C4_v[1] * xi6_v[5];

        topLeftDataL0[19] += C4_u[2] * C4_u[2] + C4_v[2] * C4_v[2];
        topLeftDataL0[20] += C4_u[2] * C4_u[3] + C4_v[2] * C4_v[3];
        topLeftDataL0[21] += C4_u[2] * xi6_u[0] + C4_v[2] * xi6_v[0];
        topLeftDataL0[22] += C4_u[2] * xi6_u[1] + C4_v[2] * xi6_v[1];
        topLeftDataL0[23] += C4_u[2] * xi6_u[2] + C4_v[2] * xi6_v[2];
        topLeftDataL0[24] += C4_u[2] * xi6_u[3] + C4_v[2] * xi6_v[3];
        topLeftDataL0[25] += C4_u[2] * xi6_u[4] + C4_v[2] * xi6_v[4];
        topLeftDataL0[26] += C4_u[2] * xi6_u[5] + C4_v[2] * xi6_v[5];

        topLeftDataL0[27] += C4_u[3] * C4_u[3] + C4_v[3] * C4_v[3];
        topLeftDataL0[28] += C4_u[3] * xi6_u[0] + C4_v[3] * xi6_v[0];
        topLeftDataL0[29] += C4_u[3] * xi6_u[1] + C4_v[3] * xi6_v[1];
        topLeftDataL0[30] += C4_u[3] * xi6_u[2] + C4_v[3] * xi6_v[2];
        topLeftDataL0[31] += C4_u[3] * xi6_u[3] + C4_v[3] * xi6_v[3];
        topLeftDataL0[32] += C4_u[3] * xi6_u[4] + C4_v[3] * xi6_v[4];
        topLeftDataL0[33] += C4_u[3] * xi6_u[5] + C4_v[3] * xi6_v[5];

        topLeftDataL0[34] += xi6_u[0] * xi6_u[0] + xi6_v[0] * xi6_v[0];
        topLeftDataL0[35] += xi6_u[0] * xi6_u[1] + xi6_v[0] * xi6_v[1];
        topLeftDataL0[36] += xi6_u[0] * xi6_u[2] + xi6_v[0] * xi6_v[2];
        topLeftDataL0[37] += xi6_u[0] * xi6_u[3] + xi6_v[0] * xi6_v[3];
        topLeftDataL0[38] += xi6_u[0] * xi6_u[4] + xi6_v[0] * xi6_v[4];
        topLeftDataL0[39] += xi6_u[0] * xi6_u[5] + xi6_v[0] * xi6_v[5];

        topLeftDataL0[40] += xi6_u[1] * xi6_u[1] + xi6_v[1] * xi6_v[1];
        topLeftDataL0[41] += xi6_u[1] * xi6_u[2] + xi6_v[1] * xi6_v[2];
        topLeftDataL0[42] += xi6_u[1] * xi6_u[3] + xi6_v[1] * xi6_v[3];
        topLeftDataL0[43] += xi6_u[1] * xi6_u[4] + xi6_v[1] * xi6_v[4];
        topLeftDataL0[44] += xi6_u[1] * xi6_u[5] + xi6_v[1] * xi6_v[5];

        topLeftDataL0[45] += xi6_u[2] * xi6_u[2] + xi6_v[2] * xi6_v[2];
        topLeftDataL0[46] += xi6_u[2] * xi6_u[3] + xi6_v[2] * xi6_v[3];
        topLeftDataL0[47] += xi6_u[2] * xi6_u[4] + xi6_v[2] * xi6_v[4];
        topLeftDataL0[48] += xi6_u[2] * xi6_u[5] + xi6_v[2] * xi6_v[5];

        topLeftDataL0[49] += xi6_u[3] * xi6_u[3] + xi6_v[3] * xi6_v[3];
        topLeftDataL0[50] += xi6_u[3] * xi6_u[4] + xi6_v[3] * xi6_v[4];
        topLeftDataL0[51] += xi6_u[3] * xi6_u[5] + xi6_v[3] * xi6_v[5];

        topLeftDataL0[52] += xi6_u[4] * xi6_u[4] + xi6_v[4] * xi6_v[4];
        topLeftDataL0[53] += xi6_u[4] * xi6_u[5] + xi6_v[4] * xi6_v[5];

        topLeftDataL0[54] += xi6_u[5] * xi6_u[5] + xi6_v[5] * xi6_v[5];

        float a = feature_res[0];
        float b = feature_res[1];

        topRightDataL0[0] += a * C4_u[0] + b * C4_v[0];
        topRightDataL0[1] += a * C4_u[1] + b * C4_v[1];
        topRightDataL0[2] += a * C4_u[2] + b * C4_v[2];
        topRightDataL0[3] += a * C4_u[3] + b * C4_v[3];
        topRightDataL0[4] += a * xi6_u[0] + b * xi6_v[0];
        topRightDataL0[5] += a * xi6_u[1] + b * xi6_v[1];
        topRightDataL0[6] += a * xi6_u[2] + b * xi6_v[2];
        topRightDataL0[7] += a * xi6_u[3] + b * xi6_v[3];
        topRightDataL0[8] += a * xi6_u[4] + b * xi6_v[4];
        topRightDataL0[9] += a * xi6_u[5] + b * xi6_v[5];

        botRightDataL0[0] += a * a + b * b;

        hNum++;
        hNumL0++;
        shiftUp(false);
    }

    void done()
    {
        H.setZero(11,11);
        shiftUp(true);
        assert(hNumL0 == hNumL1 == 0);

        int topLeftDataIdx = 0;
        for(int i = 0; i < 10; i++)
            for(int j = 0; j < 10; j++)
            {
                H(i,j) = H(j,i) = topLeftDataL2[topLeftDataIdx++];
            }

        int topRightDataIdx = 0;
        for(int i = 0; i < 10; i++)
            for(int j = 10; j < 11; j++)
            {
                H(i,j) = H(j,i) = topRightDataL2[topRightDataIdx++];
            }

        int botRightDataIdx = 0;
        for(int i = 10; i <11; i++)
            for(int j = 10; j < 11; j++)
            {
                H(i,j) = H(j,i) = botRightDataL2[botRightDataIdx++];
            }
    }

public:
//    Mat1111f H;
private:
    float topLeftDataL0[56];    // 11*11的H矩阵的左上部分,10*10大小的对称方阵，需要保存的矩阵元素个数为(10*10-10)/2+10 = 55
    float topLeftDataL1[56];        // SSE2指令集使用8个128位(16字节)寄存器，对于float数据(32位，4字节)，一个寄存器可以储存4个float数据
    float topLeftDataL2[56];

    float topRightDataL0[16];   // H右上，3*6 = 18
    float topRightDataL1[16];
    float topRightDataL2[16];

    float botRightDataL0[8];    // H右下，3*3-3/2+3 = 6
    float botRightDataL1[8];
    float botRightDataL2[8];

    void shiftUp(bool force)
    {
        if(hNumL0 > 1000 || force)
        {
             // L1 = L1 + L0, L0 = 0
            for(int i = 0; i < 56; i=i+4)   //  i=i+4 一次对4组float数据进行处理
                _mm_store_ps(topLeftDataL1 + i, _mm_add_ps(_mm_load_ps(topLeftDataL1 + i), _mm_load_ps(topLeftDataL0 + i)));
            for(int i = 0; i < 16; i = i + 4)
                _mm_store_ps(topRightDataL1 + i, _mm_add_ps(_mm_load_ps(topRightDataL1 + i), _mm_load_ps(topRightDataL0 + i)));
            for(int i = 0; i < 8; i = i + 4)
                _mm_store_ps(botRightDataL1 + i, _mm_add_ps(_mm_load_ps(botRightDataL1 + i), _mm_load_ps(botRightDataL0 + i)));

            hNumL1 += hNumL0;
            hNumL0 = 0;

            memset(topLeftDataL0, 0, sizeof(topLeftDataL0));
            memset(topRightDataL0, 0, sizeof(topRightDataL0));
            memset(botRightDataL0, 0, sizeof(botRightDataL0));
        }

        if(hNumL1 > 1000 || force)
        {
            // L2 = L2 + L1, L1 = 0
           for(int i = 0; i < 56; i=i+4)   //  i=i+4 一次对4组float数据进行处理
               _mm_store_ps(topLeftDataL2 + i, _mm_add_ps(_mm_load_ps(topLeftDataL2 + i), _mm_load_ps(topLeftDataL1 + i)));
           for(int i = 0; i < 32; i = i + 4)
               _mm_store_ps(topRightDataL2 + i, _mm_add_ps(_mm_load_ps(topRightDataL2 + i), _mm_load_ps(topRightDataL1 + i)));
           for(int i = 0; i < 8; i = i + 4)
               _mm_store_ps(botRightDataL2 + i, _mm_add_ps(_mm_load_ps(botRightDataL2 + i), _mm_load_ps(botRightDataL1 + i)));

           hNumL2 += hNumL1;
           hNumL1 = 0;

           memset(topLeftDataL1, 0, sizeof(topLeftDataL1));
           memset(topRightDataL1, 0, sizeof(topRightDataL1));
           memset(botRightDataL1, 0, sizeof(botRightDataL1));
        }
    }
};

/* 特征点法不优化相机内参 6 +1 =7 */
class AccumulateMatrix77f : public AccumulateMatrixXXf
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    AccumulateMatrix77f() = default;
    ~AccumulateMatrix77f() = default;
    AccumulateMatrix77f(const AccumulateMatrix77f& other) = delete;
    AccumulateMatrix77f& operator =(const AccumulateMatrix77f& other) = delete;
public:
    void init()
    {
        memset(topLeftDataL0, 0, sizeof(topLeftDataL0));
        memset(topLeftDataL1, 0, sizeof(topLeftDataL1));
        memset(topLeftDataL2, 0, sizeof(topLeftDataL2));

        memset(topRightDataL0, 0, sizeof(topRightDataL0));
        memset(topRightDataL1, 0, sizeof(topRightDataL1));
        memset(topRightDataL2, 0, sizeof(topRightDataL2));

        memset(botRightDataL0, 0, sizeof(botRightDataL0));
        memset(botRightDataL1, 0, sizeof(botRightDataL1));
        memset(botRightDataL2, 0, sizeof(botRightDataL2));

        hNum = hNumL0 = hNumL1 = hNumL2 = 0;
    }
    void upDate(const float* const xi6_u, const float* const xi6_v,
                const Vec2f& feature_res
                )
    {
        topLeftDataL0[0] += xi6_u[0] * xi6_u[0] + xi6_v[0] * xi6_v[0];
        topLeftDataL0[1] += xi6_u[0] * xi6_u[1] + xi6_v[0] * xi6_v[1];
        topLeftDataL0[2] += xi6_u[0] * xi6_u[2] + xi6_v[0] * xi6_v[2];
        topLeftDataL0[3] += xi6_u[0] * xi6_u[3] + xi6_v[0] * xi6_v[3];
        topLeftDataL0[4] += xi6_u[0] * xi6_u[4] + xi6_v[0] * xi6_v[4];
        topLeftDataL0[5] += xi6_u[0] * xi6_u[5] + xi6_v[0] * xi6_v[5];

        topLeftDataL0[6]  += xi6_u[1] * xi6_u[1] + xi6_v[1] * xi6_v[1];
        topLeftDataL0[7]  += xi6_u[1] * xi6_u[2] + xi6_v[1] * xi6_v[2];
        topLeftDataL0[8]  += xi6_u[1] * xi6_u[3] + xi6_v[1] * xi6_v[3];
        topLeftDataL0[9]  += xi6_u[1] * xi6_u[4] + xi6_v[1] * xi6_v[4];
        topLeftDataL0[10] += xi6_u[1] * xi6_u[5] + xi6_v[1] * xi6_v[5];

        topLeftDataL0[11] += xi6_u[2] * xi6_u[2] + xi6_v[2] * xi6_v[2];
        topLeftDataL0[12] += xi6_u[2] * xi6_u[3] + xi6_v[2] * xi6_v[3];
        topLeftDataL0[13] += xi6_u[2] * xi6_u[4] + xi6_v[2] * xi6_v[4];
        topLeftDataL0[14] += xi6_u[2] * xi6_u[5] + xi6_v[2] * xi6_v[5];

        topLeftDataL0[15] += xi6_u[3] * xi6_u[3] + xi6_v[3] * xi6_v[3];
        topLeftDataL0[16] += xi6_u[3] * xi6_u[4] + xi6_v[3] * xi6_v[4];
        topLeftDataL0[17] += xi6_u[3] * xi6_u[5] + xi6_v[3] * xi6_v[5];

        topLeftDataL0[18] += xi6_u[4] * xi6_u[4] + xi6_v[4] * xi6_v[4];
        topLeftDataL0[19] += xi6_u[4] * xi6_u[5] + xi6_v[4] * xi6_v[5];

        topLeftDataL0[20] += xi6_u[5] * xi6_u[5] + xi6_v[5] * xi6_v[5];

        float a = feature_res[0];
        float b = feature_res[1];

        topRightDataL0[0] += a * xi6_u[0] + b * xi6_v[0];
        topRightDataL0[1] += a * xi6_u[1] + b * xi6_v[1];
        topRightDataL0[2] += a * xi6_u[2] + b * xi6_v[2];
        topRightDataL0[3] += a * xi6_u[3] + b * xi6_v[3];
        topRightDataL0[4] += a * xi6_u[4] + b * xi6_v[4];
        topRightDataL0[5] += a * xi6_u[5] + b * xi6_v[5];

        botRightDataL0[0] += a * a + b * b;

        hNum++;
        hNumL0++;
        shiftUp(false);

    }
    void done()
    {
        H.setZero(7,7);
        shiftUp(true);
        assert(hNumL0 == hNumL1 == 0);

        int topLeftDataIdx = 0;
        for(int i = 0; i < 6; i++)
            for(int j = 0; j < 6; j++)
            {
                H(i,j) = H(j,i) = topLeftDataL2[topLeftDataIdx++];
            }

        int topRightDataIdx = 0;
        for(int i = 0; i < 6; i++)
            for(int j = 6; j < 7; j++)
            {
                H(i,j) = H(j,i) = topRightDataL2[topRightDataIdx++];
            }

        int botRightDataIdx = 0;
        for(int i = 6; i < 7; i++)
            for(int j = 6; j < 7; j++)
            {
                H(i,j) = H(j,i) = botRightDataL2[botRightDataIdx++];
            }
    }
public:
//    Mat77f H;    // 7*7-7 / 2 + 7 = 28
private:
    float topLeftDataL0[24]; // 6*6-6/2+6 = 21
    float topLeftDataL1[24];
    float topLeftDataL2[24];

    float topRightDataL0[8]; // 6
    float topRightDataL1[8];
    float topRightDataL2[8];

    float botRightDataL0[8]; // 1
    float botRightDataL1[8];
    float botRightDataL2[8];

    void shiftUp(bool force)
    {
        if(force || hNumL0 > 1000)
        {   // L1 = L1 + L0
            for(int i = 0; i < 24; i= i + 4)
                _mm_store_ps(topLeftDataL1 + i, _mm_add_ps(_mm_load_ps(&topLeftDataL1[i]), _mm_load_ps(&topLeftDataL0[i])));
            for(int i = 0; i < 8; i = i + 4)
                _mm_store_ps(topRightDataL1 + i, _mm_add_ps(_mm_load_ps(&topRightDataL1[i]), _mm_load_ps(&topRightDataL0[i])));
            for(int i = 0; i < 8; i = i + 4)
                _mm_store_ps(botRightDataL1 + i, _mm_add_ps(_mm_load_ps(&botRightDataL1[i]), _mm_load_ps(&botRightDataL0[i])));
            hNumL1 += hNumL0;

            // L0 = 0
            memset(topLeftDataL0, 0, sizeof(topLeftDataL0));
            memset(topRightDataL0, 0, sizeof(topRightDataL0));
            memset(botRightDataL0, 0, sizeof(botRightDataL0));
            hNumL0=0;
        }

        if(force || hNumL1 > 1000)
        {
            // L2 = L2 + L1
            for(int i = 0; i < 24; i = i + 4)
                _mm_store_ps(topLeftDataL2 + i, _mm_add_ps(_mm_load_ps(&topLeftDataL2[i]), _mm_load_ps(&topLeftDataL1[i])));
            for(int i = 0; i < 8; i = i + 4)
                _mm_store_ps(topRightDataL2 + i, _mm_add_ps(_mm_load_ps(&topRightDataL2[i]), _mm_load_ps(&topRightDataL1[i])));
            for(int i = 0; i < 8; i = i + 4)
                _mm_store_ps(botRightDataL2 + i, _mm_add_ps(_mm_load_ps(&botRightDataL2[i]), _mm_load_ps(&botRightDataL1[i])));

            // L1 = 0
            memset(topLeftDataL1, 0, sizeof(topLeftDataL1));
            memset(topRightDataL1, 0, sizeof(topRightDataL1));
            memset(botRightDataL1, 0, sizeof(botRightDataL1));
            hNumL1=0;
        }
    }

};

/* 得到全局的H,b */
/* 等到下一个版本再使用类模板，这一版本先用函数参数的形式代替模板参数 */
//template <TrackMode T1, CameraMode T2> class AccumulatedTopHessian      // finished 9.6 afternoon
//{
//public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//public:
//    AccumulatedTopHessian()
//    {
//        for(int thd_idx = 0; thd_idx < ThreadNum; thd_idx++)
//        {
//            resNum[thd_idx] = 0;
//            frameNum[thd_idx] = 0;
//            accH[thd_idx] = NULL;
//        }
//    }
//    ~AccumulatedTopHessian()
//    {
//        for(int thd_idx = 0; thd_idx < ThreadNum; thd_idx++)
//        {
//            if(accH[ThreadNum] != NULL)
//            {
//                delete accH[ThreadNum];
//                accH[ThreadNum] = NULL;
//            }
//        }
//    }

//    void setZero(int frameN, int min = 0, int max = 1, Vec10 *basicUnit = NULL, int thd_idx = 0)
//    {
//        if(frameN != frameNum[thd_idx])
//        {
//            if(accH[thd_idx] != NULL) delete [] accH[thd_idx];

//            if(T1 == Direct && T2 == Optimize)
//                accH[thd_idx] = new AccumulateMatrix1313f[frameN * frameN];
//            else if(T1 == Direct && T2 == Fix)
//                accH[thd_idx] = new AccumulateMatrix99f[frameN * frameN];
//            else if(T1 == Feature && T2 == Optimize)
//                accH[thd_idx] = new AccumulateMatrix1111f[frameN * frameN];
//            else if(T1 == Feature && T2 == Fix)
//                accH[thd_idx] = new AccumulateMatrix77f[frameN * frameN];
//            else
//                assert(false && "No such type!");
//        }
//        for(int i = 0; i < frameN * frameN; i++)
//        {
//            accH[thd_idx][i].init();
//        }
//        frameNum[thd_idx] = frameN;
//        resNum[thd_idx] = 0;
//    }

//    void stitchDouble(MatXX &H, VecX& b, const EnergyFunction* const EF, int thd_idx=0)  // finished 8.30       本函数横向得到某个线程的全局的H和b(在使用中就是用一层的H,即这里没有开启多线程)
//    {
//        assert(T1 == Direct || T1 == Feature);
//        assert(T2 == Optimize || T2 == Fix);

//        if(T1 == Direct && T2 == Optimize)  // H,b都是展开的
//        {
//            // 全局的海塞H
//            H = MatXX::Zero(frameNum[thd_idx] * 8 + CPARS, frameNum[thd_idx] * 8 + CPARS);
//            // 全局的b
//            b = VecX::Zero(frameNum[thd_idx] * 8 + CPARS);

//            for(int h = 0; h < frameNum[thd_idx]; h++)
//                for(int t = 0; t < frameNum[thd_idx]; t++)
//                {
//                    int host_idx = CPARS + h * 8;
//                    int target_idx = CPARS + t * 8;
//                    int partH_idx = h + frameNum[thd_idx] * t;

//                    // 在粘合出全局的H之前，调用done(),及时整理出一个最新的块矩阵accH
//                    accH[thd_idx][partH_idx].done();
//                    if(accH[thd_idx][partH_idx].hNum == 0) continue;
//                    // 转化为双精度
//                    Mat1313 partH = accH[thd_idx][partH_idx].H.cast<double>();
///* 全局的H与局部的H的对应关系 */
//                    // 把局部的H放到全局的H中(这部分只得到了H的左下侧，后边通过对称得到右上侧)
//                    H.block<8,8>(host_idx, host_idx).noalias() += EF->adHost[partH_idx] * partH.block<8,8>(CPARS,CPARS) * EF->adHost[partH_idx].transpose();
//                    H.block<8,8>(target_idx, target_idx).noalias() += EF->adTarget[partH_idx] * partH.block<8,8>(CPARS,CPARS) * EF->adTarget[partH_idx].transpose();
//                    H.block<8,8>(host_idx, target_idx).noalias() += EF->adHost[partH_idx] * partH.block<8,8>(CPARS,CPARS) * EF->adTarget[partH_idx].transpose();
//                    H.block<8,CPARS>(host_idx, 0).noalias() += EF->adHost[partH_idx] * partH.block<8,CPARS>(CPARS,0);
//                    H.block<8,CPARS>(target_idx, 0).noalias() += EF->adTarget[partH_idx] * partH.block<8,CPARS>(CPARS,0);
//                    H.topLeftCorner<CPARS,CPARS>().noalias() += partH.block<CPARS,CPARS>(0,0);

//                    // 全局的b由局部的H的最后一列得到
//                    b.segment<8>(host_idx).noalias() += EF->adHost[partH_idx] * partH.block<8,1>(CPARS, 8+CPARS);
//                    b.segment<8>(target_idx).noalias() += EF->adTarget[partH_idx] * partH.block<8,1>(CPARS, 8+CPARS);
//                    b.head<CPARS>().noalias() += partH.block<CPARS,1>(0, 8+CPARS);
//                }
//            // 以下两个for得到了全局H的右上侧
//            for(int t = 0; t < frameNum[thd_idx]; t++)
//            {
//                int target_idx = CPARS + t * 8;
//                H.block<CPARS,8>(0, target_idx).noalias() = H.block<8,CPARS>(target_idx,0).transpose();
//            }
//            for(int h = 0; h < frameNum[thd_idx]; h++)
//            {
//                int host_idx = CPARS + h * 8;
//                H.block<CPARS,8>(0, host_idx).noalias() = H.block<8,CPARS>(host_idx,0).transpose();

//                for(int t = h + 1; t < frameNum[thd_idx]; t++)
//                {
//                    int target_idx = CPARS + t * 8;
//                    H.block<8,8>(host_idx, target_idx).noalias() += H.block<8,8>(target_idx, host_idx).transpose();
//                    H.block<8,8>(target_idx, host_idx).noalias() = H.block<8,8>(host_idx, target_idx).transpose();
//                }
//            }
//        }

//        else if(T1 == Direct && T2 == Fix)  // 6 + 2 + 1(res)    相机内参C不优化
//        {
//            // 全局的海塞H
//            H = MatXX::Zero(frameNum[thd_idx] * 8, frameNum[thd_idx] * 8);
//            // 全局的b
//            b = VecX::Zero(frameNum[thd_idx] * 8);

//            for(int h = 0; h < frameNum[thd_idx]; h++)
//                for(int t = 0; t < frameNum[thd_idx]; t++)
//                {
//                    int host_idx = h * 8;
//                    int target_idx = t * 8;
//                    int partH_idx = h + frameNum[thd_idx] * t;

//                    // 在粘合出全局的H之前，调用done(),及时整理出一个最新的块矩阵accH
//                    accH[thd_idx][partH_idx].done();
//                    if(accH[thd_idx][partH_idx].hNum == 0) continue;
//                    // 转化为双精度
//                    Mat99 partH = accH[thd_idx][partH_idx].H.cast<double>();
///* 全局的H与局部的H的对应关系 */
//                    // 把局部的H放到全局的H中(这部分只得到了H的左下侧，后边通过对称得到右上侧)
//                    H.block<8,8>(host_idx, host_idx).noalias() += EF->adHost[partH_idx] * partH.block<8,8>(0,0) * EF->adHost[partH_idx].transpose();
//                    H.block<8,8>(target_idx, target_idx).noalias() += EF->adTarget[partH_idx] * partH.block<8,8>(0,0) * EF->adTarget[partH_idx].transpose();
//                    H.block<8,8>(host_idx, target_idx).noalias() += EF->adHost[partH_idx] * partH.block<8,8>(0,0) * EF->adTarget[partH_idx].transpose();

//                    // 全局的b由局部的H的最后一列得到
//                    b.segment<8>(host_idx).noalias() += EF->adHost[partH_idx] * partH.block<8,1>(0,8);
//                    b.segment<8>(target_idx).noalias() += EF->adTarget[partH_idx] * partH.block<8,1>(0,8);
//                }
//            // 以下for得到了全局H的右上侧
//            for(int h = 0; h < frameNum[thd_idx]; h++)
//            {
//                for(int t = h + 1; t < frameNum[thd_idx]; t++)
//                {
//                    int host_idx = h * 8;
//                    int target_idx = t * 8;
//                    H.block<8,8>(host_idx, target_idx).noalias() += H.block<8,8>(target_idx, host_idx).transpose();
//                    H.block<8,8>(target_idx, host_idx).noalias() = H.block<8,8>(host_idx, target_idx).transpose();
//                }
//            }
//        }

//        else if(T1 == Feature && T2 == Optimize)    // 4 + 6 + 1(res)   没有仿射亮度参数(a,b)
//        {
//            // 全局的海塞H
//            H = MatXX::Zero(frameNum[thd_idx] * 6 + CPARS, frameNum[thd_idx] * 6 + CPARS);
//            // 全局的b
//            b = VecX::Zero(frameNum[thd_idx] * 6 + CPARS);

//            for(int h = 0; h < frameNum[thd_idx]; h++)
//                for(int t = 0; t < frameNum[thd_idx]; t++)
//                {
//                    int host_idx = CPARS + h * 6;
//                    int target_idx = CPARS + t * 6;
//                    int partH_idx = h + frameNum[thd_idx] * t;

//                    // 在粘合出全局的H之前，调用done(),及时整理出一个最新的块矩阵accH
//                    accH[thd_idx][partH_idx].done();
//                    if(accH[thd_idx][partH_idx].hNum == 0) continue;
//                    // 转化为双精度
//                    Mat1111 partH = accH[thd_idx][partH_idx].H.cast<double>();
///* 全局的H与局部的H的对应关系 */
//                    // 把局部的H放到全局的H中(这部分只得到了H的左下侧，后边通过对称得到右上侧)
//                    H.block<6,6>(host_idx, host_idx).noalias() += EF->adHost_Feature[partH_idx] * partH.block<6,6>(CPARS,CPARS) * EF->adHost_Feature[partH_idx].transpose();
//                    H.block<6,6>(target_idx, target_idx).noalias() += EF->adTarget_Feature[partH_idx] * partH.block<6,6>(CPARS,CPARS) * EF->adTarget_Feature[partH_idx].transpose();
//                    H.block<6,6>(host_idx, target_idx).noalias() += EF->adHost_Feature[partH_idx] * partH.block<6,6>(CPARS,CPARS) * EF->adTarget_Feature[partH_idx].transpose();
//                    H.block<6,CPARS>(host_idx, 0).noalias() += EF->adHost_Feature[partH_idx] * partH.block<6,CPARS>(CPARS,0);
//                    H.block<6,CPARS>(target_idx, 0).noalias() += EF->adTarget_Feature[partH_idx] * partH.block<6,CPARS>(CPARS,0);
//                    H.topLeftCorner<CPARS,CPARS>().noalias() += partH.block<CPARS,CPARS>(0,0);

//                    // 全局的b由局部的H的最后一列得到
//                    b.segment<6>(host_idx).noalias() += EF->adHost_Feature[partH_idx] * partH.block<6,1>(CPARS, 6+CPARS);
//                    b.segment<6>(target_idx).noalias() += EF->adTarget_Feature[partH_idx] * partH.block<6,1>(CPARS, 6+CPARS);
//                    b.head<CPARS>().noalias() += partH.block<CPARS,1>(0, 6+CPARS);
//                }
//            // 以下两个for得到了全局H的右上侧
//            for(int t = 0; t < frameNum[thd_idx]; t++)
//            {
//                int target_idx = CPARS + t * 6;
//                H.block<CPARS,6>(0, target_idx).noalias() = H.block<6,CPARS>(target_idx,0).transpose();
//            }
//            for(int h = 0; h < frameNum[thd_idx]; h++)
//            {
//                int host_idx = CPARS + h * 6;
//                H.block<CPARS,6>(0, host_idx).noalias() += H.block<6,CPARS>(host_idx,0).transpose();

//                for(int t = h + 1; t < frameNum[thd_idx]; t++)
//                {
//                    int target_idx = CPARS + t * 6;
//                    H.block<6,6>(host_idx, target_idx).noalias() += H.block<6,6>(target_idx, host_idx).transpose();
//                    H.block<6,6>(target_idx, host_idx).noalias() = H.block<6,6>(host_idx, target_idx).transpose();
//                }
//            }
//        }

//        else if(T1 == Feature && T2 == Fix)     // 6 + 1 (res)
//        {
//            // 全局的海塞H
//            H = MatXX::Zero(frameNum[thd_idx] * 6, frameNum[thd_idx] * 6);
//            // 全局的b
//            b = VecX::Zero(frameNum[thd_idx] * 6);

//            for(int h = 0; h < frameNum[thd_idx]; h++)
//                for(int t = 0; t < frameNum[thd_idx]; t++)
//                {
//                    int host_idx = h * 6;
//                    int target_idx = t * 6;
//                    int partH_idx = h + frameNum[thd_idx] * t;

//                    // 在粘合出全局的H之前，调用done(),及时整理出一个最新的块矩阵accH
//                    accH[thd_idx][partH_idx].done();
//                    if(accH[thd_idx][partH_idx].hNum == 0) continue;
//                    // 转化为双精度
//                    Mat77 partH = accH[thd_idx][partH_idx].H.cast<double>();
///* 全局的H与局部的H的对应关系 */
//                    // 把局部的H放到全局的H中(这部分只得到了H的左下侧，后边通过对称得到右上侧)
//                    H.block<6,6>(host_idx, host_idx).noalias() += EF->adHost_Feature[partH_idx] * partH.block<6,6>(0,0) * EF->adHost_Feature[partH_idx].transpose();
//                    H.block<6,6>(target_idx, target_idx).noalias() += EF->adTarget_Feature[partH_idx] * partH.block<6,6>(0,0) * EF->adTarget_Feature[partH_idx].transpose();
//                    H.block<6,6>(host_idx, target_idx).noalias() += EF->adHost_Feature[partH_idx] * partH.block<6,6>(0,0) * EF->adTarget_Feature[partH_idx].transpose();

//                    // 全局的b由局部的H的最后一列得到
//                    b.segment<6>(host_idx).noalias() += EF->adHost_Feature[partH_idx] * partH.block<6,1>(0,6);
//                    b.segment<6>(target_idx).noalias() += EF->adTarget_Feature[partH_idx] * partH.block<6,1>(0,6);
//                }
//            // 以下for得到了全局H的右上侧
//            for(int h = 0; h < frameNum[thd_idx]; h++)
//            {
//                for(int t = h + 1; t < frameNum[thd_idx]; t++)
//                {
//                    int host_idx = h * 6;
//                    int target_idx = t * 6;
//                    H.block<6,6>(host_idx, target_idx).noalias() += H.block<6,6>(target_idx, host_idx).transpose();
//                    H.block<6,6>(target_idx, host_idx).noalias() = H.block<6,6>(host_idx, target_idx).transpose();
//                }
//            }
//        }

//        else
//        {
//            assert(false && "No such mode!");
//        }
//    }

//    void stitchDoubleMultiThread(MultiThread<Vec10>* reduce, MatXX &H, VecX& b, const EnergyFunction* const EF, bool usePrior, bool MT)
//    {
//        if(MT)  // 多线程得到全局的H,b
//        {
//            MatXX Hs[ThreadNum];
//            VecX bs[ThreadNum];
//            for(int i = 0; i < ThreadNum; i++)
//            {
//                assert(frameNum[0] == frameNum[i]); // 确保每一层(线程）的H和b大小相同

//                if(T1 == Direct && T2 == Optimize)
//                {
//                    Hs[i] = MatXX::Zero(frameNum[0] * 8 + CPARS, frameNum[0] * 8 + CPARS);
//                    bs[i] = VecX::Zero(frameNum[0] * 8 + CPARS);
//                } else if(T1 == Direct && T2 == Fix)
//                {
//                    Hs[i] = MatXX::Zero(frameNum[0] * 8, frameNum[0] * 8);
//                    bs[i] = VecX::Zero(frameNum[0] * 8);
//                }else if(T1 == Feature && T2 == Optimize)
//                {
//                    Hs[i] = MatXX::Zero(frameNum[0] * 6 + CPARS, frameNum[0] * 6 + CPARS);
//                    bs[i] = VecX::Zero(frameNum[0] * 6 + CPARS);
//                }else if(T1 == Feature && T2 == Fix)
//                {
//                    Hs[i] = MatXX::Zero(frameNum[0] * 6, frameNum[0] * 6);
//                    bs[i] = VecX::Zero(frameNum[0] * 6);
//                }else { assert(false && "No such mode!"); }
//            }
//            // 多线程求解Hs,bs每一层(线程)的全局的Ｈ，ｂ
//            auto newFunc = std::bind(&AccumulatedTopHessian::stitchDoubleInternal, this, Hs, bs, EF, usePrior, _1, _2, _3, _4);
//            reduce->reduce(newFunc, 0, frameNum[0] * frameNum[0]);
//            // 每一层的结果进行累加，得到最终的全局的H,b
//            H = Hs[0];
//            b = bs[0];

//            for(int i = 1; i < ThreadNum; i++)
//            {
//                H.noalias() += Hs[i];
//                b.noalias() += bs[i];
//                resNum[0] += resNum[i];
//            }
//        } else  // 单线程得到全局的H,b
//        {
//            if(T1 == Direct && T2 == Optimize)
//            {
//                H = MatXX::Zero(frameNum[0] * 8 + CPARS, frameNum[0] * 8 + CPARS);
//                b = VecX::Zero(frameNum[0] * 8 + CPARS);
//            } else if(T1 == Direct && T2 == Fix)
//            {
//                H = MatXX::Zero(frameNum[0] * 8, frameNum[0] * 8);
//                b = VecX::Zero(frameNum[0] * 8);
//            }else if(T1 == Feature && T2 == Optimize)
//            {
//                H = MatXX::Zero(frameNum[0] * 6 + CPARS, frameNum[0] * 6 + CPARS);
//                b = VecX::Zero(frameNum[0] * 6 + CPARS);
//            }else if(T1 == Feature && T2 == Fix)
//            {
//                H = MatXX::Zero(frameNum[0] * 6, frameNum[0] * 6);
//                b = VecX::Zero(frameNum[0] * 6);
//            }else { assert(false && "No such mode!"); }

//            stitchDoubleInternal(&H, &b, EF, usePrior, 0, frameNum[0] * frameNum[0], NULL, 0);
//        }
//    }
//        void stitchDoubleInternal(MatXX* H, VecX* b, const EnergyFunction* const EF, bool usePrior, int min, int max, Vec10* baiscUnit, int thd_idx)      // 我认为这样的stitch才是正确的
//        {
//            int toAggregate = ThreadNum;
//            if(thd_idx == -1) { toAggregate = 1; thd_idx = 0;}          // 不使用多线程时的处理
//            if(min == max) return;

//            if(T1 == Direct && T2 == Optimize)      // finished 9.5 afternoon
//            {
//                for(int k = min; k < max; k++)
//                {
//                    int host = k % frameNum[0];     // 帧的索引
//                    int target = k / frameNum[0];

//                    int host_idx = CPARS + host * 8;    // 帧的参数的索引
//                    int target_idx = CPARS + target * 8;
//                    int idx = host + frameNum[0] * target;  // k

//                    assert(idx == k);

//                    Mat1313 accH_double = Mat1313::Zero();

//                    for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)   // 把每个线程求解出的局部H进行对应的累加，得到总的局部H: accH_double, 它在全局H中的索引为idx
//                    {
//                        accH[thd_idx2][idx].done();     // 采用float数据进行H矩阵的求解，求出来再转换成double进行后续的使用
//                        if(accH[thd_idx2][idx].hNum == 0) continue;
//                        accH_double += accH[thd_idx2][idx].H.cast<double>();
//                    }
//                    /* 和stitchDouble中的哪个更正确呢？stitchDouble中的中间项为某层(某个线程)的局部H,而这里的中间项accH_double.block<8,8>(CPARS, CPARS)为所有层的局部H */
//                    /* 我认为都正确, 一个使用多线程，而另一个stitchDouble不使用 */
//                    H[thd_idx].block<8,8>(host_idx, host_idx).noalias() += EF->adHost[idx] * accH_double.block<8,8>(CPARS, CPARS) * EF->adHost[idx].transpose();    //
//                    H[thd_idx].block<8,8>(target_idx, target_idx).noalias() += EF->adTarget[idx] * accH_double.block<8,8>(CPARS, CPARS) * EF->adTarget[idx].transpose();
//                    H[thd_idx].block<8,8>(host_idx, target_idx).noalias() += EF->adHost[idx] * accH_double.block<8,8>(CPARS, CPARS) * EF->adTarget[idx].transpose();
//                    H[thd_idx].block<8,CPARS>(host_idx, 0).noalias() += EF->adHost[idx] * accH_double.block<8, CPARS>(CPARS, 0);
//                    H[thd_idx].block<8,CPARS>(target_idx, 0).noalias() += EF->adTarget[idx] * accH_double.block<8, CPARS>(CPARS, 0);
//                    H[thd_idx].topLeftCorner<CPARS, CPARS>().noalias() += accH_double.block<CPARS, CPARS>(0,0);

//                    b[thd_idx].segment<8>(host_idx).noalias() += EF->adHost[idx] * accH_double.block<8,1>(CPARS, 8+CPARS);
//                    b[thd_idx].segment<8>(target_idx).noalias() += EF->adTarget[idx] * accH_double.block<8,1>(CPARS, 8+CPARS);
//                    b[thd_idx].head<CPARS>().noalias() += accH_double.block<CPARS,1>(0, 8+CPARS);

//                }

//                if(min == 0 && usePrior)        // ???????????    作用是什么？Active的残差不需要usePrior，而Linearized的残差需要usePrior
//                {
//                    H[thd_idx].diagonal().head<CPARS>() += EF->cPrior;
//                    b[thd_idx].head<CPARS>() += EF->cPrior.cwiseProduct(EF->deltaF_C.cast<double>());
//                    for(int h = 0; h < frameNum[thd_idx]; h++)
//                    {
//                        H[thd_idx].diagonal().segment<8>(CPARS + h * 8) += EF->allBackKeyFrames[h]->state_prior;
//                        b[thd_idx].segment<8>(CPARS + h * 8) += EF->allBackKeyFrames[h]->state_prior.cwiseProduct(EF->allBackKeyFrames[h]->state_x);
//                    }
//                }
//                // 对称得到H的右上侧
//                for(int t = 0; t < frameNum[thd_idx]; t++)
//                {
//                    int target_idx = CPARS + t * 8;
//                    H[thd_idx].block<CPARS,8>(0, target_idx).noalias() = H[thd_idx].block<8,CPARS>(target_idx,0).transpose();
//                }
//                for(int h = 0; h < frameNum[thd_idx]; h++)
//                {
//                    int host_idx = CPARS + h * 8;
//                    H[thd_idx].block<CPARS,8>(0, host_idx).noalias() = H[thd_idx].block<8,CPARS>(host_idx,0).transpose();

//                    for(int t = h + 1; t < frameNum[thd_idx]; t++)
//                    {
//                        int target_idx = CPARS + t * 8;
//                        H[thd_idx].block<8,8>(host_idx, target_idx).noalias() += H[thd_idx].block<8,8>(target_idx, host_idx).transpose();
//                        H[thd_idx].block<8,8>(target_idx, host_idx).noalias() = H[thd_idx].block<8,8>(host_idx, target_idx).transpose();
//                    }
//                }
//            }

//            else if(T1 == Direct && T2 == Fix)      //
//            {
//                for(int k = min; k < max; k++)
//                {
//                    int host = k % frameNum[0];     // 帧的索引
//                    int target = k / frameNum[0];

//                    int host_idx = host * 8;    // 帧的参数的索引
//                    int target_idx = target * 8;
//                    int idx = host + frameNum[0] * target;  // k

//                    assert(idx == k);

//                    Mat99 accH_double = Mat99::Zero();

//                    for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)   // 把每个线程求解出的局部H进行对应的累加，得到总的局部H: accH_double, 它在全局H中的索引为idx
//                    {
//                        accH[thd_idx2][idx].done();     // 采用float数据进行H矩阵的求解，求出来再转换成double进行后续的使用
//                        if(accH[thd_idx2][idx].hNum == 0) continue;
//                        accH_double += accH[thd_idx2][idx].H.cast<double>();
//                    }
//                    /* 和stitchDouble中的哪个更正确呢？stitchDouble中的中间项为某层(某个线程)的局部H,而这里的中间项accH_double.block<8,8>(CPARS, CPARS)为所有层的局部H */
//                    /* 我认为都正确, 一个使用多线程，而另一个stitchDouble不使用 */
//                    H[thd_idx].block<8,8>(host_idx, host_idx).noalias() += EF->adHost[idx] * accH_double.block<8,8>(0, 0) * EF->adHost[idx].transpose();    //
//                    H[thd_idx].block<8,8>(target_idx, target_idx).noalias() += EF->adTarget[idx] * accH_double.block<8,8>(0, 0) * EF->adTarget[idx].transpose();
//                    H[thd_idx].block<8,8>(host_idx, target_idx).noalias() += EF->adHost[idx] * accH_double.block<8,8>(0, 0) * EF->adTarget[idx].transpose();


//                    b[thd_idx].segment<8>(host_idx).noalias() += EF->adHost[idx] * accH_double.block<8,1>(0, 8);
//                    b[thd_idx].segment<8>(target_idx).noalias() += EF->adTarget[idx] * accH_double.block<8,1>(0, 8);
//                }

//                if(min == 0 && usePrior)        // ???????????    作用是什么？Active的残差不需要usePrior，而Linearized的残差需要usePrior
//                {
//                    for(int h = 0; h < frameNum[thd_idx]; h++)
//                    {
//                        H[thd_idx].diagonal().segment<8>(h * 8) += EF->allBackKeyFrames[h]->state_prior;
//                        b[thd_idx].segment<8>(h * 8) += EF->allBackKeyFrames[h]->state_prior.cwiseProduct(EF->allBackKeyFrames[h]->state_x);
//                    }
//                }
//                // 对称得到H的右上侧
//                for(int h = 0; h < frameNum[thd_idx]; h++)
//                {
//                    int host_idx = h * 8;

//                    for(int t = h + 1; t < frameNum[thd_idx]; t++)
//                    {
//                        int target_idx = t * 8;

//                        H[thd_idx].block<8,8>(host_idx, target_idx).noalias() += H[thd_idx].block<8,8>(target_idx, host_idx).transpose();
//                        H[thd_idx].block<8,8>(target_idx, host_idx).noalias() = H[thd_idx].block<8,8>(host_idx, target_idx).transpose();
//                    }
//                }
//            }

//            else if(T1 == Feature && T2 == Optimize)
//            {
//                for(int k = min; k < max; k++)
//                {
//                    int host = k % frameNum[0];     // 帧的索引
//                    int target = k / frameNum[0];

//                    int host_idx = CPARS + host * 6;    // 帧的参数的索引
//                    int target_idx = CPARS + target * 6;
//                    int idx = host + frameNum[0] * target;  // k

//                    assert(idx == k);

//                    Mat1111 accH_double = Mat1111::Zero();

//                    for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)   // 把每个线程求解出的局部H进行对应的累加，得到总的局部H: accH_double, 它在全局H中的索引为idx
//                    {
//                        accH[thd_idx2][idx].done();     // 采用float数据进行H矩阵的求解，求出来再转换成double进行后续的使用
//                        if(accH[thd_idx2][idx].hNum == 0) continue;
//                        accH_double += accH[thd_idx2][idx].H.cast<double>();
//                    }

//                    H[thd_idx].block<6,6>(host_idx, host_idx).noalias() += EF->adHost_Feature[idx] * accH_double.block<6,6>(CPARS, CPARS) * EF->adHost_Feature[idx].transpose();    //
//                    H[thd_idx].block<6,6>(target_idx, target_idx).noalias() += EF->adTarget_Feature[idx] * accH_double.block<6,6>(CPARS, CPARS) * EF->adTarget_Feature[idx].transpose();
//                    H[thd_idx].block<6,6>(host_idx, target_idx).noalias() += EF->adHost_Feature[idx] * accH_double.block<6,6>(CPARS, CPARS) * EF->adTarget_Feature[idx].transpose();
//                    H[thd_idx].block<6,CPARS>(host_idx, 0).noalias() += EF->adHost_Feature[idx] * accH_double.block<6, CPARS>(CPARS, 0);
//                    H[thd_idx].block<6,CPARS>(target_idx, 0).noalias() += EF->adTarget_Feature[idx] * accH_double.block<6, CPARS>(CPARS, 0);
//                    H[thd_idx].topLeftCorner<CPARS, CPARS>().noalias() += accH_double.block<CPARS, CPARS>(0,0);

//                    b[thd_idx].segment<6>(host_idx).noalias() += EF->adHost_Feature[idx] * accH_double.block<6,1>(CPARS, 6+CPARS);
//                    b[thd_idx].segment<6>(target_idx).noalias() += EF->adTarget_Feature[idx] * accH_double.block<6,1>(CPARS, 6+CPARS);
//                    b[thd_idx].head<CPARS>().noalias() += accH_double.block<CPARS,1>(0, 6+CPARS);

//                }

//                if(min == 0 && usePrior)        // ???????????    作用是什么？Active的残差不需要usePrior，而Linearized的残差需要usePrior
//                {
//                    H[thd_idx].diagonal().head<CPARS>() += EF->cPrior;
//                    b[thd_idx].head<CPARS>() += EF->cPrior.cwiseProduct(EF->deltaF_C.cast<double>());
//                    for(int h = 0; h < frameNum[thd_idx]; h++)
//                    {
//                        H[thd_idx].diagonal().segment<6>(CPARS + h * 6) += EF->allBackKeyFrames[h]->state_prior_Feature;
//                        b[thd_idx].segment<6>(CPARS + h * 6) += EF->allBackKeyFrames[h]->state_prior_Feature.cwiseProduct(EF->allBackKeyFrames[h]->state_x_Feature);
//                    }
//                }
//                // 对称得到H的右上侧
//                for(int t = 0; t < frameNum[thd_idx]; t++)
//                {
//                    int target_idx = CPARS + t * 6;
//                    H[thd_idx].block<CPARS,6>(0, target_idx).noalias() = H[thd_idx].block<6,CPARS>(target_idx,0).transpose();
//                }
//                for(int h = 0; h < frameNum[thd_idx]; h++)
//                {
//                    int host_idx = CPARS + h * 6;
//                    H[thd_idx].block<CPARS,6>(0, host_idx).noalias() = H[thd_idx].block<6,CPARS>(host_idx,0).transpose();

//                    for(int t = h + 1; t < frameNum[thd_idx]; t++)
//                    {
//                        int target_idx = CPARS + t * 6;
//                        H[thd_idx].block<6,6>(host_idx, target_idx).noalias() += H[thd_idx].block<6,6>(target_idx, host_idx).transpose();
//                        H[thd_idx].block<6,6>(target_idx, host_idx).noalias() = H[thd_idx].block<6,6>(host_idx, target_idx).transpose();
//                    }
//                }
//            }

//            else if(T1 == Feature && T2 == Fix)
//            {
//                for(int k = min; k < max; k++)
//                {
//                    int host = k % frameNum[0];     // 帧的索引
//                    int target = k / frameNum[0];

//                    int host_idx = host * 6;    // 帧的参数的索引
//                    int target_idx = target * 6;
//                    int idx = host + frameNum[0] * target;  // k

//                    assert(idx == k);

//                    Mat77 accH_double = Mat77::Zero();

//                    for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)   // 把每个线程求解出的局部H进行对应的累加，得到总的局部H: accH_double, 它在全局H中的索引为idx
//                    {
//                        accH[thd_idx2][idx].done();     // 采用float数据进行H矩阵的求解，求出来再转换成double进行后续的使用
//                        if(accH[thd_idx2][idx].hNum == 0) continue;
//                        accH_double += accH[thd_idx2][idx].H.cast<double>();
//                    }

//                    H[thd_idx].block<6,6>(host_idx, host_idx).noalias() += EF->adHost_Feature[idx] * accH_double.block<6,6>(0, 0) * EF->adHost_Feature[idx].transpose();    //
//                    H[thd_idx].block<6,6>(target_idx, target_idx).noalias() += EF->adTarget_Feature[idx] * accH_double.block<6,6>(0, 0) * EF->adTarget_Feature[idx].transpose();
//                    H[thd_idx].block<6,6>(host_idx, target_idx).noalias() += EF->adHost_Feature[idx] * accH_double.block<6,6>(0, 0) * EF->adTarget_Feature[idx].transpose();


//                    b[thd_idx].segment<6>(host_idx).noalias() += EF->adHost_Feature[idx] * accH_double.block<6,1>(0, 6);
//                    b[thd_idx].segment<6>(target_idx).noalias() += EF->adTarget_Feature[idx] * accH_double.block<6,1>(0, 6);

//                }

//                if(min == 0 && usePrior)        // ???????????    作用是什么？Active的残差不需要usePrior，而Linearized的残差需要usePrior
//                {
//                    for(int h = 0; h < frameNum[thd_idx]; h++)
//                    {
//                        H[thd_idx].diagonal().segment<6>(h * 6) += EF->allBackKeyFrames[h]->state_prior_Feature;
//                        b[thd_idx].segment<6>(h * 6) += EF->allBackKeyFrames[h]->state_prior_Feature.cwiseProduct(EF->allBackKeyFrames[h]->state_x_Feature);
//                    }
//                }
//                // 对称得到H的右上侧
//                for(int h = 0; h < frameNum[thd_idx]; h++)
//                {
//                    int host_idx = h * 6;

//                    for(int t = h + 1; t < frameNum[thd_idx]; t++)
//                    {
//                        int target_idx = t * 6;

//                        H[thd_idx].block<6,6>(host_idx, target_idx).noalias() += H[thd_idx].block<6,6>(target_idx, host_idx).transpose();
//                        H[thd_idx].block<6,6>(target_idx, host_idx).noalias() = H[thd_idx].block<6,6>(host_idx, target_idx).transpose();
//                    }
//                }
//            }
//            else { assert(false && "No such mode!"); }
//        }

//    template<ResMode T> void addPoint(BackPixelPoint* point, const EnergyFunction* const EF, int thd_idx = 0)   // finished 9.4 14:41
//    {
//        assert(T == active || T == linearized || T == marginalize);

//        if(T1 == Direct && T2 == Optimize)
//        {
//            Vec4f delta_C = EF->deltaF_C;
//            float delta_idepth = point->delta_idepth;

//            float b_idepth_acc = 0;
//            float H_idepth_idepth_acc = 0;
//            Vec4f H_C_idepth_acc = Vec4f::Zero();

//            // 对于该点的每一个残差项
//            for(BackResidual* res : point->backResiduals)
//            {
//                if(T == active) { if(res->isLinearized || !res->isActiveAndIsGoodNEW) continue; }
//                else if(T == linearized) { if(!res->isLinearized || !res->isActiveAndIsGoodNEW) continue; }
//                else if(T == marginalize){ if(!res->isLinearized || !res->isActiveAndIsGoodNEW) continue; }
//                else{ assert(false && "No such mode!");}

//                ResidualJacobian* resJ = res->J;
//                int idx = res->backHostIndex + res->backTargetIndex * frameNum[thd_idx];
//                Mat18f delta = EF->adHTdeltaF[idx];

//                VecNRf resApprox;   // 直接法residual pattern点的个数
//                if(T == active)
//                    resApprox = resJ->res_direct;
//                else if(T == linearized)
//                {
//                    // 计算J_p_d_xi * delta
//                    __m128 J_p_delta_x = _mm_set1_ps(resJ->J_p_d_xi[0].dot(delta.head<6>()) + resJ->J_p_d_cam[0].dot(delta_C) + resJ->J_p_d_idepth[0] * delta_idepth);
//                    __m128 J_p_delta_y = _mm_set1_ps(resJ->J_p_d_xi[1].dot(delta.head<6>()) + resJ->J_p_d_cam[1].dot(delta_C) + resJ->J_p_d_idepth[1] * delta_idepth);
//                    __m128 delta_a = _mm_set1_ps((float)(delta[6]));
//                    __m128 delta_b = _mm_set1_ps((float)(delta[7]));

//                    for(int i = 0; i < res_pattern_pointNum; i = i+4)   // float数据占４字节
//                    {
//                        // delta_res: residual pattern中的残差变化量
//                        // delta_res = res_toZeroF_direct - [J_I_d_p*J_p_d_xi J_res_d_aff] * delta
//                        __m128 delta_res = _mm_load_ps(((float*)&res->res_toZeroF_direct)+i);
//                        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps(((float*)(resJ->J_I_d_p))+i), J_p_delta_x));
//                        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps(((float*)(resJ->J_I_d_p+1))+i), J_p_delta_y));
//                        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps(((float*)(resJ->J_res_d_aff))+i), delta_a));
//                        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps(((float*)(resJ->J_res_d_aff+1))+i), delta_b));
//                        _mm_store_ps(((float*)&resApprox)+i, delta_res);
//                    }
//                }
//                else if(T == marginalize)
//                {
//                    resApprox = res->res_toZeroF_direct;
//                }
//                else {
//                    assert(false && "No such mode!");
//                }

//                // compute J_I_d_p^T * r, and J_res_d_aff^T * r. (both are 2-vectors).
//                Vec2f J_I_d_p_J_r(0,0);
//                Vec2f J_res_d_aff_J_r(0,0);
//                float J_r_J_r;
//                for(int i = 0; i < res_pattern_pointNum; i++)
//                {
//                    J_I_d_p_J_r[0] += resApprox[i] * resJ->J_I_d_p[0][i];
//                    J_I_d_p_J_r[1] += resApprox[i] * resJ->J_I_d_p[1][i];
//                    J_res_d_aff_J_r[0] += resApprox[i] * resJ->J_res_d_aff[0][i];
//                    J_res_d_aff_J_r[1] += resApprox[i] * resJ->J_res_d_aff[1][i];
//                    J_r_J_r += resApprox[i] * resApprox[i];
//                }

//                {
//                    // H 第一部分,与深度相关的H, 将hostframe和targetframe之间的该点的所有残差项的Hessian进行求和
//                    Vec2f J_I_d_p_2_J_p_d_idepth = resJ->J_I_d_p_2 * resJ->J_p_d_idepth;
//                    b_idepth_acc += J_I_d_p_J_r.dot(resJ->J_p_d_idepth);
//                    H_idepth_idepth_acc += J_I_d_p_2_J_p_d_idepth.dot(resJ->J_p_d_idepth);
//                    H_C_idepth_acc += resJ->J_p_d_cam[0] * J_I_d_p_2_J_p_d_idepth[0] + resJ->J_p_d_cam[1] * J_I_d_p_2_J_p_d_idepth[1];

//                    // H 第二部分，与深度无关的H, 将hostframe和targetframe之间的该点的所有残差项的Hessian进行求和，得到hostframe和targetframe之间的accH
//                    accH[thd_idx][idx].upDate(resJ->J_p_d_cam[0].data(), resJ->J_p_d_xi[0].data(),
//                            resJ->J_p_d_cam[1].data(), resJ->J_p_d_xi[1].data(),
//                            resJ->J_I_d_p_2,
//                            resJ->J_res_d_aff_J_I_d_p,
//                            J_I_d_p_J_r,
//                            resJ->J_res_d_aff_2,
//                            J_res_d_aff_J_r,
//                            J_r_J_r);
//                }

//                resNum[thd_idx]++;
//            }

//            if(T == active)
//            {
//                point->H_idepth_idepth_accAF = H_idepth_idepth_acc;
//                point->b_idepth_accAF = b_idepth_acc;
//                point->H_C_idepth_accAF = H_C_idepth_acc;
//            }
//            if(T == linearized || T == marginalize)
//            {
//                point->H_idepth_idepth_accLF = H_idepth_idepth_acc;
//                point->b_idepth_accLF = b_idepth_acc;
//                point->H_C_idepth_accLF = H_C_idepth_acc;
//            }
//            if(T == marginalize)
//            {
//                point->H_idepth_idepth_accAF = 0;
//                point->H_C_idepth_accAF.setZero(4);
//                point->b_idepth_accAF = 0;
//            }
//        }
//        else if(T1 == Direct && T2 == Fix)      // 直接法固定相机内参
//        {
//            Vec4f delta_C = Vec4f::Zero();  // 直接法相机内参不优化，因此delta_C = 0
//            float delta_idepth = point->delta_idepth;

//            float b_idepth_acc = 0;
//            float H_idepth_idepth_acc = 0;
//            // Vec4f H_C_idepth_acc 这项不存在

//            // 对于该点的每一个残差项
//            for(BackResidual* res : point->backResiduals)
//            {
//                if(T == active) {if(res->isLinearized || !res->isActiveAndIsGoodNEW) continue;}
//                else if(T == linearized) {if(!res->isLinearized || !res->isActiveAndIsGoodNEW) continue;}
//                else if(T == marginalize) {if(!res->isLinearized || !res->isActiveAndIsGoodNEW) continue;}
//                else { assert(false && "No such mode!"); }

//                ResidualJacobian* resJ = res->J;
//                int idx = res->backHostIndex + res->backTargetIndex * frameNum[thd_idx];
//                Mat18f delta = EF->adHTdeltaF[idx];

//                VecNRf resApprox;   // 直接法residual pattern点的个数
//                if(T == active)
//                    resApprox = resJ->res_direct;
//                else if(T == linearized)
//                {
//                    assert(resJ->J_p_d_cam[0] == Vec4f::Zero() && resJ->J_p_d_cam[1] == Vec4f::Zero() && "直接法固定相机内参时，雅可比中J_p_d_cam应为0");
//                    // 计算J_p_d_xi * delta
//                    __m128 J_p_delta_x = _mm_set1_ps(resJ->J_p_d_xi[0].dot(delta.head<6>()) + resJ->J_p_d_cam[0].dot(delta_C) + resJ->J_p_d_idepth[0] * delta_idepth);
//                    __m128 J_p_delta_y = _mm_set1_ps(resJ->J_p_d_xi[1].dot(delta.head<6>()) + resJ->J_p_d_cam[1].dot(delta_C) + resJ->J_p_d_idepth[1] * delta_idepth);
//                    __m128 delta_a = _mm_set1_ps((float)(delta[6]));
//                    __m128 delta_b = _mm_set1_ps((float)(delta[7]));

//                    for(int i = 0; i < res_pattern_pointNum; i = i+4)   // float数据占４字节
//                    {
//                        // delta_res: residual pattern中的残差变化量
//                        // delta_res = res_toZeroF - [J_IJ_p_d_xi J_I_d_aff] * delta
//                        __m128 delta_res = _mm_load_ps(((float*)&res->res_toZeroF_direct)+i);
//                        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps(((float*)(resJ->J_I_d_p))+i), J_p_delta_x));
//                        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps(((float*)(resJ->J_I_d_p+1))+i), J_p_delta_y));
//                        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps(((float*)(resJ->J_res_d_aff))+i), delta_a));
//                        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps(((float*)(resJ->J_res_d_aff+1))+i), delta_b));
//                        _mm_store_ps(((float*)&resApprox)+i, delta_res);
//                    }

//                }
//            }

//            if(T == active)
//            {
//                point->H_idepth_idepth_accAF = H_idepth_idepth_acc;
//                point->b_idepth_accAF = b_idepth_acc;
//            }
//            if(T == linearized || T == marginalize)
//            {
//                point->H_idepth_idepth_accLF = H_idepth_idepth_acc;
//                point->b_idepth_accLF = b_idepth_acc;
//            }
//            if(T == marginalize)
//            {
//                point->H_idepth_idepth_accAF = 0;
//                point->H_C_idepth_accAF.setZero(4);
//                point->b_idepth_accAF = 0;
//            }
//        }
//        else if(T1 == Feature && T2 == Optimize)    // 特征点法优化相机内参           没有(a,b),不是residual pattern，不存在J_I_d_p, J_res_d_aff      这段需要仔细检查，很有可能有错
//        {
//            Vec4f delta_C = EF->deltaF_C;
//            float delta_idepth = point->delta_idepth;

//            float b_idepth_acc = 0;
//            float H_idepth_idepth_acc = 0;
//            Vec4f H_C_idepth_acc = Vec4f::Zero();

//            // 对于该点的每一个残差项
//            for(BackResidual* res : point->backResiduals)
//            {
//                if(T == active) { if(res->isLinearized || res->isActiveAndIsGoodNEW) continue; }
//                else if(T == linearized) { if(!res->isLinearized || !res->isActiveAndIsGoodNEW) continue; }
//                else if(T == marginalize) { if(!res->isLinearized || !res->isActiveAndIsGoodNEW) continue; }
//                else {assert(false && "No such mode!");}

//                ResidualJacobian* resJ = res->J;
//                int idx = res->backHostIndex + res->backTargetIndex * frameNum[thd_idx];
//                Mat16f delta = EF->adHTdeltaF_feature[idx];

//                Vec2f resApprox;
//                if(T == active)
//                    resApprox = resJ->res_feature_uv;
//                else if(T == linearized)
//                {
///* 重写J_p_delta_x, J_p_delta_y */
////                    // 计算J_p_d_xi * delta
////                    __m128 J_p_delta_x = _mm_set1_ps(resJ->J_p_d_xi[0].dot(delta.head<6>()) + resJ->J_p_d_cam[0].dot(delta_C) + resJ->J_p_d_idepth[0] * delta_idepth);
////                    __m128 J_p_delta_y = _mm_set1_ps(resJ->J_p_d_xi[1].dot(delta.head<6>()) + resJ->J_p_d_cam[1].dot(delta_C) + resJ->J_p_d_idepth[1] * delta_idepth);

////                    // delta_res: 单个点的残差变化量
////                    // delta_res = res_toZeroF_feature - [J_p_d_xi] * delta
////                    __m128 delta_res = _mm_load_ps((float*)(&res->res_toZeroF_feature));
////                    delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(resJ->J_p_d_xi)), J_p_delta_x));
////                    delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(resJ->J_p_d_xi+1)),J_p_delta_y));
////                    _mm_store_ps(((float*)&resApprox), delta_res);

////                    __m128 J_p_d_xi_delta_x = _mm_set1_ps(resJ->J_p_d_xi[0].dot(delta.head<6>()));
////                    __m128 J_p_d_xi_delta_y = _mm_set1_ps(resJ->J_p_d_xi[1].dot(delta.head<6>()));

//                    // delta_res = delta_res_0 + delta_res_new
////                    __m128 delta_res = _mm_load_ps((float*)(&res->res_toZeroF_feature));    // delta_res_0
////                    delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(resJ->J_p_d_xi)), _mm_load_ps((float*)(&delta))));     //delta_res_new
////                    delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(resJ->J_p_d_xi + 1)), _mm_load_ps((float*)(&delta))));
////                    delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(resJ->J_p_d_cam)), _mm_load_ps((float*)(&delta_C))));
////                    delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(resJ->J_p_d_cam + 1)), _mm_load_ps((float*)(&delta_C))));
////                    delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(&resJ->J_p_d_idepth[0])),_mm_load_ps((float*)(&delta_idepth))));
////                    delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(&resJ->J_p_d_idepth[1])),_mm_load_ps((float*)(&delta_idepth))));
////                    _mm_store_ps(((float*)&resApprox), delta_res);

//// 如果上边的SSE2指令存在语法问题，可使用以下语句直接求解resApprox
////                    resApprox = res->res_toZeroF_feature + resJ->J_p_d_xi[0].dot(delta) + resJ->J_p_d_cam[0].dot(delta_C) + resJ->J_p_d_idepth[0] * delta_idepth
////                                                         + resJ->J_p_d_xi[1].dot(delta) + resJ->J_p_d_cam[1].dot(delta_C) + resJ->J_p_d_idepth[1] * delta_idepth;

//                    __m128 J_p_delta_x = _mm_set1_ps(resJ->J_p_d_xi[0].dot(delta) + resJ->J_p_d_cam[0].dot(delta_C) + resJ->J_p_d_idepth[0] * delta_idepth);
//                    __m128 J_p_delta_y = _mm_set1_ps(resJ->J_p_d_xi[1].dot(delta) + resJ->J_p_d_cam[1].dot(delta_C) + resJ->J_p_d_idepth[1] * delta_idepth);

//                    // delta_res: 单个点的残差变化量
//                    // delta_res = res_toZeroF_feature - [J_p_d_xi] * delta
//                    __m128 delta_res_u = _mm_load_ps((float*)(&res->res_toZeroF_feature[0]));
//                    delta_res_u = _mm_add_ps(delta_res_u, _mm_mul_ps(delta_res_u, J_p_delta_x));
//                    __m128 delta_res_v = _mm_load_ps((float*)(&res->res_toZeroF_feature[1]));
//                    delta_res_v = _mm_add_ps(delta_res_v, _mm_mul_ps(delta_res_v, J_p_delta_y));
//                    _mm_store_ps(((float*)&resApprox[0]), delta_res_u);
//                    _mm_store_ps(((float*)&resApprox[1]), delta_res_v);


//                }
//                else if(T == marginalize) { resApprox = res->res_toZeroF_feature; }
//                else { assert(false && "No such mode!"); }

//                 // 不存在J_I_d_p, J_res_d_aff
////                float J_r_J_r = resApprox * resApprox;
//                {
//                    // H 第一部分,与深度相关的H, 将hostframe和targetframe之间的该点的所有残差项的Hessian进行求和

//                    b_idepth_acc += resApprox.dot(resJ->J_p_d_idepth);
//                    H_idepth_idepth_acc += resJ->J_p_d_idepth.dot(resJ->J_p_d_idepth);
//                    H_C_idepth_acc += resJ->J_p_d_cam[0] * resJ->J_p_d_idepth[0] + resJ->J_p_d_cam[1] * resJ->J_p_d_idepth[1];

//                    // H 第二部分，与深度无关的H, 将hostframe和targetframe之间的该点的所有残差项的Hessian进行求和，得到hostframe和targetframe之间的accH
//                    accH[thd_idx][idx].upDate(resJ->J_p_d_cam[0].data(), resJ->J_p_d_xi[0].data(),
//                                              resJ->J_p_d_cam[1].data(), resJ->J_p_d_xi[1].data(),
//                                              resApprox);
//                }

//                resNum[thd_idx]++;

//            }

//            if(T == active)
//            {
//                point->H_idepth_idepth_accAF = H_idepth_idepth_acc;
//                point->b_idepth_accAF = b_idepth_acc;
//                point->H_C_idepth_accAF = H_C_idepth_acc;
//            }
//            if(T == linearized || T == marginalize)
//            {
//                point->H_idepth_idepth_accLF = H_idepth_idepth_acc;
//                point->b_idepth_accLF = b_idepth_acc;
//                point->H_C_idepth_accLF = H_C_idepth_acc;
//            }
//            if(T == marginalize)
//            {
//                point->H_idepth_idepth_accAF = 0;
//                point->H_C_idepth_accAF.setZero(4);
//                point->b_idepth_accAF = 0;
//            }
//        }
//        else if(T1 == Feature && T2 == Fix)             // 没有C
//        {
//            Vec4f delta_C = EF->deltaF_C;
//            float delta_idepth = point->delta_idepth;

//            float b_idepth_acc = 0;
//            float H_idepth_idepth_acc = 0;

//            for(BackResidual* res : point->backResiduals)
//            {
//                if(T == active) { if(res->isLinearized || res->isActiveAndIsGoodNEW) continue; }
//                else if(T == linearized) { if(!res->isLinearized || !res->isActiveAndIsGoodNEW) continue; }
//                else if(T == marginalize) { if(!res->isLinearized || !res->isActiveAndIsGoodNEW) continue; }
//                else {assert(false && "No such mode!");}

//                ResidualJacobian* resJ = res->J;
//                int idx = res->backHostIndex + res->backTargetIndex * frameNum[thd_idx];
//                Mat16f delta = EF->adHTdeltaF_feature[idx];

//                Vec2f resApprox;
//                if(T == active)
//                    resApprox = resJ->res_feature_uv;
//                else if(T == linearized)
//                {
//                    __m128 J_p_delta_x = _mm_set1_ps(resJ->J_p_d_xi[0].dot(delta) + resJ->J_p_d_idepth[0] * delta_idepth);
//                    __m128 J_p_delta_y = _mm_set1_ps(resJ->J_p_d_xi[1].dot(delta) + resJ->J_p_d_idepth[1] * delta_idepth);

//                    // delta_res: 单个点的残差变化量
//                    // delta_res = res_toZeroF_feature - [J_p_d_xi] * delta
//                    __m128 delta_res_u = _mm_load_ps((float*)(&res->res_toZeroF_feature[0]));
//                    delta_res_u = _mm_add_ps(delta_res_u, _mm_mul_ps(delta_res_u, J_p_delta_x));
//                    __m128 delta_res_v = _mm_load_ps((float*)(&res->res_toZeroF_feature[1]));
//                    delta_res_v = _mm_add_ps(delta_res_v, _mm_mul_ps(delta_res_v, J_p_delta_y));
//                    _mm_store_ps(((float*)&resApprox[0]), delta_res_u);
//                    _mm_store_ps(((float*)&resApprox[1]), delta_res_v);
//                }
//                else if(T == marginalize) { resApprox = res->res_toZeroF_feature; }
//                else { assert(false && "No such mode!"); }

//                {
//                    // H 第一部分,与深度相关的H, 将hostframe和targetframe之间的该点的所有残差项的Hessian进行求和

//                    b_idepth_acc += resApprox.dot(resJ->J_p_d_idepth);
//                    H_idepth_idepth_acc += resJ->J_p_d_idepth.dot(resJ->J_p_d_idepth);

//                    // H 第二部分，与深度无关的H, 将hostframe和targetframe之间的该点的所有残差项的Hessian进行求和，得到hostframe和targetframe之间的accH
//                    accH[thd_idx][idx].upDate(resJ->J_p_d_xi[0].data(), resJ->J_p_d_xi[1].data(),
//                                              resApprox);
//                }

//                resNum[thd_idx]++;

//            }

//            if(T == active)
//            {
//                point->H_idepth_idepth_accAF = H_idepth_idepth_acc;
//                point->b_idepth_accAF = b_idepth_acc;
//                point->H_C_idepth_accAF = Vec4f::Zero();
//            }
//            if(T == linearized || T == marginalize)
//            {
//                point->H_idepth_idepth_accLF = H_idepth_idepth_acc;
//                point->b_idepth_accLF = b_idepth_acc;
//                point->H_C_idepth_accLF = Vec4f::Zero();
//            }
//            if(T == marginalize)
//            {
//                point->H_idepth_idepth_accAF = 0;
//                point->H_C_idepth_accAF.setZero(4);
//                point->b_idepth_accAF = 0;
//            }
//        }
//        else{ assert(false && "No such mode!"); }


//    }

//public:
//    int resNum[ThreadNum];
//    int frameNum[ThreadNum];
//    AccumulateMatrixXXf* accH[ThreadNum];
//};



class BackPixelPoint;
class EnergyFunction;
/* 累加的全局的H,b */
/* 第一版先用函数参数代替模板参数 */
class AccumulatedTopHessian
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    AccumulatedTopHessian(TrackMode T1, CameraMode T2)
    {
        this->T1 = T1;
        this->T2 = T2;
        for(int thd_idx = 0; thd_idx < ThreadNum; thd_idx++)
        {
            resNum[thd_idx] = 0;
            frameNum[thd_idx] = 0;
            accH[thd_idx] = NULL;
        }
    }
    ~AccumulatedTopHessian()
    {
        for(int thd_idx = 0; thd_idx < ThreadNum; thd_idx++)
        {
            if(accH[ThreadNum] != NULL)
            {
                delete accH[ThreadNum];
                accH[ThreadNum] = NULL;
            }
        }
    }

    void setZero(int frameN, int min = 0, int max = 1, Vec10 *basicUnit = NULL, int thd_idx = 0);

    void stitchDouble(MatXX &H, VecX& b, const EnergyFunction* const EF, int thd_idx=0);

    void stitchDoubleMultiThread(MultiThread<Vec10>* reduce, MatXX &H, VecX& b, const EnergyFunction* const EF, bool usePrior, bool MT);
        void stitchDoubleInternal(MatXX* H, VecX* b, const EnergyFunction* const EF, bool usePrior, int min, int max, Vec10* baiscUnit, int thd_idx);

    template<ResMode T> void addPointMultiThread(vector<BackPixelPoint*>* backPoints, const EnergyFunction* const EF, int min=0, int max=1, Vec10* basicUnit=0, int thd_idx = 0)
        {
            for(int i = min; i < max; i++) addPoint<T>((*backPoints)[i], EF, thd_idx);
        }

        template<ResMode T> void addPoint(BackPixelPoint* point, const EnergyFunction* const EF, int thd_idx = 0);


public:
    TrackMode T1;
    CameraMode T2;
    int resNum[ThreadNum];
    int frameNum[ThreadNum];
    AccumulateMatrixXXf* accH[ThreadNum];
};
}
#endif // ACCUMULATEDHESSIAN_H
