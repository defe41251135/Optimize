#ifndef HESSIAN_H
#define HESSIAN_H

#include "util/all_util_include.h"
#include "FullSystem/Settings.h"
#include "Backend/BackEnd.h"
#include "Backend/BackPart.h"

using namespace std;
using namespace world3000;

namespace SLAMSystem
{
/* schur complement 舒尔补求解，将H中的逆深度部分去掉，得到Hsc, 用于增量方程的求解(求解完再更新逆深度) */


template<int i, int j> class AccumulateMatrixPartXXf
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    void init()
    {
        matPartL0.setZero(i, j);
        matPartL1.setZero(i, j);
        matPartL2.setZero(i, j);

        num = numL0 = numL1 = numL2 = 0;
    }
    void done()
    {
        shiftUp(true);
        assert(numL0 == numL1 == 0);
        num = numL2;
    }
    void upDate(const Eigen::Matrix<float, i, 1>& L, const Eigen::Matrix<float, j, 1>& R, float w)
    {
        matPartL0 += w * L * R.transpose();
        numL0++;
        shiftUp(false);
    }
public:
    Eigen::Matrix<float, i, j> matPartL0;
    Eigen::Matrix<float, i, j> matPartL1;
    Eigen::Matrix<float, i, j> matPartL2;
    size_t num;

private:
    int numL0, numL1, numL2;

    void shiftUp(bool force)
    {
        if(numL0 > 1000 || force)
        {
            matPartL1 += matPartL0;
            matPartL0.setZero(i, j);

            numL1 += numL0;
            numL0 = 0;
        }
        if(numL1 > 1000 || force)
        {
            matPartL2 += matPartL1;
            matPartL1.setZero(i, j);

            numL2 += numL1;
            numL1 = 0;
        }
    }
};

template<int i> class AccumulateMatrixPartX1f
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    void init()
    {
        matPartL0.setZero(i, 1);
        matPartL1.setZero(i, 1);
        matPartL2.setZero(i ,1);

        num = numL0 = numL1 = numL2 = 0;
    }

    void done()
    {
        shiftUp(true);
        assert(numL0 == numL1 == 0);
        num = numL2;
    }

    void upDate(const Eigen::Matrix<float, i, 1>& L, float w)
    {
        matPartL0 += w * L;
        numL0++;
        shiftUp(false);
    }

    void upDateNoWeight(const Eigen::Matrix<float, i, 1>& L)
    {
        matPartL0 += L;
        numL0++;
        shiftUp(false);
    }

public:
    Eigen::Matrix<float, i, 1> matPartL0;
    Eigen::Matrix<float, i, 1> matPartL1;
    Eigen::Matrix<float, i, 1> matPartL2;
    size_t num;
private:
    int numL0, numL1, numL2;

    void shiftUp(bool force)
    {
        if(numL0 > 1000 || force)
        {
            matPartL1 += matPartL0;
            matPartL0.setZero(i, 1);

            numL1 += numL0;
            numL0 = 0;
        }
        if(numL1 > 1000 || force)
        {
            matPartL2 += matPartL1;
            matPartL1.setZero(i, 1);

            numL2 += numL1;
            numL1 = 0;
        }
    }
};

class AccumulatedSCHessian      /* finished 9.11 afternoon */
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    AccumulatedSCHessian(TrackMode T1, CameraMode T2)
    {
        this->T1 = T1;
        this->T2 = T2;
        for(int thd_idx = 0; thd_idx < ThreadNum; thd_idx++)
        {
            acc_HCF[thd_idx] = NULL;
            acc_HCF_Feature[thd_idx] = NULL;
            acc_HFF[thd_idx] = NULL;
            acc_bF[thd_idx] = NULL;
            frameNum[thd_idx] = 0;
        }
    }
    ~AccumulatedSCHessian()
    {
        for(int thd_idx = 0; thd_idx < ThreadNum; thd_idx++)
        {
            if(acc_HCF[thd_idx] != NULL) delete [] acc_HCF[thd_idx];      // 删除二维数组的第二维
            if(acc_HFF[thd_idx] != NULL) delete [] acc_HFF[thd_idx];
            if(acc_bF[thd_idx] != NULL) delete [] acc_bF[thd_idx];
        }
    }

    void setZero(int frameN, int min = 0, int max = 1, Vec10* basicUnit = NULL, int thd_idx = 0);

    void stitchDoubleMultiThread(MultiThread<Vec10>* reduce, MatXX &H, VecX &b, const EnergyFunction* const EF, bool MT);
    void stitchDoubleInternal(MatXX* H, VecX* b, const EnergyFunction* const EF, int min, int max, Vec10* basicUnit, int thd_idx);

    void addPointMultiThread(vector<BackPixelPoint*>* points, bool shiftPriorToZero, int min = 0, int max = 1, Vec10* basicUnit = NULL, int thd_idx = 0);
    void addPoint(BackPixelPoint* point, bool shiftPriorToZero, int thd_idx=0);
public:
    TrackMode T1;
    CameraMode T2;

    // H矩阵中相机(C)与帧(F)对应的块
    AccumulateMatrixPartXXf<8, CPARS>* acc_HCF[ThreadNum];     // 指针 ＋ 数组---->动态二维数组(行确定，列不确定)
    AccumulateMatrixPartXXf<6, CPARS>* acc_HCF_Feature[ThreadNum];
    // H矩阵中帧(F)与帧(F)对应的块
    AccumulateMatrixPartXXf<8, 8>* acc_HFF[ThreadNum];
    AccumulateMatrixPartXXf<6, 6>* acc_HFF_Feature[ThreadNum];
    // b向量中帧(F)对应的块
    AccumulateMatrixPartX1f<8>* acc_bF[ThreadNum];
    AccumulateMatrixPartX1f<6>* acc_bF_Feature[ThreadNum];

    // H矩阵中相机(C)与相机(C)对应的块(左上角4 * 4那个)
    AccumulateMatrixPartXXf<CPARS, CPARS> acc_HCC[ThreadNum];     // 只有一个这样的块，一维数组就可以表示了
    // b向量中相机(C)对应的块
    AccumulateMatrixPartX1f<CPARS> acc_bC[ThreadNum];

    int frameNum[ThreadNum];
};

}
#endif // HESSIAN_H
