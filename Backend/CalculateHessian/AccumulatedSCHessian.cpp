#include "Backend/CalculateHessian/AccumulatedSCHessian.h"

namespace SLAMSystem {

void AccumulatedSCHessian::setZero(int frameN, int min, int max, Vec10 *basicUnit, int thd_idx)            // todo 9.7
{
    if(acc_HCF[thd_idx] != NULL) delete [] acc_HCF[thd_idx];
    if(acc_HFF[thd_idx] != NULL) delete [] acc_HFF[thd_idx];
    if(acc_bF[thd_idx] != NULL) delete [] acc_bF[thd_idx];

    if(T1 == Direct && T2 == Optimize)
    {
        acc_HCF[thd_idx] = new AccumulateMatrixPartXXf<8, CPARS>[frameN * frameN];
        acc_HFF[thd_idx] = new AccumulateMatrixPartXXf<8, 8>[frameN * frameN * frameN]; // ???
        acc_bF[thd_idx] = new AccumulateMatrixPartX1f<8>[frameN * frameN];

        acc_HCF_Feature[thd_idx] = NULL;
        acc_HFF_Feature[thd_idx] = NULL;
        acc_bF_Feature[thd_idx] = NULL;

        for(int i = 0; i < frameN * frameN; i++)
        {
            acc_HCF[thd_idx][i].init();
            acc_bF[thd_idx][i].init();

            for(int j = 0; j < frameN; j++)
            {
                acc_HFF[thd_idx][i*frameN+j].init();
            }
        }

        acc_HCC[thd_idx].init();
        acc_bC[thd_idx].init();

    } else if(T1 == Direct && T2 == Fix)    // 相机内参固定不优化
    {
        acc_HCF[thd_idx] = NULL;    // 没有这一项，置零
        acc_HFF[thd_idx] = new AccumulateMatrixPartXXf<8, 8>[frameN * frameN * frameN];
        acc_bF[thd_idx] = new AccumulateMatrixPartX1f<8>[frameN * frameN];

        acc_HCF_Feature[thd_idx] = NULL;
        acc_HFF_Feature[thd_idx] = NULL;
        acc_bF_Feature[thd_idx] = NULL;

        for(int i = 0; i < frameN * frameN; i++)
        {
            acc_bF[thd_idx][i].init();
            for(int j = 0; j < frameN; j++)
            {
                acc_HFF[thd_idx][i*frameN+j].init();
            }
        }
        acc_HCC[thd_idx].init();    // 没有这一项，置零
        acc_bC[thd_idx].init();     // 没有这一项，置零
    } else if(T1 == Feature && T2 == Optimize)  // 8*8 --> 6*6
    {
        acc_HCF_Feature[thd_idx] = new AccumulateMatrixPartXXf<6, CPARS>[frameN * frameN];
        acc_HFF_Feature[thd_idx] = new AccumulateMatrixPartXXf<6, 6>[frameN * frameN * frameN];
        acc_bF_Feature[thd_idx] = new AccumulateMatrixPartX1f<6>[frameN * frameN];

        acc_HCF[thd_idx] = NULL;
        acc_HFF[thd_idx] = NULL;
        acc_bF[thd_idx] = NULL;

        for(int i = 0; i < frameN * frameN; i++)
        {
            acc_HCF[thd_idx][i].init();
            acc_bF[thd_idx][i].init();

            for(int j = 0; j < frameN; j++)
            {
                acc_HFF[thd_idx][i*frameN+j].init();
            }
        }

        acc_HCC[thd_idx].init();
        acc_bC[thd_idx].init();
    } else if(T1 == Feature && T2 == Fix)
    {
        acc_HCF_Feature[thd_idx] = NULL;
        acc_HFF_Feature[thd_idx] = new AccumulateMatrixPartXXf<6, 6>[frameN * frameN * frameN];
        acc_bF_Feature[thd_idx] = new AccumulateMatrixPartX1f<6>[frameN * frameN];

        acc_HCF[thd_idx] = NULL;
        acc_HFF[thd_idx] = NULL;
        acc_bF[thd_idx] = NULL;

        for(int i = 0; i < frameN * frameN; i++)
        {
            acc_bF[thd_idx][i].init();
            for(int j = 0; j < frameN; j++)
            {
                acc_HFF[thd_idx][i*frameN+j].init();
            }
        }
        acc_HCC[thd_idx].init();
        acc_bC[thd_idx].init();
    } else { assert(false && "No such mode!"); }
    frameNum[thd_idx] = frameN;
}

void AccumulatedSCHessian::stitchDoubleMultiThread(MultiThread<Vec10> *reduce, MatXX &H, VecX &b, const EnergyFunction * const EF, bool MT)
{

    if(MT)  // 多线程得到全局的H,b
    {
        MatXX Hs[ThreadNum];
        VecX bs[ThreadNum];
        for(int i = 0; i < ThreadNum; i++)
        {
            assert(frameNum[0] == frameNum[i]); // 确保每一层(线程）的H和b大小相同

            if(T1 == Direct && T2 == Optimize)
            {
                Hs[i] = MatXX::Zero(frameNum[0] * 8 + CPARS, frameNum[0] * 8 + CPARS);
                bs[i] = VecX::Zero(frameNum[0] * 8 + CPARS);
            } else if(T1 == Direct && T2 == Fix)
            {
                Hs[i] = MatXX::Zero(frameNum[0] * 8, frameNum[0] * 8);
                bs[i] = VecX::Zero(frameNum[0] * 8);
            }else if(T1 == Feature && T2 == Optimize)
            {
                Hs[i] = MatXX::Zero(frameNum[0] * 6 + CPARS, frameNum[0] * 6 + CPARS);
                bs[i] = VecX::Zero(frameNum[0] * 6 + CPARS);
            }else if(T1 == Feature && T2 == Fix)
            {
                Hs[i] = MatXX::Zero(frameNum[0] * 6, frameNum[0] * 6);
                bs[i] = VecX::Zero(frameNum[0] * 6);
            }else { assert(false && "No such mode!"); }
        }
        // 多线程求解Hs,bs每一层(线程)的全局的Ｈ，ｂ
        auto newFunc = std::bind(&AccumulatedSCHessian::stitchDoubleInternal, this, Hs, bs, EF, _1, _2, _3, _4);
        reduce->reduce(newFunc, 0, frameNum[0] * frameNum[0]);
        // 每一层的结果进行累加，得到最终的全局的H,b
        H = Hs[0];
        b = bs[0];

        for(int i = 1; i < ThreadNum; i++)
        {
            H.noalias() += Hs[i];
            b.noalias() += bs[i];
        }
    } else  // 单线程得到全局的H,b
    {
        if(T1 == Direct && T2 == Optimize)
        {
            H = MatXX::Zero(frameNum[0] * 8 + CPARS, frameNum[0] * 8 + CPARS);
            b = VecX::Zero(frameNum[0] * 8 + CPARS);
        } else if(T1 == Direct && T2 == Fix)
        {
            H = MatXX::Zero(frameNum[0] * 8, frameNum[0] * 8);
            b = VecX::Zero(frameNum[0] * 8);
        }else if(T1 == Feature && T2 == Optimize)
        {
            H = MatXX::Zero(frameNum[0] * 6 + CPARS, frameNum[0] * 6 + CPARS);
            b = VecX::Zero(frameNum[0] * 6 + CPARS);
        }else if(T1 == Feature && T2 == Fix)
        {
            H = MatXX::Zero(frameNum[0] * 6, frameNum[0] * 6);
            b = VecX::Zero(frameNum[0] * 6);
        }else { assert(false && "No such mode!"); }

        stitchDoubleInternal(&H, &b, EF, 0, frameNum[0] * frameNum[0], NULL, 0);
    }

}

void AccumulatedSCHessian::stitchDoubleInternal(MatXX *H, VecX *b, const EnergyFunction * const EF, int min, int max, Vec10 *basicUnit, int thd_idx)
{
    int toAggregate = ThreadNum;
    if(thd_idx == -1) { toAggregate = 1; thd_idx = 0; }         // 不使用多线程时的处理
    if(min == max) return;

    int frameN = frameNum[0];
    int frameN2 = frameN * frameN;

    if(T1 == Direct && T2 == Optimize)
    {
        for(int k = min; k < max; k++)
        {
            int host = k % frameN;  // 帧的索引
            int target = k / frameN;

            int host_idx = CPARS + host * 8;    // 帧的参数的索引
            int target_idx = CPARS + target * 8;
            int partH_idx = host + frameN * target; // 局部的H,b的索引

            assert(partH_idx == k);

            Mat84 H_CF = Mat84::Zero();
            Vec8 b_F = Vec8::Zero();

            for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)
            {
                acc_HCF[thd_idx2][partH_idx].done();
                acc_bF[thd_idx2][partH_idx].done();
                H_CF += acc_HCF[thd_idx2][partH_idx].matPartL2.cast<double>();
                b_F += acc_bF[thd_idx2][partH_idx].matPartL2.cast<double>();
            }

            H[thd_idx].block<8, CPARS>(host_idx, 0) += EF->adHost[partH_idx] * H_CF;
            H[thd_idx].block<8, CPARS>(target_idx, 0) += EF->adTarget[partH_idx] * H_CF;
            b[thd_idx].segment<8>(host_idx) += EF->adHost[partH_idx] * b_F;
            b[thd_idx].segment<8>(target_idx) += EF->adTarget[partH_idx] * b_F;

            for(int k = 0; k < frameN; k++)
            {
                int residual_idx = CPARS + k * 8;
                int host_target_residual_idx = partH_idx + k * frameN2;
                int host_residual_Idx = host + frameN * k;

                Mat88 HFF = Mat88::Zero();

                for(int thd_idx3 = 0; thd_idx3 < toAggregate; thd_idx3++)
                {
                    acc_HFF[thd_idx3][host_target_residual_idx].done();
                    if(acc_HFF[thd_idx3][host_target_residual_idx].num == 0) continue;
                    HFF += acc_HFF[thd_idx3][host_target_residual_idx].matPartL2.cast<double>();
                }

                H[thd_idx].block<8,8>(host_idx, host_idx) += EF->adHost[partH_idx] * HFF * EF->adHost[host_residual_Idx].transpose();
                H[thd_idx].block<8,8>(target_idx, residual_idx) += EF->adTarget[partH_idx] * HFF * EF->adTarget[host_residual_Idx].transpose();
                H[thd_idx].block<8,8>(target_idx, host_idx) += EF->adTarget[partH_idx] * HFF * EF->adHost[host_residual_Idx].transpose();
                H[thd_idx].block<8,8>(host_idx, residual_idx) += EF->adHost[partH_idx] * HFF * EF->adTarget[host_residual_Idx].transpose();
            }


        }

        for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)
        {
            acc_HCC[thd_idx2].done();
            acc_bC[thd_idx2].done();

            H[thd_idx].topLeftCorner<CPARS, CPARS>() += acc_HCC[thd_idx2].matPartL2.cast<double>();
            b[thd_idx].head<CPARS>() += acc_bC[thd_idx2].matPartL2.cast<double>();
        }
    }

    else if(T1 == Direct && T2 == Fix)
    {
        for(int k = min; k < max; k++)
        {
            int host  = k % frameN; // 帧的索引
            int target = k/ frameN;

            int host_idx = host * 8;
            int target_idx = target * 8;
            int partH_idx = host + frameN * target;

            assert(partH_idx == k);

            Vec8 b_F = Vec8::Zero();

            for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)
            {
                acc_bF[thd_idx2][partH_idx].done();
                b_F += acc_bF[thd_idx2][partH_idx].matPartL2.cast<double>();
            }

            b[thd_idx].segment<8>(host_idx) += EF->adHost[partH_idx] * b_F;
            b[thd_idx].segment<8>(target_idx) += EF->adTarget[partH_idx] * b_F;

            for(int k = 0; k < frameN; k++)
            {
                int residual_idx = k * 8;
                int host_target_residual_idx = partH_idx + k * frameN2;
                int host_residual_Idx = host + frameN * k;

                Mat88 HFF = Mat88::Zero();

                for(int thd_idx3 = 0; thd_idx3 < toAggregate; thd_idx3++)
                {
                    acc_HFF[thd_idx3][host_target_residual_idx].done();
                    if(acc_HFF[thd_idx3][host_target_residual_idx].num == 0) continue;
                    HFF += acc_HFF[thd_idx3][host_target_residual_idx].matPartL2.cast<double>();
                }

                H[thd_idx].block<8,8>(host_idx, host_idx) += EF->adHost[partH_idx] * HFF * EF->adHost[host_residual_Idx].transpose();
                H[thd_idx].block<8,8>(target_idx, residual_idx) += EF->adTarget[partH_idx] * HFF * EF->adTarget[host_residual_Idx].transpose();
                H[thd_idx].block<8,8>(target_idx, host_idx) += EF->adTarget[partH_idx] * HFF * EF->adHost[host_residual_Idx].transpose();
                H[thd_idx].block<8,8>(host_idx, residual_idx) += EF->adHost[partH_idx] * HFF * EF->adTarget[host_residual_Idx].transpose();
            }
        }


    }

    else if(T1 == Feature && T2 == Optimize)
    {
        for(int k = min; k < max; k++)
        {
            int host = k % frameN;  // 帧的索引
            int target = k / frameN;

            int host_idx = CPARS + host * 6;    // 帧的参数的索引
            int target_idx = CPARS + target * 6;
            int partH_idx = host + frameN * target; // 局部的H,b的索引

            assert(partH_idx == k);

            Mat64 H_CF_Feature = Mat64::Zero();
            Vec6 b_F_Feature = Vec6::Zero();

            for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)
            {
                acc_HCF_Feature[thd_idx2][partH_idx].done();
                acc_bF_Feature[thd_idx2][partH_idx].done();
                H_CF_Feature += acc_HCF_Feature[thd_idx2][partH_idx].matPartL2.cast<double>();
                b_F_Feature += acc_bF_Feature[thd_idx2][partH_idx].matPartL2.cast<double>();
            }

            H[thd_idx].block<6, CPARS>(host_idx, 0) += EF->adHost_Feature[partH_idx] * H_CF_Feature;
            H[thd_idx].block<6, CPARS>(target_idx, 0) += EF->adTarget_Feature[partH_idx] * H_CF_Feature;
            b[thd_idx].segment<6>(host_idx) += EF->adHost_Feature[partH_idx] * b_F_Feature;
            b[thd_idx].segment<6>(target_idx) += EF->adTarget_Feature[partH_idx] * b_F_Feature;

            for(int k = 0; k < frameN; k++)
            {
                int residual_idx = CPARS + k * 6;
                int host_target_residual_idx = partH_idx + k * frameN2;
                int host_residual_Idx = host + frameN * k;

                Mat66 HFF_Feature = Mat66::Zero();

                for(int thd_idx3 = 0; thd_idx3 < toAggregate; thd_idx3++)
                {
                    acc_HFF_Feature[thd_idx3][host_target_residual_idx].done();
                    if(acc_HFF_Feature[thd_idx3][host_target_residual_idx].num == 0) continue;
                    HFF_Feature += acc_HFF_Feature[thd_idx3][host_target_residual_idx].matPartL2.cast<double>();
                }

                H[thd_idx].block<6,6>(host_idx, host_idx) += EF->adHost_Feature[partH_idx] * HFF_Feature * EF->adHost_Feature[host_residual_Idx].transpose();
                H[thd_idx].block<6,6>(target_idx, residual_idx) += EF->adTarget_Feature[partH_idx] * HFF_Feature * EF->adTarget_Feature[host_residual_Idx].transpose();
                H[thd_idx].block<6,6>(target_idx, host_idx) += EF->adTarget_Feature[partH_idx] * HFF_Feature * EF->adHost_Feature[host_residual_Idx].transpose();
                H[thd_idx].block<6,6>(host_idx, residual_idx) += EF->adHost_Feature[partH_idx] * HFF_Feature * EF->adTarget_Feature[host_residual_Idx].transpose();
            }
        }

        for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)
        {
            acc_HCC[thd_idx2].done();
            acc_bC[thd_idx2].done();

            H[thd_idx].topLeftCorner<CPARS, CPARS>() += acc_HCC[thd_idx2].matPartL2.cast<double>();
            b[thd_idx].head<CPARS>() += acc_bC[thd_idx2].matPartL2.cast<double>();
        }
    }

    else if(T1 == Feature && T2 == Fix)
    {
        for(int k = min; k < max; k++)
        {
            int host = k % frameN;  // 帧的索引
            int target = k / frameN;

            int host_idx = host * 6;    // 帧的参数的索引
            int target_idx = target * 6;
            int partH_idx = host + frameN * target; // 局部的H,b的索引

            assert(partH_idx == k);


            Vec6 b_F_Feature = Vec6::Zero();

            for(int thd_idx2 = 0; thd_idx2 < toAggregate; thd_idx2++)
            {

                acc_bF_Feature[thd_idx2][partH_idx].done();
                b_F_Feature += acc_bF_Feature[thd_idx2][partH_idx].matPartL2.cast<double>();
            }

            b[thd_idx].segment<6>(host_idx) += EF->adHost_Feature[partH_idx] * b_F_Feature;
            b[thd_idx].segment<6>(target_idx) += EF->adTarget_Feature[partH_idx] * b_F_Feature;

            for(int k = 0; k < frameN; k++)
            {
                int residual_idx = k * 6;
                int host_target_residual_idx = partH_idx + k * frameN2;
                int host_residual_Idx = host + frameN * k;

                Mat66 HFF_Feature = Mat66::Zero();

                for(int thd_idx3 = 0; thd_idx3 < toAggregate; thd_idx3++)
                {
                    acc_HFF_Feature[thd_idx3][host_target_residual_idx].done();
                    if(acc_HFF_Feature[thd_idx3][host_target_residual_idx].num == 0) continue;
                    HFF_Feature += acc_HFF_Feature[thd_idx3][host_target_residual_idx].matPartL2.cast<double>();
                }

                H[thd_idx].block<6,6>(host_idx, host_idx) += EF->adHost_Feature[partH_idx] * HFF_Feature * EF->adHost_Feature[host_residual_Idx].transpose();
                H[thd_idx].block<6,6>(target_idx, residual_idx) += EF->adTarget_Feature[partH_idx] * HFF_Feature * EF->adTarget_Feature[host_residual_Idx].transpose();
                H[thd_idx].block<6,6>(target_idx, host_idx) += EF->adTarget_Feature[partH_idx] * HFF_Feature * EF->adHost_Feature[host_residual_Idx].transpose();
                H[thd_idx].block<6,6>(host_idx, residual_idx) += EF->adHost_Feature[partH_idx] * HFF_Feature * EF->adTarget_Feature[host_residual_Idx].transpose();
            }
        }
    }
    else { assert(false && "No such mode!"); }

}

void AccumulatedSCHessian::addPointMultiThread(vector<BackPixelPoint *> *points, bool shiftPriorToZero, int min, int max, Vec10 *basicUnit, int thd_idx)
{
    for(int i = min; i < max; i++) addPoint((*points)[i], shiftPriorToZero, thd_idx);
}

void AccumulatedSCHessian::addPoint(BackPixelPoint *point, bool shiftPriorToZero, int thd_idx)
{
    int goodResNum = 0;
    for(BackResidual* res : point->backResiduals) if(res->isActiveAndIsGoodNEW) goodResNum++;
    if(goodResNum == 0)
    {
        point->H_depth = 0;
        point->b_Sum = 0;
        point->idepth_hessian=0;
        point->maxRelBaseline=0;
        return;
    }

    //*******累加得到H中与点的深度相关的部分．
    float H_point_idepth = point->H_idepth_idepth_accAF + point->H_idepth_idepth_accLF + point->prior;
    if(H_point_idepth < 1e-10) H_point_idepth = 1e-10;

    point->idepth_hessian = H_point_idepth;     // 线性增量方程H * deltaX = -b 中与点深度相关的H
    point->H_depth = 1 / H_point_idepth;

    point->b_Sum = point->b_idepth_accAF + point->b_idepth_accLF;   // 线性增量方程H * deltaX = -b 中与点深度相关的b
    if(shiftPriorToZero) point->b_Sum += point->prior * point->delta_idepth;   //

    /* 更新acc_HCC，acc_bC */
    Vec4f H_C_idepth = point->H_C_idepth_accAF + point->H_C_idepth_accLF;   /* 不优化相机内参时该项为0 */
    if(T2 == Optimize)
    {
        acc_HCC[thd_idx].upDate(H_C_idepth, H_C_idepth, point->H_depth);
        acc_bC[thd_idx].upDate(H_C_idepth, point->b_Sum * point->H_depth);

        assert(std::isfinite((float)(point->H_depth)));
    } else if(T2 == Fix)
    {
        assert(H_C_idepth == Vec4f::Zero() && "不优化相机内参时H_C_idepth应为0");
    } else { assert(false && "No such mode!"); }

    int frameN2 = frameNum[thd_idx] * frameNum[thd_idx];
    for(BackResidual* res : point->backResiduals)
    {
        if(!res->isActiveAndIsGoodNEW) continue;    // 跳过坏点
        int res_h2t = res->backHostIndex + res->backTargetIndex * frameNum[thd_idx];
        /* 更新acc_HFF */
        for(BackResidual* res2 : point->backResiduals)
        {
            if(!res2->isActiveAndIsGoodNEW) continue;
            if(T1 == Direct)
                acc_HFF[thd_idx][res_h2t + res2->backTargetIndex * frameN2].upDate(res->JpJdF, res2->JpJdF, point->H_depth);
            else if( T1 == Feature)
                acc_HFF_Feature[thd_idx][res_h2t + res2->backTargetIndex * frameN2].upDate(res->JpJdF_Feature, res2->JpJdF_Feature, point->H_depth);
        }

        /* 更新acc_HCF, acc_bF */
        if(T1 == Direct && T2 == Optimize)
        {
            acc_HCF[thd_idx][res_h2t].upDate(res->JpJdF, H_C_idepth, point->H_depth);
            acc_bF[thd_idx][res_h2t].upDate(res->JpJdF, point->H_depth * point->b_Sum);
        } else if(T1 == Direct && T2 == Fix)
        {
            acc_bF[thd_idx][res_h2t].upDate(res->JpJdF, point->H_depth * point->b_Sum);
        } else if(T1 == Feature && T2 == Optimize)
        {
            acc_HCF_Feature[thd_idx][res_h2t].upDate(res->JpJdF_Feature, H_C_idepth, point->H_depth);
            acc_bF_Feature[thd_idx][res_h2t].upDate(res->JpJdF_Feature, point->H_depth * point->b_Sum);
        } else if(T1 == Feature && T2 == Fix)
        {
            acc_bF_Feature[thd_idx][res_h2t].upDate(res->JpJdF_Feature, point->H_depth * point->b_Sum);
        } else { assert(false && "No such mode!");}
    }
}









}
