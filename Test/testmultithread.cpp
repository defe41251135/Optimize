#include "util/multithread_3000world.h"
using namespace world3000;
void floorFloats(vector<int>* const optimized, const vector<float>* const toOptimize/*, int min, int max*/)
{
    cout << "GOOD!\n";
//    optimized->resize(toOptimize->size());

//    if(toOptimize->size() < max) max = toOptimize->size();
//    if(min > toOptimize->size()) cout << "arrange error!\n";

//    for(int i = min; i < max; i++)
    for(int i = 0; i < toOptimize->size(); i++)
    {
        (*optimized)[i] = floor((*toOptimize)[i]);
        (*optimized)[i]++;
    }
//    for(auto i:optimized) cout << i << " "; cout << endl;
}




int main()
{


    vector<float> v1;
    for(int i=0;i<10000;i++)
        v1.push_back(i+0.5);
//    for(auto i:v1) cout << i << " "; cout << endl;
    vector<int> v2; v2.resize(v1.size());
    for(auto i:v2) cout << i << " ";
    auto bindFunction = std::bind(&floorFloats, &v2, &v1/*, _1, _2*/);   // 此处应按指针传值而不是引用传值，引用的话值不会成功传递，原因暂时不明
    MultiThread<int> threads;
    threads.reduce(bindFunction, 0, v1.size()-1, 0);
//    floorFloats(v2,v1,1,3);

    cout << endl << endl << "...\n";
    for(auto i:v2) cout << i << " "; cout << endl;

    float wG[3];
    int w = 256;
    for (int level = 1; level < 3; ++ level)
    {
//        wG[level] = w >> level;
        wG[level] = w / pow(2,level);
        cout << wG[level] << "\t";
    }

    vector<float> vv1 = {3,5,1,6,2,11,22,0};

    auto index1 = find(vv1.begin(), vv1.end(), 2);
    cout << distance(vv1.begin(), index1) << endl;

    nth_element(vv1.begin(),vv1.begin()+4,vv1.end());
    for(auto i:vv1) cout << i << " "; cout << endl;

    MatXXf H;
    H.setZero(3,3);
    cout << H <<endl;
    cout << H.cols() << endl;

    // 语法问题
//    {
//        float res_toZeroF_feature = 0.01;
//        float resApprox;

//        Vec6f J_p_d_xi[2], delta; Vec4f J_p_d_cam[2], delta_C; Vec2f J_p_d_idepth; float delta_idepth;
//        J_p_d_xi[1] << 5,4,3,2,1,1;
//        J_p_d_cam[0] << 0.1, 0.2, 0.3, 0.4;
//        J_p_d_cam[1] << 0.5, 0.6, 0.7, 0.8;
//        delta_C << 1,2,3,4;
//        J_p_d_idepth << 2,3;
//        delta_idepth = 0.05;
//        __m128 delta_res = _mm_load_ps((float*)(&res_toZeroF_feature));    // delta_res_0
//        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(J_p_d_xi)), _mm_load_ps((float*)(&delta))));     //delta_res_new
//        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(J_p_d_xi + 1)), _mm_load_ps((float*)(&delta))));
//        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(J_p_d_cam)), _mm_load_ps((float*)(&delta_C))));
//        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(J_p_d_cam + 1)), _mm_load_ps((float*)(&delta_C))));
//        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(&J_p_d_idepth[0])),_mm_load_ps((float*)(&delta_idepth))));
//        delta_res = _mm_add_ps(delta_res, _mm_mul_ps(_mm_load_ps((float*)(&J_p_d_idepth[1])),_mm_load_ps((float*)(&delta_idepth))));
//        _mm_store_ps(((float*)&resApprox), delta_res);

//        cout << resApprox << endl;
//    }

    return 0;
}

