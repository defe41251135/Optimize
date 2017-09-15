#ifndef ENUM_3000WORLD_H
#define ENUM_3000WORLD_H

#include "util/all_util_include.h"

namespace world3000 {

enum TrackMode{Direct, Feature};     // 直接法或特征点法
enum CameraMode{Optimize, Fix};      // 相机内参优化或者固定
enum ResMode{active, linearized, marginalize};

enum RawPixelPointState {
                    RPS_GOOD=0,
                    RPS_OOB,
                    RPS_OUTLIER,
                    RPS_UNINIT};    // undetermined

enum ResState {Res_INIT,Res_IN, Res_OOB, Res_OUT};

}
#endif // ENUM_3000WORLD_H
