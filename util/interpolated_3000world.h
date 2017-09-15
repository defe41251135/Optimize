#ifndef INTERPOLATED_3000WORLD_H
#define INTERPOLATED_3000WORLD_H

#include "util/all_util_include.h"

/*
 * 数据插值
 */
namespace world3000 {

inline Vec3f getInterpolatedElement33(const Vec3f* const mat, const float x, const float y, const int width)
{
    int x_int = (int)x, y_int = (int)y;
    float delta_x = x - x_int, delta_y = y - y_int;
    float delta_xy = delta_x * delta_y;
    const Vec3f* bp = mat + x_int + y_int * width;

    return delta_xy *  *(const Vec3f*)(bp+1+width) + (delta_y - delta_xy) *  *(const Vec3f*)(bp+width) +
           (delta_x - delta_xy) *  *(const Vec3f*)(bp+1) + (1 - delta_x - delta_y + delta_xy) *  *(const Vec3f*)(bp);
}

}
#endif // INTERPOLATED_3000WORLD_H
