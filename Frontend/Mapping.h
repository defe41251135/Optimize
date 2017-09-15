#ifndef MAPPING_H
#define MAPPING_H
#include "util/all_util_include.h"

#include "Frontend/Part.h"
using namespace std;

namespace SLAMSystem {

class Mapping
{
public:
    Mapping()=default;
    void Run();
public:
    // todo
   vector<Frame*> _pallFrames;       // 系统所有的关键帧
   vector<Residual*> _pallResiduals;       // 系统所有的残差项
};
}
#endif // MAPPING_H
