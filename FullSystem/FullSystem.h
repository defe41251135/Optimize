#ifndef FULLSYSTEM_H
#define FULLSYSTEM_H
#include "Frontend/Tracking.h"
#include "Frontend/Mapping.h"
#include "Frontend/Part.h"
using namespace std;


namespace SLAMSystem {

//enum TrackMode{Direct, Feature};     // 直接法或特征点法
//enum CameraMode{Optimize, Fix};      // 相机内参优化或者固定
//enum ResMode{active, linearized, marginalize};

class EnergyFunction;
class Mapping;

class FullSystem
{
public:
    FullSystem(const string& datasetFilename);
    FullSystem()=default;
    ~FullSystem();
//    RunSystem();

public:
    Tracking* tracking;
    EnergyFunction* EF;     //

    Mapping* mapping;
    thread* thdMapping;


};

}
#endif // FULLSYSTEM_H
