#include "FullSystem/FullSystem.h"

namespace SLAMSystem
{



FullSystem::FullSystem(const string& dataFilename)
{
    tracking = new Tracking();
    EF = new EnergyFunction();
    tracking->setEnergyFunction(EF);
    EF->setTracking(tracking);

    tracking->camera = new Camera();
    tracking->camera->backData = new BackCamera(tracking->camera);

    EF->threadReduce1 = &tracking->threadReduce1;

    mapping = new Mapping();
    thdMapping = new thread(&Mapping::Run, mapping);

    // viewer todo 8.14



}

FullSystem::~FullSystem()
{
    // todo
}




}
