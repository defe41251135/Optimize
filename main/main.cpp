#include "util/all_util_include.h"
using namespace std;
using namespace world3000;


class Back
{
public:
    Back()
        : Id(2)
    {

    }

    vector<float> timeStamps;
    int Id;
};
class Front
{
public:
    Front()
        : Id(1)
    {

    }

public:
    Back* backData;
    int Id;
};
class System
{
public:
    System()
    {
        front = new Front();
//        back = new Back();
//        front->backData = back;

        front->backData = new Back();
//        try
//        {vector<Front*> fronts;
//        fronts.reserve(20);}
//        catch (exception e)
//        {
//            cout << "error!\n";
//        }
        fun1();
    }
    ~System()
    {
        delete front->backData;
        front->backData = NULL;
        delete front;
        front = NULL;

//        cout << back->Id << endl;
        cout << front->backData->Id << endl;
        cout << front->Id << endl;
    }
    void fun1()
    {
        vector<Front*> allFronts;
//        allFronts.resize(20,NULL);
        allFronts.reserve(20);
        allFronts.push_back(front);
    }
    void add1To10000(int first, int last, int result)
    {
        while(first <= last)
        {
            result += first;
            first++;
        }
    }
public:
    Front* front;
    Back* back;
};


void foo()
{
    MatXX H = Mat44::Constant(0);
    cout << H <<endl;
    const int a=6;
    const int b=6;
    H.conservativeResize(a,b);
//    H = Eigen::Matrix<double,a,b>::Constant(1);
    cout << H << endl;
}

int mainmmm()
{
//    foo();

    System* sys = new System();
//    sys->fun1();

//    float rmse = optimize(5);




}


