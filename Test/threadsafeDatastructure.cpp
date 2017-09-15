#include "util/all_util_include.h"
using namespace std;

/*
 * 使用条件变量的多线程安全的队列模板
 */
template<typename T>
class ThreadSafeQueue
{
private:
    queue<T> dataQueue;
    mutable mutex mtx;      // 为将互斥量应用于常量成员函数以及变量成员函数，互斥量的性质必须是可变的
    condition_variable condVar;
public:
    ThreadSafeQueue() = default;
    ThreadSafeQueue(const ThreadSafeQueue& another)     // 拷贝构造函数
    {
        unique_lock<mutex> lk(another.mtx);
        unique_lock<mutex> lk2(mtx);
        dataQueue = another.dataQueue;
    }
    ThreadSafeQueue& operator=(const ThreadSafeQueue& another) = default;   // 使用默认的拷贝赋值操作符
    ~ThreadSafeQueue() = default;
public:
    void push(T aData)
    {
        unique_lock<mutex> lk(mtx);
        dataQueue.push(aData);
        condVar.notify_one();
    }
    void waitAndPop(T& aData)
    {
        unique_lock<mutex> lk(mtx);
        condVar.wait(lk,[this]{return !dataQueue.empty();});
        aData = dataQueue.front();
        dataQueue.pop();
    }
    std::shared_ptr<T> waitAndPop()
    {
        unique_lock<mutex> lk(mtx);
        condVar.wait(lk,[this]{return !dataQueue.empty();});
        std::shared_ptr<T> res(std::make_shared<T>(dataQueue.front()));
        dataQueue.pop();
        return res;
    }
    bool tryPop(T& aData)
    {
        unique_lock<mutex> lk(mtx);
        if(dataQueue.empty())   return false;
        aData = dataQueue.front();
        dataQueue.pop();
        return true;
    }
    std::shared_ptr<T> tryPop()
    {
        unique_lock<mutex> lk(mtx);
        if(dataQueue.empty())   return std::shared_ptr<T>();
        std::shared_ptr<T> res(std::make_shared<T>(dataQueue.front()));
        dataQueue.pop();
        return res;
    }
    bool empty() const
    {
        unique_lock<mutex> lk(mtx);
        return dataQueue.empty();
    }
};

class empty_stack: std::exception
{

};
