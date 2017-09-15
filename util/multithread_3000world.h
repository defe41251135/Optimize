#ifndef MULTITHREAD_3000WORLD_H
#define MULTITHREAD_3000WORLD_H

#include "util/include_3000world.h"
#include "util/typedef_3000world.h"
#include "util/settings_3000world.h"




using namespace std;


namespace world3000 {

template<typename T> class MultiThread
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    MultiThread()
        : isRunning(true)
    {
        _inputFunction = std::bind(&MultiThread::_inputFunctionDefault, this, _1, _2, _3, _4);
        for(int thdIdx = 0; thdIdx < ThreadNum; thdIdx++)
        {
            _isDone[thdIdx] = false;
            _gotOne[thdIdx] = true;
            _workerThreads[thdIdx] = thread(&MultiThread::workerLoop, this, thdIdx);
        }
    }
    ~MultiThread()
    {
        isRunning = false;
        _exMutex.lock();
        _todoSignal.notify_all();
        _exMutex.unlock();

        for(int thdIdx = 0; thdIdx < ThreadNum; thdIdx++)
        {
            _workerThreads[thdIdx].join();
        }

        cout << "Destructor MultiThread!\n";
    }

    void reduce(function<void(int, int, T*, int)> inputFunction, int first, int last, int stepSize = 0) // 函数参数：要并行处理的多线程函数，函数参数的起始位，终止位，每个线程处理的步长
    {
        memset(&basicUnit, 0, sizeof(basicUnit));

        if(stepSize == 0)
            stepSize = ((last-first)+ThreadNum-1)/ThreadNum;

        unique_lock<mutex> lock(_exMutex);

        _inputFunction = inputFunction;
        _nextIdx = first;
        _maxIdx = last;
        _stepSize = stepSize;
        for(int thdIdx = 0; thdIdx < ThreadNum; thdIdx++)
        {
            _isDone[thdIdx] = false;
            _gotOne[thdIdx] = false;
        }
        _todoSignal.notify_all();

        while(1)
        {
            _doneSignal.wait(lock);

            bool allDone = true;
            for(int thdIdx = 0; thdIdx < ThreadNum; thdIdx++)
            {
                allDone = allDone && _isDone[thdIdx];
            }
            if(allDone)
            {   cout << "ALL DONE!\n";
                break;
            }
        }

        // 还原
        _nextIdx = 0;
        _maxIdx = 0;
        _inputFunction = std::bind(&MultiThread::_inputFunctionDefault, this, _1, _2, _3, _4);

    }

public:
    T basicUnit;        // 多线程处理的基本单元

private:
    void workerLoop(int idx)        // 一个工人的工作循环
    {
        unique_lock<mutex> lock(_exMutex);

        while(isRunning)
        {
            //try to get something to do.
            int todo = 0;
            bool gotSomething = false;
            if(_nextIdx < _maxIdx)
            {
                // got something!
                todo = _nextIdx;
                _nextIdx += _stepSize;
                gotSomething = true;
            }
            // if got something: do it (unlock in the meantime)
            if(gotSomething)
            {
                lock.unlock();
                assert(_inputFunction != 0);
                cout << "worker " << idx << " is working..\n";
                T t;
                memset(&t, 0, sizeof(t));
                _inputFunction(todo, min(todo+_stepSize, _maxIdx), &t, idx);  // worker idxth is working
                _gotOne[idx] = true;

                lock.lock();

                basicUnit += t;
            }
            // otherwise wait on signal, releasing lock in the meantime.
            else
            {
                _isDone[idx] = true;
                cout << "worker " << idx << " is waiting..\n";
                _doneSignal.notify_all();
                _todoSignal.wait(lock);
            }
        }
    }

private:
    thread _workerThreads[ThreadNum];      // n个线程
    bool _isDone[ThreadNum];
    bool _gotOne[ThreadNum];

    mutex _exMutex;
    condition_variable _todoSignal;
    condition_variable _doneSignal;

    int _nextIdx;
    int _maxIdx;
    int _stepSize;

    bool isRunning;

    function<void(int,int,T*,int)> _inputFunction;
    void _inputFunctionDefault(int i, int j, T* t, int tid)
    {
//        assert(false);
    }
};
}
#endif // MULTITHREAD_3000WORLD_H
