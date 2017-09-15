#include "util/all_util_include.h"

// 多线程的生产者－－消费者模型
// 实时不稳定
class ProducterConsumer
{
public:
    ProducterConsumer()
        : produceDone(false), consumeDone(false)
    {}
private:
    queue<int> dataQueue;                       // 在两个线程之间传递数据的队列
    mutex mtx;                                  // 保护该队列的互斥量
    condition_variable condVarConsumer;                 // 条件变量,用于唤醒消费者线程
    condition_variable condVarProductor;                // 条件变量,用于唤醒生产者线程

    bool produceDone;                           // 生产者线程结束标志
    bool consumeDone;                           // 消费者线程结束标志

public:
    void produce()
    {

        for(unsigned i = 0; i < 10; i++)
        {
            this_thread::sleep_for(chrono::seconds(1));
            int data = rand()%10;                       // 准备数据
            unique_lock<mutex> lk(mtx);                 // 数据准备好后，对互斥量上锁
            dataQueue.push(data);                       // 操作被保护数据
            cout << "Productor: " << data << endl;
            condVarConsumer.notify_one();                       // 唤醒一个等待condVar的线程
            condVarProductor.wait(lk, [this]{return true;});
        }
        produceDone = true;                             // 生产者线程结束
    }
    void consume()
    {
        while(!produceDone)
        {
            unique_lock<mutex> lk2(mtx);                // 对互斥量上锁
//            while(!notified)
//            {
//                condVar.wait(lk2);                      // wait()会检查条件(通过调用所提供的lambda函数)，如果条件不满足，那么wait()将解锁互斥量，并将该线程设置为睡眠状态．
//                                                        // 当准备数据的线程调用notify_one()通知条件变量时，处理数据的线程被唤醒，重新对互斥量上锁，并对条件再次检查，
//                                                        // 如果条件满足，那么从wait()返回，继续执行处理数据的线程;如果条件不满足，继续对互斥量解锁，线程进入睡眠状态
//            }
            condVarConsumer.wait(lk2, [this]{return !dataQueue.empty();});
            while(!dataQueue.empty())
            {
                int data = dataQueue.front();
                dataQueue.pop();
                cout << "Consumer: " << data+1 << endl; // 执行消费者线程的操作

            }
            condVarProductor.notify_one();
        consumeDone = true;
        }
    }
};

int useProducterConsumer()
{
    ProducterConsumer* pAndC = new ProducterConsumer();
    thread t1(&ProducterConsumer::produce, pAndC);
    thread t2(&ProducterConsumer::consume, pAndC);
    t1.join();
    t2.join();
}

//int main()
//{

//    useProducterConsumer();

//}

void null(){
//class ConditionVariableExample
//{
//private:
//    queue<int> dataQueue;
//    mutex mtx;
//    condition_variable dataCond;

//    bool notified;
//public:
//    bool moreDataToPrepare(vector<int>& vec)
//    {
//        return !vec.empty();
//    }

//    int prepareData(vector<int>& vec)
//    {
//        int rt = vec.back();
//        vec.pop_back();
//        return rt;
//    }

//    void process(int data)
//    {
//        cout << data+1 << " ";
//    }

//    /* 准备数据的线程 */
//    void dataPreparationThread(vector<int>& vec)
//    {
//        while(moreDataToPrepare(vec))
//        {
//            const int data = prepareData(vec);     // 准备数据
//            std::unique_lock<mutex> lk(mtx);
//            dataQueue.push(data);
//            notified = true;
//            dataCond.notify_one();              // 对等待的线程(如果有等待线程)进行通知

//        }
//    }
//    /* 处理数据的线程 */
//    void dataProcessingThread()
//    {
//        while(true)
//        {
//            std::unique_lock<mutex> lk(mtx);    // 对互斥量上锁

//            while (!notified) {
//                dataCond.wait(lk)
//            }
//            dataCond.wait(lk, tempbool/*[]{return !dataQueue.empty();}*/);  // 调用std::condition_variable的成员函数wait(),传递了一个锁lk和一个lambda表达式
//            // wait()会检查条件(通过调用所提供的lambda函数)，如果条件不满足，那么wait()将解锁互斥量，并将该线程设置为睡眠状态．
//            // 当准备数据的线程调用notify_one()通知条件变量时，处理数据的线程被唤醒，重新对互斥量上锁，并对条件再次检查，
//            // 如果条件满足，那么从wait()返回，继续执行处理数据的线程;如果条件不满足，继续对互斥量解锁，线程进入睡眠状态
//            int data = dataQueue.front();       // 取出队列中的第一个数据
//            dataQueue.pop();
//            lk.unlock();                        // 互斥量解锁
//            process(data);                      // 处理数据
////            if(isLastInt(data)) break; //

//        }
//    }
//};

//int main()
//{
//    vector<int> vec;
//    for(unsigned i = 0; i < 100; i++) vec.push_back(i);
//    ConditionVariableExample* example = new ConditionVariableExample();
//    thread t1(&ConditionVariableExample::dataPreparationThread, example, vec);
//    thread t2(&ConditionVariableExample::dataProcessingThread, example);
//    t1.join();
//    t2.join();
//}

}



