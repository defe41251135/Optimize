#include "util/all_util_include.h"
using namespace std;
/* parallel accumulate
 * 实现一个并行版的std::accumulate
 */


/*
 * 负责一个块的累加任务
 */
template<typename Iterator, typename T>
void accumulateBlock(Iterator first, Iterator last, T& blockResult)
{
    blockResult += accumulate(first, last, blockResult);
}

/* 多线程并行累加计算,并行多个块的累加任务，将每个块的累加结果累加，得到最终的累加结果
 * 采用了多线程分层级的处理方式.
 */
template<typename Iterator, typename T>
T parallelAccumulate(Iterator first, Iterator last, T init)
{
    const unsigned long length = std::distance(first, last);
    if(length == 0) return init;    // 输入范围为空，则返回初始值init

    const unsigned long minTaskPerThread = 25;
    const unsigned long maxThreadsNum = (length + minTaskPerThread - 1) / minTaskPerThread; // 范围内元素的总数量除以线程最小任务数，确定启动线程的最大数量
    const unsigned long hardwareThreadsNum = thread::hardware_concurrency();

    const unsigned long threadsNum = min(hardwareThreadsNum == 0 ? 2 : hardwareThreadsNum, maxThreadsNum);
    const unsigned long blockSize = length/threadsNum;

    vector<T> results(threadsNum);
    vector<thread> threads(threadsNum - 1);

    Iterator blockBegin = first;
    for(unsigned long i = 0; i < threadsNum-1; i++)
    {
        Iterator blockEnd = blockBegin;
        std::advance(blockEnd,blockSize);

        threads[i] = thread(&accumulateBlock<Iterator, T>, blockBegin, blockEnd, ref(results[i]));

        blockBegin = blockEnd;
    }
    accumulateBlock<Iterator, T>(blockBegin, last, results[threadsNum - 1]);

    //sum them all
    for_each(threads.begin(), threads.end(), mem_fn(&thread::join));

    T finalResult = accumulate(results.begin(), results.end(), init);
    return finalResult;
}


int main221()
{
    vector<int> vec;
    for(unsigned i = 0; i<200; i++) vec.push_back(i+1);
    vector<int>::iterator begin = vec.begin();
    vector<int>::iterator end = vec.end();
    int result = parallelAccumulate<vector<int>::iterator, int>(begin, end, 0);
    cout <<result << endl;

}
