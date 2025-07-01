#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include "Defines.h"
struct CompareTasks {
    bool operator()(const std::pair<ui, std::function<void()>>& lhs, const std::pair<ui, std::function<void()>>& rhs) {
        return lhs.first > rhs.first;
    }
};
class ThreadPool {
public:
    ThreadPool(size_t threads);
    ~ThreadPool();

    template<class F, class... Args>
    auto enqueue(ui p, F&& f, Args&&... args) -> std::future<typename std::invoke_result<F, Args...>::type>;
    // auto enqueue(F&& f, Args&&... args) -> std::future<typename std::invoke_result<F, Args...>::type>;

private:
    std::vector<std::thread> workers;
    // std::queue<std::function<void()>> tasks;
    std::priority_queue<std::pair<ui, std::function<void()>>, std::vector<std::pair<ui, std::function<void()>>>, CompareTasks>  tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

inline ThreadPool::ThreadPool(size_t threads) : stop(false) {
    for(size_t i = 0; i < threads; ++i){
        workers.emplace_back(
            [this] {
                for(;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                            [this]{ return this->stop || !this->tasks.empty(); });
                        if(this->stop && this->tasks.empty())
                            return;
                        // task = std::move(this->tasks.front());
                        task = std::move(this->tasks.top().second);
                        this->tasks.pop();
                    }

                    task();
                }
            }
        );
        // std::cout << "Thread " << i << " started.\n";
    }
}

template<class F, class... Args>
auto ThreadPool::enqueue(ui gra_pri,F&& f, Args&&... args) -> std::future<typename std::invoke_result<F, Args...>::type> {
    using return_type = typename std::invoke_result<F, Args...>::type;

    auto task = std::make_shared< std::packaged_task<return_type()> >(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace(gra_pri,[task](){ 
            // std::cout << "Task started by thread: " << std::this_thread::get_id() << "\n";
            (*task)(); });
        // tasks.emplace([task](){ 
        //     // std::cout << "Task started by thread: " << std::this_thread::get_id() << "\n";
        //     (*task)(); });
    }
    condition.notify_one();
    return res;
}

inline ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(std::thread &worker: workers)
        worker.join();
}

#endif

