#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>
#include <atomic>

/**
 * ThreadPool - A C++11 thread pool implementation for parallel task execution
 *
 * Features:
 * - Configurable number of worker threads
 * - Task queue with automatic load balancing
 * - Support for task futures and return values
 * - Graceful shutdown with task completion
 * - Thread-safe task submission
 */
class ThreadPool {
private:
    // Worker threads
    std::vector<std::thread> workers;

    // Task queue
    std::queue<std::function<void()>> tasks;

    // Synchronization primitives
    std::mutex queue_mutex;
    std::condition_variable condition;

    // Pool state
    std::atomic<bool> stop;
    std::atomic<size_t> active_tasks;

public:
    /**
     * Constructor - Creates thread pool with specified number of threads
     * @param num_threads Number of worker threads (default: hardware concurrency)
     */
    explicit ThreadPool(size_t num_threads = std::thread::hardware_concurrency())
        : stop(false), active_tasks(0)
    {
        // Ensure at least one thread
        if (num_threads == 0) num_threads = 1;

        // Create worker threads
        for (size_t i = 0; i < num_threads; ++i) {
            workers.emplace_back([this] {
                this->worker_thread();
            });
        }
    }

    /**
     * Destructor - Waits for all tasks to complete and joins threads
     */
    ~ThreadPool() {
        shutdown();
    }

    /**
     * Enqueue a task for execution
     * @param f Function to execute
     * @param args Arguments to pass to function
     * @return Future that will contain the result
     */
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::result_of<F(Args...)>::type>
    {
        using return_type = typename std::result_of<F(Args...)>::type;

        // Create packaged task
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<return_type> result = task->get_future();

        // Add to queue
        {
            std::unique_lock<std::mutex> lock(queue_mutex);

            // Don't allow enqueueing after stopping the pool
            if (stop) {
                throw std::runtime_error("enqueue on stopped ThreadPool");
            }

            tasks.emplace([task]() { (*task)(); });
        }

        // Notify one waiting thread
        condition.notify_one();

        return result;
    }

    /**
     * Wait for all tasks to complete
     */
    void wait() {
        std::unique_lock<std::mutex> lock(queue_mutex);
        condition.wait(lock, [this] {
            return tasks.empty() && active_tasks == 0;
        });
    }

    /**
     * Get number of worker threads
     */
    size_t size() const {
        return workers.size();
    }

    /**
     * Get number of pending tasks
     */
    size_t pending() {
        std::unique_lock<std::mutex> lock(queue_mutex);
        return tasks.size();
    }

    /**
     * Shutdown the thread pool and wait for all threads to finish
     */
    void shutdown() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if (stop) return;  // Already stopped
            stop = true;
        }

        // Wake up all threads
        condition.notify_all();

        // Join all threads
        for (std::thread& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }

private:
    /**
     * Worker thread function - processes tasks from the queue
     */
    void worker_thread() {
        while (true) {
            std::function<void()> task;

            // Wait for a task or stop signal
            {
                std::unique_lock<std::mutex> lock(queue_mutex);

                condition.wait(lock, [this] {
                    return stop || !tasks.empty();
                });

                // Exit if stopping and no tasks remain
                if (stop && tasks.empty()) {
                    return;
                }

                // Get next task
                if (!tasks.empty()) {
                    task = std::move(tasks.front());
                    tasks.pop();
                    active_tasks++;
                }
            }

            // Execute task
            if (task) {
                task();
                active_tasks--;

                // Notify waiting threads that a task completed
                condition.notify_all();
            }
        }
    }
};

/**
 * ResultAggregator - Thread-safe accumulator for gene analysis results
 *
 * Collects results from multiple threads and merges them into global sets/maps
 */
class ResultAggregator {
private:
    std::mutex mutex;
    std::set<std::string>& AllGeneBirth;
    std::set<std::string>& AllGeneDuplication;
    std::map<int, int>& AllGeneLoss;

public:
    ResultAggregator(std::set<std::string>& birth,
                    std::set<std::string>& duplication,
                    std::map<int, int>& loss)
        : AllGeneBirth(birth), AllGeneDuplication(duplication), AllGeneLoss(loss)
    {}

    /**
     * Aggregate results from a single family processing
     * Thread-safe: Uses mutex to protect shared data structures
     */
    void aggregate(const std::set<std::string>& GeneBirth_local,
                   const std::set<std::string>& GeneDuplication_local,
                   const std::map<int, int>& GeneLoss_local)
    {
        std::lock_guard<std::mutex> lock(mutex);

        // Merge birth events
        for (const auto& gene : GeneBirth_local) {
            AllGeneBirth.insert(gene);
        }

        // Merge duplication events
        for (const auto& gene : GeneDuplication_local) {
            AllGeneDuplication.insert(gene);
        }

        // Merge loss events (accumulate counts)
        for (const auto& pair : GeneLoss_local) {
            AllGeneLoss[pair.first] += pair.second;
        }
    }
};

/**
 * FamilyProcessingContext - Contains all data needed to process one gene family
 *
 * Isolates thread-local state to prevent race conditions
 */
struct FamilyProcessingContext {
    // Input data
    int family_id;
    std::set<std::string> genes;

    // Thread-local working variables (prevent race conditions)
    std::map<std::string, int> visited_local;
    std::vector<std::string> group_local;

    // Thread-local result accumulators
    std::vector<std::string> AllTrees_local;
    std::vector<std::vector<std::string>> AllTreeGeneName_local;
    std::set<std::string> GeneBirth_local;
    std::set<std::string> GeneDuplication_local;
    std::map<int, int> GeneLoss_local;

    // Shared read-only data (safe to share)
    const std::string* speciesTree;
    const std::map<std::string, int>* species;
    const std::map<std::string, std::vector<std::string>>* adjacency;
    const std::map<std::pair<std::string, std::string>, double>* edges;
    int S;

    FamilyProcessingContext()
        : family_id(-1), speciesTree(nullptr), species(nullptr),
          adjacency(nullptr), edges(nullptr), S(0)
    {}
};

/**
 * IOTask - Contains data for parallel file I/O operations
 */
struct IOTask {
    int i, j;  // Species indices
    std::string filename;

    // Thread-local results
    std::vector<std::tuple<std::string, std::string, double>> ortholog_pairs;

    IOTask(int si, int sj, const std::string& fname)
        : i(si), j(sj), filename(fname)
    {}
};

#endif // THREADPOOL_H
