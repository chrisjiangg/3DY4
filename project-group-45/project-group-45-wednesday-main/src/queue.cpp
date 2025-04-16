
#include <sstream>
#include <queue>
#include <mutex>
#include "dy4.h"
#include <condition_variable>

class Queue{
    private:
        std::queue<std::vector<real>> my_queue;
        mutable std::mutex m;
        std::condition_variable condition;
    
    public:
        //constructor
        Queue(): my_queue(), m(), condition() {}

        //check if the queue is empty
        bool is_empty() const{
            std::lock_guard<std::mutex> lock(m);
            return my_queue.empty();
        }

        //add an item to the queue
        void enqueue(std::vector<real> item){
            std::lock_guard<std::mutex> lock(m);
            my_queue.push(item);
            condition.notify_one();
        }

        //remove and return the front item
        void dequeue(std::vector<real>& data){
            std::unique_lock<std::mutex> lock(m);
            while(my_queue.empty()){
                condition.wait(lock);
            }
            data = my_queue.front();
            my_queue.pop();
        }

        //return the number of items in the queue
        int size() const{
            std::lock_guard<std::mutex> lock(m);
            return my_queue.size();
        }
};

