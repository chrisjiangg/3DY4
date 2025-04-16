#ifndef DY4_QUEUE_H
#define DY4_QUEUE_H

#include <stdexcept>
#include <vector>
#include <string>
#include "dy4.h"

class Queue {
    private:
        std::vector<std::vector<real>> items; //uses vector to store items

    public:
        Queue(); //constructor
        bool is_empty() const; //check if the queue is empty
        void enqueue(std::vector<real>  item); //add item to queue
        std::vector<real>  dequeue(); //remove and return the front item
        size_t size() const; //return the number of items in the queue
};

#endif