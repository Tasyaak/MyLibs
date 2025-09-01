#pragma once
#include <cassert>
#include "core/config.hpp"

template <typename T>
class RingBuffer
{
private:
    T* buffer;
    size_t capacity;
    size_t start;
    size_t end;
    bool full;

public:
    RingBuffer() : buffer(nullptr), capacity(0), start(0), end(0), full(true) {}
    RingBuffer(size_t size) : buffer(new T[size]), capacity(size), start(0), end(0), full(false) {}
    RingBuffer(T* a, size_t n) : capacity(n), start(0), end(0), full(true)
    {
        assert(a != nullptr && "Pointer cannot be nullptr");

        buffer = new T[n];
        std::copy(a, a + n, buffer);
    }

    void add(T value)
    {
        assert(!empty() && "RingBuffer is empty");

        buffer[end] = value;
        if (full)
            start = (start + 1) % capacity;
        end = (end + 1) % capacity;
        full = (end == start);
    }

    size_t size() const
    {
        assert(!empty() && "RingBuffer is empty");

        if (full)
            return capacity;
        if (end >= start)
            return end - start;
        return capacity - start + end;
    }

    bool empty() const { return buffer == nullptr; }

    T operator[](size_t index) const
    {
        assert(index < size() && "RingBuffer index out of range");

        return buffer[(start + index) % capacity];
    }
};