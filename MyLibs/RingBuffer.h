#pragma once

template <typename T>
class RingBuffer
{
private:
    T* buffer;
    int capacity;
    int start;
    int end;
    bool full;

public:
    inline RingBuffer() : buffer(nullptr), capacity(0), start(0), end(0), full(true) {}
    inline RingBuffer(int size) : buffer(new T[size]), capacity(size), start(0), end(0), full(false) {}
    inline RingBuffer(T* a, int n) : capacity(n), start(0), end(0), full(true)
    {
        buffer = new T[n];
        for (int i = 0; i < n; ++i)
            buffer[i] = a[i];
    }

    inline void add(T value)
    {
        buffer[end] = value;
        if (full)
            start = (start + 1) % capacity;
        end = (end + 1) % capacity;
        full = (end == start);
    }

    inline int size() const
    {
        if (full)
            return capacity;
        if (end >= start)
            return end - start;
        return capacity - start + end;
    }

    inline bool empty() const { return buffer == nullptr; }

    inline T operator[](int index) const
    {
        if (index < 0 || index >= size())
            throw std::out_of_range("Index out of range");
        return buffer[(start + index) % capacity];
    }
};