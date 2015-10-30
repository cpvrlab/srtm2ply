
#ifndef UTILITIES_H
#define UTILITIES_H

#if defined(_MSC_VER)
#define NOMINMAX
#endif

#include <string>
#include <chrono>
#include <string>
#include <iostream>

class ScopedTimer
{
public:
    ScopedTimer(std::string name):
        _name(name),
        _start()
    {
        std::cout << "Starting task \"" << _name << "\"..." << "\n";
        _start = std::chrono::high_resolution_clock::now();
    }

    ~ScopedTimer()
    {
        auto stop = std::chrono::high_resolution_clock::now();
        std::cout << "Task \"" << _name << "\" ended. Elapsed Time: "
                  << (std::chrono::duration_cast<std::chrono::milliseconds>(stop-_start).count()) << "ms."
                  << std::endl;
    }

private:
    std::string _name;
    std::chrono::high_resolution_clock::time_point _start;
};

template<class F>
void time(std::string name, F callback)
{
    ScopedTimer timer(name);
    callback();
}

std::string getCwd();

std::string absolutePath(const std::string &path);

#endif // UTILITIES_H
