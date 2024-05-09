#include "timer.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath> 


using namespace tools;

#ifndef NO_TIMER

template<typename T>
std::string to_pretty_str( T n )
{
    std::stringstream ss;
    ss.precision(3);
    if( n<1e4 )
    {
        ss << n;
    }
    else if( n<1e6 )
    {
        int thousands = std::floor( n/1000 );
        ss << thousands << ","
           << std::setw(3) << std::setfill('0')
           << n-1000*thousands;
    }
    else if( n<1e9 )
    {
        ss << n/1.e6 << " million";
    }
    else if( n<1e12 )
    {
        ss << n/1.e9 << " billion";
    }
    else
    {
        ss << n/1.e12 << " trillion";
    }
    return ss.str();
}

bool timer::silent = false;

timer::timer(const std::string& label)
    : label(label), running(false), start_time(), elapsed_time()
{}

timer & timer::start() {
    if (!running) {
        start_time = high_resolution_clock::now(); 
        running = true;
    } else if (!silent) {
        std::cout << "Timer is already running." << std::endl;
    }
    return *this;
}


timer & timer::stop() {
    if (!running && !silent) {
        std::cout << "Timer not running." << std::endl;
    }
    elapsed_time += high_resolution_clock::now() - start_time;
    running = false;
    return *this;
}

const double timer::elapsed() {
    return elapsed_time.count();
}

timer & timer::reset() {
    elapsed_time = duration<double>::zero();
    running = false;
    return *this;
}

timer & timer::reset_and_print() {
    print();
    reset();
    return *this;
}

const timer & timer::print() const {
    std::cout << to_pretty_str(elapsed_time.count()) << " seconds elapsed for " << label << std::endl;
    return *this;
}

void timer::set_label(std::string new_label) {
    label = new_label;
}

void timer::silence() {
    silent = !silent;
}

timer::~timer() {
    if (running) {
        stop();
    }
    if (!silent) {
        print();
    }
}

#else

// dummy class
timer::timer(const std::string& label) {}
timer & timer::start() { return *this; }
timer & timer::stop() { return *this; }
duration<double> & timer::elapsed() { return elapsed_time; }
timer & timer::reset() { return *this; }
timer & timer::reset_and_print() { return *this; }
const timer & timer::print() const { return *this; }
void timer::set_label(std::string new_label) {}
void timer::silence() {}
timer::~timer() {}

#endif 