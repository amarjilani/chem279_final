#pragma once
#include <chrono>
#include <string>

namespace tools {
    using namespace std::chrono;
    
    class timer {
        public:
            explicit timer(const std::string & label);
            ~timer();
            static void silence();
            void set_label(std::string label);
            const timer & print() const;
            timer & start();
            timer & stop();
            const double elapsed();
            timer & reset();
            timer & reset_and_print();
        private:
            std::string label;
            bool running;
            high_resolution_clock::time_point start_time;
            duration<double> elapsed_time;
            static bool silent;
    };
}