#ifndef TIMER
#define TIMER

#include "Tools.h"

using namespace std::chrono;

class Clock {
 public:
  Clock() : Clock("default_clock") {}
  explicit Clock(const std::string &name) { name_ = name; }

  const char *Start() {
    startTime_ = duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count();
    sprintf(str_, "[\033[1m\033[44;30m%-12s\033[0m] Start...", name_.c_str());
    return str_;
  }

  const char *Count(const char *fmt = "", ...) {
    va_list args;
    char str2[1000];
    va_start(args, fmt);
    vsprintf(str2, fmt, args);
    va_end(args);
    uint64_t end_time = duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count();
    double t = double(end_time - startTime_) / 1e6;
    sprintf(str_, "[\033[1m\033[44;31m%-12s\033[0m] %.6lfs   %s", name_.c_str(), t, str2);
    return str_;
  }

 private:
  char str_[1000]{};
  std::string name_{};
  uint64_t startTime_{};
};

class timeCount{
private:
    uint64_t _startTime;
    uint64_t _endTime;
    uint64_t _lastTime;
    bool _timeState;

public:
//explicit transform
    explicit timeCount(){
        _startTime = _endTime = _lastTime = 0;
        _timeState = false;
    }
    uint64_t GetSysTime(){
        return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    }
    void StartTime(){
        _timeState = true;
        _startTime = _lastTime = GetSysTime();
    }
    void EndTime(){
        _timeState = false;
        _endTime = GetSysTime();
    }
    void UpdateLastTime()
    {
        _lastTime = GetSysTime();
    }
    void Gaptime(const char* str){
        uint64_t curTime = GetSysTime();
        std::cout << CYAN << "timeCount Info" << RESET;
        std::cout << std::setfill('-') << std::setw(35) << str << " : "<< \
        std::fixed << std::setprecision(4) <<  static_cast<double>(curTime - _lastTime)/ 1e6 << "s";
        std::cout << " Total Time Usage: " << static_cast<double>(curTime - _startTime) / 1e6 << "s" << std::endl;
        _lastTime = curTime;
        std::cout << std::setfill(' ');
    }
    double QueryTime(){
        return static_cast<double>(_endTime - _startTime) / 1e6;
    }
    double QueryTimeMilliseSecond(){
        return static_cast<double>(_endTime - _startTime) / 1e3;
    }

};
#endif