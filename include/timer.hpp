#ifndef TIMER_H
#define TIMER_H
#include <chrono>
#include <iostream>
#include <string>

/** A simple microsecond timer for convenience*/
class timer{

    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    std::string name="";
    bool is_running=false;

    void pretty_print_me(long duration){
      if(name != ""){std::cout << "Time taken by " << name << " is ";
      }else{std::cout << "Time taken is ";}
      std::cout << (float)duration/1.0e6 << " seconds\n";
    }
 
  public:
  
    void begin(std::string name = ""){
      this->name = name;
      is_running=true;
      start = std::chrono::high_resolution_clock::now();
    }
    void pause(){
      stop = std::chrono::high_resolution_clock::now();
      is_running=false;
    }
    void end(){
      stop = std::chrono::high_resolution_clock::now();
      is_running=false;
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      pretty_print_me(duration.count());
    }
   long get_current_time(){
      // If running, get current count, else get count when stopped. Count in microseconds
      std::chrono::time_point<std::chrono::high_resolution_clock> now;
      if(is_running){now = std::chrono::high_resolution_clock::now();
      }else{now = stop;}

      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now - start);
      return duration.count();
    }
    void print_current_time(){pretty_print_me(get_current_time());}
 
};

#endif