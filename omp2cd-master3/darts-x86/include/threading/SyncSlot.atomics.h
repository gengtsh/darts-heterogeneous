#pragma once
//#include <cstdint>
#include <iostream>
#include <stdint.h>
//#include "Atomics.h"
#include <atomic>

namespace darts
{

  /*
	 * Class: SyncSlot
   * This class contains to ints which are used as a counter and reset. This is to allow
   * a codelet to run.
	 * 
	 * See Also:
	 * <Atomics>
	*/
  
  class SyncSlot
  {
  private:
      std::atomic<unsigned int> counter_;
      std::atomic<unsigned int> reset_;
  public:
    SyncSlot(uint32_t dep, uint32_t res):
    counter_(dep),
    reset_(res) { }
    
    void
    initSyncSlot(uint32_t dep, uint32_t res)
    {
        counter_ = dep;
        reset_ = res;
    }
    
    //dec the counter
    bool 
    decCounter(void)
    {
        return counter_--;
        //return (1==Atomics::fetchSub(counter_, 1U));
    }
    
    //inc the counter
    void
    incCounter(void)
    {
        counter_ += 1;
        reset_   += 1;
        //Atomics::fetchAdd(counter_, 1U);
        //Atomics::fetchAdd(reset_, 1U);
    }
    
    //resets the counter
    void 
    resetCounter(void)
    {
        counter_ = reset_;
        //return Atomics::boolcompareAndSwap(counter_,0U,reset_);
    }
    
    uint32_t getCounter(void) const{
        return counter_;
    }
    
    //returns the if the counter has reached zero
    bool ready(void) const{
        return (counter_ == 0);
    }  
  };

} // namespace darts
