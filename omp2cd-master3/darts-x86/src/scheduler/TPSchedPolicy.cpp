#include "TPSchedPolicy.h"
#include "Codelet.h"
#include "MicroScheduler.h"
#include "TPScheduler.h"
#include "tpClosure.h"
#include <cstdlib>
#ifdef TRACE
#include "getClock.h"
#endif

uint64_t TPSleepTime;

#define MUST_USE_GPU(field) ((field)& 0x4)

namespace darts {
inline void cleanupOrSleep(bool scheduled, TPScheduler* myTPS,
    std::queue<tpClosure*>& pendingClosuresToDel,
    std::queue<ThreadedProcedure*>& pendingTPsToDel)
{
    /* If nothing was scheduled or there are too many closures or TPs
* to delete, do it. If there is nothing to delete, just go to
* sleep.*/
    if (pendingClosuresToDel.size() > PENDINGTHRESHOLD) {
        while (pendingClosuresToDel.empty() == false) {
            tpClosure* closureToDel = pendingClosuresToDel.front();
            pendingClosuresToDel.pop();
            delete closureToDel;
        }
    } else if (pendingTPsToDel.size() > PENDINGTHRESHOLD) {
        while (pendingTPsToDel.empty() == false) {
            ThreadedProcedure* TPToDel = pendingTPsToDel.front();
            pendingTPsToDel.pop();
            delete TPToDel;
        }
    } else if (scheduled == false) {
        usleep(myTPS->getSleepTime());
        if (myTPS->getSleepTime() <= TPSleepTime) {
            myTPS->setSleepTime(myTPS->getSleepTime() << 2);
        }
    } else {
        myTPS->setSleepTime(1);
    }
}

bool TPPushFull(TPScheduler* myTPS)
{
/*TPs and closures are deleted on the run */
    while (myTPS->alive()) {
        tpClosure* tempClosure = myTPS->popReadyHead();
        if (!tempClosure) {
            size_t numPeers = myTPS->getNumPeers();
            if (numPeers)
                tempClosure = myTPS->getPeer(rand() % numPeers)->popReadyTail();
        }
        if (tempClosure) {
            tempClosure->factory(tempClosure);
            delete tempClosure;
        }
        Codelet* tempCodelet = myTPS->popCodelet();
        while (tempCodelet && myTPS->alive()) {
            bool fail = true;
            for (size_t i = 0; i < myTPS->getNumSub(); i++) {
                MScheduler* myCDS = static_cast<MScheduler*>(myTPS->getSubScheduler(tempCodelet->getID() % myTPS->getNumSub()));
                if (myCDS->placeCodelet(tempCodelet)) {
                    fail = false;
                    break;
                }
            }
            if (fail == false)
                tempCodelet = myTPS->popCodelet();

            if (tempCodelet)
            {
                ThreadedProcedure* checkTP = tempCodelet->getTP();
                bool deleteTP = (checkTP) ? checkTP->checkParent() : false;
                tempCodelet->fire();
                if (deleteTP)
                    if (checkTP->decRef())
                        delete checkTP;
            }
            tempCodelet = myTPS->popCodelet();
        }
    }
    return true;
}

bool TPRoundRobin(TPScheduler* myTPS)
{
	while (myTPS->alive()) {
		//Check if we have any work in our deque
		tpClosure* tempClosure = myTPS->popReadyHead();

		if (!tempClosure) {
			size_t numPeers = myTPS->getNumPeers();
			if (numPeers)
				tempClosure = myTPS->getPeer(rand() % numPeers)->popReadyTail();
		}

		if (tempClosure) {
			tempClosure->factory(tempClosure);
			delete tempClosure;
			//Get the work ready!
		}
		//Lets do the work!
		Codelet* tempCodelet = myTPS->popCodelet();

		while (tempCodelet && myTPS->alive()) {
			MScheduler* myCDS = static_cast<MScheduler*>(myTPS->getSubScheduler(myTPS->getSubIndexInc()));
			while (!myCDS->placeCodelet(tempCodelet)) {
				myCDS = static_cast<MScheduler*>(myTPS->getSubScheduler(myTPS->getSubIndexInc()));
			}
			tempCodelet = myTPS->popCodelet();
		}
	}
    return true;
}

bool TPSingleLevel(TPScheduler* myTPS)
{
    /*TPs and closures are deleted on the run */
    while (myTPS->alive()) {
        tpClosure* tempClosure = myTPS->popReadyHead();
        if (!tempClosure) {
            size_t numPeers = myTPS->getNumPeers();
            if (numPeers)
                tempClosure = myTPS->getPeer(rand() % numPeers)->popReadyTail();
        }
        if (tempClosure) {
            tempClosure->factory(tempClosure);
            delete tempClosure;
        }
        Codelet* tempCodelet = myTPS->popCodelet();
        while (tempCodelet && myTPS->alive()) {
            ThreadedProcedure* checkTP = tempCodelet->getTP();
            tempCodelet->fire();
            if (checkTP->decRef()) {
                bool isParent = checkTP->checkParent();
                if (isParent)
                    delete checkTP;
                else
                    checkTP->parentTP_ = (ThreadedProcedure*)1;
            }
            tempCodelet = myTPS->popCodelet();
        }
    }
    return true;
}

bool TPTwoLevelSteal(TPScheduler* myTPS)
{
    /*TPs and closures are deleted on the run */
    while (myTPS->alive()) {
        tpClosure* tempClosure = myTPS->popReadyHead();

        if (!tempClosure) {
            size_t numPeers = myTPS->getNumPeers();
            if (numPeers)
                tempClosure = myTPS->getPeer(rand())->popReadyTail();
        }

        if (tempClosure) {
            tempClosure->factory(tempClosure);
            delete tempClosure;
        }
        Codelet* tempCodelet = myTPS->popCodelet();
        while (tempCodelet && myTPS->alive()) {
            tempCodelet->fire();

            ThreadedProcedure* checkTP = tempCodelet->getTP();
            if (checkTP && checkTP->checkParent())
                if (checkTP->decRef())
                    delete checkTP;

            tempCodelet = myTPS->popCodelet();
        }
    }
    return true;
}

bool TPTwoLevelSteal_Sleep(TPScheduler* myTPS)
{
    /*TPs and closures are deleted when scheduler is idle */
    std::queue<ThreadedProcedure*> pendingTPsToDel;
    std::queue<tpClosure*> pendingClosuresToDel;

    while (myTPS->alive()) {
        bool scheduled = false;
        tpClosure* tempClosure = myTPS->popReadyHead();
        if (!tempClosure) {
            size_t numPeers = myTPS->getNumPeers();
            if (numPeers)
                tempClosure = myTPS->getPeer(rand())->popReadyTail();
        }
        if (tempClosure) {
            scheduled = true;
            tempClosure->factory(tempClosure);
            pendingClosuresToDel.push(tempClosure);
        }

        Codelet* tempCodelet = myTPS->popCodelet();
        while (tempCodelet && myTPS->alive()) {
            tempCodelet->fire();
            ThreadedProcedure* checkTP = tempCodelet->getTP();
            bool deleteTP = (checkTP) ? checkTP->checkParent() : false;
            if (deleteTP)
                if (checkTP->decRef())
                    pendingTPsToDel.push(checkTP);

            tempCodelet = myTPS->popCodelet();
        }
        cleanupOrSleep(scheduled, myTPS, pendingClosuresToDel, pendingTPsToDel);
    }
    cleanupOrSleep(false, myTPS, pendingClosuresToDel, pendingTPsToDel);
    return true;
}

bool TPTwoLevelSteal_DeleteAtFinal(TPScheduler* myTPS)
{
    /*TPs and closures are deleted when scheduler is idle */
    std::queue<ThreadedProcedure*> pendingTPsToDel;
    std::queue<tpClosure*> pendingClosuresToDel;

    while (myTPS->alive()) {
        tpClosure* tempClosure = myTPS->popReadyHead();
        if (!tempClosure) {
            size_t numPeers = myTPS->getNumPeers();
            if (numPeers)
                tempClosure = myTPS->getPeer(rand())->popReadyTail();
        }
        if (tempClosure) {
            tempClosure->factory(tempClosure);
            pendingClosuresToDel.push(tempClosure);
        }
        Codelet* tempCodelet = myTPS->popCodelet();
        while (tempCodelet && myTPS->alive()) {
            tempCodelet->fire();
            ThreadedProcedure* checkTP = tempCodelet->getTP();
            bool deleteTP = (checkTP) ? checkTP->checkParent() : false;
            if (deleteTP)
                if (checkTP->decRef())
                    pendingTPsToDel.push(checkTP);

            tempCodelet = myTPS->popCodelet();
        }
    }
    cleanupOrSleep(false, myTPS, pendingClosuresToDel, pendingTPsToDel);
    return true;
}


/* TPGPU scheduler 
 * if metadata & 01 == true, push codelet to GPU MC 
 * else push codelet to CPU MC
 * if GPU/other CPU cores are busy, run codelet directily
 */
bool TPGpuPushFull(TPScheduler* myTPS)
{
/*TPs and closures are deleted on the run */
    while (myTPS->alive()) {
        tpClosure* tempClosure = myTPS->popReadyHead();
        if (!tempClosure) {
            size_t numPeers = myTPS->getNumPeers();
            if (numPeers)
                tempClosure = myTPS->getPeer(rand() % numPeers)->popReadyTail();
        }
        if (tempClosure) {
            tempClosure->factory(tempClosure);
            delete tempClosure;
        }
        Codelet* tempCodelet = myTPS->popCodelet();
		while (tempCodelet && myTPS->alive()) {
		
			bool fail = true;
		
			//std::cout<<"check Codelet is for GPU or not! \n"<<std::endl;	 
			uint32_t stat = tempCodelet->getStatus();
			if(MUST_USE_GPU(stat)) {
            
				MScheduler* myCDS = static_cast<MScheduler*>(myTPS->getSubScheduler(GPU_MC_ID)); //only 1 GPU and GPUID=0
				
				//std::cout<<"check GPU is "<<myCDS->isGpuCapable()<<" (0: not available, 1 available) \n"<<std::endl;	 
			
				if (myCDS->isGpuCapable()){
					myCDS->setUsingGpu(true);
				
					//std::cout<<"GPU is available ! \n"<<std::endl;	 
					myCDS->placeGpuCodelet(tempCodelet); 
					tempCodelet = myTPS->popCodelet();
					fail = false;		
					continue;
				}
			}
			
			//std::cout<<"CPU core prepare!\n"<<std::endl;	
			uint32_t Id=tempCodelet->getID();
			std::size_t numSub = myTPS->getNumSub()-1;// -1 : 1 GPU
			uint32_t subSchedId = Id%numSub +1;
			//std::size_t numSub = myTPS->getNumSub();
			//uint32_t subSchedId = Id%numSub;
			//std::cout<<"numSub:"<<numSub<<",subSchedId:"<<subSchedId<<"\n"<<std::endl;			
			MScheduler* myCDS = static_cast<MScheduler*>(myTPS->getSubScheduler(subSchedId ));
			myCDS->setUsingGpu(false);
			//std::cout<<"push codelet to CPU core !\n"<<std::endl;	
			if (myCDS->placeCodelet(tempCodelet)) {
				fail = false;
				tempCodelet = myTPS->popCodelet();
				continue;
            }
			
            if (fail == true)
                tempCodelet = myTPS->popCodelet();

			//std::cout<<"TPSche run himself!"<<std::endl;	
            if (tempCodelet)
            {
                ThreadedProcedure* checkTP = tempCodelet->getTP();
                bool deleteTP = (checkTP) ? checkTP->checkParent() : false;
                tempCodelet->fire();
                if (deleteTP)
                    if (checkTP->decRef())
                        delete checkTP;
            }
            tempCodelet = myTPS->popCodelet();
        }
    }
    return true;
}

} // Namespace darts
