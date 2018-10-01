
#include <unistd.h>
#include <climits>
#include <cstdio>
#include <cerrno>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <vector>
#include "DARTS.h"
#include "tmp1.h"
using std::vector;

int main(void){

    Runtime *rt;
    
	ThreadAffinity affin(4, 1, COMPACT, TPROUNDROBIN, MCSTANDARD);

	if (affin.generateMask()){

        rt = new Runtime(&affin);
        
        rt->run(launch<TMP_TP>(2,4,&Runtime::finalSignal));

    }
    std::cout<<"finish!"<<std::endl;
    
    delete rt;
    return 0;

}
