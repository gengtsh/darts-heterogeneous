#include "Spmv.h"
#include "conf.h"


// ****************************************************************************
// Function: addBenchmarkSpecOptions
//
// Purpose:
//   Add benchmark specific options parsing.
//
// Arguments:
//   op: the options parser / parameter database
//
// Programmer: Lukasz Wesolowski
// Creation: June 21, 2010
// Returns:  nothing
//
// ****************************************************************************
void addBenchmarkSpecOptions(OptionParser &op)
{
    
    //op.addOption("iterations", OPT_INT, "100", "Number of SpMV iterations "
    //             "per pass");
    op.addOption("iterations", OPT_INT, "1", "Number of SpMV iterations "
                 "per pass");
    //op.addOption("mm_filename", OPT_STRING, "random", "Name of file "
    //             "which stores the matrix in Matrix Market format");
    //op.addOption("mm_filename", OPT_STRING, "1138_bus.mtx", "1138_bus.mtx"
    //             "which stores the matrix in Matrix Market format");
    
    //op.addOption("mm_filename", OPT_STRING, "FullChip.mtx", "FullChip.mtx"
    //             "which stores the matrix in Matrix Market format");
    
    op.addOption("maxval", OPT_FLOAT, "10", "Maximum value for random "
                 "matrices");
}

bool sortbyseclg(const pair<int,int>&a,const pair<int,int>&b){
    return (a.second > b.second);
}

bool sortbyfstle(const pair<int,int>&a,const pair<int,int>&b){
    return (a.first < b.first);
}

struct Thd{
    int first;
    int second;
    int third;
};


bool sortbythdlg(const Thd &a,const Thd &b){
    return (a.third > b.third);
}


bool sortle(const int a, const int b){
    return (a<b);
}


void preProcessNonZeroes(vector<pair<int,int>> &dst,int *src, int st, int ed){

    for(int i=st; i<ed;++i){
        dst.push_back(make_pair(i,src[i+1]-src[i]));
    }
    sort(dst.begin(),dst.end(),sortbyseclg);

}


void preProcessNonZeroes(int *dst,int *src, int st, int ed){

    for(int i=st; i<ed;++i){
        dst[i]=src[i+1]-src[i];
    }
}


void preProcessNonZeroes(vector<pair<int,int>> &dstVec,int *dstArray,int *src, int st, int ed){

    for(int i=st; i<ed;++i){
        int tt=src[i+1]-src[i];
        dstArray = tt;
        dstVec.push_back(make_pair(i,tt));
    }
    sort(dstVec.begin(),dstVec.end(),sortbyseclg);
    //int pt = (ed-st)*0.995;
    //dstVec.erase(dstVec.begin()+pt,dstVec.end());

}


//calc groups , the total number is less than num; in one group, the diffs is less than pt (percent), then cluster in each group up to the total value is thd since the overhead of create a codelet 
//groups<int,int>: start, end
void calcGroups(int *refArray,vector<vector<pair<int,int>>> &groups,vector<pair<int,int>> &srcVec, int num, double pt,int thd =0){
    int idx=0;
    int refIdx = 0;
    int nextIdx = 0;
    int ref = 0;
    int next = 0;
   
    vector<pair<int,int>> all;
    while (idx<num){
        vector<int> init;
        vector<pair<int,int>> group;
#ifdef DARTS_DEBUG
        //std::cout<<"get initial group <idx, value> and sort init from less idx to bigger idx"<<std::endl;
#endif
        refIdx = srcVec[idx].first;
        ref    = srcVec[idx].second;
        init.push_back(refIdx);
        ++idx;
        nextIdx = srcVec[idx].first;
        next    = srcVec[idx].second;

        while (next > pt*ref && idx < num){
            init.push_back(nextIdx);
            ++idx;
            nextIdx = srcVec[idx].first;
            next    = srcVec[idx].second;
        } 
        sort(init.begin(),init.end(),sortle); 
       
#ifdef DARTS_DEBUG
        //std::cout<<"further process group: cluster them if possible"<<std::endl;
#endif
        int val = 0;
        int start = init[0];
        int end = 0;
        for (int i=0; i<init.size();++i){
            end  = init[i] +1;
            int tmp  = refArray[end]-refArray[start];
            
            //check where priority value is located in the lower interval
            if(thd > 0){
                for(int j=0; j<all.size();++j){
                    int pStart = all[j].first ;
                    int pEnd   = all[j].second;
                    if((pStart > start) && (pEnd < end)){
                        tmp -= refArray[pEnd]-refArray[pStart];
                    }
                }
            }
            val += tmp; 
            if(i == init.size()-1){
                group.push_back(make_pair(start,end));
                all.push_back(make_pair(start,end));
            }else if(val > thd){
                group.push_back(make_pair(start,end));
                all.push_back(make_pair(start,end));
                start = init[i+1];
                val = 0;
            }

        }
        
        groups.push_back(group);
    }

#ifdef DARTS_DEBUG 
    for(int i = 0; i< groups.size(); ++i){
       for(int j=0; j<groups[i].size();++j){
            std::cout<<groups[i][j].first<<","<<groups[i][j].second<<std::endl;
       } 
       std::cout<<std::endl;
    }
#endif
}

void memsetVal(int *src, int val, int st, int ed){
    int sz = ed - st;
//    memset(src+st, val, sz*sizeof(int));
    std::fill(src+st,src+ed,val);

}


int buildGroups(vector<GroupAttr> &groups,int *src, int st, int ed, int topn, double pt, int thd=0){
    //get nonZeroesPerRowVec = valVec
    vector<pair<int,int>> valVec;
    
    for(int i=st; i<ed;++i){
        int tt=src[i+1]-src[i];
        valVec.push_back(make_pair(i,tt));
    }
    sort(valVec.begin(),valVec.end(),sortbyseclg);
    //int pt = (ed-st)*0.995;
    //valVec.erase(valVec.begin()+pt,valVec.end());

    //calc groups
    int idx=0;
    int refIdx = 0;
    int nextIdx = 0;
    int ref = 0;
    int next = 0;
    int gId = 0; 
    vector<pair<int,int>> all;
    while (idx<topn){
        vector<int> init;
        
        //"get initial group <idx, value> and sort init from less idx to bigger idx";
        refIdx = valVec[idx].first;
        ref    = valVec[idx].second;
        init.push_back(refIdx);
        ++idx;
        nextIdx = valVec[idx].first;
        next    = valVec[idx].second;

        while (next > pt*ref && idx < topn){
            init.push_back(nextIdx);
            ++idx;
            nextIdx = valVec[idx].first;
            next    = valVec[idx].second;
        } 
        sort(init.begin(),init.end(),sortle); 
       
        //"further process group: cluster them if possible";
        //memset src from start to end to zero when group is built. it is because once the src were process by high level group, the lower group will no need processagain.
        
        int val = 0;
        int start = init[0];
        int end = 0;
        for (int i=0; i<init.size();++i){
            end  = init[i] +1;
            int tmp  = src[end]-src[start];
            
            //check where priority value is located in the lower interval
            if(thd > 0){
                for(int j=0; j<all.size();++j){
                    int pStart = all[j].first ;
                    int pEnd   = all[j].second;
                    if((pStart > start) && (pEnd < end)){
                        tmp -= src[pEnd]-src[pStart];
                    }
                }
            }
            
            val += tmp; 
            if(i == init.size()-1){
                groups.push_back({start,end,gId});
                all.push_back(make_pair(start,end));
                //memsetVal(src,0,start,end);
            }else if(val > thd){
                groups.push_back({start,end,gId});
                all.push_back(make_pair(start,end));
                //memsetVal(src,0,start,end);
                start = init[i+1];
                val = 0;
            }

        }
        
        ++gId ;
    }

#ifdef DARTS_DEBUG 
       for(int j=0; j<groups.size();++j){
            std::cout<<groups[j].start<<","<<groups[j].end<<", "<<groups[j].groupId<<std::endl;
       } 
       std::cout<<std::endl;
#endif
    return gId;
}


int buildGroups_v2(vector<GroupAttr> &groups,int *src, int st, int ed, int topn, double pt, int thd=0){
    //get nonZeroesPerRowVec = valVec
    vector<pair<int,int>> valVec;
    //int *tA = new int [ed-st];   
    for(int i=st; i<ed;++i){
        int tt=src[i+1]-src[i];
        valVec.push_back(make_pair(i,tt));
    //    tA[i]=tt;
    }
    sort(valVec.begin(),valVec.end(),sortbyseclg);
    //int pt = (ed-st)*0.995;
    //valVec.erase(valVec.begin()+pt,valVec.end());
    //string csv = "xx.csv";
    //ofstream out;
    //out.open(csv, std::ofstream::out);
    //for(int i=0; i<valVec.size();++i){
    //    out<<valVec[i].second<<",";
    //    if(i%10 == 0)
    //        out<<"\n";
    //}
    //out.close();
    //double mean = calc_mean(tA,st,ed);
    //double sd   = calc_sd(tA,st,ed);
    //double cv   = calc_cv(tA,st,ed);
    //std::cout<<"st: "<<st<<",ed: "<<ed<<",mean: "<<mean<<",sd: "<<sd<<",cv: "<<cv<<", max: "<<valVec[0].second<<std::endl;
    //delete []tA;

   //valVec only keep topn  
    valVec.erase(valVec.begin()+topn, valVec.end());
    
    sort(valVec.begin(),valVec.end(), sortbyfstle);
  
#ifdef DARTS_DEBUG 
    for(int i=0; i<valVec.size();++i){
        std::cout<<"valVec["<<i<<"].first = "<< valVec[i].first<<std::endl;
    }

#endif
    vector<Thd> gp;
    int start =0;
    int end   = 0;
    for(int i=0;i<valVec.size();++i){

            start = valVec[i].first;
            end   = valVec[i].first+1;

            while(i+1<valVec.size() && valVec[i+1].first-valVec[i].first<thd ){
                end = valVec[i+1].first+1;
                ++i;
            }

            gp.push_back({start,end,src[end]-src[start]});
    }

#ifdef DARTS_DEBUG 
    for(int i=0; i<gp.size();++i){

        std::cout<<"gp["<<i<<"].first = "<< gp[i].first<<", gp["<<i<<"].second = "<<gp[i].second<<std::endl;
    }

#endif
    // further calc groups
    
    sort(gp.begin(),gp.end(),sortbythdlg);


#ifdef DARTS_DEBUG 
    for(int i=0; i<gp.size();++i){

        std::cout<<"gp["<<i<<"].first = "<< gp[i].first<<", gp["<<i<<"].second = "<<gp[i].second<<", gp["<<i<<"].third = "<<gp[i].third<<std::endl;
    }

#endif

    int idx=0;
    int ref = 0;
    int next = 0;
    int gId = 0; 
    int sz = gp.size();
    while(idx<sz){
        
        groups.push_back({gp[idx].first,gp[idx].second,gId});

        ref  = gp[idx].third;
        ++idx;
        next = gp[idx].third;
        while(next > pt*ref && idx<sz){
            groups.push_back({gp[idx].first,gp[idx].second,gId});
            ++idx;
            next = gp[idx].third;
        }
        ++gId;
    }
    


#ifdef DARTS_DEBUG 
       for(int j=0; j<groups.size();++j){
            std::cout<<groups[j].start<<","<<groups[j].end<<", "<<groups[j].groupId<<std::endl;
       } 
       std::cout<<std::endl;
#endif
    return gId;
}
