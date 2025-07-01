#include "testing/DensityDecom.h"


int main(int argc, char **argv){
    
    MultilayerGraph mg;
    char *path = argv[1];
    string path_str(path);

    MultilayerGraph::GetMulGraphClock().Start();
    mg.LoadFromFile(path_str);
    #ifdef VSHT
    DD_VSHT(mg);
    #endif

    #ifdef HPDD
    DD_HPDD(mg);
    #endif

    #ifdef HPDDPLUS
    DD_HPDD_Plus(mg);
    #endif


    mg.PrintStatistics();
    return 0;
}