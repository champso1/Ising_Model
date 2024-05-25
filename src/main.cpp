#include "utils.hpp"
#include "IsingModel.hpp"


int main(void) {


    IsingModel *ismdl = new IsingModel();
    ismdl->Run();
    delete ismdl;
    
    return 0;
}