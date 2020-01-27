/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-27 17:56:36
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-27 18:24:29
 */
#include "SM_cxSM.h"
#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
    CXSM model;
    model.SetInput(10,50,0,0.1);
    model.PrintPotentialParameter();
    model.PrintLocalMinima();
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
    return 0;
}
