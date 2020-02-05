/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-02-04 13:29:13
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-05 12:13:37
 */
#include "Relaxation.h"


Relaxation::Relaxation()
{
    SetDOF(1,1);
    SetGrid();
    SetMaxIteration();
    SetConvergeCriterion();
}
Relaxation::Relaxation(double NField, int NB_Left)
{
    SetDOF(NField,NB_Left);
    SetGrid();
    SetMaxIteration();
    SetConvergeCriterion();
}
void Relaxation::SetGrid(int MeshPoints)
{
    _MeshPoints = MeshPoints;
}
void Relaxation::SetMaxIteration(int itemax)
{
    _ITE_MAX = itemax;
}
void Relaxation::SetDOF(int NField, int NB_Left)
{
    _NField = NField;
    _NE = 2*_NField;
    _NB_Left = NB_Left;
}
void Relaxation::SetConvergeCriterion(double conv)
{
    _conv = conv;
}