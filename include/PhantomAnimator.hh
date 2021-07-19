#ifndef PhantomAnimator_class
#define PhantomAnimator_class

#include <iostream>
#include <ctime>
#include <functions.h>
#include <functional>

#include <igl/readTGF.h>
#include <igl/writeTGF.h>
#include <igl/writeDMAT.h>
#include <igl/writeMESH.h>
#include <igl/writePLY.h>
#include <igl/readDMAT.h>
#include <igl/readMESH.h>
#include <igl/readPLY.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/boundary_conditions.h>
#include <igl/mat_max.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/SparseCore>

class Viewer;
class PhantomAnimator
{
    //functions
public:
    PhantomAnimator();
    PhantomAnimator(string prefix);
    ~PhantomAnimator();

    bool ReadFiles(string prefix);
    bool CalculateWeights();

    VectorXd GetWeight(int col) { return W.col(col); }

    void GetMeshes(MatrixXd &_V, MatrixXi &_F, MatrixXd &_C, MatrixXi &_BE)
    {
        _V = V;
        _F = F;
        _C = C;
        _BE = BE;
    }
    MatrixXd GetV() { return V; }
    MatrixXi GetF() { return F; }
    MatrixXd GetC() { return C; }
    MatrixXi GetBE() { return BE; }

    //cutting arms
    bool InitArms();
    void CutArms(double threshold);
private:
    string name;
    MatrixXd ReadNode(string fileName);
    MatrixXi ReadEle(string fileName);
    MatrixXd V0;
    MatrixXi T0;
    VectorXd armW;

    //variables
    MatrixXd C, V;
    MatrixXi F, T;
    MatrixXi BE;

    //LBS
    MatrixXd W;
};

#endif