#ifndef PhantomAnimator_class
#define PhantomAnimator_class

#include <iostream>
#include <ctime>
#include <functions.h>
#include <functional>

//#include "bodytracking.hh"

#include <igl/readTGF.h>
#include <igl/writeTGF.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/writeMESH.h>
#include <igl/readMESH.h>
#include <igl/readPLY.h>
#include <igl/directed_edge_parents.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/boundary_conditions.h>
#include <igl/directed_edge_parents.h>
#include <igl/directed_edge_orientations.h>
#include <igl/deform_skeleton.h>
#include <igl/forward_kinematics.h>
#include <igl/doublearea.h>
#include <igl/dqs.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/Timer.h>
#include <igl/mat_max.h>
#include <igl/arap.h>
#include <igl/shortest_edge_and_midpoint.h>

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
    bool ReadOBJ(string fileName);
    bool ReadOBJW(string fileName);
    bool SetSimpleF(int, bool onlyV = false);
    VectorXd GetWeight(int i, int col);
    void Animate(RotationList vQ, MatrixXd &CT, MatrixXi &BET);

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

    MatrixXd GetV_simple() { return V_simple; }
    MatrixXi GetF_simple() { return F_simple; }

    void PreComputeSmooth()
    {
        smoothingM.clear();
        double epsilon(1e-5);
        bw = VectorXd::Zero(extractedVertices[125].size());
        for (int i = 0; i < extractedVertices[125].size(); i++)
        {
            int v = extractedVertices[125][i];
            if (cleanWeights[v].find(3)->second > epsilon || cleanWeights[v].find(4)->second > epsilon)
            {
                smoothingM[i] = {};
                bw(i) = 1;
            }
        }

        for (int i = 0; i < F_simple.rows(); i++)
        {
            if (smoothingM.find(F_simple(i, 0)) != smoothingM.end())
            {
                smoothingM[F_simple(i, 0)].push_back(F_simple(i, 1));
                smoothingM[F_simple(i, 0)].push_back(F_simple(i, 2));
            }
            if (smoothingM.find(F_simple(i, 1)) != smoothingM.end())
            {
                smoothingM[F_simple(i, 1)].push_back(F_simple(i, 0));
                smoothingM[F_simple(i, 1)].push_back(F_simple(i, 2));
            }
            if (smoothingM.find(F_simple(i, 2)) != smoothingM.end())
            {
                smoothingM[F_simple(i, 2)].push_back(F_simple(i, 1));
                smoothingM[F_simple(i, 2)].push_back(F_simple(i, 0));
            }
        }
    }

    void LaplacianSmooth(int iterNum, double w)
    {
        for (int i = 0; i < iterNum; i++)
        {
            for (auto iter : smoothingM)
            {
                Vector3d target(0,0,0);
                for(int v:iter.second)
                    target += V_simple.row(v);
                target /= (double) iter.second.size();
                V_simple.row(iter.first) = V_simple.row(iter.first).transpose()*(1-w) + target*w;
            }
        }
    }

    //variables
private:
    MatrixXd C, V, V_simple, V0;
    MatrixXi BE, F, F_simple;
    VectorXi P;
    vector<map<int, double>> cleanWeights;
    map<int, pair<int, int>> shells;
    map<int, vector<int>> extractedVertices;

    //LBS
    MatrixXd M;

    //smoothing
    map<int, vector<int>> smoothingM;

public:
    VectorXd bw;
};

#endif