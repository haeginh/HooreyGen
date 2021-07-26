#ifndef PhantomAnimator_class
#define PhantomAnimator_class

#include <iostream>
#include <ctime>
#include <functions.h>
#include <functional>

//#include "bodytracking.hh"

#include <igl/readTGF.h>
#include <igl/writeTGF.h>
#include <igl/writeDMAT.h>
#include <igl/writeMESH.h>
#include <igl/writePLY.h>
#include <igl/readDMAT.h>
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

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/SparseCore>

class Viewer;
class PhantomAnimator
{
public:
    enum JOINT
    {
        SHOULDER,
        ELBOW,
        WRIST
    };

    PhantomAnimator();
    PhantomAnimator(string prefix);
    ~PhantomAnimator();

    bool ReadFiles(string prefix);
    bool CalculateWeights(double w_smooth, bool boneOnly = false);
    bool WriteWeights(string prefix)
    {
        igl::writeDMAT(prefix + ".W", W, false);
        igl::writeDMAT(prefix + ".S", smoothMM, false);
        return true;
    }
    bool ReadW(string prefix);

    VectorXd GetWeight(int col) { return W.col(col); }
    void Animate(RotationList vQ, bool fixOther = false);
    void AnimateDQS(RotationList vQ);

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
    VectorXi GetP() { return P; }
    VectorXd GetSmoothW(JOINT j)
    {
        return smoothMM.col((int)j);
    }
    double GetVolume(bool initial = false){
        if(initial) return volume0;
        return CalculateVolume(V, F);
    }

    void PreComputeAdjacent()
    {
        for (int i = 0; i < F.rows(); i++)
        {
            adjacent[F(i, 0)].push_back(F(i, 1));
            adjacent[F(i, 0)].push_back(F(i, 2));
            adjacent[F(i, 1)].push_back(F(i, 0));
            adjacent[F(i, 1)].push_back(F(i, 2));
            adjacent[F(i, 2)].push_back(F(i, 1));
            adjacent[F(i, 2)].push_back(F(i, 0));
        }
        for (auto &iter : adjacent)
        {
            sort(iter.second.begin(), iter.second.end());
            iter.second.erase(unique(iter.second.begin(), iter.second.end()), iter.second.end());
        }
    }


    void LaplacianSmooth(JOINT j, int iterNum) //all
    {
        for (int i = 0; i < iterNum; i++)
        {
            for (auto iter : smoothM[(int)j])
            {
                Vector3d v(0, 0, 0);
                for (auto iter2 : iter.second)
                {
                    v += V.row(iter2.first) * iter2.second;
                }
                V.row(iter.first) = v;
            }
        }
    }

    void ReleaseRest(double w = 1.)
    {
        MatrixXd U0 = V0.array().colwise() * W.col(0).array() * w;
        MatrixXd U = V.array().colwise() * (1 - W.col(0).array() * w);
        V = U0 + U;
    }

    void ArmOffset(MatrixXd normalsV, double w = 0.1)
    {
        VectorXd arms = W.col(2) + W.col(6);
        //arms = arms.array().pow(3);
        double e(1e-3);
        for (int i = 0; i < V.rows(); i++)
        {
            if (arms(i) < e)
                continue;
            V.row(i) += normalsV.row(i) * arms(i) * w;
        }
    }

    void ShoulderOffset(MatrixXd normalsV, double w = 0.1)
    {
        VectorXd up = smoothMM.col(JOINT::SHOULDER);
        up /= up.array().maxCoeff();
        //V.col(2) += up*w;
        double e(1e-3);
        for (int i = 0; i < V.rows(); i++)
        {
            if (up(i) < e)
                continue;
            V.row(i) += normalsV.row(i) * up(i) * w;
        }
    }

    void InitializeV()
    {
        V = V0;
        C = C0;
    }


private:
    double CalculateVolume(const MatrixXd &V1, const MatrixXi &F1)
    {
        Eigen::MatrixXd V2(V1.rows() + 1, V1.cols());
        V2.topRows(V1.rows()) = V;
        V2.bottomRows(1).setZero();
        Eigen::MatrixXi T(F.rows(), 4);
        T.leftCols(3) = F;
        T.rightCols(1).setConstant(V.rows());
        Eigen::VectorXd vol;
        igl::volume(V2, T, vol);
        return std::abs(vol.sum());
    }

    //variables
    MatrixXd C, C0, V, V0;
    MatrixXi BE, BE0, F;
    VectorXi P;
    double volume0;

    //dqs
    vector<map<int, double>> cleanWeights;

    //LBS
    MatrixXd W;

    //smoothing
    map<int, vector<int>> adjacent;
    vector<map<int, map<int, double>>> smoothM;
    MatrixXd smoothMM;

    map<int, vector<int>> adjacent_T;
};

#endif