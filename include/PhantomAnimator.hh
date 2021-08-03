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
#include <igl/copyleft/cgal/intersect_other.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include "igl/opengl/glfw/Viewer.h"

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/SparseCore>

#define PI 3.14159265358979323846

class Viewer;
class PhantomAnimator
{
    //functions
public:
    PhantomAnimator();
    PhantomAnimator(string prefix);
    ~PhantomAnimator();

    bool ReadFiles(string prefix);
    void Animate(RotationList vQ);

    MatrixXd GetV() { return V; }
    MatrixXi GetF() { return F; }
    MatrixXd GetVI() { return VI; }
    MatrixXi GetFI() { return FI; }
    MatrixXd GetC() { return C; }
    MatrixXi GetBE() { return BE; }
    VectorXi GetP() { return P; }

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

    void LaplacianSmooth(int j, int iterNum, double fix0, double fix1) //all
    {
        for (int i = 0; i < iterNum; i++)
        {
            for (int n = Vnum[j].first; n < Vnum[j].first + Vnum[j].second; n++)
            {
                Vector3d v(0, 0, 0);
                for (int a : adjacent[n])
                {
                    v += V.row(a);
                }
                v /= (double)adjacent[n].size();
                double w = min(abs(W1(n, j) - fix0), abs(W1(n, j) - fix1));
                if (w < 0.001)
                    continue;
                w *= 0.01;
                V.row(n) = V.row(n) * (1 - w) + v.transpose() * w;
            }
        }
    }

    void AdjustVolume(int j, double fix0, double fix1, const MatrixXd &normal)
    {
        double diff = GetVolume(j) / GetVolume(-j) - 1.f;
        while (abs(diff) > 0.0001)
        {
            if (diff > 0)
                LaplacianSmooth(j, 3, fix0, fix1);
            else
            {
                for (int n = Vnum[j].first; n < Vnum[j].first + Vnum[j].second; n++)
                {
                    double w = min(abs(W1(n, j) - fix0), abs(W1(n, j) - fix1));
                    if (w < 0.001)
                        continue;
                    V.row(n) += normal.row(n) * w * 0.01;
                }
            }
            diff = GetVolume(j) / GetVolume(-j) - 1.f;
        }
    }

    void WriteOBJ(string fileN)
    {
        ofstream ofs(fileN);
        ofs << "#result file for arm bones" << endl
            << endl;
        ofs << "mtllib " + fileN.substr(0, fileN.length() - 3) + "mtl" << endl
            << endl;
        ofs.setf(ios::fixed);
        ofs.setf(ios::showpoint);
        ofs.precision(5);
        for (int i = 0; i < V.rows(); i++)
            ofs << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << endl;
        for (int s = 0; s <  BE.rows(); s++)
        {
            ofs << endl
                << "g " << s << endl;
            ofs << "usemtl " << s << endl;
            ofs << "s" << endl
                << endl;
            for (int i = 0; i < Fnum[s].second; i++)
                ofs << "f " << F(Fnum[s].first + i, 0) + 1 << " " << F(Fnum[s].first + i, 1) + 1 << " " << F(Fnum[s].first + i, 2) + 1 << endl;
        }
    }

    //intersection detection
    int DetectIntersections(MatrixXd &Color, int a, int b = -1);
    int DetectIntersections(MatrixXd &Color);
    int SelfIntersection(MatrixXd &Color)
    {
        igl::copyleft::cgal::RemeshSelfIntersectionsParam param;
        param.detect_only = true;
        MatrixXd VV;
        MatrixXi FF, IF, TT;
        VectorXi J, IM;
        igl::copyleft::cgal::remesh_self_intersections(V, F, param, VV, FF, IF, J, IM);
        if (IF.rows() > 0)
        {
            vector<int> Kvec(IF.data(), IF.data() + IF.size());
            sort(Kvec.begin(), Kvec.end());
            Kvec.erase(unique(Kvec.begin(), Kvec.end()), Kvec.end());
            MatrixXi K = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(Kvec.data(), Kvec.size());
            MatrixXd R = RowVector3d(1., 0.3, 0.3).replicate(K.rows(), 1);
            igl::slice_into(R, K, 1, Color);
        }
        return IF.rows();
    }

    //giving weights
    void CalculateElbowW(int a, Vector3d axis, double inner2 = 1.f, double outer2 = 30.f);
    double CalculateShoulderW(int a, Vector3d axis, double inner2 = 1.f, double outer2 = 25.f);
    void CalculateCleanWeights();
    void InitializeW()
    {
        W = MatrixXd::Zero(V.rows(), BE.rows());
        for (int i = 0; i < BE.rows(); i++)
            W.block(Vnum[i].first, i, Vnum[i].second, 1) = VectorXd::Ones(Vnum[i].second);
    }

    VectorXd GetWeight(int i) { return W.col(i); }
    void LeaveOnlyShoulderWeights()
    {
        for (int i = 0; i < BE.rows(); i++)
        {
            if (i == 1 || i == 5)
                continue;
            W.block(Vnum[i].first, 0, Vnum[i].second, BE.rows()) = MatrixXd::Zero(Vnum[i].second, BE.rows());
            W.block(Vnum[i].first, 0, Vnum[i].second, 1) = VectorXd::Ones(Vnum[i].second);
        }
    }

    //set W1 for vol. adjustment
    void SetW1(){W1=W;}

    void InitializeV()
    {
        V = V0;
        C = C0;
    }
    vector<int> GetShellVec(){
        vector<int> shells;
        for(auto iter:Fnum) shells.push_back(iter.first);
        return shells;
    }

    //volume check
    double GetVolume(int id)
    {
        if (id < 0)
            return volumes[-id];
        else if (id<100)
            return CalculateVolume(V, F.block(Fnum[id].first, 0, Fnum[id].second, 3));
        else
            return CalculateVolume(VI, FI.block(Fnum[id].first, 0, Fnum[id].second, 3));
    }
private:
    double CalculateVolume(const MatrixXd &V1, const MatrixXi &F1)
    {
        Eigen::MatrixXd V2(V1.rows() + 1, V1.cols());
        V2.topRows(V1.rows()) = V1;
        V2.bottomRows(1).setZero();
        Eigen::MatrixXi T(F1.rows(), 4);
        T.leftCols(3) = F1;
        T.rightCols(1).setConstant(V1.rows());
        Eigen::VectorXd vol;
        igl::volume(V2, T, vol);
        return std::abs(vol.sum());
    }
    void ReadOBJ(string fName);

    //variables
    string name;
    MatrixXd C, C0, V, V0;
    MatrixXi BE, BE0, F;
    VectorXi P;

    //dqs
    vector<map<int, double>> cleanWeights;

    //LBS
    MatrixXd W, W1; // W1 is for volume adjustment

    //smoothing
    map<int, vector<int>> adjacent;

    //bone ply
    map<int, pair<int, int>> Vnum, Fnum; //#V, #F
    //map<int, int> boneID;

    //volumes
    map<int, double> volumes;

    //inner structures
    MatrixXd VI, VI0, WI;
    MatrixXi FI;
};

#endif