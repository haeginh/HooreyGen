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
        for (int i = 0; i < FI.rows(); i++)
        {
            adjacentI[FI(i, 0)].push_back(FI(i, 1));
            adjacentI[FI(i, 0)].push_back(FI(i, 2));
            adjacentI[FI(i, 1)].push_back(FI(i, 0));
            adjacentI[FI(i, 1)].push_back(FI(i, 2));
            adjacentI[FI(i, 2)].push_back(FI(i, 1));
            adjacentI[FI(i, 2)].push_back(FI(i, 0));
        }
        for (auto &iter : adjacentI)
        {
            sort(iter.second.begin(), iter.second.end());
            iter.second.erase(unique(iter.second.begin(), iter.second.end()), iter.second.end());
        }
    }

    void LaplacianSmooth(int j, int iterNum, double fix0, double fix1) //all
    {
        for (int i = 0; i < iterNum; i++)
        {
            for (int n = Vnum[j].first; n < Vnum[j].first + surface[j]; n++)
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
                for (int n = Vnum[j].first; n < Vnum[j].first + surface[j]; n++)
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

        int count(1);
        for (auto iter : Vnum)
        {
            MatrixXd V1;
            MatrixXi F1;
            if (iter.first < 100)
            {
                V1 = V.block(iter.second.first, 0, iter.second.second, 3);
                F1 = F.block(Fnum[iter.first].first, 0, Fnum[iter.first].second, 3).array() - iter.second.first + count;
            }
            else
            {
                V1 = VI.block(iter.second.first, 0, iter.second.second, 3);
                F1 = FI.block(Fnum[iter.first].first, 0, Fnum[iter.first].second, 3).array() - iter.second.first + count;
            }
            count += iter.second.second;

            for (int i = 0; i < V1.rows(); i++)
                ofs << "v " << V1(i, 0) << " " << V1(i, 1) << " " << V1(i, 2) << endl;

            ofs << endl
                << "g " << iter.first << endl;
            ofs << "usemtl " << iter.first << endl;
            ofs << "s" << endl
                << endl;

            for (int i = 0; i < F1.rows(); i++)
                ofs << "f " << F1(i, 0) << " " << F1(i, 1) << " " << F1(i, 2) << endl;
        }
    }

    //intersection detection
    int DetectIntersections(MatrixXd &Color, int a, int b = -1);
    int DetectIntersections(MatrixXd &Color);
    int DetectIntersections(MatrixXd &color, MatrixXd &colorI, bool fix=false);
    enum SHELL
    {
        INNER,
        OUTER
    };
    int SelfIntersection(MatrixXd &Color, SHELL shell, double fix = 0)
    {
        igl::copyleft::cgal::RemeshSelfIntersectionsParam param;
        param.detect_only = true;
        MatrixXd VV, V1;
        MatrixXi FF, IF, TT, F1, FD1;
        VectorXi J, IM;
        switch (shell)
        {
        case INNER:
            igl::copyleft::cgal::remesh_self_intersections(VI, FI, param, VV, FF, IF, J, IM);
            break;
        case OUTER:
            FD1 = igl::slice(F, FD, 1);
            igl::copyleft::cgal::remesh_self_intersections(V, FD1, param, VV, FF, IF, J, IM);
            break;
        }
        if (IF.rows() > 0)
        {
            vector<int> Kvec(IF.data(), IF.data() + IF.size());
            sort(Kvec.begin(), Kvec.end());
            Kvec.erase(unique(Kvec.begin(), Kvec.end()), Kvec.end());
            VectorXi K = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(Kvec.data(), Kvec.size());
            if(shell==OUTER) 
                for(int i=0;i<K.rows();i++) K(i) = FD(K(i));
            MatrixXd R = RowVector3d(1., 0.3, 0.3).replicate(K.rows(), 1);
            igl::slice_into(R, K, 1, Color);
            if (fix > 0)
            {
                MatrixXi face;
                if (shell == INNER)
                    face = igl::slice(FI, K, 1);
                else if (shell == OUTER)
                    face = igl::slice(F, K, 1);
                vector<int> Vvec(face.data(), face.data() + face.size());
                sort(Vvec.begin(), Vvec.end());
                Vvec.erase(unique(Vvec.begin(), Vvec.end()), Vvec.end());
                if (shell == INNER)
                {
                    for (int v : Vvec)
                    {
                        Vector3d c(0.,0.,0.);
                        for (int n : adjacentI[v]) c += VI.row(n);
                        c /= (double) adjacentI[v].size();
                        VI.row(v) = VI.row(v)*(1-fix) + c.transpose()*fix;
                    }
                }
                else if (shell == OUTER)
                {
                    for (int v : Vvec)
                    {
                        Vector3d c(0.,0.,0.);
                        for (int n : adjacent[v]) c += V.row(n);
                        c /= (double) adjacent[v].size();
                        V.row(v) = V.row(v)*(1-fix) + c.transpose()*fix;
                    }
                }               
            }
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

    //inner structure
    void SetInnerVertices()
    {
        VI = bary * V;
    }

    void AdjustInnerVolume(int id, int j, double fix0, double fix1, MatrixXd &normal)
    {
        MatrixXd W2(V.rows(), 2);
        W2.col(0) = W1.col(j);
        W2.col(1) = W1.col(j + 4);
        W2 = bary * W2;
        double diff = GetVolume(id) / GetVolume(-id) - 1.f;
        cout<<"[vol. adj #"<<id<<"]"<<flush;
        while (abs(diff) > 0.0001)
        {
            //diff *= 1000.f;

            for (int n = Vnum[id].first; n < Vnum[id].first + Vnum[id].second; n++)
            {
                int col(0);
                if (VI(n, 0) < 0)
                    col = 1;
                //cout<<W2(n, col)<<endl;
                double w = min(abs(W2(n, col) - fix0), abs(W2(n, col) - fix1));
                // if (w < 0.001)
                //     continue;
                // else if (diff > 0)
                //w *= -diff;
                VI.row(n) -= normal.row(n) * w * diff;
            }
            diff = GetVolume(id) / GetVolume(-id) - 1.f;
            cout<<" -> "<<diff*100.f<<flush;
        }
        cout<<" -> "<<diff*100.f<<endl;
    }

    void AdjustInnerVolume2(int id, int j, double fix0, double fix1, MatrixXd &normal)
    {
        MatrixXd W2(V.rows(), 2);
        W2.col(0) = W.col(j);
        W2.col(1) = W.col(j + 4);
        W2 = bary * W2;
        double diff = GetVolume(id) / GetVolume(-id) - 1.f;
        while (abs(diff) > 0.0001)
        {
            for (int n = Vnum[id].first; n < Vnum[id].first + Vnum[id].second; n++)
            {
                int col(0);
                if (VI(n, 0) < 0)
                    col = 1;
                //cout<<W2(n, col)<<endl;
                double w = min(abs(W2(n, col) - fix0), abs(W2(n, col) - fix1));
                if (w < 0.001)
                    continue;
                else if (diff > 0)
                    w = -w;
                VI.row(n) += normal.row(n) * w * 0.001;
            }
            diff = GetVolume(id) / GetVolume(-id) - 1.f;
        }
    }

    //set W1 for vol. adjustment
    void SetW1() { W1 = W; }

    void InitializeV()
    {
        V = V0;
        C = C0;
    }
    vector<int> GetShellVec()
    {
        vector<int> shells;
        for (auto iter : Fnum)
            shells.push_back(iter.first);
        return shells;
    }

    //volume check
    double GetVolume(int id)
    {
        if (id < 0)
            return volumes[-id];
        else if (id < 100)
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
    MatrixXd GetNormal(MatrixXd _V, MatrixXi _F)
    {
        MatrixXd F_normals, V_normals;
        igl::per_face_normals(_V, _F, F_normals);
        igl::per_vertex_normals(_V, _F, F_normals, V_normals);
        return V_normals;
    }

    //variables
    string name;
    MatrixXd C, C0, V, V0;
    MatrixXi BE, BE0, F, T;
    VectorXi P;
    vector<int> surface;

    //dqs
    vector<map<int, double>> cleanWeights;

    //LBS
    MatrixXd W, W1; // W1 is for volume adjustment

    //smoothing
    map<int, vector<int>> adjacent, adjacentI;

    //bone ply
    map<int, pair<int, int>> Vnum, Fnum, Tnum; //#V, #F
    //map<int, int> boneID;

    //volumes
    map<int, double> volumes;

    //inner structures
    MatrixXd VI, VI0, WI;
    MatrixXi FI;
    SparseMatrix<double> bary;
    VectorXi FD;
};

#endif