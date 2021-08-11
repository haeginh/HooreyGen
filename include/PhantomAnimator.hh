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
#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/intersect_other.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/ray_mesh_intersect.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/SparseCore>

class PhantomAnimator
{
public:
    enum JOINT
    {
        SHOULDER,
        ELBOW,
        WRIST,
        CLAVICLE
    };

    PhantomAnimator();
    PhantomAnimator(string prefix);
    ~PhantomAnimator();

    bool ReadTGF();
    bool CalculateWeights(double w_smooth, bool innerOnly = false);
    bool WriteWeights()
    {
        igl::writeDMAT(prefix + ".W", W, false);
        igl::writeDMAT(prefix + ".S", smoothMM, false);
        return true;
    }
    bool ReadW();

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
    MatrixXd GetVo() { return Vo; }
    MatrixXi GetFo() { return Fo; }
    MatrixXd GetVbone() { return Vbone; }
    MatrixXi GetFbone() { return Fbone; }
    MatrixXd GetC() { return C; }
    MatrixXi GetBE() { return BE; }
    VectorXi GetP() { return P; }
    VectorXd GetSmoothW(JOINT j)
    {
        return smoothMM.col((int)j);
    }
    double GetVolume(bool initial = false)
    {
        if (initial)
            return volume0;
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
        for (int i = 0; i < Fo.rows(); i++)
        {
            adjacentO[Fo(i, 0)].push_back(Fo(i, 1));
            adjacentO[Fo(i, 0)].push_back(Fo(i, 2));
            adjacentO[Fo(i, 1)].push_back(Fo(i, 0));
            adjacentO[Fo(i, 1)].push_back(Fo(i, 2));
            adjacentO[Fo(i, 2)].push_back(Fo(i, 1));
            adjacentO[Fo(i, 2)].push_back(Fo(i, 0));
            adjacentOF[Fo(i, 0)].push_back(i);
            adjacentOF[Fo(i, 1)].push_back(i);
            adjacentOF[Fo(i, 2)].push_back(i);            
        }
        for (auto &iter : adjacentO)
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

    VectorXd ArmOffset(MatrixXd normalsV, double w = 0.1)
    {
        Vector3d refR = C.row(9) - C.row(8);
        refR = -refR.cross(Vector3d(0, 0, -1)).cross(refR);
        refR.normalize();
        Vector3d refL = C.row(4) - C.row(3);
        refL = -refL.cross(Vector3d(0, 0, -1)).cross(refL);
        refL.normalize();
        VectorXd normL = (normalsV.array().rowwise() * refL.transpose().array()).rowwise().sum();
        VectorXd normR = (normalsV.array().rowwise() * refR.transpose().array()).rowwise().sum();

        VectorXd armL = W.col(2);
        VectorXd armR = W.col(6);
        armL /= armL.array().maxCoeff();
        armR /= armR.array().maxCoeff();
        VectorXd shoul = smoothMM.col(JOINT::SHOULDER);
        shoul /= shoul.array().maxCoeff();
        armL = armL.array() * shoul.array() * normL.array();
        armR = armR.array() * shoul.array() * normR.array();

        armL += W.col(1) * 0.2;
        armR += W.col(5) * 0.2;
        VectorXd tot = (armL + armR);
        double e(1e-3);
        for (int i = 0; i < V.rows(); i++)
        {
            //double norm = ref.dot(Vector3d())
            if (tot(i) > 0.5)
                tot(i) = 0.5;
            if (tot(i) > e)
                V.row(i) += normalsV.row(i) * tot(i) * w;
            else
                tot(i) = 0;
        }
        return tot;
    }

    VectorXd LowerArmOffset(MatrixXd normalsV, double w = 0.1)
    {
        VectorXd tot = smoothMM.col(0);
        double e(1e-3);
        for (int i = 0; i < V.rows(); i++)
        {
            //double norm = ref.dot(Vector3d())
            // if (tot(i) > 0.5)
            //     tot(i) = 0.5;
            if (tot(i) > e)
                V.row(i) += normalsV.row(i) * tot(i) * w;
            else
                tot(i) = 0;
        }

        return tot;
    }

    VectorXd ShoulderShrink(const MatrixXd &normalsV, double w = 0.1)
    {
        Vector3d refR = C.row(9) - C.row(8);
        refR = refR.cross(Vector3d(0, 0, -1)).cross(refR);
        refR.normalize();
        Vector3d refL = C.row(4) - C.row(3);
        refL = refL.cross(Vector3d(0, 0, -1)).cross(refL);
        refL.normalize();
        VectorXd normL = (normalsV.array().rowwise() * refL.transpose().array()).rowwise().sum();
        VectorXd normR = (normalsV.array().rowwise() * refR.transpose().array()).rowwise().sum();

        VectorXd shoul = W.col(1).array() * normL.array() + W.col(5).array() * normR.array();
        shoul /= shoul.array().maxCoeff();
        double e(1e-3);
        for (int i = 0; i < V.rows(); i++)
        {
            if (shoul(i) > e)
                V.row(i) -= normalsV.row(i) * shoul(i);
            else
                shoul(i) = 0;
        }
        return shoul;
    }

    int DetectSelfInterSection(MatrixXd &color, double tryFix = 0)
    {
        igl::copyleft::cgal::RemeshSelfIntersectionsParam param;
        param.detect_only = true;
        MatrixXd Vtmp;
        MatrixXi Ftmp, IF;
        VectorXi J, IM;
        igl::copyleft::cgal::remesh_self_intersections(V, F, param, Vtmp, Ftmp, IF, J, IM);
        if (IF.rows())
        {
            vector<int> Kvec(IF.data(), IF.data() + IF.size());
            sort(Kvec.begin(), Kvec.end());
            Kvec.erase(unique(Kvec.begin(), Kvec.end()), Kvec.end());
            VectorXi K = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(Kvec.data(), Kvec.size());
            MatrixXd R = RowVector3d(1., 0.3, 0.3).replicate(K.rows(), 1);
            igl::slice_into(R, K, 1, color);
            MatrixXi FF = igl::slice(F, K, 1);
            vector<int> Fvec(FF.data(), FF.data() + FF.size());
            if (tryFix > 0)
            {
                for (int i : Fvec)
                {
                    Vector3d c(0., 0., 0.);
                    for (int n : adjacent[i])
                        c += V.row(n);
                    c /= (double)adjacent[i].size();
                    V.row(i) = V.row(i) * (1 - tryFix) + c.transpose() * tryFix;
                }
            }
        }
        return IF.rows();
    }

    void InitializeV()
    {
        V = V0;
        C = C0;
    }

    void MoveIntersected(vector<int> vIdx, double w = 0.1)
    {
        sort(vIdx.begin(), vIdx.end());
        vIdx.erase(unique(vIdx.begin(), vIdx.end()), vIdx.end());

        ComputeNormal();
        vector<int> near;
        set<int> chk(vIdx.begin(), vIdx.end());
        for (int i : vIdx)
        {
            // chk.insert(F(i,0)); chk.insert(F(i,1)); chk.insert(F(i,2));
            // V.row(F(i,0)) += F_normals.row(i) * w;
            // V.row(F(i,1)) += F_normals.row(i) * w;
            // V.row(F(i,2)) += F_normals.row(i) * w;
            // near.insert(near.end(), adjacent[F(i,0)].begin(), adjacent[F(i,0)].end());
            // near.insert(near.end(), adjacent[F(i,1)].begin(), adjacent[F(i,1)].end());
            // near.insert(near.end(), adjacent[F(i,2)].begin(), adjacent[F(i,2)].end());
            V.row(i) += V_normals.row(i) * w;
            near.insert(near.end(), adjacent[i].begin(), adjacent[i].end());
        }
        sort(near.begin(), near.end());
        near.erase(unique(near.begin(), near.end()), near.end());
        for (int i : near)
        {
            if(chk.find(i)!=chk.end()) continue;
            Vector3d c(0, 0, 0);
            for (int j : adjacent[i])
                c += V.row(j);
            c /= (double)adjacent[i].size();
            // V.row(i) = V.row(i)*0.8 + c.transpose()*0.2;
            V.row(i) = (V.row(i) + V_normals.row(i) * w * 0.5) * 0.8 + c.transpose()*0.2;
        }
    }
    
    // outer skin
    void MergeOBJ();
    void WriteMergedOBJ();
    void SetSkinLayers(MatrixXd &color)
    {
        cout << "<Skin layer generation>" << endl;
        //initial offset
        MatrixXd DV = baryO * V;
        ComputeNormal();
        MatrixXd V_normalsO = baryO * V_normals;
        // for (auto iter : o2rst)
        //     Vo.row(iter.first) = V.row(iter.second) + V_normals.row(iter.second) * 0.2;
        for (int i = 0; i < D.rows(); i++)
            DV.row(i) = DV.row(i) + V_normalsO.row(i) * 0.2;
        igl::slice_into(DV, D, 1, Vo);
        function<bool()> AdjustVolumeO = [&]() -> bool
        {
            bool fixed = false;
            double diff = CalculateVolume(Vo, Fo) / volumeO0 - 1.f;
            V_normalsO = GetNormal(Vo, Fo) * 0.001;
            cout << "[vol. adjustment]" << flush;
            DV = igl::slice(Vo, D, 1);
            while (abs(diff) > 0.0001)
            {
                fixed = true;
                cout << " -> " << diff * 100. << flush;
                if (diff < 0)
                    for (int i = 0; i < D.rows(); i++)
                        DV.row(i) = DV.row(i) + V_normalsO.row(D(i));
                // for (auto iter : o2rst)
                //     Vo.row(iter.first) = Vo.row(iter.first) + V_normals.row(iter.second);
                else
                    for (int i = 0; i < D.rows(); i++)
                        DV.row(i) = DV.row(i) - V_normalsO.row(D(i));
                // for (auto iter : o2rst)
                //     Vo.row(iter.first) = Vo.row(iter.first) - V_normals.row(iter.second);
                igl::slice_into(DV, D, 1, Vo);
                diff = CalculateVolume(Vo, Fo) / volumeO0 - 1.f;
            }
            cout << " -> " << diff * 100 << endl;
            return fixed;
        };

        igl::copyleft::cgal::RemeshSelfIntersectionsParam param;
        param.detect_only = true;
        MatrixXd Vtmp;
        MatrixXi Ftmp, IF;
        VectorXi J, IM;
        function<bool()> FixSelfIntersectionO = [&]() -> bool
        {
            bool fixed(false);
            cout << "[self-intersection]" << flush;
            while (1)
            {
                igl::copyleft::cgal::remesh_self_intersections(Vo, Fo, param, Vtmp, Ftmp, IF, J, IM);
                cout << " -> " << IF.rows() << flush;
                if (IF.rows() == 0)
                    break;
                fixed = true;
                vector<int> Kvec(IF.data(), IF.data() + IF.size());
                sort(Kvec.begin(), Kvec.end());
                Kvec.erase(unique(Kvec.begin(), Kvec.end()), Kvec.end());
                VectorXi K = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(Kvec.data(), Kvec.size());
                MatrixXi FF = igl::slice(Fo, K, 1);
                vector<int> Fvec(FF.data(), FF.data() + FF.size());
                sort(Fvec.begin(), Fvec.end());
                Fvec.erase(unique(Fvec.begin(), Fvec.end()), Fvec.end());

                for (int i : Fvec)
                {
                    Vector3d c(0., 0., 0.);
                    for (int n : adjacentO[i])
                        c += Vo.row(n);
                    c /= (double)adjacentO[i].size();
                    Vo.row(i) = c;//Vo.row(i) * 0.5 + c.transpose() * 0.5;
                }
            }
            cout << endl;
            return fixed;
        };
        //if(FixSelfIntersectionO()) AdjustVolumeO();

        function<bool()> FixProtruded = [&]()->bool
        {
            bool fixed(false);
            MatrixXd F_normals;
            igl::per_face_normals(Vo, Fo, F_normals);

            for(int i=0;i<D.rows();i++)
            {
                Vector3d normal = F_normals.row(adjacentOF[D(i)][0]);
                for(int n:adjacentOF[D(i)])
                {
                    if(normal.dot(F_normals.row(n))<0)
                    {
                        Vector3d c(0.,0.,0.);
                        for(int m:adjacentO[D(i)])
                            c += Vo.row(m);
                        Vo.row(D(i)) = c / (double) adjacentO[D(i)].size();
                        fixed = true;
                        break;
                    }
                }
            }
            return fixed;
        };

        function<bool()> FixDegen = [&]()->bool
        {
            bool fixed(false);
            set<int> chk;
            VectorXd WV = GetSmoothW(SHOULDER);
            for(int i=0;i<D.rows();i++) chk.insert(D(i));
            for(int i=0;i<Fo.rows();i++)
            {
                if(chk.find(Fo(i, 0))==chk.end()) continue;
                if(chk.find(Fo(i, 1))==chk.end()) continue;
                if(chk.find(Fo(i, 2))==chk.end()) continue;
                if(WV(Fo(i, 0))+WV(Fo(i, 1))+WV(Fo(i, 2))<0.1) continue;
                Vector3d a = (Vo.row(Fo(i,0))-Vo.row(Fo(i, 1))).normalized();
                Vector3d b = (Vo.row(Fo(i,1))-Vo.row(Fo(i, 2))).normalized();
                Vector3d c = (Vo.row(Fo(i,2))-Vo.row(Fo(i, 0))).normalized();
                if(c.dot(-a)<0.95 && a.dot(-b)<0.95 && b.dot(-c)<0.95)  continue;
                for(int j=0;j<3; j++)
                {
                    Vector3d cc(0.,0.,0.);
                    for(int n:adjacentO[Fo(i, j)]) cc += Vo.row(n);
                    Vo.row(Fo(i, j)) = Vo.row(Fo(i, j))*0.3 + cc.transpose()/(double)adjacentO[Fo(i, j)].size()*0.7;
                }
                fixed = true;
            }
            return fixed;
        };

        function<bool()> CheckRSTIntersection = [&]() -> bool
        {
            bool fixed = false;
            cout << "[fix RST intersection]" << flush;
            while (1)
            {
                V_normalsO = GetNormal(Vo, Fo);
                MatrixXd V100 = Vo - V_normalsO * 0.01;
                if (!igl::copyleft::cgal::intersect_other(V100, Fo, V, F, false, IF))
                    break;
                cout << " -> " << IF.rows() << flush;
                fixed = true;
                MatrixXi FF = igl::slice(Fo, IF.col(0), 1);
                VectorXi faceVI(Map<VectorXi>(FF.data(), FF.cols() * FF.rows()));
                vector<int> Fvec(faceVI.data(), faceVI.data() + faceVI.size());
                sort(Fvec.begin(), Fvec.end());
                Fvec.erase(unique(Fvec.begin(), Fvec.end()), Fvec.end());
                //vector<int> sVec;
                for (int i : Fvec)
                {
                    Vo.row(i) += V_normalsO.row(i) * 0.01;
                //    sVec.insert(sVec.end(), adjacentO[i].begin(), adjacentO[i].end());
                }
                // sort(sVec.begin(), sVec.end());
                // sVec.erase(unique(sVec.begin(), sVec.end()), sVec.end());
                // for (int i : sVec)
                // {
                //     Vector3d c(0., 0., 0.);
                //     for (int n : adjacentO[i])
                //         c += Vo.row(n);
                //     c /= (double)adjacentO[i].size();
                //     Vo.row(i) = Vo.row(i) * 0.7 + c.transpose() * 0.3;
                // }
            }
            cout << " -> 0" << endl;
            return fixed;
        };
        //CheckRSTIntersection();

        function<bool()> CheckLayerIntersection = [&]() -> bool
        {
            bool fixed = false;
            cout << "[layer-intersection]" << flush;
            while (1)
            {
                V_normalsO = GetNormal(Vo, Fo);
                MatrixXd V100 = Vo - V_normalsO * 0.01;
                MatrixXd V50 = Vo - V_normalsO * 0.005;
                bool chk1 = igl::copyleft::cgal::intersect_other(V100, Fo, V50, Fo, false, IF);
                MatrixXi FF = igl::slice(Fo, IF.col(0), 1);
                VectorXi faceVI(Map<VectorXi>(FF.data(), FF.cols() * FF.rows()));
                vector<int> Fvec(faceVI.data(), faceVI.data() + faceVI.size());
                bool chk2 = igl::copyleft::cgal::intersect_other(V50, Fo, Vo, Fo, false, IF);
                if (!chk1 && !chk2)
                    break;
                FF = igl::slice(Fo, IF.col(0), 1);
                faceVI = VectorXi(Map<VectorXi>(FF.data(), FF.cols() * FF.rows()));
                Fvec.insert(Fvec.end(), faceVI.data(), faceVI.data() + faceVI.size());
                fixed = true;
                sort(Fvec.begin(), Fvec.end());
                Fvec.erase(unique(Fvec.begin(), Fvec.end()), Fvec.end());
                cout << " -> " << Fvec.size() << flush;
                vector<int> sVec = Fvec;
                for (int i : Fvec)
                {
                    for (int n : adjacentO[i])
                        sVec.insert(sVec.end(), adjacentO[i].begin(), adjacentO[i].end());
                }
                sort(sVec.begin(), sVec.end());
                sVec.erase(unique(sVec.begin(), sVec.end()), sVec.end());
                for (int i : sVec)
                {
                    Vector3d c(0., 0., 0.);
                    for (int n : adjacentO[i])
                        c += Vo.row(n);
                    c /= (double)adjacentO[i].size();
                    Vo.row(i) = Vo.row(i) * 0.9 + c.transpose() * 0.1;
                }
            }
            cout << endl;
            return fixed;
        };
        CheckLayerIntersection();
        FixDegen();
        FixProtruded();
        FixSelfIntersectionO();
        FixProtruded();
        AdjustVolumeO();
        while (1)
        {
            if(!FixSelfIntersectionO() && !CheckLayerIntersection() && !CheckRSTIntersection()) break;
            AdjustVolumeO();
        }

        // if (IF.rows())
        // {
        //     cout << "There are " << IF.rows() << " intersections between RST and skin-100um (blue)" << endl;
        //     MatrixXd B = RowVector3d(0.3, 0.3, 1.0).replicate(IF.rows(), 1);
        //     igl::slice_into(B, IF.col(0), 1, color);
        // }
        // else
        //     cout << "There is no intersection between RST and skin-100um" << endl;
        // igl::copyleft::cgal::intersect_other(V100, Fo, V50, Fo, false, IF);
        // if (IF.rows())
        // {
        //     cout << "There are " << IF.rows() << " intersections between skin-50um and skin-100um (green)" << endl;
        //     MatrixXd G = RowVector3d(0.3, 1.0, 0.3).replicate(IF.rows(), 1);
        //     igl::slice_into(G, IF.col(0), 1, color);
        // }
        // else
        //     cout << "There is no intersection between skin-50um and skin-100um" << endl;
        // igl::copyleft::cgal::intersect_other(Vo, Fo, V50, Fo, false, IF);
        // if (IF.rows())
        // {
        //     cout << "There are " << IF.rows() << " intersections between skin-50um and skin-outer (green)" << endl;
        //     MatrixXd G = RowVector3d(0.3, 1.0, 0.3).replicate(IF.rows(), 1);
        //     igl::slice_into(G, IF.col(0), 1, color);
        // }
        // else
        //     cout << "There is no intersection between skin-50um and skin-outer" << endl;
    }
    void ConfirmSkinLayers();

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
    MatrixXd GetNormal(MatrixXd _V, MatrixXi _F)
    {
        MatrixXd _F_normals, _V_normals;
        igl::per_face_normals(_V, _F, _F_normals);
        igl::per_vertex_normals(_V, _F, _F_normals, _V_normals);
        return _V_normals;
    }
    void ComputeNormal(){
        MatrixXd _F_normals;
        igl::per_face_normals(V, F, _F_normals);
        igl::per_vertex_normals(V, F, _F_normals, V_normals);
    }

    //variables
    MatrixXd C, C0, V, V0;
    MatrixXi BE, BE0, F;
    VectorXi P;
    double volume0;
    string prefix;
    MatrixXd V_normals;

    //dqs
    vector<map<int, double>> cleanWeights;

    //LBS
    MatrixXd W;

    //smoothing
    map<int, vector<int>> adjacent, adjacentO, adjacentOF;
    vector<map<int, map<int, double>>> smoothM;
    MatrixXd smoothMM;
    map<int, vector<int>> adjacent_T;

    //merging
    map<int, vector<Vector3d>> VV;
    map<int, vector<Vector3i>> FF;
    map<int, string> shellName;
    MatrixXd Vo, Vbone;
    MatrixXi Fo, Fbone;
    map<int, int> o2rst;
    VectorXi D;
    SparseMatrix<double> baryO;
    map<int, vector<int>> handsO;
    double volumeO0;

public:
    MatrixXd GetW() { return W; }
    MatrixXd GetsmoothMM() { return smoothMM; }
};

#endif