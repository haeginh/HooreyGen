#include "PhantomAnimator.hh"
#include "igl/opengl/glfw/Viewer.h"
PhantomAnimator::PhantomAnimator() {}

PhantomAnimator::PhantomAnimator(string prefix)
{
    ReadFiles(prefix);
    PreComputeAdjacent();
    volume0 = CalculateVolume(V0, F);
    //CalculateWeights(0.1);
}

bool PhantomAnimator::ReadFiles(string prefix)
{
    cout << "Read " + prefix + ".tgf" << endl;
    igl::readTGF(prefix + ".tgf", C, BE0);
    BE = BE0.block(0, 0, 9, 2);
    C0 = C;
    igl::directed_edge_parents(BE, P);
    cout << "Read " + prefix + ".ply" << endl;
    igl::readPLY(prefix + ".ply", V, F);
    V0 = V;

    //SetSimpleF();

    return true;
}

bool PhantomAnimator::CalculateWeights(double w_smooth, bool boneOnly)
{
    //perform BBW
    MatrixXd boneP = GenerateBonePoints(C, BE0, 1.);
    MatrixXd V1(V.rows() + boneP.rows(), 3);
    V1 << V, boneP;
    cout << "<Tetrahedralization>" << endl;
    MatrixXd VT;
    MatrixXi TT, FT;
    igl::copyleft::tetgen::tetrahedralize(V1, F, "pYq", VT, TT, FT);

    cout << "<Calculate Bone Weights>" << endl;
    MatrixXd bc;
    VectorXi b;
    igl::boundary_conditions(VT, TT, C, VectorXi(), BE0, MatrixXi(), b, bc);
    for (int i = 9; i < BE0.rows(); i++)
    {
        bc.col(0) += bc.col(i);
    }
    bc.conservativeResize(b.rows(), 9);

    cout << bc.rows() << " X " << bc.cols() << endl;
    igl::BBWData bbw_data;
    bbw_data.active_set_params.max_iter = 10;
    bbw_data.verbosity = 2;
    MatrixXd WT;
    if (!igl::bbw(VT, TT, b, bc, bbw_data, WT))
        return EXIT_FAILURE;
    
    cout << "<Calculate Joint Weights>" << endl;
    vector<int> jn(C.rows(), 0);
    jn[4] = 1;
    jn[5] = 2;
    jn[6] = 3;
    jn[9] = 4;
    jn[10] = 5;
    jn[11] = 6;
    b.resize(C.rows());
    for (int i = 0; i < C.rows(); i++)
        b(i) = V.rows() + i;
    bc = MatrixXd::Zero(C.rows(), 7);
    for (int i = 0; i < C.rows(); i++)
        bc(i, jn[i]) = 1;
    cout << bc.rows() << " X " << bc.cols() << endl;
    MatrixXd WTJ;
    if (!igl::bbw(VT, TT, b, bc, bbw_data, WTJ))
        return EXIT_FAILURE;

    //matching between tetra & ply
    cout << "matching between tetra & ply.." << flush;
    auto grid = GenerateGrid(VT);
    int count(0);
    W.resize(V.rows(), BE.rows());
    MatrixXd WJ(V.rows(), 7);
    for (int i = 0; i < V.rows(); i++)
    {
        int x = floor(V(i, 0) + 0.5);
        int y = floor(V(i, 1) + 0.5);
        int z = floor(V(i, 2) + 0.5);
        auto key = make_tuple(x, y, z);
        for (int n : grid[key])
        {
            if (fabs(V(n, 0) - VT(i, 0)) > 0.01)
                continue;
            if (fabs(V(n, 1) - VT(i, 1)) > 0.01)
                continue;
            if (fabs(V(n, 2) - VT(i, 2)) > 0.01)
                continue;
            W.row(i) = WT.row(n);
            WJ.row(i) = WTJ.row(n);
            cout << "\rmatching between tetra & ply.." << ++count << "/" << V.rows() << flush;
            break;
        }
    }
    igl::normalize_row_sums(W, W);
    igl::normalize_row_sums(WJ, WJ);
    smoothMM.resize(V.rows(), 3);
    smoothMM.col(0) = WJ.col(1).array() + WJ.col(4).array();
    smoothMM.col(1) = WJ.col(2).array().pow(3) + WJ.col(5).array().pow(3);
    smoothMM.col(2) = WJ.col(3).array().pow(3) + WJ.col(6).array().pow(3);
    cout << endl;

    double epsilon = 1e-5;
    smoothMM *= w_smooth;
    cleanWeights.clear();
    smoothM.clear();
    smoothM.resize(3);

    for (int i = 0; i < V.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (smoothMM(i, j) < epsilon)
                continue;
            map<int, double> coeff;
            coeff[i] = 1 - smoothMM(i, j);
            double f = smoothMM(i, j) / (double)adjacent[i].size();
            for (int j : adjacent[i])
                coeff[j] = f;
            smoothM[j][i] = coeff;
        }

        map<int, double> weight;
        for (int j = 0; j < W.cols(); j++)
        {
            if (W(i, j) < 1e-4)
                continue;
            weight[j] = W(i, j);
        }
        cleanWeights.push_back(weight);
    }
    return true;
}
bool PhantomAnimator::ReadW(string prefix)
{
    if (!igl::readDMAT(prefix + ".W", W))
        return false;
    if (!igl::readDMAT(prefix + ".S", smoothMM))
        return false;

    cleanWeights.clear();
    smoothM.clear();
    smoothM.resize(3);

    double epsilon = 1e-5;
    for (int i = 0; i < V.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (smoothMM(i, j) < epsilon)
                continue;
            map<int, double> coeff;
            coeff[i] = 1 - smoothMM(i, j);
            double f = smoothMM(i, j) / (double)adjacent[i].size();
            for (int j : adjacent[i])
                coeff[j] = f;
            smoothM[j][i] = coeff;
        }

        map<int, double> weight;
        for (int j = 0; j < W.cols(); j++)
        {
            if (W(i, j) < 1e-4)
                continue;
            weight[j] = W(i, j);
        }
        cleanWeights.push_back(weight);
    }

    return true;
}

#include "igl/forward_kinematics.h"
#include "igl/lbs_matrix.h"
void PhantomAnimator::Animate(RotationList vQ, bool fixOthers)
{
    MatrixXd M, T;
    if (fixOthers)
    {
        T.resize(BE.rows() * 4, 3);
        for (int i = 0; i < BE.rows(); i++)
        {
            Affine3d a = Translation3d(Vector3d(C.row(BE(i, 0)))) * vQ[i] * Translation3d(-Vector3d(C.row(BE(i, 0))));
            T.block(i * 4, 0, 4, 3) = a.matrix().transpose().block(0, 0, 4, 3);
        }
    }
    else
    {
        igl::forward_kinematics(C, BE, P, vQ, T);
        for (int e = 0; e < BE.rows(); e++)
        {
            C.row(BE(e, 1)) = C.row(BE(e, 1)).homogeneous() * T.block(e * 4, 0, 4, 3); // a * Vector3d(C.row(BE(e, 1)));
        }
    }

    igl::lbs_matrix(V, W, M);
    V = M * T;
}

void PhantomAnimator::AnimateDQS(RotationList vQ)
{
    MatrixXd T;
    igl::forward_kinematics(C, BE, P, vQ, T);

    vector<Vector3d> vT;
    igl::forward_kinematics(C, BE, P, vQ, vQ, vT);
    myDqs(V, cleanWeights, vQ, vT, V);

    for (int e = 0; e < BE.rows(); e++)
    {
        C.row(BE(e, 1)) = C.row(BE(e, 1)).homogeneous() * T.block(e * 4, 0, 4, 3);
    }
}
