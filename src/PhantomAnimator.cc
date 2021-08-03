#include "PhantomAnimator.hh"
#include "igl/opengl/glfw/Viewer.h"
PhantomAnimator::PhantomAnimator() {}

PhantomAnimator::PhantomAnimator(string prefix)
    : name(prefix)
{
    ReadFiles(prefix);
    ReadOBJ(prefix);

    cout<<"<Initial bone volumes>"<<endl;
    for(auto iter:Fnum)
    {
        double vol;
        if(iter.first<100) vol = CalculateVolume(V0, F.block(iter.second.first, 0, iter.second.second, 3));
        else vol = CalculateVolume(VI0, FI.block(iter.second.first, 0, iter.second.second, 3));
        volumes[iter.first] = vol;
        cout<<iter.first<<" : "<<vol<<endl;
    }
    PreComputeAdjacent();
    //CalculateWeights(0.1);
}

bool PhantomAnimator::ReadFiles(string prefix)
{
    cout << "Read " + prefix + ".tgf" << endl;
    igl::readTGF(prefix + ".tgf", C, BE0);
    BE = BE0.block(0, 0, 9, 2);
    C0 = C;
    igl::directed_edge_parents(BE, P);
    //cout << "Read " + prefix + ".ply" << endl;
    // igl::readPLY(prefix + ".ply", V, F);
    // V0 = V;

    //SetSimpleF();

    return true;
}

// bool PhantomAnimator::CalculateWeights(double w_smooth)
// {
//     //perform BBW
//     MatrixXd boneP = GenerateBonePoints(C, BE0, 1.);
//     MatrixXd V1(V.rows() + boneP.rows(), 3);
//     V1 << V, boneP;
//     cout << "<Tetrahedralization>" << endl;
//     MatrixXd VT;
//     MatrixXi TT, FT;
//     igl::copyleft::tetgen::tetrahedralize(V1, F, "pYq", VT, TT, FT);

//     cout << "<Calculate Bone Weights>" << endl;
//     MatrixXd bc;
//     VectorXi b;
//     igl::boundary_conditions(VT, TT, C, VectorXi(), BE0, MatrixXi(), b, bc);
//     for (int i = 9; i < BE0.rows(); i++)
//     {
//         bc.col(0) += bc.col(i);
//     }
//     bc.conservativeResize(b.rows(), 9);

//     cout << bc.rows() << " X " << bc.cols() << endl;
//     igl::BBWData bbw_data;
//     bbw_data.active_set_params.max_iter = 10;
//     bbw_data.verbosity = 2;
//     MatrixXd WT;
//     if (!igl::bbw(VT, TT, b, bc, bbw_data, WT))
//         return EXIT_FAILURE;

//     cout << "<Calculate Joint Weights>" << endl;
//     vector<int> jn(C.rows(), 0);
//     jn[4] = 1;
//     jn[5] = 2;
//     jn[6] = 3;
//     jn[9] = 4;
//     jn[10] = 5;
//     jn[11] = 6;
//     b.resize(C.rows());
//     for (int i = 0; i < C.rows(); i++)
//         b(i) = V.rows() + i;
//     bc = MatrixXd::Zero(C.rows(), 7);
//     for (int i = 0; i < C.rows(); i++)
//         bc(i, jn[i]) = 1;
//     cout << bc.rows() << " X " << bc.cols() << endl;
//     MatrixXd WTJ;
//     if (!igl::bbw(VT, TT, b, bc, bbw_data, WTJ))
//         return EXIT_FAILURE;

//     //matching between tetra & ply
//     cout << "matching between tetra & ply.." << flush;
//     auto grid = GenerateGrid(VT);
//     int count(0);
//     W.resize(V.rows(), BE.rows());
//     MatrixXd WJ(V.rows(), 7);
//     for (int i = 0; i < V.rows(); i++)
//     {
//         int x = floor(V(i, 0) + 0.5);
//         int y = floor(V(i, 1) + 0.5);
//         int z = floor(V(i, 2) + 0.5);
//         auto key = make_tuple(x, y, z);
//         for (int n : grid[key])
//         {
//             if (fabs(V(n, 0) - VT(i, 0)) > 0.01)
//                 continue;
//             if (fabs(V(n, 1) - VT(i, 1)) > 0.01)
//                 continue;
//             if (fabs(V(n, 2) - VT(i, 2)) > 0.01)
//                 continue;
//             W.row(i) = WT.row(n);
//             WJ.row(i) = WTJ.row(n);
//             cout << "\rmatching between tetra & ply.." << ++count << "/" << V.rows() << flush;
//             break;
//         }
//     }
//     igl::normalize_row_sums(W, W);
//     igl::normalize_row_sums(WJ, WJ);
//     smoothMM.resize(V.rows(), 3);
//     smoothMM.col(0) = WJ.col(1).array() + WJ.col(4).array();
//     smoothMM.col(1) = WJ.col(2).array().pow(3) + WJ.col(5).array().pow(3);
//     smoothMM.col(2) = WJ.col(3).array().pow(3) + WJ.col(6).array().pow(3);
//     cout << endl;

//     int selected = 0;
//     Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
//     igl::opengl::glfw::Viewer viewer;
//     viewer.data().set_mesh(V, F);
//     viewer.data().set_data(smoothMM.col(0));
//     viewer.data().set_edges(C, BE, sea_green);
//     viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods) -> bool
//     {
//         switch (key)
//         {
//         case '.':
//             selected++;
//             selected = std::min(std::max(selected, 0), (int)W.cols() - 1);
//             viewer.data().set_data(W.col(selected));
//             break;
//         case ',':
//             selected--;
//             selected = std::min(std::max(selected, 0), (int)W.cols() - 1);
//             viewer.data().set_data(W.col(selected));
//             break;
//         case '[':
//             viewer.data().set_data(smoothMM.col(1));
//             break;
//         }
//         return true;
//     };
//     viewer.launch();

//     double epsilon = 1e-5;
//     smoothMM *= w_smooth;
//     cleanWeights.clear();
//     smoothM.clear();
//     smoothM.resize(3);

//     for (int i = 0; i < V.rows(); i++)
//     {
//         for (int j = 0; j < 3; j++)
//         {
//             if (smoothMM(i, j) < epsilon)
//                 continue;
//             map<int, double> coeff;
//             coeff[i] = 1 - smoothMM(i, j);
//             double f = smoothMM(i, j) / (double)adjacent[i].size();
//             for (int j : adjacent[i])
//                 coeff[j] = f;
//             smoothM[j][i] = coeff;
//         }

//         map<int, double> weight;
//         for (int j = 0; j < W.cols(); j++)
//         {
//             if (W(i, j) < 1e-4)
//                 continue;
//             weight[j] = W(i, j);
//         }
//         cleanWeights.push_back(weight);
//     }
//     return true;
// }
// bool PhantomAnimator::ReadW(string prefix)
// {
//     if (!igl::readDMAT(prefix + ".W", W))
//         return false;
//     if (!igl::readDMAT(prefix + ".S", smoothMM))
//         return false;

//     cleanWeights.clear();
//     smoothM.clear();
//     smoothM.resize(3);

//     double epsilon = 1e-5;
//     for (int i = 0; i < V.rows(); i++)
//     {
//         for (int j = 0; j < 3; j++)
//         {
//             if (smoothMM(i, j) < epsilon)
//                 continue;
//             map<int, double> coeff;
//             coeff[i] = 1 - smoothMM(i, j);
//             double f = smoothMM(i, j) / (double)adjacent[i].size();
//             for (int j : adjacent[i])
//                 coeff[j] = f;
//             smoothM[j][i] = coeff;
//         }

//         map<int, double> weight;
//         for (int j = 0; j < W.cols(); j++)
//         {
//             if (W(i, j) < 1e-4)
//                 continue;
//             weight[j] = W(i, j);
//         }
//         cleanWeights.push_back(weight);
//     }
//     // ifstream ifs(prefix + ".W");
//     // if (!ifs.is_open())
//     //     return false;
//     // string aLine;
//     // cleanWeights.clear();
//     // while (getline(ifs, aLine))
//     // {
//     //     if (aLine.length() == 0)
//     //         continue;
//     //     stringstream ss(aLine);
//     //     int num;
//     //     ss >> num;
//     //     map<int, double> weight;
//     //     for (int i = 0; i < num; i++)
//     //     {
//     //         int idx;
//     //         double w;
//     //         ss >> idx >> w;
//     //         weight[idx] = w;
//     //     }        igl::writeDMAT(prefix + ".W", W, false);

//     //     cleanWeights.push_back(weight);
//     // }
//     // ifs.close();
//     return true;
// }

#include "igl/forward_kinematics.h"
#include "igl/lbs_matrix.h"
void PhantomAnimator::Animate(RotationList vQ)
{
    MatrixXd M, T;
    igl::forward_kinematics(C, BE, P, vQ, T);
    // igl::lbs_matrix(V, W, M);
    // V = M * T;

    vector<Vector3d> vT;
    igl::forward_kinematics(C, BE, P, vQ, vQ, vT);

    myDqs(V, cleanWeights, vQ, vT, V);

    for (int e = 1; e < BE.rows(); e++)
    {
        C.row(BE(e, 1)) = C.row(BE(e, 1)).homogeneous() * T.block(e * 4, 0, 4, 3);
    }
}

// void PhantomAnimator::AnimateDQS(RotationList vQ)
// {
//     MatrixXd T;
//     igl::forward_kinematics(C, BE, P, vQ, T);

//     vector<Vector3d> vT;
//     igl::forward_kinematics(C, BE, P, vQ, vQ, vT);
//     myDqs(V, cleanWeights, vQ, vT, V);

//     for (int e = 0; e < BE.rows(); e++)
//     {
//         C.row(BE(e, 1)) = C.row(BE(e, 1)).homogeneous() * T.block(e * 4, 0, 4, 3);
//         //        if(e>0 && e<9)
//         // {
//         //     plyV[e-1] = plyV[e-1].rowwise().homogeneous() * T.block(e * 4, 0, 4, 3);
//         // }
//     }
// }

int PhantomAnimator::DetectIntersections(MatrixXd &Color, int a, int b)
{
    if (b < 0)
        b = P(a);
    MatrixXi IF;
    MatrixXi Fa = F.block(Fnum[a].first, 0, Fnum[a].second, 3);
    MatrixXi Fb = F.block(Fnum[b].first, 0, Fnum[b].second, 3);
    if (igl::copyleft::cgal::intersect_other(V, Fa, V, Fb, false, IF))
    {
        IF.col(0) = IF.col(0).array() + Fnum[a].first;
        IF.col(1) = IF.col(1).array() + Fnum[b].first;
        vector<int> Kvec(IF.data(), IF.data() + IF.size());
        sort(Kvec.begin(), Kvec.end());
        Kvec.erase(unique(Kvec.begin(), Kvec.end()), Kvec.end());
        MatrixXi K = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(Kvec.data(), Kvec.size());
        MatrixXd R = RowVector3d(1., 0.3, 0.3).replicate(K.rows(), 1);
        igl::slice_into(R, K, 1, Color);
        return IF.rows();
    }
    else
        return 0;
}

int PhantomAnimator::DetectIntersections(MatrixXd &Color)
{
    int count(0);
    for (int i = 1; i < BE.rows(); i++)
    {
        if (P(i) < 0)
            count += DetectIntersections(Color, i, 0);
        else
            count += DetectIntersections(Color, i);
    }
    return count;
}

void PhantomAnimator::CalculateElbowW(int a, Vector3d axis, double inner2, double outer2)
{
    //check intersection
    int b = P(a);
    MatrixXd Va = V.block(Vnum[a].first, 0, Vnum[a].second, 3);
    MatrixXd Vb = V.block(Vnum[b].first, 0, Vnum[b].second, 3);
    MatrixXi Fa = F.block(Fnum[a].first, 0, Fnum[a].second, 3).array() - Vnum[a].first;
    MatrixXi Fb = F.block(Fnum[b].first, 0, Fnum[b].second, 3).array() - Vnum[b].first;
    vector<int> vecA, vecB;
    Vector3d center = C.row(BE(a, 0));
    for (double theta = 0;; theta += 5)
    {
        Affine3d a = Translation3d(center) * AngleAxisd(theta / 180. * PI, axis) * Translation3d(-center);
        MatrixXd T = a.matrix().transpose().block(0, 0, 4, 3);
        MatrixXd Va1 = Va.rowwise().homogeneous() * T;
        MatrixXi IF;
        if (!igl::copyleft::cgal::intersect_other(Va1, Fa, Vb, Fb, false, IF))
            continue;

        Va = Va1;
        MatrixXi IFa = igl::slice(Fa, IF.col(0));
        MatrixXi IFb = igl::slice(Fb, IF.col(1));
        vector<int> vecA1(IFa.data(), IFa.data() + IFa.size());
        vector<int> vecB1(IFb.data(), IFb.data() + IFb.size());
        sort(vecA1.begin(), vecA1.end());
        sort(vecB1.begin(), vecB1.end());
        vecA1.erase(unique(vecA1.begin(), vecA1.end()), vecA1.end());
        vecB1.erase(unique(vecB1.begin(), vecB1.end()), vecB1.end());
        vecA = vecA1; vecB = vecB1;
        // center = Vector3d::Zero();
        // for (int i : vecA)
        //     center += Va.row(i);
        // for (int i : vecB)
        //     center += Vb.row(i);
        // center = center / (double)(vecA.size() + vecB.size());
        break;
    }

    MatrixXd VT;
    MatrixXi TT, FT;
    igl::copyleft::tetgen::tetrahedralize(Va, Fa, "pYq", VT, TT, FT);
   
    vector<int> half, rigid;
    for (int i = 0; i < Va.rows(); i++)
    {
        double dist2 = (Va.row(vecA[0]).transpose() - Vector3d(VT.row(i))).squaredNorm();
        for(int n:vecA) dist2 = min(dist2, (Va.row(n).transpose() - Vector3d(VT.row(i))).squaredNorm());
        if (dist2 < inner2)
            half.push_back(i);
        else if (dist2 > outer2)
            rigid.push_back(i);
    }

    VectorXi B(half.size() + rigid.size());
    MatrixXd BC = MatrixXd::Zero(half.size() + rigid.size(), 2);
    for (int i = 0; i < half.size(); i++)
    {
        B(i) = half[i];
        BC(i, 0) = 0.5;
        BC(i, 1) = 0.5;
    }
    for (int i = 0; i < rigid.size(); i++)
    {
        B(half.size() + i) = rigid[i];
        BC(half.size() + i, 1) = 1;
    }
    igl::BBWData bbw_data;
    bbw_data.active_set_params.max_iter = 10;
    bbw_data.verbosity = 2;
    MatrixXd WT;
    igl::bbw(VT, TT, B, BC, bbw_data, WT);
    W.block(Vnum[a].first, a, Vnum[a].second, 1) = WT.block(0, 1, Vnum[a].second, 1);
    W.block(Vnum[a].first, b, Vnum[a].second, 1) = WT.block(0, 0, Vnum[a].second, 1);

    igl::copyleft::tetgen::tetrahedralize(Vb, Fb, "pYq", VT, TT, FT);
    half.clear();
    rigid.clear();
    for (int i = 0; i < Vb.rows(); i++)
    {
        double dist2 = (Vb.row(vecB[0]).transpose() - Vector3d(VT.row(i))).squaredNorm();
        for(int n:vecB) dist2 = min(dist2, (Vb.row(n).transpose() - Vector3d(VT.row(i))).squaredNorm());
        if (dist2 < inner2)
            half.push_back(i);
        else if (dist2 > outer2)
            rigid.push_back(i);
    }

    B.resize(half.size() + rigid.size());
    BC = MatrixXd::Zero(half.size() + rigid.size(), 2);
    for (int i = 0; i < half.size(); i++)
    {
        B(i) = half[i];
        BC(i, 0) = 0.5;
        BC(i, 1) = 0.5;
    }
    for (int i = 0; i < rigid.size(); i++)
    {
        B(half.size() + i) = rigid[i];
        BC(half.size() + i, 1) = 1;
    }
    igl::bbw(VT, TT, B, BC, bbw_data, WT);
    W.block(Vnum[b].first, b, Vnum[b].second, 1) = WT.block(0, 1, Vnum[b].second, 1);
    W.block(Vnum[b].first, a, Vnum[b].second, 1) = WT.block(0, 0, Vnum[b].second, 1);
}

double PhantomAnimator::CalculateShoulderW(int a, Vector3d axis, double inner2, double outer2)
{
    //check intersection
    int b = P(a);
    if (b < 0)
        b = 0;
    MatrixXd Va = V.block(Vnum[a].first, 0, Vnum[a].second, 3);
    MatrixXd Vb = V.block(Vnum[b].first, 0, Vnum[b].second, 3);
    MatrixXi Fa = F.block(Fnum[a].first, 0, Fnum[a].second, 3).array() - Vnum[a].first;
    MatrixXi Fb = F.block(Fnum[b].first, 0, Fnum[b].second, 3).array() - Vnum[b].first;
    vector<int> vecA;
    Vector3d center = C.row(BE(a, 0));
    double theta;
    for (theta = 0;; theta += 5)
    {
        Affine3d a = Translation3d(center) * AngleAxisd(theta / 180. * PI, axis) * Translation3d(-center);
        MatrixXd T = a.matrix().transpose().block(0, 0, 4, 3);
        MatrixXd Va1 = Va.rowwise().homogeneous() * T;
        MatrixXi IF;
        if (!igl::copyleft::cgal::intersect_other(Va1, Fa, Vb, Fb, false, IF))
            continue;

        Va = Va1;
        MatrixXi IFa = igl::slice(Fa, IF.col(0));
        vector<int> vecA1(IFa.data(), IFa.data() + IFa.size());
        sort(vecA1.begin(), vecA1.end());
        vecA1.erase(unique(vecA1.begin(), vecA1.end()), vecA1.end());
        vecA = vecA1;
        break;
    }
    MatrixXd VT;
    MatrixXi TT, FT;
    igl::copyleft::tetgen::tetrahedralize(Va, Fa, "pYq", VT, TT, FT);

    vector<int> half, rigid;
    for (int i = 0; i < Va.rows(); i++)
    {
        double dist2 = (Va.row(vecA[0]).transpose() - Vector3d(VT.row(i))).squaredNorm();
        for(int n:vecA) dist2 = min(dist2, (Va.row(n).transpose() - Vector3d(VT.row(i))).squaredNorm());
        if (dist2 < inner2)
            half.push_back(i);
        else if (dist2 > outer2)
            rigid.push_back(i);
    }

    VectorXi B(half.size() + rigid.size());
    MatrixXd BC = MatrixXd::Zero(half.size() + rigid.size(), 2);
    for (int i = 0; i < half.size(); i++)
    {
        B(i) = half[i];
        BC(i, 0) = 1;
    }
    for (int i = 0; i < rigid.size(); i++)
    {
        B(half.size() + i) = rigid[i];
        BC(half.size() + i, 1) = 1;
    }
    igl::BBWData bbw_data;
    bbw_data.active_set_params.max_iter = 10;
    bbw_data.verbosity = 2;
    MatrixXd WT;
    igl::bbw(VT, TT, B, BC, bbw_data, WT);
    W.block(Vnum[a].first, 0, Vnum[a].second, BE.rows()) = MatrixXd::Zero(Vnum[a].second, BE.rows());
    W.block(Vnum[a].first, a, Vnum[a].second, 1) = WT.block(0, 1, Vnum[a].second, 1);
    W.block(Vnum[a].first, b, Vnum[a].second, 1) = WT.block(0, 0, Vnum[a].second, 1);
    return theta;
}

void PhantomAnimator::CalculateCleanWeights(){
    cleanWeights.clear();
    double e(1e-5);
    for(int i=0;i<W.rows();i++)
    {
        double sum(0.);
        map<int, double> weight;
        for(int j=0;j<W.cols();j++)
        {
            if(W(i,j)<e) continue;
            sum += W(i,j);
            weight[j] = W(i,j);
        }
        for(auto &iter:weight) iter.second/=sum;
        cleanWeights.push_back(weight);
    }
}
void PhantomAnimator::ReadOBJ(string fName)
{
    ifstream ifs(fName + "_bone.obj");
    vector<Vector3d> objV_vec, objV_vec1;
    vector<Vector3i> objF_vec, objF_vec1;
    int prevV(0), prevF(0);
    int shell(0);
    int outerV(0), outerF(0);
    while (!ifs.eof())
    {
        string aLine;
        getline(ifs, aLine);
        stringstream ss(aLine);
        string first;
        ss >> first;
        if (first == "v")
        {
            double x, y, z;
            ss >> x >> y >> z;
            objV_vec.push_back(Vector3d(x, y, z));
        }
        else if (first == "g")
        {
            if (objF_vec.size())
            {   
                Fnum[shell] = make_pair(prevF-outerF, objF_vec.size() - prevF);
                prevF = objF_vec.size();
            }
            ss >> shell;
            if(outerV==0 && shell>100) 
            {
                outerF = prevF;
                outerV = prevV;
            }
            Vnum[shell] = make_pair(prevV-outerV, objV_vec.size() - prevV);
            prevV = objV_vec.size();
        }
        else if (first == "f")
        {
            int a, b, c;
            ss >> a >> b >> c;
            objF_vec.push_back(Vector3i(a - 1, b - 1, c - 1));
        }
    }
    Fnum[shell] = make_pair(prevF-outerF, objF_vec.size() - prevF);

    //outer shells should come first
    V.resize(outerV, 3);
    for (int i = 0; i < V.rows(); i++)
        V.row(i) = objV_vec[i];
    V0 = V;
    F.resize(outerF, 3);
    for (int i = 0; i < F.rows(); i++)
        F.row(i) = objF_vec[i];

    //inner structure
    VI.resize(objV_vec.size()-outerV,3);
    for(int i=0;i<VI.rows();i++)
        VI.row(i) = objV_vec[outerV+i];
    VI0 = VI;
    FI.resize(objF_vec.size()-outerF,3);
    for(int i=0;i<FI.rows();i++)
        FI.row(i) = objF_vec[outerF+i];
    FI = FI.array() - outerV;        
}