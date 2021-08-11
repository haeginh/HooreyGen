#include "PhantomAnimator.hh"
PhantomAnimator::PhantomAnimator() {}

PhantomAnimator::PhantomAnimator(string _prefix)
    : prefix(_prefix)
{
    ReadTGF();
    MergeOBJ();
    PreComputeAdjacent();
    volume0 = CalculateVolume(V0, F);

    //CalculateWeights(0.1);
}

PhantomAnimator::~PhantomAnimator()
{
}

bool PhantomAnimator::ReadTGF()
{
    cout << "Read " + prefix + ".tgf" << endl;
    igl::readTGF(prefix + ".tgf", C, BE0);
    BE = BE0.block(0, 0, 9, 2);
    C0 = C;
    igl::directed_edge_parents(BE, P);

    //SetSimpleF();

    return true;
}

bool PhantomAnimator::CalculateWeights(double w_smooth, bool innerOnly)
{
    //perform BBW
    MatrixXd boneP = GenerateBonePoints(C, BE0, 1.);

    cout << "<Tetrahedralization>" << endl;
    MatrixXd VT, V1;
    MatrixXi TT, FT, F1;
    if (innerOnly)
    {
        V1.resize(V.rows() + boneP.rows(), 3);
        V1 << V, boneP;
        F1 = F;
    }
    else
    {
        V1.resize(Vo.rows() + boneP.rows(), 3);
        V1 << Vo, boneP;
        F1 = Fo;
    }

    if (igl::copyleft::tetgen::tetrahedralize(V1, F1, "p/0.001YqO/3", VT, TT, FT) != 0)
        return EXIT_FAILURE;

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

    if (!innerOnly)
    {
        cout << endl
             << "<Extract deformable vertices of outer skin>" << endl;
        MatrixXd F_normals, V_normals;
        igl::per_face_normals(V, F, F_normals);
        igl::per_vertex_normals(V, F, F_normals, V_normals);
        MatrixXd offset = V + V_normals * 0.02;
        //vector<bool> chkRST(offset.rows(), false);
        vector<int> deformable;
        for (int i = 0; i < Vo.rows(); i++)
        {
            if (WT(i, 4) > 0.9)
            {
                handsO[4].push_back(i);
                continue;
            }
            else if (WT(i, 8) > 0.9)
            {
                handsO[8].push_back(i);
                continue;
            }
            else if (WT(i, 0) > 0.9)
                continue;
            deformable.push_back(i);
            cout << "\rDeformable..." << deformable.size() << flush;
            // double dist2(1);
            // int id(-1);
            // for (int j = 0; j < offset.rows(); j++)
            // {
            //     if (chkRST[j])
            //         continue;
            //     double _dist2 = (offset.row(j) - Vo.row(i)).squaredNorm();
            //     if (_dist2 > dist2)
            //         continue;
            //     id = j;
            //     dist2 = _dist2;
            // }
            // chkRST[id] = true;
            // o2rst[i] = id;
            // cout << "\rExtracting..." << o2rst.size() << flush;
        }
        cout << " (hands: " << handsO[4].size() << " / " << handsO[8].size() << ")" << endl;
        D.resize(deformable.size());
        for (int i = 0; i< D.rows(); i++)
            D(i) = deformable[i];
        MatrixXd DV = igl::slice(Vo, D, 1);
        VectorXd dist;
        VectorXi I; MatrixXd T;
        igl::point_mesh_squared_distance(DV, offset, F, dist, I, T);
        typedef Eigen::Triplet<double> TRI;
        std::vector<TRI> coeff;
        for(int i= 0; i<DV.rows();i++)
        {
            Vector3d s0 = offset.row(F(i,0))-T.row(i);
            Vector3d s1 = offset.row(F(i,1))-T.row(i);
            Vector3d s2 = offset.row(F(i,2))-T.row(i);
            double a0 = s1.cross(s2).norm();
            double a1 = s0.cross(s2).norm();
            double a2 = s1.cross(s0).norm();
            double sum = a0 + a1 + a2;
            coeff.push_back(TRI(i, F(i, 0), a0/sum));
            coeff.push_back(TRI(i, F(i, 1), a1/sum));
            coeff.push_back(TRI(i, F(i, 2), a2/sum));
        }
        baryO.resize(DV.rows(), V.rows());
        baryO.setFromTriplets(coeff.begin(), coeff.end());
    }

    cout << endl
         << "<Calculate Joint Weights>" << endl;
    vector<int> jn(C.rows(), 0);
    jn[3] = 7;
    jn[4] = 1;
    jn[5] = 2;
    jn[6] = 3;
    jn[8] = 8;
    jn[9] = 4;
    jn[10] = 5;
    jn[11] = 6;
    b.resize(C.rows());
    for (int i = 0; i < C.rows(); i++)
        if(innerOnly) b(i) = V.rows() + i;
        else b(i) = Vo.rows() + i;
    bc = MatrixXd::Zero(C.rows(), 9);
    for (int i = 0; i < C.rows(); i++)
        bc(i, jn[i]) = 1;
    cout << bc.rows() << " X " << bc.cols() << endl;
    MatrixXd WTJ;
    if (!igl::bbw(VT, TT, b, bc, bbw_data, WTJ))
        return EXIT_FAILURE;

    cout << endl
         << "Organizing and recording weights.." << flush;
    igl::normalize_row_sums(WT, WT);
    igl::normalize_row_sums(WTJ, WTJ);
    smoothMM.resize(VT.rows(), 4);
    smoothMM.col(0) = WTJ.col(1).array() + WTJ.col(4).array();
    smoothMM.col(1) = WTJ.col(2).array().pow(3) + WTJ.col(5).array().pow(3);
    smoothMM.col(2) = WTJ.col(3).array().pow(3) + WTJ.col(6).array().pow(3);
    smoothMM.col(3) = WTJ.col(7).array().pow(3) + WTJ.col(8).array().pow(3);
    smoothMM *= w_smooth;
    cout << "done" << endl;

    if (innerOnly)
    {
        W = WT.block(0, 0, V.rows(), WT.cols());
        smoothMM = smoothMM.block(0, 0, V.rows(), smoothMM.cols());
    }
    else
    {
        //write
        igl::writeMESH(prefix + ".mesh", VT, TT, FT);
        igl::writeDMAT(prefix + ".W", WT, false);
        igl::writeDMAT(prefix + ".S", smoothMM, false);

        //set inner
        cout << endl
             << "<Assignemt of weights for inner vertices>" << endl;
        auto baryCoord = GenerateBarycentricCoord(VT, TT, V);
        SparseMatrix<double> bary = GenerateBarySparse(baryCoord, VT.rows());
        W = bary * WT;
        smoothMM = bary * smoothMM;
        cout << endl;
    }
    double epsilon = 1e-5;
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
    return EXIT_SUCCESS;
}

bool PhantomAnimator::ReadW()
{
    if (!igl::readDMAT(prefix + ".W", W))
        return false;
    if (!igl::readDMAT(prefix + ".S", smoothMM))
        return false;
    MatrixXd VT;
    MatrixXi TT, FT;
    if (!igl::readMESH(prefix + ".mesh", VT, TT, FT))
        return false;

    cout << endl
         << "<Extract deformable vertices of outer skin>" << endl;
    // MatrixXd F_normals, V_normals;
    // igl::per_face_normals(V, F, F_normals);
    // igl::per_vertex_normals(V, F, F_normals, V_normals);
    //vector<bool> chkRST(offset.rows(), false);
    vector<int> deformable;
        for (int i = 0; i < Vo.rows(); i++)
        {
            if (W(i, 4) > 0.9)
            {
                handsO[4].push_back(i);
                continue;
            }
            else if (W(i, 8) > 0.9)
            {
                handsO[8].push_back(i);
                continue;
            }
            else if (W(i, 0) > 0.9)
                continue;
            deformable.push_back(i);
            cout << "\rDeformable..." << deformable.size() << flush;
            // double dist2(1);
            // int id(-1);
            // for (int j = 0; j < offset.rows(); j++)
            // {
            //     if (chkRST[j])
            //         continue;
            //     double _dist2 = (offset.row(j) - Vo.row(i)).squaredNorm();
            //     if (_dist2 > dist2)
            //         continue;
            //     id = j;
            //     dist2 = _dist2;
            // }
            // chkRST[id] = true;
            // o2rst[i] = id;
            // cout << "\rExtracting..." << o2rst.size() << flush;
        }
        cout << " (hands: " << handsO[4].size() << " / " << handsO[8].size() << ")" << endl;
        D.resize(deformable.size());
        for (int i = 0; i< D.rows(); i++)
            D(i) = deformable[i];
        MatrixXd offset = Vo - GetNormal(Vo, Fo) * 0.2;
        MatrixXd DV = igl::slice(offset, D, 1);
        VectorXd dist;
        VectorXi I; MatrixXd T;
        igl::point_mesh_squared_distance(DV, V, F, dist, I, T);
        typedef Eigen::Triplet<double> TRI;
        std::vector<TRI> coeff;
        for(int i= 0; i<DV.rows();i++)
        {
            Vector3d s0 = offset.row(F(I(i),0))-T.row(i);
            Vector3d s1 = offset.row(F(I(i),1))-T.row(i);
            Vector3d s2 = offset.row(F(I(i),2))-T.row(i);
            double a0 = s1.cross(s2).norm();
            double a1 = s0.cross(s2).norm();
            double a2 = s1.cross(s0).norm();
            double sum = 1.f/(a0 + a1 + a2);
            if(a0*a1*a2<0) {
                if(a0>a1 && a0>a2) coeff.push_back(TRI(i, F(I(i), 0), 1)); continue;
                if(a1>a0 && a1>a2) coeff.push_back(TRI(i, F(I(i), 1), 1)); continue;
                if(a2>a1 && a2>a0) coeff.push_back(TRI(i, F(I(i), 2), 1)); continue;
            }
            coeff.push_back(TRI(i, F(I(i), 0), a0*sum));
            coeff.push_back(TRI(i, F(I(i), 1), a1*sum));
            coeff.push_back(TRI(i, F(I(i), 2), a2*sum));
            //cout<<dist(i)<<" "<<offset.row(F(I(i), 0))<<" / "<<offset.row(F(I(i), 1))<<" / "<<offset.row(F(I(i), 2))<<" // "<<a0/sum<<" "<<a1/sum<<" "<<a2/sum<<flush; getchar();
        }
        baryO.resize(DV.rows(), V.rows());
        baryO.setFromTriplets(coeff.begin(), coeff.end());
    
    // for (int i = 0; i < Vo.rows(); i++)
    // {
    //     if (W(i, 4) > 0.9)
    //     {
    //         handsO[4].push_back(i);
    //         continue;
    //     }
    //     else if (W(i, 8) > 0.9)
    //     {
    //         handsO[8].push_back(i);
    //         continue;
    //     }
    //     else if (W(i, 0) > 0.9)
    //         continue;
    //     double dist2(1);
    //     int id(-1);
    //     for (int j = 0; j < offset.rows(); j++)
    //     {
    //         if (chkRST[j])
    //             continue;
    //         double _dist2 = (offset.row(j) - Vo.row(i)).squaredNorm();
    //         if (_dist2 > dist2)
    //             continue;
    //         id = j;
    //         dist2 = _dist2;
    //     }
    //     chkRST[id] = true;
    //     o2rst[i] = id;
    //     cout << "\rExtracting..." << o2rst.size() << flush;
    // }
    // cout << " (hands: " << handsO[4].size() << " / " << handsO[8].size() << ")" << endl;

    //set inner
    cout << endl
         << "<Assignemt of weights for inner vertices>" << endl;
    auto baryCoord = GenerateBarycentricCoord(VT, TT, V);
    SparseMatrix<double> bary = GenerateBarySparse(baryCoord, VT.rows());
    W = bary * W;
    smoothMM = bary * smoothMM;
    cout << endl;

    double epsilon = 1e-5;
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
        for (int v : handsO[e])
            Vo.row(v) = Vo.row(v).homogeneous() * T.block(e * 4, 0, 4, 3);
    }
}

void PhantomAnimator::MergeOBJ()
{
    //whole
    cout << "Reading " + prefix + ".obj..." << flush;
    int prevV(0), id;
    int count(1);
    VV.clear();
    FF.clear();
    shellName.clear();
    ifstream ifs(prefix + ".obj");
    while (!ifs.eof())
    {
        string aLine, first;
        getline(ifs, aLine);
        stringstream ss(aLine);
        ss >> first;
        if (first == "v")
        {
            double x, y, z;
            ss >> x >> y >> z;
            VV[-1].push_back(Vector3d(x, y, z));
        }
        else if (first == "g")
        {
            string name;
            ss >> name;
            id = atoi(name.substr(0, name.find("_")).c_str());
            if (id == 10600 || id == 900 || id == 910 || id == 10400)
                id = -1;
            else
            {
                while (shellName.find(id) != shellName.end())
                    id += 10;
                shellName[id] = name;
            }
            count += prevV;
            prevV = VV[-1].size();
            VV[id] = VV[-1];
            VV[-1].clear();
            cout << "\rReading " + prefix + ".obj..." << name.substr(0, 10) << flush;
        }
        else if (first == "f")
        {
            if (id < 0)
                continue;
            int a, b, c;
            ss >> a >> b >> c;
            FF[id].push_back(Vector3i(a - count, b - count, c - count));
        }
    }
    ifs.close();
    VV.erase(-1);
    cout << "\rReading " + prefix + ".obj...done            " << endl;

    //set RST (11600) vertices
    V.resize(VV[11600].size(), 3);
    for (int i = 0; i < VV[11600].size(); i++)
        V.row(i) = VV[11600][i];
    V0 = V;
    F.resize(FF[11600].size(), 3);
    for (int i = 0; i < FF[11600].size(); i++)
        F.row(i) = FF[11600][i];

    //set outer skin (12200) vertices
    Vo.resize(VV[12200].size(), 3);
    for (int i = 0; i < VV[12200].size(); i++)
        Vo.row(i) = VV[12200][i];
    Fo.resize(FF[12200].size(), 3);
    for (int i = 0; i < FF[12200].size(); i++)
        Fo.row(i) = FF[12200][i];
    volumeO0 = CalculateVolume(Vo, Fo);

    //change bone
    ifstream ifs2(prefix + "_result.obj");
    map<int, vector<Vector3d>> v_tmp;
    map<int, vector<Vector3i>> f_tmp;
    prevV = 0;
    count = 1;
    while (!ifs2.eof())
    {
        string aLine, first;
        getline(ifs2, aLine);
        stringstream ss(aLine);
        ss >> first;
        if (first == "v")
        {
            double x, y, z;
            ss >> x >> y >> z;
            v_tmp[-1].push_back(Vector3d(x, y, z));
        }
        else if (first == "g")
        {
            count += prevV;
            ss >> id;
            v_tmp[id] = v_tmp[-1];
            prevV = v_tmp[-1].size();
            v_tmp[-1].clear();
        }
        else if (first == "f")
        {
            int a, b, c;
            ss >> a >> b >> c;
            f_tmp[id].push_back(Vector3i(a - count, b - count, c - count));
        }
    }
    ifs2.close();

    //organize bone obj
    int vNum(0), fNum(0);
    for (int i = 0; i < 9; i++)
    {
        vNum += v_tmp[i].size();
        fNum += f_tmp[i].size();
    }
    Vbone.resize(vNum, 3);
    Fbone.resize(fNum, 3);
    for (int i = 0, vNum = 0, fNum = 0; i < 9; i++)
    {
        for (Vector3i f : f_tmp[i])
            Fbone.row(fNum++) = f.array() + vNum;
        for (Vector3d v : v_tmp[i])
            Vbone.row(vNum++) = v;
    }

    //separate 1 and 5
    VV.erase(4500);
    FF.erase(4500);
    for (int i = 1; i < 6; i += 4)
    {
        vector<bool> chk(v_tmp[i].size(), false);
        vector<bool> chkF(f_tmp[i].size(), false);
        chk[0] = true;
        int num(1);
        while (num)
        {
            num = 0;
            for (int j = 0; j < f_tmp[i].size(); j++)
            {
                if (chkF[j])
                    continue;
                Vector3i f = f_tmp[i][j];
                if (!chk[f(0)] && !chk[f(1)] && !chk[f(2)])
                    continue;
                num++;
                chk[f(0)] = true;
                chk[f(1)] = true;
                chk[f(2)] = true;
                chkF[j] = true;
            }
        }
        double y0(-__DBL_MAX__), y1(-__DBL_MAX__); // max
        vector<Vector3d> v0, v1;
        vector<int> o2n;
        for (int j = 0; j < chk.size(); j++)
        {
            if (chk[j])
            {
                y0 = max(y0, v_tmp[i][j](1));
                o2n.push_back(v0.size());
                v0.push_back(v_tmp[i][j]);
            }
            else
            {
                y1 = max(y1, v_tmp[i][j](1));
                o2n.push_back(v1.size());
                v1.push_back(v_tmp[i][j]);
            }
        }
        int start = VV[4500].size();
        if (y0 > y1) // y0 is scapulae(0)
        {
            VV[4500].insert(VV[4500].end(), v0.begin(), v0.end());
            v_tmp[i] = v1;
            vector<Vector3i> f1;
            for (int j = 0; j < chkF.size(); j++)
            {
                Vector3i f = f_tmp[i][j];
                Vector3i f_new = Vector3i(o2n[f(0)], o2n[f(1)], o2n[f(2)]);
                if (chkF[j])
                    FF[4500].push_back(f_new.array() + start);
                else
                    f1.push_back(f_new);
            }
            f_tmp[i] = f1;
        }
        else
        {
            VV[4500].insert(VV[4500].end(), v1.begin(), v1.end());
            v_tmp[i] = v0;
            vector<Vector3i> f0;
            for (int j = 0; j < chkF.size(); j++)
            {
                Vector3i f = f_tmp[i][j];
                Vector3i f_new = Vector3i(o2n[f(0)], o2n[f(1)], o2n[f(2)]);
                if (!chkF[j])
                    FF[4500].push_back(f_new.array() + start);
                else
                    f0.push_back(f_new);
            }
            f_tmp[i] = f0;
        }
    }

    for (int i = 1; i < 5; i++)
    {
        int n = v_tmp[i].size();
        v_tmp[i].insert(v_tmp[i].end(), v_tmp[i + 4].begin(), v_tmp[i + 4].end());
        v_tmp.erase(i + 4);
        for (int j = 0; j < f_tmp[i + 4].size(); j++)
            f_tmp[i].push_back(f_tmp[i + 4][j].array() + n);
    }

    v_tmp[2400] = v_tmp[1];
    v_tmp.erase(1);
    f_tmp[2400] = f_tmp[1];
    f_tmp.erase(1);
    v_tmp[1300] = v_tmp[2];
    v_tmp.erase(2);
    f_tmp[1300] = f_tmp[2];
    f_tmp.erase(2);
    v_tmp[1900] = v_tmp[3];
    v_tmp.erase(3);
    f_tmp[1900] = f_tmp[3];
    f_tmp.erase(3);
    v_tmp[2200] = v_tmp[4];
    v_tmp.erase(4);
    f_tmp[2200] = f_tmp[4];
    f_tmp.erase(4);

    //replace!
    for (auto &iter : VV)
    {
        if (v_tmp.find(iter.first) == v_tmp.end())
            continue;
        iter.second = v_tmp[iter.first];
        FF[iter.first] = f_tmp[iter.first];
    }
}

void PhantomAnimator::WriteMergedOBJ()
{
    ConfirmSkinLayers();
    //output
    ofstream ofs(prefix + "_CT.obj");
    ofs << "#CT poture phantom for " + prefix << endl
        << endl;
    ofs << "mtllib " + prefix + ".mtl" << endl
        << endl;
    ofs.setf(ios::fixed);
    ofs.setf(ios::showpoint);
    ofs.precision(5);

    int count = 1;
    int aaa(0);
    for (auto iter : VV)
    {
        for (Vector3d v : iter.second)
            ofs << "v " << v.transpose() << endl;

        ofs << endl
            << "g " << shellName[iter.first] << endl;
        ofs << "usemtl " << shellName[iter.first] << endl;
        ofs << "s" << endl
            << endl;

        for (Vector3i f : FF[iter.first])
            ofs << "f " << f.transpose().array() + count << endl;

        count += iter.second.size();
    }
    ofs.close();
    cout << prefix << "_CT.obj was written" << endl;
}
void PhantomAnimator::ConfirmSkinLayers()
{
    MatrixXd F_normals, V_normals;
    igl::per_face_normals(Vo, Fo, F_normals);
    igl::per_vertex_normals(Vo, Fo, F_normals, V_normals);
    FF[12210] = FF[12200];
    FF[12201] = FF[12200];
    for (int i = 0; i < V.rows(); i++)
        VV[11600][i] = V.row(i);
    for (int i = 0; i < Vo.rows(); i++)
    {
        VV[12200][i] = Vo.row(i);
        VV[12201][i] = Vo.row(i) - V_normals.row(i) * 0.005;
        VV[12210][i] = Vo.row(i) - V_normals.row(i) * 0.01;
    }
}