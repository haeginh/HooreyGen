#include "PhantomAnimator.hh"
#include "igl/opengl/glfw/Viewer.h"
PhantomAnimator::PhantomAnimator() {}

PhantomAnimator::PhantomAnimator(string prefix)
{
    ReadFiles(prefix);
}

bool PhantomAnimator::ReadFiles(string prefix)
{
    cout << "Read " + prefix + ".tgf" << endl;
    igl::readTGF(prefix + ".tgf", C, BE);
    igl::directed_edge_parents(BE, P);
    cout << "Read " + prefix + ".ply" << endl;
    MatrixXd V_ply;
    MatrixXi F_ply;
    igl::readPLY(prefix + ".ply", V_ply, F_ply);

    string objName = prefix + "_simple.obj";
    cout << "Read " + objName << endl;
    ReadOBJ(objName);
    V0=V;

    cout << "Read " + objName + "_W" << endl;
    if (!ReadOBJW(objName + "_W"))
    {
        cout << "-> not exists" << endl;
        //perform BBW
        MatrixXd VT, WT;
        MatrixXi TT, FT;
        cout << "Read " + prefix + ".W and " + prefix + ".mesh" << endl;
        if (!igl::readDMAT(prefix + ".W", WT) || (!igl::readMESH(prefix + ".mesh", VT, TT, FT)))
        {
            cout << "-> not exists" << endl;
            MatrixXd W;
            W.resize(V_ply.rows(), BE.rows());
            MatrixXd boneP = GenerateBonePoints(C, BE, 1.);
            MatrixXd V1(V_ply.rows() + boneP.rows(), 3);
            V1 << V_ply, boneP;
            cout << "<Tetrahedralization>" << endl;
            igl::copyleft::tetgen::tetrahedralize(V1, F_ply, "pYq", VT, TT, FT);

            cout << "<Calculate Bone Weights>" << endl;
            MatrixXd bc;
            VectorXi b;
            igl::boundary_conditions(VT, TT, C, VectorXi(), BE, MatrixXi(), b, bc);
            cout << bc.rows() << " " << bc.cols() << endl;
            igl::BBWData bbw_data;
            bbw_data.active_set_params.max_iter = 10;
            bbw_data.verbosity = 2;
            if (!igl::bbw(VT, TT, b, bc, bbw_data, WT))
                return EXIT_FAILURE;
            igl::normalize_row_sums(WT, WT);

            //matching between tetra & ply
            cout << "matching between tetra & ply.." << flush;
            auto grid = GenerateGrid(VT);
            int count(0);
            for (int i = 0; i < V_ply.rows(); i++)
            {
                int x = floor(V_ply(i, 0) + 0.5);
                int y = floor(V_ply(i, 1) + 0.5);
                int z = floor(V_ply(i, 2) + 0.5);
                auto key = make_tuple(x, y, z);
                for (int n : grid[key])
                {
                    if (fabs(V_ply(n, 0) - VT(i, 0)) > 0.01)
                        continue;
                    if (fabs(V_ply(n, 1) - VT(i, 1)) > 0.01)
                        continue;
                    if (fabs(V_ply(n, 2) - VT(i, 2)) > 0.01)
                        continue;
                    W.row(i) = WT.row(n);
                    cout << "\rmatching between tetra & ply.." << ++count << "/" << V_ply.rows() << flush;
                    break;
                }
            }
            cout << endl;

            int selected = 0;
            Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
            igl::opengl::glfw::Viewer viewer;
            viewer.data().set_mesh(V_ply, F_ply);
            viewer.data().set_data(W.col(selected));
            viewer.data().set_edges(C, BE, sea_green);
            viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods) -> bool
            {
                switch (key)
                {
                case '.':
                    selected++;
                    selected = std::min(std::max(selected, 0), (int)W.cols() - 1);
                    viewer.data().set_data(W.col(selected));
                    break;
                case ',':
                    selected--;
                    selected = std::min(std::max(selected, 0), (int)W.cols() - 1);
                    viewer.data().set_data(W.col(selected));
                    break;
                }
                return true;
            };
            viewer.launch();
            igl::writeMESH(prefix + ".mesh", VT, TT, FT);
            igl::writeDMAT(prefix + ".W", WT, false);
        }

        map<int, map<int, double>> baryCoords = GenerateBarycentricCoord(VT, TT, V);
        cleanWeights.clear();
        double epsilon(1e-5);
        for (size_t i = 0; i < baryCoords.size(); i++)
        {
            if (baryCoords.find(i) == baryCoords.end())
            {
                cerr << "missing coordinate for vertex #" << i << endl;
                exit(0);
            }
            VectorXd w = VectorXd::Zero(BE.rows());
            for (auto iter : baryCoords[i])
                w += iter.second*WT.row(iter.first);
            
            map<int, double> weight;
            double sum(0);
            for (int j=0;j<w.size();j++)
            {
                if(w(j)<epsilon) continue;
                weight[j] = w(j);
                sum += w(j);
            }
            for (auto &w : weight)
                w.second /= sum;
            cleanWeights.push_back(weight);
        }
        ofstream ofs(objName + "_W");
        for (auto iter : cleanWeights)
        {
            ofs << iter.size() << " ";
            for (auto w : iter)
                ofs << w.first << " " << w.second << " ";
            ofs << endl;
        }
        ofs.close();
    }
    //SetSimpleF();

    return true;
}

bool PhantomAnimator::ReadOBJ(string fileName)
{
    ifstream ifs(fileName);
    if (!ifs.is_open())
        return false;
    string aLine;
    vector<Vector3d> vertices;
    vector<Vector3i> faces;
    int idx;
    while (getline(ifs, aLine))
    {
        stringstream ss(aLine);
        string buf;
        ss >> buf;
        if (buf == "v")
        {
            double x, y, z;
            ss >> x >> y >> z;
            vertices.push_back(Vector3d(x, y, z));
        }
        else if (buf == "g")
        {
            string name;
            ss >> name;
            idx = atoi(name.substr(2, name.find('_', 2) - 2).c_str());
            if (shells.find(idx) != shells.end())
            {
                cerr << fileName << " has multiple id number: " << idx << endl;
                exit(0);
            }
            shells[idx] = make_pair(faces.size(), 0);
        }
        else if (buf == "f")
        {
            int a, b, c;
            ss >> a >> b >> c;
            faces.push_back(Vector3i(a - 1, b - 1, c - 1));
            shells[idx].second++;
        }
    }
    ifs.close();

    V.resize(vertices.size(), 3);
    F.resize(faces.size(), 3);
    for (size_t i = 0; i < vertices.size(); i++)
        V.row(i) = vertices[i];
    for (size_t i = 0; i < faces.size(); i++)
        F.row(i) = faces[i];
    return true;
}

bool PhantomAnimator::ReadOBJW(string fileName)
{
    ifstream ifs(fileName);
    if (!ifs.is_open())
        return false;
    string aLine;
    cleanWeights.clear();
    while (getline(ifs, aLine))
    {
        if (aLine.length() == 0)
            continue;
        stringstream ss(aLine);
        int num;
        ss >> num;
        map<int, double> weight;
        for (int i = 0; i < num; i++)
        {
            int idx;
            double w;
            ss >> idx >> w;
            weight[idx] = w;
        }
        cleanWeights.push_back(weight);
    }
    ifs.close();
    return true;
}

bool PhantomAnimator::SetSimpleF(int idx, bool onlyV)
{
    if (shells.find(idx) == shells.end())
        return false;
    F_simple.resize(shells[idx].second, 3);
    vector<int> extractedV;
    if (extractedVertices.find(idx) != extractedVertices.end())
        extractedV = extractedVertices[idx];
    else
    {
        for (int j = 0; j < shells[idx].second; j++)
        {
            F_simple.row(j) = F.row(shells[idx].first + j);
            extractedV.push_back(F_simple(j, 0));
            extractedV.push_back(F_simple(j, 1));
            extractedV.push_back(F_simple(j, 2));
        }

        sort(extractedV.begin(), extractedV.end());
        extractedV.erase(unique(extractedV.begin(), extractedV.end()), extractedV.end());
        extractedVertices[idx] = extractedV;
    }

    V_simple.resize(extractedV.size(), 3);
    map<int, int> w2e;
    int num = 0;
    for (int v : extractedV)
    {
        V_simple.row(num) = V.row(v);
        w2e[v] = num++;
    }
    
    if(onlyV) return true;
    for (int r = 0; r < F_simple.rows(); r++)
    {
        F_simple(r, 0) = w2e[F_simple(r, 0)];
        F_simple(r, 1) = w2e[F_simple(r, 1)];
        F_simple(r, 2) = w2e[F_simple(r, 2)];
    }
    return true;
}

VectorXd PhantomAnimator::GetWeight(int idx, int col)
{
    if (extractedVertices.find(idx) == extractedVertices.end())
    {
        cout << "there is no data for " << idx << endl;
        return VectorXd::Zero(1);
    }
    VectorXd wVec = VectorXd::Zero(extractedVertices[idx].size());
    for (int i = 0; i < extractedVertices[idx].size(); i++)
    {
        int v = extractedVertices[idx][i];
        auto iter = cleanWeights[v].find(col);
        if (iter == cleanWeights[v].end())
            continue;
        wVec(i) = iter->second;
    }
    return wVec;
}

#include "igl/forward_kinematics.h"
#include "igl/lbs_matrix.h"
void PhantomAnimator::Animate(RotationList vQ, MatrixXd &CT, MatrixXi &BET)
{
    if(M.rows()==0){
        MatrixXd W = MatrixXd::Zero(V.rows(), BE.rows());
        for(int i=0;i<cleanWeights.size();i++) 
            for(auto iter:cleanWeights[i]) W(i, iter.first) = iter.second;
        igl::lbs_matrix(V0, W, M);
    }
     MatrixXd T;
     igl::forward_kinematics(C, BE, P, vQ, T);
     igl::deform_skeleton(C, BE, T, CT, BET);
     V = M * T;
    //     vector<Vector3d> vT;
    // igl::forward_kinematics(C,BE, P, vQ, vQ, vT);
    // myDqs(V0, cleanWeights, vQ, vT, V);
}
