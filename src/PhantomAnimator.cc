#include "PhantomAnimator.hh"
#include "igl/opengl/glfw/Viewer.h"
PhantomAnimator::PhantomAnimator() {}

PhantomAnimator::PhantomAnimator(string prefix)
    : name(prefix)
{
    ReadFiles(prefix);
    //CalculateWeights(0.1);
}

bool PhantomAnimator::ReadFiles(string prefix)
{
    cout << "Read " + prefix + ".tgf" << endl;
    igl::readTGF(prefix + ".tgf", C, BE);
    cout << "Read " + prefix + ".mesh" << endl;
    if (!igl::readMESH(prefix + ".mesh", V, T, F) | !igl::readDMAT(prefix + ".W", W))
    {
        cout << "Read " + prefix + ".ply" << endl;
        igl::readPLY(prefix + ".ply", V, F);
        CalculateWeights();
        igl::writeDMAT(prefix + ".W", W, false);
        igl::writeMESH(prefix + ".mesh", V, T, F);
    }

    //SetSimpleF();

    return true;
}

bool PhantomAnimator::CalculateWeights()
{
    //perform BBW
    MatrixXd boneP = GenerateBonePoints(C, BE, 1.);
    MatrixXd V1(V.rows() + boneP.rows(), 3);
    V1 << V, boneP;
    cout << "<Tetrahedralization>" << endl;
    MatrixXd VT;
    MatrixXi FT;
    igl::copyleft::tetgen::tetrahedralize(V1, F, "pYq", VT, T, FT);

    cout << "<Calculate Bone Weights>" << endl;
    MatrixXd bc;
    VectorXi b;
    igl::boundary_conditions(VT, T, C, VectorXi(), BE, MatrixXi(), b, bc);
    cout << bc.rows() << " X " << bc.cols() << endl;
    igl::BBWData bbw_data;
    bbw_data.active_set_params.max_iter = 10;
    bbw_data.verbosity = 2;
    if (!igl::bbw(VT, T, b, bc, bbw_data, W))
        return EXIT_FAILURE;

    igl::normalize_row_sums(W, W);

    int selected = 0;
    Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(VT, FT);
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
    V = VT;
    F = FT;
    return true;
}

bool PhantomAnimator::InitArms()
{
    cout << "read tet files..." << flush;
    V0 = ReadNode(name + ".node");
    if (V0.cols() != 3)
        return false;
    T0 = ReadEle(name + ".ele");
    if (T0.cols() != 5)
        return false;
    cout << "done" << endl;

    if (!igl::readDMAT(name + ".armW", armW))
    {
        cout << "calculate centers..." << flush;
        MatrixXd centers(T0.rows(), 3);
        for (int i = 0; i < T0.rows(); i++)
            centers.row(i) = (V0.row(T0(i, 0)) + V0.row(T0(i, 1)) + V0.row(T0(i, 2)) + V0.row(T0(i, 3))) * 0.25;
        cout << "done" << endl;
        map<int, map<int, double>> bary = GenerateBarycentricCoord(V, T, centers);
        cout << "setting arm weights..." << flush;
        VectorXd armW0 = W.col(2) + W.col(3) + W.col(4) + W.col(6) + W.col(7) + W.col(8);
        armW = VectorXd::Zero(T0.rows());
        for (int i; i < T0.rows(); i++)
        {
            if (bary.find(i) == bary.end()) //eye vertices
            {
                armW(i) = 0;
                continue;
            }
            for (auto iter : bary[i])
                armW(i) += armW0(iter.first) * iter.second;
        }
        igl::writeDMAT(name + ".armW", armW, false);
    }
    else cout<<"read armW file"<<endl;
    return true;
}

void PhantomAnimator::CutArms(double threshold)
{
    vector<int> extV, extT;
    for(int i=0;i<armW.rows();i++)
    {
        if(armW(i)>threshold) continue;
        extT.push_back(i);
        extV.push_back(T0(i,0));
        extV.push_back(T0(i,1));
        extV.push_back(T0(i,2));
        extV.push_back(T0(i,3));
    }
    sort(extV.begin(), extV.end());
    extV.erase(unique(extV.begin(), extV.end()), extV.end());
    
    map<int, int> w2e;
    string name1 = name + "_" + to_string(threshold);
    ofstream ofsNode(name1 + ".node");
    ofsNode<<extV.size()<<"  3  0  0"<<endl;
    ofsNode.setf(ios::fixed); ofsNode.setf(ios::showpoint); ofsNode.precision(10);
    for(int i = 0; i<extV.size(); i++) 
    {
        w2e[extV[i]] = i;
        ofsNode<<i<<" "<<V0(extV[i], 0)<<" "<<V0(extV[i], 1)<<" "<<V0(extV[i], 2)<<endl;
    }
    ofsNode.close();
  
    ofstream ofsEle(name1 + ".ele");
    ofsEle<<extT.size()<<"  4  1"<<endl;
    int n(0);
    for(int i:extT)
    {
        ofsEle<<n++<<" "<<w2e[T0(i, 0)]<<" "<<w2e[T0(i, 1)]<<" "<<w2e[T0(i, 2)]<<" "<<w2e[T0(i, 3)]<<" "<<T0(i, 4)<<endl;
    }
    ofsEle.close();
}

MatrixXd PhantomAnimator::ReadNode(string fileName)
{
    ifstream ifs(fileName);
    if (!ifs.is_open())
        return MatrixXd::Zero(1, 1);
    int num, tmp;
    ifs >> num >> tmp >> tmp >> tmp;
    MatrixXd V0(num, 3);
    for (int i = 0; i < num; i++)
        ifs >> tmp >> V0(i, 0) >> V0(i, 1) >> V0(i, 2);
    return V0;
}

MatrixXi PhantomAnimator::ReadEle(string fileName)
{
    ifstream ifs(fileName);
    if (!ifs.is_open())
        return MatrixXi::Zero(1, 1);
    int num, tmp;
    ifs >> num >> tmp >> tmp;
    MatrixXi T0(num, 5);
    for (int i = 0; i < num; i++)
        ifs >> tmp >> T0(i, 0) >> T0(i, 1) >> T0(i, 2) >> T0(i, 3) >> T0(i, 4);
    return T0;
}