#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include "igl/readTGF.h"
#include "igl/writeTGF.h"
#include "igl/writePLY.h"
using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    map<int, int> idCnvt;
    ifstream ifs2("idConvert.txt");
    while(!ifs2.eof())
    {
        string aLine;
        getline(ifs2, aLine);
        if(aLine.empty()) continue;
        stringstream ss(aLine);
        int a, b;
        ss>>a>>b;
        idCnvt[a] = b;
    }
    ifs2.close();

    string fileN(argv[1]);
    cout<<"\rReading "<<fileN<<"...0/"<<idCnvt.size()+1<<flush;
    ifstream ifs(fileN);
    map<int, vector<Vector3d>> V;
    map<int, vector<Vector3i>> F;
    map<int, pair<Vector3d, Vector3d>> boundingBox;
    int startNext(1), start;
    int id(-1);
    while(!ifs.eof())
    {
        string aLine;
        getline(ifs, aLine);
        if(aLine.empty()) continue;

        stringstream ss(aLine);
        string first;
        ss>>first;
        if(first=="v")
        {
            double x, y, z;
            ss>>x>>y>>z;
            V[-1].push_back(Vector3d(x,y,z));
        }
        else if(first=="g")
        {
            start = startNext;
            startNext += V[-1].size();
            string name;
            ss>>name;
            id = atoi(name.substr(0, name.find("_")).c_str());

            //set bbox
            if(boundingBox.find(id)==boundingBox.end())
            {
                double dmax(__DBL_MAX__);
                boundingBox[id] = make_pair(Vector3d(dmax, dmax, dmax), Vector3d(-dmax, -dmax, -dmax));  // left  (x>0)
                boundingBox[-id] = make_pair(Vector3d(dmax, dmax, dmax), Vector3d(-dmax, -dmax, -dmax)); // right (x<0)
            }
            for(Vector3d v:V[-1])
            {
                int i = id;
                if(v(0)<0) i = -id;
                for(int j=0;j<3;j++)
                {
                    if(boundingBox[i].first(j)>v(j)) boundingBox[i].first(j) = v(j);
                    else if(boundingBox[i].second(j)<v(j)) boundingBox[i].second(j) = v(j);                   
                }
            }

            if(idCnvt.find(id)!=idCnvt.end()){// || id==11600){ //idCnvt or RST
                if(V.find(id)!= V.end())
                {
                    cerr<<name<<" is the second shell of "<<id<<endl;
                    exit(1);
                }                
                V[id] = V[-1];
                cout<<"\rReading "<<fileN<<"..."<<V.size()-1<<"/"<<idCnvt.size()<<flush;
            }
            else id = -1;
            V[-1].clear();
        }
        else if (first == "f")
        {
            if(id<0) continue;
            int a, b, c;
            ss>>a>>b>>c;
            F[id].push_back(Vector3i(a-start, b-start, c-start));
        }
    }    
    ifs.close();
    cout<<endl;

    //set C
    MatrixXd C;
    MatrixXi BE;
    igl::readTGF("mrcp_am.tgf", C, BE);
    //less important
    function<double(int j, double w, int dim)> intrpl = [&](int j, double w, int dim)->double
    {
        return boundingBox[j].first(dim) * (1-w) + boundingBox[j].second(dim) * w;
    };

    //MRCPAM
    C.row(0) = Vector3d(boundingBox[-4100].first(0)*0.5+boundingBox[4100].second(0)*0.5,0,intrpl(4100, 0.25, 2));
    C.row(1) = Vector3d(C(0,0),C(0,1),intrpl(5100, 0.5, 2));
    C.row(2) = Vector3d(C(0,0),C(0,1),intrpl(4900, 0.3, 2));
    C.row(7) = Vector3d(C(0,0),C(0,1),boundingBox[4900].second(2));
    C.row(22) = Vector3d(C(0,0),C(0,1),boundingBox[3900].second(2));
    C.row(12) = Vector3d(C(0,0)*0.5 + boundingBox[4100].second(0)*0.5,0,C(0,2));
    C.row(15) = Vector3d(C(0,0)*0.5 + boundingBox[-4100].first(0)*0.5,0,C(0,2));
    C.row(13) = Vector3d(intrpl(3400, 0.5, 0),intrpl(2800, 0.5, 1),boundingBox[2800].first(2)*0.8 + boundingBox[3400].second(2)*0.2);
    C.row(16) = Vector3d(intrpl(-3400, 0.5, 0),intrpl(-2800, 0.5, 1),boundingBox[-2800].first(2)*0.8 + boundingBox[-3400].second(2)*0.2);
    C.row(14) = Vector3d(C(13, 0),C(13, 1), boundingBox[3400].first(2));
    C.row(17) = Vector3d(C(16, 0),C(16, 1), boundingBox[3400].first(2));
    C.row(20) = Vector3d(intrpl(3700, 0.5, 0),intrpl(3700, 0.5, 1), intrpl(3700, 0.2, 2));
    C.row(21) = Vector3d(intrpl(-3700, 0.5, 0),intrpl(-3700, 0.5, 1), intrpl(-3700, 0.2, 2));
    //important
    C.row(3) = Vector3d(intrpl(2400, 0.3, 0),intrpl(2400, 0.5, 1),intrpl(2400, 0.5, 2));
    C.row(8) = Vector3d(intrpl(-2400, 0.7, 0),intrpl(-2400, 0.5, 1),intrpl(-2400, 0.5, 2));
    C.row(4) = Vector3d(intrpl(1400, 0.65, 0),intrpl(1400, 0.5, 1),intrpl(1400, 0.7, 2));
    C.row(9) = Vector3d(intrpl(-1400, 0.35, 0),intrpl(-1400, 0.5, 1),intrpl(-1400, 0.7, 2));
    C.row(5) = Vector3d(intrpl(1700, 0.5, 0),intrpl(1700, 0.5, 1), boundingBox[1700].first(2));
    C.row(10) = Vector3d(intrpl(-1700, 0.5, 0),intrpl(-1700, 0.5, 1),boundingBox[-1700].first(2));
    C.row(6) = Vector3d(intrpl(1900, 0.75, 0),intrpl(1900, 0.5, 1), boundingBox[1900].first(2));
    C.row(11) = Vector3d(intrpl(-1900, 0.25, 0),intrpl(-1900, 0.5, 1),boundingBox[-1900].first(2));
    C.row(18) = Vector3d(intrpl(2200, 0.8, 0),intrpl(2200, 0.5, 1), intrpl(2200, 0.5, 2));
    C.row(19) = Vector3d(intrpl(-2200, 0.2, 0),intrpl(-2200, 0.5, 1),intrpl(2200, 0.5, 2));

    //MRCPAF
    // C.row(0) = Vector3d(boundingBox[-4100].first(0)*0.5+boundingBox[4100].second(0)*0.5,0,intrpl(4100, 0.3, 2));
    // C.row(1) = Vector3d(C(0,0),C(0,1),intrpl(5100, 0.5, 2));
    // C.row(2) = Vector3d(C(0,0),C(0,1),intrpl(4900, 0.3, 2));
    // C.row(7) = Vector3d(C(0,0),C(0,1),intrpl(4700, 0.25, 2));
    // C.row(22) = Vector3d(C(0,0),C(0,1),boundingBox[3900].second(2));
    // C.row(12) = Vector3d(C(0,0)*0.5 + boundingBox[4100].second(0)*0.6,0,C(0,2));
    // C.row(15) = Vector3d(C(0,0)*0.5 + boundingBox[-4100].first(0)*0.4,0,C(0,2));
    // C.row(13) = Vector3d(intrpl(3400, 0.5, 0),intrpl(3400, 0.5, 1),boundingBox[2800].first(2));
    // C.row(16) = Vector3d(intrpl(-3400, 0.5, 0),intrpl(-3400, 0.5, 1),boundingBox[-2800].first(2));
    // C.row(14) = Vector3d(intrpl(3400, 0.3, 0),intrpl(3400, 0.75, 1), boundingBox[3400].first(2));
    // C.row(17) = Vector3d(intrpl(-3400, 0.7, 0),intrpl(-3400, 0.75, 1), boundingBox[3400].first(2));
    // C.row(20) = Vector3d(intrpl(3700, 0.5, 0),intrpl(3700, 0.5, 1), intrpl(3700, 0.2, 2));
    // C.row(21) = Vector3d(intrpl(-3700, 0.5, 0),intrpl(-3700, 0.5, 1), intrpl(-3700, 0.2, 2));
    // //important
    // C.row(3) = Vector3d(intrpl(2400, 0.3, 0),intrpl(2400, 0.7, 1),intrpl(2400, 0.6, 2));
    // C.row(8) = Vector3d(intrpl(-2400, 0.7, 0),C(3, 1),C(3, 2));
    // C.row(4) = Vector3d(intrpl(1400, 0.55, 0),intrpl(1400, 0.5, 1),intrpl(1400, 0.7, 2));
    // C.row(9) = Vector3d(intrpl(-1400, 0.45, 0),intrpl(-1400, 0.5, 1),intrpl(-1400, 0.7, 2));
    // C.row(5) = Vector3d(intrpl(1700, 0.5, 0),intrpl(1700, 0.3, 1), boundingBox[1700].first(2));
    // C.row(10) = Vector3d(intrpl(-1700, 0.5, 0),intrpl(-1700, 0.3, 1),boundingBox[-1700].first(2));
    // C.row(6) = Vector3d(intrpl(1900, 0.7, 0),intrpl(1900, 0.5, 1), boundingBox[1900].first(2));
    // C.row(11) = Vector3d(intrpl(-1900, 0.45, 0),intrpl(-1900, 0.5, 1),boundingBox[-1900].first(2));
    // C.row(18) = Vector3d(intrpl(2200, 0.8, 0),intrpl(2200, 0.5, 1), intrpl(2200, 0.5, 2));
    // C.row(19) = Vector3d(intrpl(-2200, 0.2, 0),intrpl(-2200, 0.5, 1),intrpl(2200, 0.5, 2));

    ifstream ifsS("sphere.ply");
    string str;
    int vNum, fNum;
    while(ifsS>>str)
    {
        if(str=="vertex") ifsS>>vNum;
        else if(str=="face") ifsS>>fNum;
        else if(str=="end_header") break;
    }
    vector<Vector3d> vSph;
    for(int i=0;i<vNum;i++)
    {
        double x, y, z;
        ifsS>>x>>y>>z;
        vSph.push_back(Vector3d(x,y,z));
    }
    vector<Vector3i> fSph;
    for(int i=0;i<fNum;i++)
    {
        int a, b, c, tmp;
        ifsS>>tmp>>a>>b>>c;
        fSph.push_back(Vector3i(a,b,c));
    }
    ifsS.close();

    ofstream ofs2(fileN.substr(0,fileN.size()-4)+"_joint.obj");
    ofs2<<"# joint OBJ file for "+fileN<<endl<<endl;
    start = 1;
    for(int i=0;i<C.rows();i++)
    {
        for(Vector3d v:vSph) ofs2<<"v "<<v.transpose() + C.row(i)<<endl;
        ofs2<<endl<<"g "<<i<<endl<<endl;
        for(Vector3i f:fSph) ofs2<<"f "<<(f.array() + start).transpose()<<endl;
        start += vSph.size();
    }
    ofs2.close();

    igl::writeTGF(fileN.substr(0,fileN.size()-3) + "tgf", C, BE);

    vector<int> idList;
    for(auto iter:idCnvt)
    {
        idList.push_back(iter.second);
        if(iter.second>100) continue;
        else if(iter.second>0) //right/left
        {
            idList.push_back(iter.second+4);
            vector<Vector3i> right, left;
            vector<bool> leftChk(V[iter.first].size(), false);
            for(Vector3i f:F[iter.first])
            {
                double x = V[iter.first][f(0)](0);
                if(x<0)
                {
                    left.push_back(f);
                    leftChk[f(0)] = true;
                    leftChk[f(1)] = true;
                    leftChk[f(2)] = true;
                }
                else right.push_back(f);
            }

            map<int, int> cnvt;
            int l(0), r(0);
            if(V.find(iter.second)!=V.end())
            {
                l = V[iter.second+4].size();
                r = V[iter.second].size();
            }
            for(int i=0;i<leftChk.size();i++)
            {
                if(leftChk[i])
                {
                    cnvt[i] = l++;
                    V[iter.second+4].push_back(V[iter.first][i]);
                }
                else
                {
                    cnvt[i] = r++;
                    V[iter.second].push_back(V[iter.first][i]);
                }
            }
            
            for(Vector3i f:left)
                F[iter.second+4].push_back(Vector3i(cnvt[f(0)], cnvt[f(1)], cnvt[f(2)]));
            for(Vector3i f:right)
                F[iter.second].push_back(Vector3i(cnvt[f(0)], cnvt[f(1)], cnvt[f(2)]));
            
            V[iter.first].clear();
            F[iter.first].clear();
        }
        else if(iter.second==0)
        {
            if(V.find(0)==V.end())
            {
                V[0]=V[iter.first];
                F[0]=F[iter.first];
                V[iter.first].clear();
                F[iter.first].clear();
            }
            else
            {
                int start = V[0].size();
                V[0].insert(V[0].end(), V[iter.first].begin(), V[iter.first].end());
                for(Vector3i f:F[iter.first])
                    F[0].push_back(f.array() + start);
            }
        }

    }
    sort(idList.begin(), idList.end());
    idList.erase(unique(idList.begin(), idList.end()), idList.end());

    ofstream ofs(fileN.substr(0, fileN.size()-4)+"_bone.obj");
    start = 1;
    ofs<<"# bone file for "+fileN<<endl<<endl;
    ofs.setf(ios::fixed);
    ofs.setf(ios::showpoint);
    ofs.precision(5);
    for(int id:idList)
    {
        for(Vector3d v:V[id])
            ofs<<"v "<<v(0)<<" "<<v(1)<<" "<<v(2)<<endl;
        ofs<<endl<<"g "<<id<<endl<<endl;
        for(Vector3i f:F[id])
            ofs<<"f "<<f(0)+start<<" "<<f(1)+start<<" "<<f(2)+start<<endl;
        start += V[id].size();
    }
    ofs.close();

    // MatrixXd V_rst(V[11600].size(), 3);
    // MatrixXi F_rst(F[11600].size(), 3);
    // for(int i=0;i<V_rst.rows();i++) V_rst.row(i) = V[11600][i];
    // for(int i=0;i<F_rst.rows();i++) F_rst.row(i) = F[11600][i];
    // igl::writePLY(fileN.substr(0, fileN.size()-3) + "ply", V_rst, F_rst, false);
    return true;
}
