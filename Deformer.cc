#include <iostream>
#include <ctime>
#include <functional>

#include "PhantomAnimator.hh"
#include "igl/opengl/glfw/Viewer.h"
#include "igl/writePLY.h"
#define PI 3.14159265358979323846

void PrintUsage()
{
    cout << "<Usage>" << endl;
    cout << "./DCIR_RDC [phantom prefix (ply, tgf)]" << endl;
    exit(1);
}

using namespace Eigen;
using namespace std;
typedef Triplet<double> T;

int main(int argc, char **argv)
{
    if (argc != 2)
        PrintUsage();

    //phantom animator
    PhantomAnimator *phantom = new PhantomAnimator(string(argv[1]));

    Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
    igl::opengl::glfw::Viewer viewer;
    int selected(0);
    vector<int> appended = {125};
    phantom->SetSimpleF(appended[0]);
    phantom->PreComputeSmooth();
    viewer.data().set_mesh(phantom->GetV_simple(), phantom->GetF_simple());
    viewer.data().clear_labels();
    viewer.data().set_data(phantom->bw);
    //  viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
    RotationList vQ(22, Quaterniond::Identity());
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods) -> bool
    {
        int idx, itNum(-1);
        double degree, weight(-1.);
        MatrixXd CT;
        MatrixXi BET;
        switch (key)
        {
        case '.':
            selected++;
            selected = std::min(std::max(selected, 0), (int)phantom->GetBE().rows() - 1);
            for (int i = 0; i < appended.size(); i++)
                viewer.data_list[i].set_data(phantom->GetWeight(appended[i], selected));
            break;
        case ',':
            selected--;
            selected = std::min(std::max(selected, 0), (int)phantom->GetBE().rows() - 1);
            for (int i = 0; i < appended.size(); i++)
            {
                viewer.data_list[i].set_data(phantom->GetWeight(appended[i], selected));
            }
            break;
        case 'S':
            cout << "set id: " << flush;
            cin >> idx;
            phantom->SetSimpleF(idx);
            for (auto &data : viewer.data_list)
                data.clear();
            viewer.data_list.resize(1);
            viewer.selected_data_index = 0;
            appended = {idx};
            viewer.data().set_mesh(phantom->GetV_simple(), phantom->GetF_simple());
            //  viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
            viewer.data().set_data(phantom->GetWeight(idx, selected));
            viewer.data().set_visible(true);
            viewer.data().clear_labels();
            break;
        case 'A':
            cout << "append id: " << flush;
            cin >> idx;
            phantom->SetSimpleF(idx);
            viewer.append_mesh();
            viewer.data().set_mesh(phantom->GetV_simple(), phantom->GetF_simple());
            viewer.data().set_data(phantom->GetWeight(idx, selected));
            viewer.data().clear_labels();
            appended.push_back(idx);
            break;
        case 'R':
            cout << "arm degree: " << flush;
            cin >> degree;
            vQ[11] = Quaterniond(AngleAxisd(degree * PI / 180., Vector3d(0, 1, 0))); //Right shoulder
            vQ[4] = Quaterniond(AngleAxisd(-degree * PI / 180., Vector3d(0, 1, 0))); //Left shoulder
            vQ[10] = Quaterniond(AngleAxisd(10 * PI / 180., Vector3d(0, 1, 0)));     //Right shoulder
            vQ[3] = Quaterniond(AngleAxisd(-10 * PI / 180., Vector3d(0, 1, 0)));     //Left shoulder

            vQ[11] = vQ[11] * Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(0, 0, 1))) * Quaterniond(AngleAxisd(-20 * PI / 180., Vector3d(1, 0, 0)));
            vQ[4] = vQ[4] * Quaterniond(AngleAxisd(30 * PI / 180., Vector3d(0, 0, 1))) * Quaterniond(AngleAxisd(-20 * PI / 180., Vector3d(1, 0, 0)));
            vQ[12] = Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0))) * Quaterniond(AngleAxisd(70 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(0, 0, 1)));
            vQ[5] = Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0))) * Quaterniond(AngleAxisd(-70 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(30 * PI / 180., Vector3d(0, 0, 1)));
            vQ[13] = Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(0, 0, 1)));
            vQ[6] = Quaterniond(AngleAxisd(30 * PI / 180., Vector3d(0, 0, 1)));
            cout << "deforming..." << flush;
            phantom->Animate(vQ, CT, BET);
            cout << "done" << endl;

            phantom->SetSimpleF(125, true);
            viewer.data().set_vertices(phantom->GetV_simple());
            // for (int i = 0; i < appended.size(); i++)
            // {
            //     phantom->SetSimpleF(appended[i], true);
            //     viewer.data_list[i].set_vertices(phantom->GetV_simple());
            // }
            //viewer.data().show_overlay_depth = true;
            //  viewer.data_list.back().set_edges(CT, BET, sea_green);
            break;
        case 'W':
            cout << "smooth mesh (1) iterating num: " << flush; cin>>itNum;
            cout << "smooth mesh (1) weight: " << flush; cin>>weight;
            cout << "performing smoothing process.."<<flush;
            if(itNum<0) itNum = 10;
            if(weight<0) weight = 0.01;
            phantom->LaplacianSmooth(itNum, weight);
            cout << "done" << endl;
            viewer.data().set_vertices(phantom->GetV_simple());
            //viewer.data().show_overlay_depth = true;
            //  viewer.data_list.back().set_edges(CT, BET, sea_green);
            break;
        case 'P':
            phantom->SetSimpleF(126);
            igl::writePLY("skin.ply", phantom->GetV_simple(), phantom->GetF_simple());
            break;
        }
        return true;
    };
    viewer.launch();

    return EXIT_SUCCESS;
}
