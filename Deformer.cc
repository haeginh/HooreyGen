#include <iostream>
#include <ctime>
#include <functional>

#include "PhantomAnimator.hh"
#include "igl/opengl/glfw/Viewer.h"

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
    phantom->InitializeW();
    phantom->CalculateCleanWeights();

    Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
    igl::opengl::glfw::Viewer viewer;

    viewer.data().set_mesh(phantom->GetVI(), phantom->GetFI());

    viewer.append_mesh();
    viewer.data().set_mesh(phantom->GetV(), phantom->GetF());
    viewer.data().clear_labels();
    viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
    MatrixXd C = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetF().rows(), 1);
    viewer.data().set_colors(C);

    viewer.data().show_overlay_depth = false;
    RotationList vQ(9, Quaterniond::Identity());
    vQ[5] = Quaterniond(AngleAxisd(20 * PI / 180., Vector3d(0, 1, 0)));  //Right shoulder
    vQ[1] = Quaterniond(AngleAxisd(-20 * PI / 180., Vector3d(0, 1, 0))); //Left shoulder
    vQ[6] = Quaterniond(AngleAxisd(130 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(-35 * PI / 180., Vector3d(0, 0, 1))) * Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0)));
    vQ[2] = Quaterniond(AngleAxisd(-130 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(35 * PI / 180., Vector3d(0, 0, 1))) * Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0)));
    vQ[7] = Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0))) * Quaterniond(AngleAxisd(60 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(-20 * PI / 180., Vector3d(0, 0, 1)));
    vQ[3] = Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0))) * Quaterniond(AngleAxisd(-60 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(20 * PI / 180., Vector3d(0, 0, 1)));
    vQ[8] = Quaterniond(AngleAxisd(-35 * PI / 180., Vector3d(0, 0, 1)));
    vQ[4] = Quaterniond(AngleAxisd(35 * PI / 180., Vector3d(0, 0, 1)));

    AngleAxisd aa7 = AngleAxisd(vQ[7]);
    AngleAxisd aa3 = AngleAxisd(vQ[3]);

    AngleAxisd aa6 = AngleAxisd(vQ[6]);
    AngleAxisd aa2 = AngleAxisd(vQ[2]);

    // ofstream ofs("hoorey_AF.txt");
    // for (int i = 0; i < vQ.size(); i++)
    //     ofs << vQ[i].w() <<" "<< vQ[i].x() <<" "<< vQ[i].y() <<" "<< vQ[i].z() << endl;
    // ofs.close();

    int selected(0);
    // if (!phantom->ReadW(string(argv[1])))
    // {
    //     phantom->CalculateWeights(0.1);
    //     phantom->WriteWeights(string(argv[1]));
    // }
    int a(0), b(1);
    vector<int> shells = phantom->GetShellVec();
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods) -> bool
    {
        int idx, itNum, intTmp(0);
        double f;
        RotationList vQ1(9, Quaterniond::Identity());
        MatrixXd C1;
        switch (key)
        {
        case 'D':
            // intTmp = phantom->DetectIntersections(C1, a);
            // intTmp += phantom->DetectIntersections(C1, b);
            C1 = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetF().rows(), 1);
            intTmp = phantom->DetectIntersections(C1);
            if (intTmp)
            {
                cout << intTmp << " intersection was detected!" << endl;
            }
            else
                cout << "no intersection was detected!" << endl;
            viewer.data().set_colors(C1);
            break;
        case 'C':
        {
            MatrixXd VV;
            MatrixXi FF, TT;
            igl::copyleft::tetgen::tetrahedralize(phantom->GetV(), phantom->GetF(), "p/0.001YqO/3", VV, TT, FF);
        }
            C1 = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetF().rows(), 1);
            intTmp = phantom->SelfIntersection(C1);
            if (intTmp)
            {
                cout << intTmp << " intersection was detected!" << endl;
            }
            else
                cout << "no intersection was detected!" << endl;
            viewer.data().set_colors(C1);
            // phantom->CalculateElbowW(3, aa3.axis());//, 9, 100);
            // phantom->CalculateElbowW(7, aa7.axis());//, 9, 100);
            // phantom->CalculateShoulderW(1, -Vector3d(0, 1, 0));
            // phantom->CalculateShoulderW(5, Vector3d(0, 1, 0));
            // phantom->CalculateCleanWeights();
            break;
        case 'A':
            phantom->AdjustVolume(1, 0, 1, viewer.data().V_normals);
            phantom->AdjustVolume(5, 0, 1, viewer.data().V_normals);
            phantom->AdjustVolume(2, 1, 0.5, viewer.data().V_normals);
            phantom->AdjustVolume(6, 1, 0.5, viewer.data().V_normals);
            phantom->AdjustVolume(3, 1, 0.5, viewer.data().V_normals);
            phantom->AdjustVolume(7, 1, 0.5, viewer.data().V_normals);
            //phantom->LaplacianSmooth(1, 10, 0, 1);
            viewer.data().set_vertices(phantom->GetV());
            break;
        case 'X':
            viewer.data().is_visible = !viewer.data().is_visible;
            break;
        case 'V':
            for (int i : shells)
            {
                double vol = phantom->GetVolume(i);
                cout << i << " : " << vol << " / " << phantom->GetVolume(-i) << " (" << vol / phantom->GetVolume(-i) * 100.f - 100.f << "%)" << endl;
            }
            break;
        case '1':                                    //wrist & elbow
            phantom->CalculateElbowW(3, aa3.axis()); //, 9, 100);
            phantom->CalculateElbowW(7, aa7.axis()); //, 9, 100);
            phantom->CalculateShoulderW(1, -Vector3d(0, 1, 0));
            phantom->CalculateShoulderW(5, Vector3d(0, 1, 0));
            phantom->CalculateCleanWeights();
            phantom->SetW1();
            vQ1 = vQ;
            vQ1[6] = Quaterniond(AngleAxisd(80 * PI / 180., aa6.axis()));
            vQ1[2] = Quaterniond(AngleAxisd(80 * PI / 180., aa2.axis()));
            cout << "deforming..." << flush;
            phantom->Animate(vQ1);
            cout << "done" << endl;
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
            break;
        case '2': //shoulder
        {
            double theta6 = phantom->CalculateShoulderW(6, aa6.axis()) - 5.;
            double theta2 = phantom->CalculateShoulderW(2, aa2.axis()) - 5.;
            if (theta6 > 53.)
                theta6 = 53.;
            if (theta2 > 53.)
                theta2 = 53.;
            vQ1[6] = Quaterniond(AngleAxisd(theta6 * PI / 180., aa6.axis()));
            vQ1[2] = Quaterniond(AngleAxisd(theta2 * PI / 180., aa2.axis()));
            phantom->Animate(vQ1);
            phantom->CalculateCleanWeights();
            //viewer.data().set_data(phantom->GetSmoothW());

            vQ1[6] = Quaterniond(AngleAxisd((53. - theta6) * PI / 180., aa6.axis()));
            vQ1[2] = Quaterniond(AngleAxisd((53. - theta2) * PI / 180., aa2.axis()));
            phantom->Animate(vQ1);
            cout << "Deforamtion is done (" << theta6 << "/" << theta2 << ")" << endl;
        }
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
            break;

        case 'I':
            phantom->InitializeV();
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            break;

        case 'P':
            phantom->WriteOBJ(string(argv[1]) + "_result.obj");
            break;
        case '[':
            selected = std::max(0, std::min(--selected, (int)phantom->GetBE().rows() - 1));
            viewer.data().set_data(phantom->GetWeight(selected));
            break;
        case ']':
            selected = std::max(0, std::min(++selected, (int)phantom->GetBE().rows() - 1));
            viewer.data().set_data(phantom->GetWeight(selected));
            break;
        }
        return true;
    };
    viewer.launch();

    return EXIT_SUCCESS;
}
