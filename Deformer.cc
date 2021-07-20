#include <iostream>
#include <ctime>
#include <functional>

#include "PhantomAnimator.hh"
#include "igl/opengl/glfw/Viewer.h"

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

    viewer.data().set_mesh(phantom->GetV(), phantom->GetF());
    viewer.data().clear_labels();
    viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);

    viewer.data().show_overlay_depth = false;
    RotationList vQ(9, Quaterniond::Identity());
    vQ[5] = Quaterniond(AngleAxisd(20 * PI / 180., Vector3d(0, 1, 0)));  //Right shoulder
    vQ[1] = Quaterniond(AngleAxisd(-20 * PI / 180., Vector3d(0, 1, 0))); //Left shoulder
    vQ[6] = Quaterniond(AngleAxisd(130 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(-35 * PI / 180., Vector3d(0, 0, 1))) * Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0)));
    vQ[2] = Quaterniond(AngleAxisd(-130 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(35 * PI / 180., Vector3d(0, 0, 1))) * Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0)));
    vQ[7] = Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0))) * Quaterniond(AngleAxisd(60 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(-20 * PI / 180., Vector3d(0, 0, 1)));
    vQ[3] = Quaterniond(AngleAxisd(-30 * PI / 180., Vector3d(1, 0, 0))) * Quaterniond(AngleAxisd(-60 * PI / 180., Vector3d(0, 1, 0))) * Quaterniond(AngleAxisd(20 * PI / 180., Vector3d(0, 0, 1)));
    vQ[8] = Quaterniond(AngleAxisd(-45 * PI / 180., Vector3d(0, 0, 1)));
    vQ[4] = Quaterniond(AngleAxisd(45 * PI / 180., Vector3d(0, 0, 1)));

    AngleAxisd aa6 = AngleAxisd(vQ[6]);
    AngleAxisd aa2 = AngleAxisd(vQ[2]);

    int selected(0);
    if (!phantom->ReadW(string(argv[1])))
    {
        phantom->CalculateWeights(0.1);
        phantom->WriteWeights(string(argv[1]));
    }

    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods) -> bool
    {
        int idx, itNum;
        double f;
        RotationList vQ1(9, Quaterniond::Identity());
        switch (key)
        {
        case 'C':
            phantom->CalculateWeights(0.1);
        case '0': //elbow & wrist
            vQ1[3] = vQ[3];
            vQ1[4] = vQ[4];
            vQ1[7] = vQ[7];
            vQ1[8] = vQ[8];
            cout << "deforming..." << flush;
            phantom->AnimateDQS(vQ1);
            cout << "done" << endl;
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
            break;
        case '1': //shoulder
            cout << "degree: " << flush;
            cin >> f;
            vQ1[6] = Quaterniond(AngleAxisd(f * PI / 180., aa6.axis()));
            vQ1[2] = Quaterniond(AngleAxisd(f * PI / 180., aa2.axis()));
            cout << "deforming..." << flush;
            //viewer.data().set_data(phantom->GetSmoothW());
            phantom->Animate(vQ1);
            cout << "done" << endl;
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
            break;
        case '2': //scapulae
            vQ1[5] = vQ[5];
            vQ1[1] = vQ[1];
            cout << "deforming..." << flush;
            //viewer.data().set_data(phantom->GetSmoothW());
            phantom->Animate(vQ1);
            cout << "done" << endl;
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
            break;
        case '3': //shoulder front (optional, only slight move)
            cout << "degree: " << flush;
            cin >> f;
            vQ1[5] = Quaterniond(AngleAxisd(f * PI / 180., Vector3d(1., 0., 0.)));
            vQ1[1] = Quaterniond(AngleAxisd(f * PI / 180., Vector3d(1., 0., 0.)));
            cout << "deforming..." << flush;
            //viewer.data().set_data(phantom->GetSmoothW());
            phantom->AnimateDQS(vQ1);
            cout << "done" << endl;
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
            break;
        case '4': //shoulder front
            phantom->ArmOffSet(viewer.data().V_normals, 0.1);
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            break;
        case 'Q':
            cout << "iter num: " << flush;
            cin >> idx;
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::WRIST, idx);
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            break;
        case 'W':
            cout << "iter num: " << flush;
            cin >> idx;
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::ELBOW, idx);
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            break;
        case 'E':
            cout << "iter num: " << flush;
            cin >> idx;
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, idx);
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            break;
        case 'I':
            phantom->InitializeV();
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            break;
        case 'S':
            phantom->ReleaseRest();
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            break;
        case 'P':
            igl::writePLY("skin.ply", phantom->GetV(), phantom->GetF());
            break;
        case ',':
            selected = min(max(--selected, 0), 8);
            viewer.data().set_data(phantom->GetWeight(selected));
            break;
        case '.':
            selected = min(max(++selected, 0), 8);
            viewer.data().set_data(phantom->GetWeight(selected));
            break;
        case '[':
            selected = min(max(--selected, 0), 3);
            viewer.data().set_data(phantom->GetSmoothW(PhantomAnimator::JOINT(selected)));
            break;
        case ']':
            selected = min(max(++selected, 0), 3);
            viewer.data().set_data(phantom->GetSmoothW(PhantomAnimator::JOINT(selected)));
            break;
        }
        return true;
    };
    viewer.launch();

    return EXIT_SUCCESS;
}
