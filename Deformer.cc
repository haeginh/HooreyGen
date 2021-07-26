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
    ifstream ifs("hoorey.txt");
    for(int i=0;i<vQ.size();i++)
        ifs>>vQ[i].w()>>vQ[i].x()>>vQ[i].y()>>vQ[i].z();
    ifs.close();

    int selected(0);
    if (!phantom->ReadW(string(argv[1])))
    {
        phantom->CalculateWeights(0.1);
        phantom->WriteWeights(string(argv[1]));
    }

    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods) -> bool
    {
        double vol;
        RotationList vQ1(9, Quaterniond::Identity());
        bool update(false);
        switch (key)
        {
        case '1': //elbow & wrist
            cout << "Deformation + Smoothing (1)..." << flush;
            vQ1 = vQ;
            vQ1[6] = Quaterniond(AngleAxisd(80. * PI / 180., AngleAxisd(vQ[6]).axis()));
            vQ1[2] = Quaterniond(AngleAxisd(80. * PI / 180., AngleAxisd(vQ[2]).axis()));
            phantom->AnimateDQS(vQ1);
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::WRIST, 5);
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::ELBOW, 15);
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 100);
            cout << "done" << endl;
            update = true;
            break;
        case '2': //elbow & wrist
            cout << "Deformation + Smoothing + Recovering (2)..." << flush;
            vQ1[6] = Quaterniond(AngleAxisd(50. * PI / 180., AngleAxisd(vQ[6]).axis()));
            vQ1[2] = Quaterniond(AngleAxisd(50. * PI / 180., AngleAxisd(vQ[2]).axis()));
            phantom->AnimateDQS(vQ1);
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 100);
            phantom->ReleaseRest();
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 100);
            phantom->ReleaseRest();
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 200);
            cout << "done" << endl;
            update = true;
            break;
        case '3':
            phantom->CalculateWeights(0.1);
            break;
        case '4': //shoulder up
            cout << "deforming..." << flush;
            vQ1[5] = Quaterniond(AngleAxisd(5 * PI / 180., Vector3d(0, 1, 0)));
            vQ1[1] = Quaterniond(AngleAxisd(-5 * PI / 180., Vector3d(0, 1, 0)));
            phantom->Animate(vQ1, true);
            phantom->Animate(vQ1, true);                                     // up 10
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            phantom->ShoulderOffset(viewer.data().V_normals, 0.3);           // shoulder offset 0.3
            phantom->ReleaseRest();                                          // relax
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 100); // smoothing 100
            phantom->Animate(vQ1, true);                                     // up 5
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            phantom->ShoulderOffset(viewer.data().V_normals, 0.3);           // shoulder offset 0.3
            phantom->ReleaseRest();                                          // relax
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 100); // smoothing 100
            phantom->Animate(vQ1, true);                                     // up 5
            phantom->ReleaseRest();                                          // relax
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 150); // smoothing 150
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            phantom->ShoulderOffset(viewer.data().V_normals, 0.2);            // shoulder offset 0.2
            phantom->ReleaseRest();                                           // relax
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 50);   // smoothing 50
            phantom->ArmOffset(viewer.data().V_normals, 0.3);                 // arm offset 0.3
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 50);   // smoothing 50
            cout << "done" << endl;
            update = true;
            break;
        case 'S': //ignore the message "^tetrahedralize: Tetgen failed to create tets"
            {
                MatrixXd VT; MatrixXi TT, FT; // dummies
                igl::copyleft::tetgen::tetrahedralize(phantom->GetV(), phantom->GetF(), "d",VT, TT, FT);
            }
            break;
        case 'V':
            vol = phantom->GetVolume();
            cout<<"current volume : "<<vol<<" / "<<phantom->GetVolume(true)<<" ("<<vol/phantom->GetVolume(true) * 100.f<<"%)"<<endl;
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
        if(update)
        {
            viewer.data().set_vertices(phantom->GetV());
            viewer.data().compute_normals();
            viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
        }
        return true;
    };
    viewer.launch();

    return EXIT_SUCCESS;
}
