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
    phantom->InitArms();

    Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
    igl::opengl::glfw::Viewer viewer;

    viewer.data().set_mesh(phantom->GetV(), phantom->GetF());
    viewer.data().clear_labels();
    viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);

    viewer.data().show_overlay_depth = false;
    viewer.data().double_sided = true;

    int selected(0);
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods) -> bool
    {
        double thres;
        switch (key)
        {
        case ' ': //0.5 is appropriate
            cout<<"threshold: "<<flush; cin>>thres;
            cout<<"cutting..."<<flush;
            phantom->CutArms(thres);
            cout<<"done"<<endl;
            break;
        case ',':
            selected = min(max(--selected, 0), 21);
            viewer.data().set_data(phantom->GetWeight(selected));
            break;
        case '.':
            selected = min(max(++selected, 0), 21);
            viewer.data().set_data(phantom->GetWeight(selected));
            break;
        }
        return true;
    };
    viewer.launch();

    return EXIT_SUCCESS;
}
