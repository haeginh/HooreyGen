#include <iostream>
#include <ctime>
#include <functional>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>

#include "PhantomAnimator.hh"

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

    //pre-generated bone file
    string name = string(argv[argc - 1]);

    //phantom animator
    PhantomAnimator *phantom = new PhantomAnimator(name);

    Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
    igl::opengl::glfw::Viewer viewer;

    viewer.data().set_mesh(phantom->GetV(), phantom->GetF());
    viewer.data().clear_labels();
    viewer.data().set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
    viewer.data().show_overlay_depth = false;
    viewer.append_mesh();
    viewer.data().set_mesh(phantom->GetVbone(), phantom->GetFbone());
    viewer.data().set_colors(RowVector3d(1., 1., 1.));
    viewer.append_mesh();
    viewer.data().set_mesh(phantom->GetVo(), phantom->GetFo());
    viewer.data().is_visible = false;
    Eigen::MatrixXd CM;
    igl::parula(Eigen::VectorXd::LinSpaced(100,0,1).eval(),false,CM);
    viewer.data(0).set_colormap(CM);
    viewer.selected_data_index = 0;
    viewer.core().background_color = Vector4f(1., 1., 1., 1.);

    RotationList vQ(9, Quaterniond::Identity());
    ifstream ifs("hoorey.txt");
    for (int i = 0; i < vQ.size(); i++)
        ifs >> vQ[i].w() >> vQ[i].x() >> vQ[i].y() >> vQ[i].z();
    ifs.close();

    if (!phantom->ReadW())
    {
        phantom->CalculateWeights(0.1);
    }

    enum DEFORM
    {
        RELAX,
        SMOOTH1,
        SMOOTH10,
        SMOOTH50,
        SMOOTH100,
        ARMPIT,
        INNERARM,
        INNERARM_SMALL,
        INNERARM_VERYSMALL,
        LOWERARM,
        LOWERARM_OUTVERYSMALL,
        LOWERARM_INVERYSMALL,
        COVER,
        SHOULDERUP
    };
    function<void(DEFORM)> Deformation = [&](DEFORM d)
    {
        MatrixXi IF;
        MatrixXd C0, C1;
        RotationList vQ1(9, Quaterniond::Identity());
        int iteration(0);
        switch (d)
        {
        case RELAX:
            phantom->ReleaseRest(0.5);
            viewer.data(0).set_data(phantom->GetWeight(0));
            cout << "relax" << endl;
            break;
        case SMOOTH1:
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 1);
            viewer.data(0).set_data(phantom->GetSmoothW(PhantomAnimator::JOINT::SHOULDER));
            cout << "smooth 1" << endl;
            break;
        case SMOOTH10:
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 10);
            viewer.data(0).set_data(phantom->GetSmoothW(PhantomAnimator::JOINT::SHOULDER));
            cout << "smooth 10" << endl;
            break;
        case SMOOTH50:
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 50);
            viewer.data(0).set_data(phantom->GetSmoothW(PhantomAnimator::JOINT::SHOULDER));
            cout << "smooth 50" << endl;
            break;
        case SMOOTH100:
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 100);
            viewer.data(0).set_data(phantom->GetSmoothW(PhantomAnimator::JOINT::SHOULDER));
            cout << "smooth 100" << endl;
            break;
        case ARMPIT:
            viewer.data(0).set_data(phantom->ShoulderShrink(viewer.data(0).V_normals, 0.02));
            cout << "armpit 0.02" << endl;
            break;
        case INNERARM:
            viewer.data(0).set_data(phantom->ArmOffset(viewer.data(0).V_normals, 1));
            cout << "inner arm 1" << endl;
            break;
        case INNERARM_SMALL:
            viewer.data(0).set_data(phantom->ArmOffset(viewer.data(0).V_normals, 0.5));
            cout << "inner arm 0.5" << endl;
            break;
        case LOWERARM:
            viewer.data(0).set_data(phantom->LowerArmOffset(viewer.data(0).V_normals, 0.5));
            cout << "inner arm -0.01" << endl;
            break;
        case LOWERARM_OUTVERYSMALL:
            viewer.data(0).set_data(phantom->LowerArmOffset(viewer.data(0).V_normals, 0.01));
            cout << "inner arm 0.01" << endl;
            break;
        case LOWERARM_INVERYSMALL:
            viewer.data(0).set_data(phantom->LowerArmOffset(viewer.data(0).V_normals, -0.01));
            cout << "inner arm -0.01" << endl;
            break;
        case COVER:
            cout << "fixing intersection" << flush;
            C0 = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetF().rows(), 1);
            C1 = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetFbone().rows(), 1);
            while (igl::copyleft::cgal::intersect_other(phantom->GetV(), phantom->GetF(), phantom->GetVbone(), phantom->GetFbone(), false, IF))
            {
                cout << " -> " << IF.rows() << flush;
                MatrixXi F0 = igl::slice(phantom->GetF(), IF.col(0), 1);
                vector<int> vec(F0.data(), F0.data() + F0.cols()*F0.rows());
                // VectorXi F0 = IF.col(0);
                // vector<int> vec(F0.data(), F0.data() + F0.cols()*F0.rows());
                // sort(vec.begin(), vec.end());
                // vec.erase(unique(vec.begin(), vec.end()), vec.end());
                phantom->MoveIntersected(vec, 0.02);
                MatrixXd R = RowVector3d(1., 0.3, 0.3).replicate(IF.rows(), 1);
                igl::slice_into(R, IF.col(0), 1, C0);
                igl::slice_into(R, IF.col(1), 1, C1);
            }
            cout << " -> " << IF.rows() << endl;
            viewer.data(1).set_colors(C1);
            viewer.data(0).set_colors(C0);
            break;
        case SHOULDERUP:
            vQ1[5] = Quaterniond(AngleAxisd(15. * PI / 180., Vector3d(0, 1, 0))); // male15
            vQ1[1] = Quaterniond(AngleAxisd(-15. * PI / 180., Vector3d(0, 1, 0)));
            phantom->Animate(vQ1, true);
            viewer.data(0).set_data(phantom->GetWeight(1) + phantom->GetWeight(5));
            cout << "shoulder up 15 deg." << endl;
            break;
        }
    };

    int selected(0), dataID(0);
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
            viewer.data(2).set_vertices(phantom->GetVo());
            viewer.data(2).compute_normals();

            phantom->LaplacianSmooth(PhantomAnimator::JOINT::WRIST, 5);
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::ELBOW, 5);
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 100);
            cout << "done" << endl;
            update = true;
            break;
        case '2': //elbow & wrist
            cout << "Deformation + Smoothing + Recovering (2)..." << flush;
            //male 50, female 55
            vQ1[6] = Quaterniond(AngleAxisd(53. * PI / 180., AngleAxisd(vQ[6]).axis()));
            vQ1[2] = Quaterniond(AngleAxisd(53. * PI / 180., AngleAxisd(vQ[2]).axis()));
            phantom->AnimateDQS(vQ1);
            viewer.data(2).set_vertices(phantom->GetVo());
            viewer.data(2).compute_normals();

            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 100);
            phantom->ReleaseRest();
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 100);
            phantom->ReleaseRest();
            phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 200);
            cout << "done" << endl;
            update = true;
            break;
        case '3':
            phantom->CalculateWeights(0.1, true);
            break;
        case '4':
        {
            MatrixXd color = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetFo().rows(), 1);
            if(!phantom->SetSkinLayers(color)) cout<<"fail!"<<endl;
            viewer.data(2).set_vertices(phantom->GetVo());
            viewer.data(2).set_colors(color);
            viewer.data(2).is_visible = true;
            viewer.data(2).compute_normals();
        }
        break;
        case 'R':
            Deformation(RELAX);
            update = true;
            break;
        case 'S':
            Deformation(SMOOTH50);
            update = true;
            break;
        case 'A':
            Deformation(INNERARM_SMALL);
            update = true;
            break;
        case 'X':
            Deformation(ARMPIT);
            update = true;
            break;
        case 'C':
            Deformation(COVER);
            update = true;
            break;
        case '/': //ignore the message "^tetrahedralize: Tetgen failed to create tets"
        {
            MatrixXd VT;
            MatrixXi TT, FT; // dummies
            igl::copyleft::tetgen::tetrahedralize(phantom->GetV(), phantom->GetF(), "d", VT, TT, FT);
        }
        break;
        case 'E':
        {
            MatrixXi IF;
            igl::copyleft::cgal::intersect_other(phantom->GetV(), phantom->GetF(), phantom->GetVbone(), phantom->GetFbone(), false, IF);
            MatrixXd C0 = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetF().rows(), 1);
            MatrixXd C1 = RowVector3d(1., 1., 1.).replicate(phantom->GetFbone().rows(), 1);
            MatrixXd R = RowVector3d(1., 0.3, 0.3).replicate(IF.rows(), 1);
            igl::slice_into(R, IF.col(0), 1, C0);
            igl::slice_into(R, IF.col(1), 1, C1);
            viewer.data(1).set_colors(C1);
            viewer.data(0).set_colors(C0);
            cout << "There are " << IF.rows() << " intersections between skin and bone" << endl;
        }
        break;
        case 'V':
            vol = phantom->GetVolume();
            cout << "current volume : " << vol << " / " << phantom->GetVolume(true) << " (" << vol / phantom->GetVolume(true) * 100.f << "%)" << endl;
            break;
        case 'P':
            igl::writePLY(string(argv[argc - 1]) + "_result.ply", phantom->GetV(), phantom->GetF());
            phantom->WriteMergedOBJ();
            break;
        case 'Z':
            viewer.data(2).is_visible = !viewer.data(2).is_visible;
            break;
        case 'B':
            viewer.data(0).is_visible = !viewer.data(0).is_visible;
            break;
        case ',':
            selected = min(max(--selected, 0), 8);
            viewer.data(0).set_data(phantom->GetWeight(selected));
            break;
        case '.':
            selected = min(max(++selected, 0), 8);
            viewer.data(0).set_data(phantom->GetWeight(selected));
            break;
        case '[':
            selected = min(max(--selected, 0), 3);
            viewer.data(0).set_data(phantom->GetSmoothW(PhantomAnimator::JOINT(selected)));
            break;
        case ']':
            selected = min(max(++selected, 0), 3);
            viewer.data(0).set_data(phantom->GetSmoothW(PhantomAnimator::JOINT(selected)));
            break;
        case ' ':
            viewer.core().is_animating = !viewer.core().is_animating;
            if (viewer.core().is_animating)
                viewer.core().animation_max_fps = 3;
            else
                viewer.core().animation_max_fps = 60;
            break;
        case '9':
            viewer.data(0).set_mesh(phantom->GetV(), phantom->GetF());
            viewer.data(0).compute_normals();
            break;
        case 'M':
            viewer.selected_data_index=dataID++%3;
            cout<<"current data ID: "<<viewer.selected_data_index<<endl;
        }
        if (update)
        {
            viewer.data(0).set_vertices(phantom->GetV());
            viewer.data(0).compute_normals();
            viewer.data(0).set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
        }
        return true;
    };

    vector<DEFORM> recipe = {RELAX, SMOOTH100, RELAX, SMOOTH100, INNERARM, INNERARM, INNERARM, ARMPIT, SMOOTH100, RELAX,
                             ARMPIT, SMOOTH100, RELAX, ARMPIT, SMOOTH50, RELAX, SMOOTH50, RELAX, INNERARM, INNERARM, RELAX, SMOOTH50,
                             SHOULDERUP, RELAX, RELAX, RELAX, INNERARM, SMOOTH50, COVER};
    int count(0);
    bool prestep(true);
    int interval(recipe.size());
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &viewer) -> bool
    {
        if (viewer.core().is_animating && count < recipe.size())
        {
            cout << count << " - ";
            Deformation(recipe[count++]);
            viewer.data(0).set_vertices(phantom->GetV());
            viewer.data(0).compute_normals();
            //if(prestep) viewer.core().is_animating = false;
            //else if(count == recipe.size()) viewer.core().is_animating = false;
        }
        if (viewer.core().is_animating && count == recipe.size())
        {
            MatrixXd Vtemp;
            MatrixXi Ttemp, Ftemp;
            if (prestep)
            {
                cout << "Start final adjustment (" << phantom->GetVolume(true) << ")" << endl;
                prestep = false;
            }
            double vol = phantom->GetVolume();
            cout << "[volume check] " << vol << " (" << vol / phantom->GetVolume(true) * 100.f - 100.f << "%)" << endl;
            double diff = vol / phantom->GetVolume(true) - 1.f;
            if (count > interval)
            {
                MatrixXd color = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetF().rows(), 1);
                cout << "[self-intersection]" << flush;
                int num = 1;
                while (num)
                {
                    num = phantom->DetectSelfInterSection(color, 0.5);
                    cout << " -> " << num << flush;
                }cout<<endl;

                // while (igl::copyleft::tetgen::tetrahedralize(phantom->GetV(), phantom->GetF(), "p/0.001YqO/3", Vtemp, Ttemp, Ftemp) != 0)
                // {
                //     phantom->LaplacianSmooth(PhantomAnimator::JOINT::CLAVICLE, 5);
                //     phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 5); //for fat
                //     Deformation(COVER);
                // }
                interval += 20;
            }
            if (abs(diff) < 0.0001)
            {
                if (igl::copyleft::tetgen::tetrahedralize(phantom->GetV(), phantom->GetF(), "p/0.001YqO/3", Vtemp, Ttemp, Ftemp) == 0)
                {
                    count++;
                    cout << "deformation finished!" << endl;
                    viewer.core().is_animating = false;
                    viewer.core().animation_max_fps = 60;
                }
                else
                {
                    MatrixXd color = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetF().rows(), 1);
                    cout << "[self-intersection]" << flush;
                    int num = 1;
                    while (num)
                    {
                        num = phantom->DetectSelfInterSection(color, 0.1);
                        cout << " -> " << num << flush;
                    }
                    // do
                    // {
                    //     phantom->LaplacianSmooth(PhantomAnimator::JOINT::CLAVICLE, 5);
                    //     phantom->LaplacianSmooth(PhantomAnimator::JOINT::SHOULDER, 5); //for fat
                    // }while (igl::copyleft::tetgen::tetrahedralize(phantom->GetV(), phantom->GetF(), "p/0.001YqO/3", Vtemp, Ttemp, Ftemp) != 0);
                    recipe.push_back(COVER);
                }
            }
            else if (diff > 0.005)
            {
                recipe.push_back(SMOOTH100);
                recipe.push_back(SMOOTH100);
                recipe.push_back(RELAX);
                recipe.push_back(COVER);
            }
            else if (diff > 0.0005)
            {
                recipe.push_back(SMOOTH50);
                recipe.push_back(RELAX);
                recipe.push_back(COVER);
            }
            else if (diff < -0.0005)
            {
                recipe.push_back(LOWERARM);
                //recipe.push_back(SMOOTH10);
                //recipe.push_back(SMOOTH1);
                recipe.push_back(RELAX);
                recipe.push_back(COVER);
            }
            else if (diff < 0)
            {
                recipe.push_back(LOWERARM_OUTVERYSMALL);
                recipe.push_back(COVER);
            }
            else if (diff > 0)
            {
                recipe.push_back(LOWERARM_INVERYSMALL);
                recipe.push_back(COVER);
            }
        }
        viewer.data(0).set_colormap(CM);
        return false;
    };
    viewer.data(0).is_visible = true;
    viewer.data(1).is_visible = true;
    viewer.core().is_animating = false;
    viewer.launch();
    return EXIT_SUCCESS;
}
