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
cout<<3<<flush;

    AngleAxisd aa7 = AngleAxisd(vQ[7]);
    AngleAxisd aa3 = AngleAxisd(vQ[3]);

    AngleAxisd aa6 = AngleAxisd(vQ[6]);
    AngleAxisd aa2 = AngleAxisd(vQ[2]);
cout<<4<<flush;

    // ofstream ofs("hoorey.txt");
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
        MatrixXd C0,C1;
        switch (key)
        {
        case 'D':
            // intTmp = phantom->DetectIntersections(C1, a);
            // intTmp += phantom->DetectIntersections(C1, b);
            C0 = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetFI().rows(), 1);
            C1 = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetF().rows(), 1);
            intTmp = 1;
            cout<<"[outer-inner] intersection"<<flush;
            {int cnt(0);
            while(intTmp)
            {
                intTmp = phantom->DetectIntersections(C1, C0, true);
                cout<<" -> "<<intTmp<<flush; cnt++;
                if(cnt>20) break;
            }cout<<endl;}
            // intTmp = phantom->DetectIntersections(C1);
            // if (intTmp)
            // {
            //     cout << intTmp << " intersection was detected!" << endl;
            // }
            // else
            //     cout << "no intersection was detected!" << endl;
            viewer.data(0).set_vertices(phantom->GetVI());
            viewer.data(1).set_colors(C1);
            viewer.data(0).set_colors(C0);
            break;
        case 'C':
        // {
        //     MatrixXd VV;
        //     MatrixXi FF, TT;
        //     MatrixXd V(phantom->GetV().rows() + phantom->GetVI().rows(), 3);
        //     V.topLeftCorner(phantom->GetV().rows(), 3) = phantom->GetV();
        //     V.bottomLeftCorner(phantom->GetVI().rows(), 3) = phantom->GetVI();
        //     MatrixXi F(phantom->GetF().rows() + phantom->GetFI().rows(), 3);
        //     F.topLeftCorner(phantom->GetF().rows(), 3) = phantom->GetF();
        //     F.bottomLeftCorner(phantom->GetFI().rows(), 3) = phantom->GetFI().array() + phantom->GetV().rows();
        //     igl::copyleft::tetgen::tetrahedralize(V, F, "p/0.001YqO/3", VV, TT, FF);
        // }
            C1 = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetF().rows(), 1);
            cout <<"[outer] self-intersection"<<flush;
            intTmp = 1;
            while(intTmp)
            {
                intTmp = phantom->SelfIntersection(C1, PhantomAnimator::SHELL::OUTER, 0.1);
                cout<<" -> "<<intTmp <<flush;
            }
            viewer.data(1).set_colors(C1);
            cout << endl<<"[inner] self-intersection"<<flush;
            C1 = RowVector3d(255. / 255, 246. / 255., 51. / 255.).replicate(phantom->GetFI().rows(), 1);
            intTmp = 1;
            while(intTmp)
            {
                intTmp = phantom->SelfIntersection(C1, PhantomAnimator::SHELL::INNER, 0.1);
                cout<<" -> "<<intTmp <<flush;
            }cout<<endl;
            viewer.data(0).set_vertices(phantom->GetVI());
            viewer.data(1).set_vertices(phantom->GetV());
            viewer.data(0).set_colors(C1);
            // intTmp = phantom->SelfIntersection(C1, PhantomAnimator::SHELL::WHOLE, viewer.data(1).V_normals, viewer.data(0).V_normals, true);
            // if (intTmp)
            //     cout <<"[whole] "<< intTmp << " intersection was detected!" << endl;
            // else
            //     cout << "[whole] no intersection was detected!" << endl;
            
            // phantom->CalculateElbowW(3, aa3.axis());//, 9, 100);
            // phantom->CalculateElbowW(7, aa7.axis());//, 9, 100);
            // phantom->CalculateShoulderW(1, -Vector3d(0, 1, 0));
            // phantom->CalculateShoulderW(5, Vector3d(0, 1, 0));
            // phantom->CalculateCleanWeights();
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
            phantom->CalculateElbowW(3, aa3.axis(), 1, 30); //, 9, 100);
            phantom->CalculateElbowW(7, aa7.axis(), 1, 30); //, 9, 100);
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
            viewer.data(1).set_vertices(phantom->GetV());
            viewer.data(1).compute_normals();
            viewer.data(1).set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
            break;
        case '2': //shoulder
        {
            double theta6 = max(min(53., phantom->CalculateShoulderW(6, aa6.axis(), 1, 30) - 5.),0.);
            double theta2 = max(min(53., phantom->CalculateShoulderW(2, aa2.axis(), 1, 30) - 5.),0.);  
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
            viewer.data(1).set_vertices(phantom->GetV());
            viewer.data(1).compute_normals();
            viewer.data(1).set_edges(phantom->GetC(), phantom->GetBE(), sea_green);
            break;
        case '3':
            phantom->AdjustVolume(1, 0, 1, viewer.data().V_normals);
            phantom->AdjustVolume(5, 0, 1, viewer.data().V_normals);
            phantom->AdjustVolume(2, 1, 0.5, viewer.data().V_normals);
            phantom->AdjustVolume(6, 1, 0.5, viewer.data().V_normals);
            phantom->AdjustVolume(3, 1, 0.5, viewer.data().V_normals);
            phantom->AdjustVolume(7, 1, 0.5, viewer.data().V_normals);
            //phantom->LaplacianSmooth(1, 10, 0, 1);
            viewer.data().set_vertices(phantom->GetV());
            break;
        case '4': //inner
            phantom->SetInnerVertices();
            viewer.data(0).set_vertices(phantom->GetVI());
            viewer.data(0).compute_normals();
            break;
        case '5':
            phantom->AdjustInnerVolume(2500, 1, 0, 1.1, viewer.data(0).V_normals); 
            phantom->AdjustInnerVolume2(1400, 2, 1.1, 0.5, viewer.data(0).V_normals);
            phantom->AdjustInnerVolume(1500, 2, 1.2, 0.3, viewer.data(0).V_normals); 
            phantom->AdjustInnerVolume(1700, 2, 1.1, 0.5, viewer.data(0).V_normals); 
            phantom->AdjustInnerVolume(2000, 3, 1.1, 0.5, viewer.data(0).V_normals); 
            phantom->AdjustInnerVolume(2100, 3, 1.1, 0.5, viewer.data(0).V_normals); 
            viewer.data(0).set_vertices(phantom->GetVI());
            break;
        case 'I':
            phantom->InitializeV();
            viewer.data(1).set_vertices(phantom->GetV());
            viewer.data(1).compute_normals();
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
