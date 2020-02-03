#include<fstream>
#include<sstream>
#include<cstring>
#include<cstdlib>
#include<ctime>
#include<float.h>
#include<algorithm>

#include "field_comp.h"
#include "point_cloud.h"
#include "skel_mesh.h"

using namespace std;

int main(int argc, char *argv[]) {

    Tree mesh;
    double size;
    clock_t timer = clock();
    unsigned int nPts;

    if (argc < 5 || argc > 5) {
        cerr << "usage: " << argv[0] << " filename initial_resolution refinement_levels div_threshold\n";
        exit(2);
    }

    double res = atof(argv[2]);
    if (log2(res)-(int) log2(res) != 0) {
        cerr << "ERROR: the resolution value must be a power of 2\n";
        exit(2);
    }

    stringstream ssname;
    ssname << argv[1];
    string str = ssname.str();
    string filename=str;
    unsigned short int zoom = atoi(argv[3]);
    double threshold=-atof(argv[4]);
    // no smoothing
    double smooth=0;
    double alpha=0.0;

    unsigned short int density_corr = 1;

    cerr << "Initializing Data Structures: ";
    clock_t start = clock();

    size = OFFReader(argv[1], mesh, res);
    node root;
    root.size = ceil(size) + 2;
    root.x = root.size / 2.0;
    root.y = root.size / 2.0;
    root.z = root.size / 2.0;

    octree oct(log2(res), &root);

    vector<node*> leaves, interior;
    oct.getLeaves(leaves, &root);
    cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

    cerr << "Distance Transform computation: ";
    start = clock();
    dt(mesh, leaves, interior, NORMAL);
    cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

    cvlab::cloud3d input;

    for (vector<node*>::iterator it = leaves.begin(); it != leaves.end(); it++) {
        if((*it)->interior)
            input.add_point(cvlab::point3d((*it)->x, (*it)->y, (*it)->z), cvlab::red);
        else
            input.add_point(cvlab::point3d((*it)->x, (*it)->y, (*it)->z), cvlab::yellow);
    }
    
    cerr << "Smoothing: ";
    start = clock();
    smooth_dm(leaves,smooth,alpha);
    cerr << "completed in " << ((double)clock() - start) / CLOCKS_PER_SEC << " sec." << endl;
    
    if (zoom > 0) {
        cerr << "Detecting the largest connected component: ";
        start = clock();
        keep_largest_cc(interior);
        cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;
    }
    cerr << "Calculating the divergence and extracting the skeleton:\n\n";
    clock_t div_start = clock();

    vector<node*> skel_input, new_leaves, skeleton, dborder;

    cerr << "Skeletonization #0...\n";
    sort(interior.begin(), interior.end(), distOrd);
    sort(leaves.begin(), leaves.end(), distOrd);

    cerr << "Grad: ";
    start = clock();
    grad(leaves);
    cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

    cerr << "Laplacian: ";
    start = clock();
    lap_from_grad(leaves, density_corr);
    cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

    if (density_corr) {
        cerr << "Density: ";
        start = clock();
        density(interior);
        cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

        cerr << "Divergence: ";
        start = clock();
        field(interior);
        cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;
    }

    cerr << "Thinning: ";
    start = clock();
    
    skeleton = skeletonize(interior, threshold, density_corr);
    while (skeleton.size() != interior.size()) {
        interior = skeleton;
        sort(interior.begin(), interior.end(), distOrd);
        skeleton = skeletonize(interior, threshold, density_corr);
    }
    
    if (zoom == 0) {
        skeleton = ensure_thinness(interior);
        while (skeleton.size() != interior.size()) {
            interior = skeleton;
            sort(interior.begin(), interior.end(), distOrd);
            skeleton = ensure_thinness(interior);
        }
    }
    cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

    double contrast = 0.4;
    double max_ = -inf;
    double min = inf;

    cerr << endl;

    for (int i = 0; i < zoom; i++) {

        stringstream ss;
        ss << i;
		smooth*=2;
		alpha+=(1.0-alpha)/zoom;
                
        cerr << "Skeletonization #" << i + 1 << "...\n";

        new_leaves.clear();
        skel_input.clear();

        for (vector<node*>::iterator it = skeleton.begin(); it != skeleton.end(); ++it) {
            (*it)->splitNode();
            for (vector<node>::iterator j = (*it)->children.begin(); j != (*it)->children.end(); ++j) {
                new_leaves.push_back(&(*j));
            }
        }
        
        cerr << "Distance Transform computation: ";
        start = clock();
        dt(mesh, new_leaves, skel_input, EXPLOSION);
        cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;
        
        vector<node*> dilated_leaves;
        cerr << "Dilation: ";
        start = clock();
        dilation(new_leaves, skel_input, dilated_leaves, &oct);
        cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

        cerr << "Distance Transform computation: ";
        start = clock();
        dt(mesh, dilated_leaves, skel_input, DILATION);
        cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;
        
        new_leaves.insert(new_leaves.end(), dilated_leaves.begin(), dilated_leaves.end());

        cerr << "Smoothing: ";
        start = clock();
        smooth_dm(new_leaves,smooth,alpha);
        cerr << "completed in " << ((double)clock() - start) / CLOCKS_PER_SEC << " sec." << endl;
            
        set_dborder(skel_input, dborder);

        sort(new_leaves.begin(), new_leaves.end(), distOrd);
        sort(skel_input.begin(), skel_input.end(), distOrd);

        cerr << "Grad: ";
        start = clock();
        grad(new_leaves);
        cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

        cerr << "Lap: ";
        start = clock();
        lap_from_grad(new_leaves, density_corr);
        cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

        if (density_corr) {
            cerr << "Density: ";
            start = clock();
            density(skel_input);
            cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;

            cerr << "Divergence: ";
            start = clock();
            field(skel_input);
            cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;
        }

        cerr << "Thinning: ";
        start = clock();
        
        skeleton = skeletonize(skel_input, threshold, density_corr);
        while (skeleton.size() != skel_input.size()) {
            skel_input = skeleton;
            sort(skel_input.begin(), skel_input.end(), distOrd);
            skeleton = skeletonize(skel_input, threshold, density_corr);
        }
        
        cerr << "completed in " << ((double) clock() - start) / CLOCKS_PER_SEC << " sec." << endl;


        reset_dborder(dborder);

        vector<node*> new_skel;
        int num = 0;

        cerr << "Skeleton dilation...";
        clock_t sd_start = clock();

        while (skel_dilation(skeleton, new_skel, &oct)) {
            
            skel_input.clear();
            dt(mesh, new_skel, skel_input, DILATION);
            set_dborder(skel_input, dborder);

            sort(new_skel.begin(), new_skel.end(), distOrd);
            sort(skel_input.begin(), skel_input.end(), distOrd);

            grad(new_skel);
            lap_from_grad(new_skel, density_corr);

            if (density_corr) {
                density(skel_input);
                field(skel_input);
            }

            skeleton.clear();

            skeleton = skeletonize(skel_input, threshold, density_corr);
            while (skeleton.size() != skel_input.size()) {
                skel_input = skeleton;
                sort(skel_input.begin(), skel_input.end(), distOrd);
                skeleton = skeletonize(skel_input, threshold, density_corr);
            }

            reset_dborder(dborder);

            num += skeleton.size();
            new_skel.clear();
        }

        cerr << num << " voxels recovered in ";
        cerr << ((double) clock() - sd_start) / CLOCKS_PER_SEC << " sec." << endl;

        new_leaves.clear();
        oct.getLeaves(new_leaves, &root);

        skeleton.clear();
        skel_input.clear();
        for (vector<node*>::iterator it = new_leaves.begin(); it != new_leaves.end(); it++) {
            (*it)->dborder = false;
            (*it)->candidate = false;

            if ((*it)->skel)
                skeleton.push_back(*it);
            if ((*it)->interior)
                skel_input.push_back(*it);
        }
        
        if (i == zoom - 1) {
            vector<node*> old_skeleton = skeleton;          
            sort(old_skeleton.begin(), old_skeleton.end(), distOrd);
            skeleton = ensure_thinness(old_skeleton);
            
            while (skeleton.size() != old_skeleton.size()) {
                old_skeleton = skeleton;
                sort(old_skeleton.begin(), old_skeleton.end(), distOrd);
                skeleton = ensure_thinness(old_skeleton);
            }
        }

        skeleton.clear();
        skel_input.clear();
        for (vector<node*>::iterator it = new_leaves.begin(); it != new_leaves.end(); it++) {
            (*it)->dborder = false;
            (*it)->candidate = false;

            if ((*it)->skel)
                skeleton.push_back(*it);
            if ((*it)->interior)
                skel_input.push_back(*it);
        }
        cerr << endl;
    }

    cerr << "\ncompleted in " << ((double) clock() - div_start) / CLOCKS_PER_SEC << " sec." << endl;

    leaves.clear();
    oct.getLeaves(leaves, &root);

    interior.clear();
    for (vector<node*>::iterator it = leaves.begin(); it != leaves.end(); ++it) {
        if ((*it)->interior)
            interior.push_back(*it);
    }

	skeleton=refine_skeleton(leaves, skeleton, mesh);
    
    if(density_corr)
        filename=filename+"_hier";
    else
        filename=filename+"_ham";
    distmap2vtk(skeleton, filename);
    cerr << "\nTOTAL time " << ((double) clock() - timer) / CLOCKS_PER_SEC << " sec." << endl;

    return 0;
}