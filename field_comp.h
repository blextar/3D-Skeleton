#ifndef FIELD_COMP_H
#define	FIELD_COMP_H

#include<algorithm>
#include<math.h>
#include<vector>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<queue>

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_eigen.h>

// miscellanea
#include"misc.h"
// octree definition
#include"octree.h"
// distance transform
#include"distTrans.h"

#include"hull.h"
#include"lap_hull.h"

static double eps=1e-4;
using namespace std;

void neighborhood(vector<node*>& neigh, node* n);
void fast_neighborhood(vector<node*>& neigh, node* n);
int CStar(node *n);
int CMinus(node *n);
void grad(vector<node*>& leaves);
void lap_from_grad(vector<node*>& input, unsigned short int density_corr);
void density(vector<node*>& input);
void dilation(vector<node*>& leaves, vector<node*>& input,vector<node*>& dleaves,octree *oct);
bool skel_dilation(vector<node*>& skel,vector<node*>& new_skel,octree *oct);
void set_dborder(vector<node*>& skel_input,vector<node*>& dborder);
void reset_dborder(vector<node*>& dborder);
vector<node*> skeletonize(vector<node*>& input, double treshold, unsigned short int density_corr);
vector<node*> ensure_thinness(vector<node*>& input);
void field(vector<node*>& input);
bool endpoint(node* n);
void keep_largest_cc(vector<node*>& interior);
void smooth_dm(vector<node*>& leaves,unsigned short int smooth,double alpha);

#endif