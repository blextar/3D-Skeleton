#include <vector>
#include "octree.h"
#include "misc.h"

using namespace std;

class h_vertex{
public:

	double normal[3];
	double w,nl;
	// pointer to the real node
	node* n;

	h_vertex();
	h_vertex(node *m, node *center);
};

class h_face{
public:
	// vertices belonging to the face
	h_vertex *vts[3];

	h_face();
	h_face(h_vertex* v1, h_vertex* v2, h_vertex* v3);

	bool operator == (const h_face &f) const{
		return (f.vts[0]->n==vts[0]->n && f.vts[1]->n==vts[1]->n && f.vts[2]->n==vts[2]->n)
		|| (f.vts[0]->n==vts[0]->n && f.vts[1]->n==vts[2]->n && f.vts[2]->n==vts[1]->n)
		|| (f.vts[0]->n==vts[1]->n && f.vts[1]->n==vts[0]->n && f.vts[2]->n==vts[2]->n)
		|| (f.vts[0]->n==vts[1]->n && f.vts[1]->n==vts[2]->n && f.vts[2]->n==vts[0]->n)
		|| (f.vts[0]->n==vts[2]->n && f.vts[1]->n==vts[0]->n && f.vts[2]->n==vts[1]->n)
		|| (f.vts[0]->n==vts[2]->n && f.vts[1]->n==vts[1]->n && f.vts[2]->n==vts[0]->n);
	}
};

class hull{
public:

	vector<h_face> faces;
	vector<h_vertex> vts;

	hull(node *n,vector<node*> &neighbors);
	hull(double x, double y, double z, vector<node*> &neighbors);
};
