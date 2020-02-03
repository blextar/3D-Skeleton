#ifndef SKEL_MESH_H
#define	SKEL_MESH_H

#include"field_comp.h"

using namespace std;

vector<node*> refine_skeleton(vector<node*>& leaves, vector<node*>& skeleton, Tree& mesh);

#endif	/* SKEL_MESH_H */

