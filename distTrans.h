#ifndef DISTTRANS_H
#define	DISTTRANS_H

#include <vector>
#include <iomanip>
#include <limits>
#include "octree.h"
#include "field_comp.h"

#define NORMAL 0
#define EXPLOSION 1
#define DILATION 2

double OFFReader(char* inFile, Tree& mesh, double res);
void dt(Tree& mesh, vector<node*>& leaves, vector<node*>& skel_input, int type);
void distmap2vtk(vector<node*> leaves, string filename);

#endif