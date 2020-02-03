#ifndef OCTREE_H
#define	OCTREE_H

#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <list>

#include"cgal.h"
    
using namespace std;

static double inf = numeric_limits<double>::infinity();

class node {
public:

    // pointers to parent and children
    node *parent;
    vector<node> children;
    // 6-neighbors
    vector<node*> neighbors;

    // background connected components of a surface voxel
    vector<node*> B, C;

    // voxel information
    double x, y, z, size, w;
    // nearest point on the mesh
    Point np;

    // values used for the skeleton computation
    double dist, lap, nlap, dens, div; //, sdist;
    // distance from the center of the voxel's faces to the border
    // used for gradient computation
    double grad[3];

    bool skel, interior, dborder, flag, candidate;
    vector<char> id;

    node() {
        parent = NULL;
        x = 0.0;
        y = 0.0;
        z = 0.0;
        size = 0.0;
        dist = inf;
        lap = 0.0;
        nlap = 0.0;
        dens = 0.0;
        div = 0.0;
        //sdist=0.0;

        flag = false;
        skel = false;
        interior = false;
        candidate = false;
        dborder = false;
    }

    node(node *p, char idx) {

        id = p->id;
        id.push_back(idx);
        parent = p;
        
        flag = false;
        skel = false;
        interior = false;
        candidate = false;
        dborder = false;

        size = p->size / 2.0;

        double offset = size / 2.0;

        if (idx == 0) {
            x = p->x - offset;
            y = p->y - offset;
            z = p->z - offset;
        } else if (idx == 1) {
            x = p->x + offset;
            y = p->y - offset;
            z = p->z - offset;
        } else if (idx == 2) {
            x = p->x - offset;
            y = p->y + offset;
            z = p->z - offset;
        } else if (idx == 3) {
            x = p->x + offset;
            y = p->y + offset;
            z = p->z - offset;
        } else if (idx == 4) {
            x = p->x - offset;
            y = p->y - offset;
            z = p->z + offset;
        } else if (idx == 5) {
            x = p->x + offset;
            y = p->y - offset;
            z = p->z + offset;
        } else if (idx == 6) {
            x = p->x - offset;
            y = p->y + offset;
            z = p->z + offset;
        } else if (idx == 7) {
            x = p->x + offset;
            y = p->y + offset;
            z = p->z + offset;
        }

        dist = inf;
        lap = 0.0;
        nlap = 0.0;
        dens = p->dens;
        div = p->div;

        // the father is no more a leaf
        p->skel = false;
    }

    bool operator ==(const node &n) const {
        if (n.id.size() != id.size())
            return false;
        else {
            for (int i = 0; i < id.size(); ++i) {
                if (id[i] != n.id[i])
                    return false;
            }
        }
        return true;
    }

    bool operator !=(const node &n) const {
        if (n.id.size() != id.size())
            return true;
        else {
            for (int i = 0; i < id.size(); ++i) {
                if (id[i] != n.id[i])
                    return true;
            }
        }
        return false;
    }

    int level();
    bool border();
    bool leaf();
    void populateNeigh();
    void recPopNeigh(node *c);
    void remoteNeigh(node *pn, unsigned int id);
    void recAdd(int level);
    void splitNode();
};

inline ostream & operator <<(ostream &out, vector<char> id) {

    for (vector<char>::iterator it = id.begin(); it != id.end(); ++it) {
        out << (int) *it;
    }
    return out;
}

inline ostream & operator <<(ostream &out, node &n) {

    double lap_sum = 0.0;
    double dens_sum = 0.0;
    double div_sum = 0.0;
    out << '(' << n.x << ',' << n.y << ',' << n.z << ')' << " LEVEL=" << n.level()
            << " ID=" << n.id << " SKEL=" << n.skel << " INTERIOR=" << n.interior
            << " DIST=" << n.dist
            << " LAP=" << n.lap;

    if (n.children.size() != 0) {
        for (vector<node>::iterator it = n.children.begin(); it != n.children.end(); ++it) {
            lap_sum += it->lap;
            dens_sum += it->dens;
            div_sum += it->div;
        }
        lap_sum /= 8.0;
        dens_sum /= 8.0;
        div_sum /= 8.0;
    }
    out << " ALAP=" << lap_sum;
    out << " DENS=" << n.dens;
    out << " ADENS=" << dens_sum;
    out << " DIV=" << n.div;
    out << " ADIV=" << div_sum;

    out << endl;
    return out;
}

class octree {
public:
    node *root;

    octree(unsigned int size, node *r) {
        root = r;
        r->recAdd(size);
    }

    void print();
    void printRec(node *n);
    int countNodes(node *n, int& count);
    void getLeaves(vector<node*>& leaves, node *n);
};

bool lapOrd(const node* n1, const node* n2);
bool divOrd(const node* n1, const node* n2);
bool distOrd(const node* n1, const node* n2);
bool posOrd(const node* n1, const node* n2);

double dist(double x1, double y1, double z1, double x2, double y2, double z2);
double dist(node *n, double x, double y, double z);
double dist(node *n, node *m);

node* findClosestNode(node *n, node* center);
bool aligned(node *n, int level, int id1, int id2, int id3, int id4);
bool isConsistent(node *n, node *c);
bool isUberConsistent(node *n, node *c, vector<node*>& neigh);

#endif