#include "octree.h"

static double eps = 1e-7;

#define DOUBLE_EQ(a,b) (fabs(a-b)<eps)

void node::recAdd(int level) {
    if (level > 0) {
        for (char i = 0; i < 8; i++) {
            node n(this, i);
            children.push_back(n);
        }
        for (unsigned int i = 0; i < 8; i++) {
            children[i].populateNeigh();
        }
        for (unsigned int i = 0; i < 8; i++)
            children[i].recAdd(level - 1);
    }
}

void node::splitNode() {
    for (unsigned int i = 0; i < 8; i++) {
        node n(this, i);
        children.push_back(n);
    }
    for (unsigned int i = 0; i < 8; i++) {
        children[i].populateNeigh();
    }
}

int node::level() {
    return id.size();
}

bool node::border() {

    return neighbors[0] == NULL || neighbors[1] == NULL || neighbors[2] == NULL
            || neighbors[3] == NULL || neighbors[4] == NULL || neighbors[5] == NULL;
}

bool node::leaf() {
    return children.size() == 0;
}

void node::populateNeigh() {

    node *p = parent;

    if (*(id.end() - 1) == 0) {

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[0], 2);
        neighbors.push_back(&(p->children[1]));
        neighbors.push_back(&(p->children[2]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[3], 1);

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[4], 4);
        neighbors.push_back(&(p->children[4]));
    } else if (*(id.end() - 1) == 1) {
        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[0], 3);

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[1], 0);

        neighbors.push_back(&(p->children[3]));
        neighbors.push_back(&(p->children[0]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[4], 5);

        neighbors.push_back(&(p->children[5]));
    } else if (*(id.end() - 1) == 2) {
        neighbors.push_back(&(p->children[0]));
        neighbors.push_back(&(p->children[3]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[2], 0);

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[3], 3);

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[4], 6);

        neighbors.push_back(&(p->children[6]));
    } else if (*(id.end() - 1) == 3) {
        neighbors.push_back(&(p->children[1]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[1], 2);

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[2], 1);

        neighbors.push_back(&(p->children[2]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[4], 7);

        neighbors.push_back(&(p->children[7]));
    } else if (*(id.end() - 1) == 4) {
        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[0], 6);

        neighbors.push_back(&(p->children[5]));
        neighbors.push_back(&(p->children[6]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[3], 5);

        neighbors.push_back(&(p->children[0]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[5], 0);
    } else if (*(id.end() - 1) == 5) {
        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[0], 7);

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[1], 4);

        neighbors.push_back(&(p->children[7]));
        neighbors.push_back(&(p->children[4]));
        neighbors.push_back(&(p->children[1]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[5], 1);
    } else if (*(id.end() - 1) == 6) {
        neighbors.push_back(&(p->children[4]));
        neighbors.push_back(&(p->children[7]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[2], 4);

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[3], 7);

        neighbors.push_back(&(p->children[2]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[5], 2);
    } else if (*(id.end() - 1) == 7) {
        neighbors.push_back(&(p->children[5]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[1], 6);

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[2], 5);

        neighbors.push_back(&(p->children[6]));
        neighbors.push_back(&(p->children[3]));

        if (level() == 1)
            neighbors.push_back(NULL);
        else
            remoteNeigh(p->neighbors[5], 3);
    }
}

// recursively check the neighbor's children in order to update their neighborhood.
// if they listed n->parent as a neighbor, and their level is >= n->level
// then it should be updated to n.

void node::recPopNeigh(node *c) {
    if (c->level() >= level()) {
        for (int i = 0; i < c->neighbors.size(); ++i) {
            if (c->neighbors[i] != NULL) {
                if (*(c->neighbors[i]) == *(parent)) {
                    c->neighbors[i] = this;
                    break;
                }
            }
        }
        if (!c->leaf()) {
            for (int i = 0; i < c->children.size(); ++i) {
                recPopNeigh(&(c->children[i]));
            }
        }
    }
}

void node::remoteNeigh(node *pn, unsigned int id) {
    if (pn == NULL)
        neighbors.push_back(NULL);
    else if (pn->leaf())
        neighbors.push_back(pn);
    else {
        neighbors.push_back(&(pn->children[id]));
        recPopNeigh(&(pn->children[id]));
    }
}

void octree::getLeaves(vector<node*>& leaves, node *n) {
    if (n->leaf())
        leaves.push_back(n);
    else {
        for (unsigned int i = 0; i < 8; i++) {
            getLeaves(leaves, &(n->children[i]));
        }
    }
}

void octree::print() {
    cerr << endl;
    printRec(root);
}

void octree::printRec(node *n) {
    for (int i = 0; i < n->level(); ++i)
        cerr << " ";
    cerr << *n << endl;
    if (!n->leaf()) {
        for (unsigned int i = 0; i < 8; i++) {
            printRec(&(n->children[i]));
        }
    }
}

int octree::countNodes(node *n, int& count) {
    count++;
    if (!n->leaf()) {
        for (unsigned int i = 0; i < 8; i++) {
            countNodes(&(n->children[i]), count);
        }
    }
    return count;
}

double dist(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt((x1 - x2)*(x1 - x2)+(y1 - y2)*(y1 - y2)+(z1 - z2)*(z1 - z2));
}

double dist(node *n, double x, double y, double z) {
    return sqrt((n->x - x)*(n->x - x)+(n->y - y)*(n->y - y)+(n->z - z)*(n->z - z));
}

double dist(node *n, node *m) {
    return sqrt((n->x - m->x)*(n->x - m->x)+(n->y - m->y)*(n->y - m->y)+(n->z - m->z)*(n->z - m->z));
}

bool divOrd(const node* n1, const node* n2) {
    return n1->div > n2->div;
}

bool lapOrd(const node* n1, const node* n2) {
    return n1->lap > n2->lap;
}

bool distOrd(const node* n1, const node* n2) {
    if (fabs(n1->dist - n2->dist) > 0.3 * n1->size) {
        return n1->dist < n2->dist;
    } else {
        return posOrd(n1, n2);
    }
}

bool posOrd(const node* n1, const node* n2) {
    return (n1->x < n2->x) || ((n1->x == n2->x) && (n1->y < n2->y))
            || ((n1->x == n2->x) && (n1->y == n2->y) && (n1->z < n2->z));
}

// look for the closest node to "center" in the tree
// whose root is "n"

node* findClosestNode(node *n, node *center) {

    if (n == NULL)
        return NULL;
    if (n->leaf() || n->level() == center->level())
        return n;

    node* closest = &(n->children[0]);

    double d_min = dist(closest, center);
    double d;

    for (int i = 0; i < n->children.size(); ++i) {
        d = dist(&(n->children[i]), center);
        if (d < d_min) {
            closest = &(n->children[i]);
            d_min = d;
        }
    }
    return findClosestNode(closest, center);
}

bool aligned(node *n, int level, int id1, int id2, int id3, int id4) {

    for (int i = n->level() - 2; i > level - 1; --i) {
        if (n->id[i] == id1 || n->id[i] == id2 || n->id[i] == id3 || n->id[i] == id4) {
            if (i == level)
                return true;
        } else
            return false;
    }
    cerr << "error\n";
}

// returns true if n is not a 6-neighbor of c

bool isConsistent(node *n, node *c) {
    for (vector<node*>::iterator it = c->neighbors.begin(); it != c->neighbors.end(); ++it) {
        if (*n == **it)
            return false;
    }
    return true;
}

// returns true if n is not a 18-neighbor of c

bool isUberConsistent(node *n, node *c, vector<node*>& neigh) {

    for (int i = 0; i < neigh.size(); i++)
        if (neigh[i] != NULL)
            if (i != 0 && i != 2 && i != 6 && i != 8 && i != 18 && i != 20 && i != 24 && i != 26)
                if (*n == *neigh[i])
                    return false;
    return true;
}
