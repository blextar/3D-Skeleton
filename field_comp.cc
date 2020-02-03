#include"field_comp.h"

#define CMINUS 0
#define CSTAR 1

static inline double nnan(double x) {
    return (isnan(x)) ? 1.0 : x;
}

void fit_hyperplane(node *n) {

    double w;
    int s;

    gsl_matrix* A = gsl_matrix_calloc(4, 4);
    gsl_vector* b = gsl_vector_calloc(4);
    gsl_vector* x = gsl_vector_calloc(4);
    gsl_permutation* p = gsl_permutation_alloc(4);

    vector<node*> neighbors;
    neighborhood(neighbors, n);
    for (vector<node*>::iterator jt = neighbors.begin(); jt != neighbors.end(); ++jt) {

        node *m = *jt;

        if (m != NULL) {
            if (*n != *m) {
                w = exp(-pow(dist(n, m), 2) / (1 / n->size * 1 / n->size));

                gsl_matrix_set(A, 0, 0, gsl_matrix_get(A, 0, 0) + w);
                gsl_matrix_set(A, 0, 1, gsl_matrix_get(A, 0, 1) + w * (m->x - n->x));
                gsl_matrix_set(A, 0, 2, gsl_matrix_get(A, 0, 2) + w * (m->y - n->y));
                gsl_matrix_set(A, 0, 3, gsl_matrix_get(A, 0, 3) + w * (m->z - n->z));

                gsl_matrix_set(A, 1, 0, gsl_matrix_get(A, 1, 0) + w * (m->x - n->x));
                gsl_matrix_set(A, 1, 1, gsl_matrix_get(A, 1, 1) + w * (m->x - n->x)*(m->x - n->x));
                gsl_matrix_set(A, 1, 2, gsl_matrix_get(A, 1, 2) + w * (m->y - n->y)*(m->x - n->x));
                gsl_matrix_set(A, 1, 3, gsl_matrix_get(A, 1, 3) + w * (m->z - n->z)*(m->x - n->x));

                gsl_matrix_set(A, 2, 0, gsl_matrix_get(A, 2, 0) + w * (m->y - n->y));
                gsl_matrix_set(A, 2, 1, gsl_matrix_get(A, 2, 1) + w * (m->x - n->x)*(m->y - n->y));
                gsl_matrix_set(A, 2, 2, gsl_matrix_get(A, 2, 2) + w * (m->y - n->y)*(m->y - n->y));
                gsl_matrix_set(A, 2, 3, gsl_matrix_get(A, 2, 3) + w * (m->z - n->z)*(m->y - n->y));

                gsl_matrix_set(A, 3, 0, gsl_matrix_get(A, 3, 0) + w * (m->z - n->z));
                gsl_matrix_set(A, 3, 1, gsl_matrix_get(A, 3, 1) + w * (m->x - n->x)*(m->z - n->z));
                gsl_matrix_set(A, 3, 2, gsl_matrix_get(A, 3, 2) + w * (m->y - n->y)*(m->z - n->z));
                gsl_matrix_set(A, 3, 3, gsl_matrix_get(A, 3, 3) + w * (m->z - n->z)*(m->z - n->z));

                gsl_vector_set(b, 0, gsl_vector_get(b, 0) + w * m->dist);
                gsl_vector_set(b, 1, gsl_vector_get(b, 1) + w * m->dist * (m->x - n->x));
                gsl_vector_set(b, 2, gsl_vector_get(b, 2) + w * m->dist * (m->y - n->y));
                gsl_vector_set(b, 3, gsl_vector_get(b, 3) + w * m->dist * (m->z - n->z));
            }
        }
    }

    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_solve(A, p, b, x);

    for (int i = 1; i < 4; i++) {
        n->grad[i - 1] = gsl_vector_get(x, i);
    }
    normalize(n->grad);

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_permutation_free(p);
}

void grad(vector<node*>& leaves) {

    for (vector<node*>::iterator it = leaves.begin(); it != leaves.end(); ++it) {
        node *n = *it;
        if (!n->border()) {
            fit_hyperplane(n);
        }
    }
}

double tetra_vol(node* n1, node* n2, node* n3, node* n4) {
    double v1[3] = {n1->x - n4->x, n1->y - n4->y, n1->z - n4->z};
    double v2[3] = {n2->x - n4->x, n2->y - n4->y, n2->z - n4->z};
    double v3[3] = {n3->x - n4->x, n3->y - n4->y, n3->z - n4->z};

    double cp[3];
    crossProd(v2, v3, cp);
    return fabs(dotProd(v1, cp)) / 6.0;
}

h_face* compute_barcoord(hull *hull){

	const double eps=1e-14;

	double l[3],theta[3],c[3],s[3],temp[3],m[3][3];
	double h,dm;

	bool skip;

	for(vector<h_face>::iterator f=hull->faces.begin(); f!=hull->faces.end(); ++f){
		skip=false;
		h=0.0;
		for(int i=0; i<3; ++i){
			for(int j=0; j<3; ++j){
				temp[j]=f->vts[(i+1)%3]->normal[j]-f->vts[(i+2)%3]->normal[j];
			}
			l[i]=length(temp);
			theta[i]=2.0*asin(l[i]/2.0);
		}
		for(int i=0; i<3; ++i)
			h+=theta[i];
		h/=2.0;

		// if pi-h<eps, then N lies on this triangle, so
		// we only need to compute the 2d barycentric coordinates

		if(M_PI-h<eps){
			for(int i=0; i<3; ++i){
				f->vts[i]->w=fabs(sin(theta[i])*f->vts[(i+2)%3]->nl*
						f->vts[(i+1)%3]->nl);
			}
			return &(*f);
		}

		// compute 3d barycentric coordinates
		for(int i=0; i<3; ++i)
			for(int j=0; j<3; ++j)
				m[i][j]=f->vts[j]->normal[i];
		dm=det(m,3);

		for(int i=0; i<3; ++i){

			c[i]=(2.0*sin(h)*sin(h-theta[i]))/(sin(theta[(i+1)%3])*sin(theta[(i+2)%3]))-1;
			s[i]=sqrt(1-c[i]*c[i]);
			s[i]*=(dm>0.0)?1:-1;

			if(fabs(s[i])<eps){
				skip=true;
			}
		}
		if(!skip){
			for(int i=0; i<3; ++i){
				f->vts[i]->w+=fabs((theta[i]-c[(i+1)%3]*theta[(i+2)%3]-c[(i+2)%3]*theta[(i+1)%3])/
				(f->vts[i]->nl*sin(theta[(i+1)%3])*s[(i+2)%3]));
			}
		}
	}
	return NULL;
}

void smooth_dm(vector<node*>& leaves, unsigned short int smooth, double alpha){

	vector<vector<node*> > ptr_hulls;

	int num_leaves=0;
	for(unsigned short int i=0; i<smooth; ++i){

		vector<node*> neigh;

		if(i==0){
			for(vector<node*>::iterator it=leaves.begin(); it!=leaves.end(); ++it){
				node *n=*it;
				if(!n->border()){
					double wtot=0.0;
					double new_dist=0.0;

					neighborhood(neigh,n);
					hull h=hull(n,neigh);
					compute_barcoord(&h);

					vector<node*> ptr_hull;
					for(vector<h_vertex>::iterator v=h.vts.begin(); v!=h.vts.end(); ++v){
						v->n->w=v->w;
						ptr_hull.push_back(v->n);
						new_dist+=v->w*v->n->dist;
						wtot+=v->w;
					}
					n->dist=(1.0-alpha)*n->dist+alpha*(new_dist/wtot);
					ptr_hulls.push_back(ptr_hull);
					num_leaves++;
				}
			}
		}
		else{

			alpha*=1/smooth;
			num_leaves=0;
			for(vector<node*>::iterator it=leaves.begin(); it!=leaves.end(); ++it){
				node *n=*it;
				if(!n->border()){

					double wtot=0.0;
					double new_dist=0.0;

					for(vector<node*>::iterator v=ptr_hulls[num_leaves].begin(); v!=ptr_hulls[num_leaves].end(); ++v){
						new_dist+=(*v)->w*(*v)->dist;
						wtot+=(*v)->w;
					}
					n->dist=(1.0-alpha)*n->dist+alpha*(new_dist/wtot);
					num_leaves++;
				}
			}
		}
	}
}

lap_hull::face::face(node *n1, node *n2, node *n3) {

    double l;
    double e1[3] = {n2->x - n1->x, n2->y - n1->y, n2->z - n1->z};
    double e2[3] = {n3->x - n1->x, n3->y - n1->y, n3->z - n1->z};

    crossProd(e1, e2, normal);
    l = length(normal);
    area = 0.5 * l;
    normal[0] *= -1;
    normal[1] *= -1;
    normal[2] *= -1;
    normalize(normal);
}

lap_hull::vertex::vertex(face *f1, face *f2, face *f3, face *f4) {

    f[0] = f1;
    f[1] = f2;
    f[2] = f3;
    f[3] = f4;
}

void lap_from_grad(vector<node*>& leaves, unsigned short int density_corr) {

    double w, vol;

    for (vector<node*>::iterator it = leaves.begin(); it != leaves.end(); ++it) {
        node *n = *it;
        if (!n->border()) {
            n->lap = 0.0;

            lap_hull::face fs[8] = {lap_hull::face(n->neighbors[0], n->neighbors[1], n->neighbors[4]),
                lap_hull::face(n->neighbors[0], n->neighbors[4], n->neighbors[3]),
                lap_hull::face(n->neighbors[0], n->neighbors[3], n->neighbors[5]),
                lap_hull::face(n->neighbors[0], n->neighbors[5], n->neighbors[1]),
                lap_hull::face(n->neighbors[2], n->neighbors[4], n->neighbors[1]),
                lap_hull::face(n->neighbors[2], n->neighbors[1], n->neighbors[5]),
                lap_hull::face(n->neighbors[2], n->neighbors[5], n->neighbors[3]),
                lap_hull::face(n->neighbors[2], n->neighbors[3], n->neighbors[4])};

            lap_hull::vertex vs[6] = {lap_hull::vertex(&fs[0], &fs[1], &fs[2], &fs[3]), 
                lap_hull::vertex(&fs[0], &fs[3], &fs[4], &fs[5]),
                lap_hull::vertex(&fs[4], &fs[5], &fs[6], &fs[7]), 
                lap_hull::vertex(&fs[1], &fs[2], &fs[6], &fs[7]),
                lap_hull::vertex(&fs[0], &fs[1], &fs[4], &fs[7]), 
                lap_hull::vertex(&fs[2], &fs[3], &fs[5], &fs[6])};

            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 4; ++j) {
                    n->lap += vs[i].f[j]->area / 3.0 * dotProd(vs[i].f[j]->normal, n->neighbors[i]->grad);
                }
            }

            double surf = 0.0;
            for (int i = 0; i < 8; ++i) {
                surf += fs[i].area;
            }
            n->nlap = n->lap / surf;

            vol = (tetra_vol(n->neighbors[0], n->neighbors[1], n->neighbors[3], n->neighbors[4]) +
                    tetra_vol(n->neighbors[1], n->neighbors[2], n->neighbors[3], n->neighbors[4]) +
                    tetra_vol(n->neighbors[0], n->neighbors[1], n->neighbors[3], n->neighbors[5]) +
                    tetra_vol(n->neighbors[1], n->neighbors[2], n->neighbors[3], n->neighbors[5]));
            n->lap /= vol;

            if (!density_corr) {
                n->lap = n->nlap;
            }
        }
    }
}

void density(vector<node*>& input) {

    double dx, dy, dz, delta_x, delta_y, delta_z;
    int xoffset, yoffset, zoffset;
    int nx, ny, nz, nxy, nxz, nyz, nxyz;

    for (vector<node*>::iterator it = input.begin(); it != input.end(); ++it) {

        node *n = *it;

        dx = n->grad[0];
        dy = n->grad[1];
        dz = n->grad[2];

        vector<node*> neighbors;
        fast_neighborhood(neighbors, n);

        xoffset = (dx < 0.0) ? 1 : -1;
        yoffset = (dy < 0.0) ? 1 : -1;
        zoffset = (dz < 0.0) ? 1 : -1;

        nx = 13 + xoffset * 3;
        ny = 13 + yoffset;
        nz = 13 + zoffset * 9;
        nxy = 13 + xoffset * 3 + yoffset;
        nxz = 13 + xoffset * 3 + zoffset * 9;
        nyz = 13 + yoffset + zoffset * 9;
        nxyz = 13 + xoffset * 3 + yoffset + zoffset * 9;

        const double delta_x = fabs(dx);
        const double delta_y = fabs(dy);
        const double delta_z = fabs(dz);

        n->dens = (
                (delta_x * (1 - delta_y)*(1 - delta_z)) * neighbors[nx]->dens +
                (delta_y * (1 - delta_x)*(1 - delta_z)) * neighbors[ny]->dens +
                (delta_z * (1 - delta_x)*(1 - delta_y)) * neighbors[nz]->dens +
                (delta_x * delta_y * (1 - delta_z)) * neighbors[nxy]->dens +
                (delta_x * delta_z * (1 - delta_y)) * neighbors[nxz]->dens +
                (delta_y * delta_z * (1 - delta_x)) * neighbors[nyz]->dens +
                (delta_x * delta_y * delta_z) * neighbors[nxyz]->dens

                - 0.5 * n->size * (

                (1 + (1 - delta_x)*(1 - delta_y)*(1 - delta_z)) * n->lap +
                (delta_x * (1 - delta_y)*(1 - delta_z)) * neighbors[nx]->lap +
                (delta_y * (1 - delta_x)*(1 - delta_z)) * neighbors[ny]->lap +
                (delta_z * (1 - delta_x)*(1 - delta_y)) * neighbors[nz]->lap +
                (delta_x * delta_y * (1 - delta_z)) * neighbors[nxy]->lap +
                (delta_x * delta_z * (1 - delta_y)) * neighbors[nxz]->lap +
                (delta_y * delta_z * (1 - delta_x)) * neighbors[nyz]->lap +
                (delta_x * delta_y * delta_z) * neighbors[nxyz]->lap

                )) / (1 - (1 - delta_x)*(1 - delta_y)*(1 - delta_z));

        if (isnan(n->dens))
            n->dens = 0.0;
    }
}

void fix_div(node *n, double densIntr, double lapIntr, double drl) {

    if (isnan(n->div)) {

        double exp1, exp2, exp3;
        double max = (n->dens - 0.5 * drl) > densIntr ? n->dens - 0.5 * drl : densIntr;

        max = n->dens > max ? n->dens : max;

        exp1 = exp(n->dens - 0.5 * drl - max);
        exp2 = exp(densIntr - max);
        exp3 = exp(n->dens - max);
        n->div = log(exp((drl * exp1 + 0.5 * n->size * (lapIntr * exp2 + n->lap * exp3)) - exp(max)));
    }

    if (!isfinite(n->div))
        n->div = 0.0;
}

void field(vector<node*>& input) {

    double dx, dy, dz, delta_x, delta_y, delta_z;
    int xoffset, yoffset, zoffset;
    int nx, ny, nz, nxy, nxz, nyz, nxyz;

    for (vector<node*>::iterator it = input.begin(); it != input.end(); ++it) {

        node *n = *it;

        dx = n->grad[0];
        dy = n->grad[1];
        dz = n->grad[2];

        vector<node*> neighbors;
        fast_neighborhood(neighbors, n);

        xoffset = (dx < 0.0) ? 1 : -1;
        yoffset = (dy < 0.0) ? 1 : -1;
        zoffset = (dz < 0.0) ? 1 : -1;

        nx = 13 + xoffset * 3;
        ny = 13 + yoffset;
        nz = 13 + zoffset * 9;
        nxy = 13 + xoffset * 3 + yoffset;
        nxz = 13 + xoffset * 3 + zoffset * 9;
        nyz = 13 + yoffset + zoffset * 9;
        nxyz = 13 + xoffset * 3 + yoffset + zoffset * 9;

        const double delta_x = fabs(dx);
        const double delta_y = fabs(dy);
        const double delta_z = fabs(dz);

        const double densIntr = (
                (delta_x * (1 - delta_y)*(1 - delta_z)) * neighbors[nx]->dens +
                (delta_y * (1 - delta_x)*(1 - delta_z)) * neighbors[ny]->dens +
                (delta_z * (1 - delta_x)*(1 - delta_y)) * neighbors[nz]->dens +
                (delta_x * delta_y * (1 - delta_z)) * neighbors[nxy]->dens +
                (delta_x * delta_z * (1 - delta_y)) * neighbors[nxz]->dens +
                (delta_y * delta_z * (1 - delta_x)) * neighbors[nyz]->dens +
                (delta_x * delta_y * delta_z) * neighbors[nxyz]->dens +
                (1 - delta_x)*(1 - delta_y)*(1 - delta_z) * n->dens
                );

        const double drl = n->dens - densIntr;

        const double lapIntr = (
                (delta_x * (1 - delta_y)*(1 - delta_z)) * neighbors[nx]->lap +
                (delta_y * (1 - delta_x)*(1 - delta_z)) * neighbors[ny]->lap +
                (delta_z * (1 - delta_x)*(1 - delta_y)) * neighbors[nz]->lap +
                (delta_x * delta_y * (1 - delta_z)) * neighbors[nxy]->lap +
                (delta_x * delta_z * (1 - delta_y)) * neighbors[nxz]->lap +
                (delta_y * delta_z * (1 - delta_x)) * neighbors[nyz]->lap +
                (delta_x * delta_y * delta_z) * neighbors[nxyz]->lap +
                (1 - delta_x)*(1 - delta_y)*(1 - delta_z) * n->lap
                );

        n->div = drl * exp(n->dens - 0.5 * drl)
                + 0.5 * n->size * (lapIntr * exp(densIntr) + n->lap * exp(n->dens));

        fix_div(n, densIntr, lapIntr, drl);
    }
}

void pop18neigh(vector<node*>& neigh, node* n, int id, int ids[8], int x, int y) {

    neigh[id] = NULL;

    if (*(n->id.end() - 1) == ids[0] || *(n->id.end() - 1) == ids[1]) {

        if (n->neighbors[x] != NULL) {
            if (n->neighbors[x]->level() == n->level() || n->neighbors[x]->level() == n->level() - 1) {
                if (isConsistent(n->neighbors[x]->neighbors[y], n)) {
                    neigh[id] = findClosestNode(n->neighbors[x]->neighbors[y], n);
                }
            } else if (n->neighbors[x]->level() < n->level() - 1) {
                if (isConsistent(n->neighbors[x]->neighbors[y], n) && aligned(n, n->neighbors[x]->level(),
                        ids[0], ids[1], ids[2], ids[3]))
                    neigh[id] = findClosestNode(n->neighbors[x]->neighbors[y], n);
            }
        }
    } else if (*(n->id.end() - 1) == ids[2] || *(n->id.end() - 1) == ids[3]) {
        if (n->neighbors[x]->neighbors[y] != NULL) {
            if (n->neighbors[x]->neighbors[y]->level() == n->level())
                neigh[id] = n->neighbors[x]->neighbors[y];
        }
    } else if (*(n->id.end() - 1) == ids[4] || *(n->id.end() - 1) == ids[5]) {
        if (n->neighbors[y]->neighbors[x] != NULL) {
            if (n->neighbors[y]->neighbors[x]->level() == n->level()) {
                neigh[id] = n->neighbors[y]->neighbors[x];
            }
        }
    } else if (*(n->id.end() - 1) == ids[6] || *(n->id.end() - 1) == ids[7])
        neigh[id] = n->neighbors[x]->neighbors[y];
}

void pop26neigh(vector<node*>& neigh, node *n, int id, int ids[8], int x, int y, int z, int id1, int id2, int id3) {

    neigh[id] = NULL;

    if (*(n->id.end() - 1) == ids[0]) {
        if (neigh[id1] != NULL) {
            if (neigh[id1]->level() == n->level() || neigh[id1]->level() == n->level() - 1) {
                if (isUberConsistent(neigh[id1]->neighbors[x], n, neigh)) {
                    neigh[id] = findClosestNode(neigh[id1]->neighbors[x], n);
                }
            } else if (neigh[id1]->level() < n->level() - 1 && aligned(n, neigh[id1]->level(),
                    ids[0], ids[2], ids[4], ids[6])) {
                if (isUberConsistent(neigh[id1]->neighbors[x], n, neigh)) {
                    neigh[id] = findClosestNode(neigh[id1]->neighbors[x], n);
                }
            }
        } else if (neigh[id2] != NULL) {
            if (neigh[id2]->level() == n->level() || neigh[id2]->level() == n->level() - 1) {
                if (isUberConsistent(neigh[id2]->neighbors[z], n, neigh)) {
                    neigh[id] = findClosestNode(neigh[id2]->neighbors[z], n);
                }
            } else if (neigh[id2]->level() < n->level() - 1 && aligned(n, neigh[id2]->level(),
                    ids[0], ids[1], ids[4], ids[5])) {
                if (isUberConsistent(neigh[id2]->neighbors[z], n, neigh)) {
                    neigh[id] = findClosestNode(neigh[id2]->neighbors[z], n);
                }
            }
        } else if (neigh[id3] != NULL) {
            // if neigh[id3] is certainly aligned with n and its y-neighbor is not
            // one of n's 18-neighbors (might happen)
            if (neigh[id3]->level() == n->level() || neigh[id3]->level() == n->level() - 1) {
                if (isUberConsistent(neigh[id3]->neighbors[y], n, neigh)) {
                    // look for the closest node in the tree rootedin neigh[id3]->neighbors[y]
                    neigh[id] = findClosestNode(neigh[id3]->neighbors[y], n);
                }
            }// else check if they're aligned, otherwise every other computation is useless
            else if (neigh[id3]->level() < n->level() - 1 && aligned(n, neigh[id3]->level(),
                    ids[0], ids[1], ids[2], ids[3])) {
                if (isUberConsistent(neigh[id3]->neighbors[y], n, neigh)) {
                    neigh[id] = findClosestNode(neigh[id3]->neighbors[y], n);
                }
            }
        }
    } else if (*(n->id.end() - 1) == ids[1]) {
        if (neigh[id1] != NULL) {
            if (neigh[id1]->level() == n->level()) {
                neigh[id] = neigh[id1]->neighbors[x];
            }
        } else if (neigh[id2] != NULL) {
            if (neigh[id2]->neighbors[z]->level() == n->level())
                neigh[id] = neigh[id2]->neighbors[z];
        } else if (neigh[id3] != NULL) {
            if (isUberConsistent(neigh[id3]->neighbors[y], n, neigh))
                neigh[id] = neigh[id3]->neighbors[y];
        }
    } else if (*(n->id.end() - 1) == ids[2]) {
        if (neigh[id2] != NULL) {
            if (neigh[id2]->level() == n->level()) {
                neigh[id] = neigh[id2]->neighbors[z];
            }
        } else if (neigh[id1] != NULL) {
            if (neigh[id1]->neighbors[x]->level() == n->level())
                neigh[id] = neigh[id1]->neighbors[x];
        } else if (neigh[id3] != NULL) {
            if (isUberConsistent(neigh[id3]->neighbors[y], n, neigh))
                neigh[id] = neigh[id3]->neighbors[y];
        }
    } else if (*(n->id.end() - 1) == ids[3]) {
        if (neigh[id3]->neighbors[y]->level() == n->level())
            neigh[id] = neigh[id3]->neighbors[y];
    } else if (*(n->id.end() - 1) == ids[4]) {
        if (neigh[id1] != NULL) {
            if (neigh[id1]->neighbors[x]->level() == n->level()) {
                neigh[id] = neigh[id1]->neighbors[x];
            }
        } else if (neigh[id2] != NULL) {
            if (neigh[id2]->neighbors[z]->level() == n->level()) {
                neigh[id] = neigh[id2]->neighbors[z];
            }
        } else if (neigh[id3] != NULL) {
            if (neigh[id3]->neighbors[y]->level() == n->level()) {
                neigh[id] = neigh[id3]->neighbors[y];
            }
        }
    } else if (*(n->id.end() - 1) == ids[5]) {
        if (neigh[id2]->neighbors[z]->level() == n->level()) {
            neigh[id] = neigh[id2]->neighbors[z];
        }
    } else if (*(n->id.end() - 1) == ids[6]) {
        if (neigh[id1]->neighbors[x]->level() == n->level()) {
            neigh[id] = neigh[id1]->neighbors[x];
        }
    } else if (*(n->id.end() - 1) == ids[7]) {
        neigh[id] = neigh[id1]->neighbors[x];
    }
}

void neighborhood(vector<node*>& neigh, node* n) {

    int ids[8];
    neigh.resize(27);
    // the node itself
    neigh[13] = n;
    // its 6-neighbors
    neigh[12] = n->neighbors[0];
    neigh[16] = n->neighbors[1];
    neigh[14] = n->neighbors[2];
    neigh[10] = n->neighbors[3];
    neigh[4] = n->neighbors[4];
    neigh[22] = n->neighbors[5];

    ids[0] = 2;
    ids[1] = 0;
    ids[2] = 3;
    ids[3] = 1;
    ids[4] = 6;
    ids[5] = 4;
    ids[6] = 7;
    ids[7] = 5;
    pop18neigh(neigh, n, 1, ids, 3, 4);

    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
    ids[4] = 4;
    ids[5] = 5;
    ids[6] = 6;
    ids[7] = 7;
    pop18neigh(neigh, n, 3, ids, 0, 4);

    ids[0] = 3;
    ids[1] = 2;
    ids[2] = 1;
    ids[3] = 0;
    ids[4] = 7;
    ids[5] = 6;
    ids[6] = 5;
    ids[7] = 4;
    pop18neigh(neigh, n, 5, ids, 2, 4);

    ids[0] = 1;
    ids[1] = 3;
    ids[2] = 0;
    ids[3] = 2;
    ids[4] = 5;
    ids[5] = 7;
    ids[6] = 4;
    ids[7] = 6;
    pop18neigh(neigh, n, 7, ids, 1, 4);

    /****/

    ids[0] = 4;
    ids[1] = 0;
    ids[2] = 6;
    ids[3] = 2;
    ids[4] = 5;
    ids[5] = 1;
    ids[6] = 7;
    ids[7] = 3;
    pop18neigh(neigh, n, 9, ids, 0, 3);

    ids[0] = 2;
    ids[1] = 6;
    ids[2] = 0;
    ids[3] = 4;
    ids[4] = 3;
    ids[5] = 7;
    ids[6] = 1;
    ids[7] = 5;
    pop18neigh(neigh, n, 11, ids, 2, 3);

    ids[0] = 1;
    ids[1] = 5;
    ids[2] = 3;
    ids[3] = 7;
    ids[4] = 0;
    ids[5] = 4;
    ids[6] = 2;
    ids[7] = 6;
    pop18neigh(neigh, n, 15, ids, 0, 1);

    ids[0] = 7;
    ids[1] = 3;
    ids[2] = 5;
    ids[3] = 1;
    ids[4] = 6;
    ids[5] = 2;
    ids[6] = 4;
    ids[7] = 0;
    pop18neigh(neigh, n, 17, ids, 2, 1);

    /****/

    ids[0] = 4;
    ids[1] = 6;
    ids[2] = 5;
    ids[3] = 7;
    ids[4] = 0;
    ids[5] = 2;
    ids[6] = 1;
    ids[7] = 3;
    pop18neigh(neigh, n, 19, ids, 3, 5);

    ids[0] = 5;
    ids[1] = 4;
    ids[2] = 7;
    ids[3] = 6;
    ids[4] = 1;
    ids[5] = 0;
    ids[6] = 3;
    ids[7] = 2;
    pop18neigh(neigh, n, 21, ids, 0, 5);

    ids[0] = 7;
    ids[1] = 5;
    ids[2] = 6;
    ids[3] = 4;
    ids[4] = 3;
    ids[5] = 1;
    ids[6] = 2;
    ids[7] = 0;
    pop18neigh(neigh, n, 25, ids, 1, 5);

    ids[0] = 6;
    ids[1] = 7;
    ids[2] = 4;
    ids[3] = 5;
    ids[4] = 2;
    ids[5] = 3;
    ids[6] = 0;
    ids[7] = 1;
    pop18neigh(neigh, n, 23, ids, 2, 5);

    /****/

    //pop26neigh(neigh,n,n_neigh,ids,x,y,z,id1,id2,id3)
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
    ids[4] = 4;
    ids[5] = 5;
    ids[6] = 6;
    ids[7] = 7;
    pop26neigh(neigh, n, 0, ids, 3, 4, 0, 3, 1, 9);

    ids[0] = 1;
    ids[1] = 3;
    ids[2] = 0;
    ids[3] = 2;
    ids[4] = 5;
    ids[5] = 7;
    ids[6] = 4;
    ids[7] = 6;
    pop26neigh(neigh, n, 6, ids, 0, 4, 1, 7, 3, 15);

    ids[0] = 3;
    ids[1] = 2;
    ids[2] = 1;
    ids[3] = 0;
    ids[4] = 7;
    ids[5] = 6;
    ids[6] = 5;
    ids[7] = 4;
    pop26neigh(neigh, n, 8, ids, 1, 4, 2, 5, 7, 17);

    ids[0] = 2;
    ids[1] = 0;
    ids[2] = 3;
    ids[3] = 1;
    ids[4] = 6;
    ids[5] = 4;
    ids[6] = 7;
    ids[7] = 5;
    pop26neigh(neigh, n, 2, ids, 2, 4, 3, 1, 5, 11);

    /****/

    ids[0] = 5;
    ids[1] = 4;
    ids[2] = 7;
    ids[3] = 6;
    ids[4] = 1;
    ids[5] = 0;
    ids[6] = 3;
    ids[7] = 2;
    pop26neigh(neigh, n, 24, ids, 1, 5, 0, 21, 25, 15);

    ids[0] = 7;
    ids[1] = 5;
    ids[2] = 6;
    ids[3] = 4;
    ids[4] = 3;
    ids[5] = 1;
    ids[6] = 2;
    ids[7] = 0;
    pop26neigh(neigh, n, 26, ids, 2, 5, 1, 25, 23, 17);

    ids[0] = 6;
    ids[1] = 7;
    ids[2] = 4;
    ids[3] = 5;
    ids[4] = 2;
    ids[5] = 3;
    ids[6] = 0;
    ids[7] = 1;
    pop26neigh(neigh, n, 20, ids, 3, 5, 2, 23, 19, 11);

    ids[0] = 4;
    ids[1] = 6;
    ids[2] = 5;
    ids[3] = 7;
    ids[4] = 0;
    ids[5] = 2;
    ids[6] = 1;
    ids[7] = 3;
    pop26neigh(neigh, n, 18, ids, 0, 5, 3, 19, 21, 9);
}

void fast_neighborhood(vector<node*>& neigh, node* n) {
    neigh.resize(27);
    neigh[0] = n->neighbors[4]->neighbors[3]->neighbors[0];
    neigh[1] = n->neighbors[4]->neighbors[3];
    neigh[2] = n->neighbors[4]->neighbors[3]->neighbors[2];
    neigh[3] = n->neighbors[4]->neighbors[0];
    neigh[4] = n->neighbors[4];
    neigh[5] = n->neighbors[4]->neighbors[2];
    neigh[6] = n->neighbors[4]->neighbors[0]->neighbors[1];
    neigh[7] = n->neighbors[4]->neighbors[1];
    neigh[8] = n->neighbors[4]->neighbors[1]->neighbors[2];
    neigh[9] = n->neighbors[3]->neighbors[0];
    neigh[10] = n->neighbors[3];
    neigh[11] = n->neighbors[3]->neighbors[2];
    neigh[12] = n->neighbors[0];
    neigh[13] = n;
    neigh[14] = n->neighbors[2];
    neigh[15] = n->neighbors[1]->neighbors[0];
    neigh[16] = n->neighbors[1];
    neigh[17] = n->neighbors[1]->neighbors[2];
    neigh[18] = n->neighbors[5]->neighbors[3]->neighbors[0];
    neigh[19] = n->neighbors[5]->neighbors[3];
    neigh[20] = n->neighbors[5]->neighbors[3]->neighbors[2];
    neigh[21] = n->neighbors[5]->neighbors[0];
    neigh[22] = n->neighbors[5];
    neigh[23] = n->neighbors[5]->neighbors[2];
    neigh[24] = n->neighbors[5]->neighbors[1]->neighbors[0];
    neigh[25] = n->neighbors[5]->neighbors[1];
    neigh[26] = n->neighbors[5]->neighbors[1]->neighbors[2];
}

/*
 * Def:
 *
 * 1. C* : the number of 26-connected components 26-adjacent to x in B ^ N*26
 * 2. C-: the number of 6-connected components 6-adjacent to x in W ^ N18
 *
 * Thm:
 *
 * If C*(x)=1 and C-(x)=1 then x is a simple point
 */

void virtual_neighborhood(node *n, vector<node*>& neigh, bool bmask[3][3][3]) {

    for (int i = 0; i < neigh.size(); ++i) {
        bmask[i % 3][(i / 3) % 3][i / 9] = true;
        if (neigh[i] == NULL)
            bmask[i % 3][(i / 3) % 3][i / 9] = false;
        else if (neigh[i]->level() > n->level() || !neigh[i]->skel)
            bmask[i % 3][(i / 3) % 3][i / 9] = false;
    }

}

void recursive_check(int id, bool bmask[3][3][3], bool visited[3][3][3], int type) {

    int x = id % 3;
    int y = (id / 3) % 3;
    int z = id / 9;

    int a, b, c, i, j, k;

    if (type == CSTAR) {
        for (i = -1; i <= 1; i++) {
            for (j = -1; j <= 1; j++) {
                for (k = -1; k <= 1; k++) {
                    a = x + i;
                    b = y + j;
                    c = z + k;
                    if (a <= 2 && a >= 0 && b <= 2 && b >= 0 && c <= 2 && c >= 0) {
                        if (!visited[a][b][c] && bmask[a][b][c]) {
                            visited[a][b][c] = true;
                            recursive_check(a + b * 3 + c * 9, bmask, visited, CSTAR);
                        }
                    }
                }
            }
        }
    }

    if (type == CMINUS) {
        if (id == 0 || id == 2 || id == 6 || id == 8 || id == 18 || id == 20 || id == 24 || id == 26)
            return;

        for (i = -1; i <= 1; i++) {
            for (j = -1; j <= 1; j++) {
                for (k = -1; k <= 1; k++) {
                    a = x + i;
                    b = y + j;
                    c = z + k;
                    if ((i == 0 && j == 0) || (i == 0 && k == 0) || (j == 0 && k == 0)) {
                        if (a <= 2 && a >= 0 && b <= 2 && b >= 0 && c <= 2 && c >= 0) {
                            if (!visited[a][b][c] && !bmask[a][b][c]) {
                                visited[a][b][c] = true;
                                recursive_check(a + b * 3 + c * 9, bmask, visited, CMINUS);
                            }
                        }
                    }
                }
            }
        }

    }
}

int CStar(node* n) {

    int a, b, c;

    vector<node*> neighbors;
    neighborhood(neighbors, n);

    bool bmask[3][3][3] = {
        {
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}
        },
        {
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}
        },
        {
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}
        }
    };

    virtual_neighborhood(n, neighbors, bmask);

    bool visited[3][3][3];
    for (int i = 0; i < 27; i++)
        visited[i % 3][(i / 3) % 3][i / 9] = false;
    visited[1][1][1] = true;

    int num = 0;
    for (int i = 0; i < 27; i++) {
        a = i % 3;
        b = (i / 3) % 3;
        c = i / 9;
        if (bmask[a][b][i / 9] && !visited[a][b][c]) {
            visited[a][b][c] = true;
            num++;
            recursive_check(i, bmask, visited, CSTAR);
        }
    }
    return num;
}

int CMinus(node *n) {

    int a, b, c;

    vector<node*> neighbors;
    neighborhood(neighbors, n);

    bool bmask[3][3][3] = {
        {
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}
        },
        {
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}
        },
        {
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}
        }
    };

    virtual_neighborhood(n, neighbors, bmask);

    bool visited[3][3][3];
    for (int i = 0; i < 27; i++)
        visited[i % 3][(i / 3) % 3][i / 9] = false;
    visited[1][1][1] = true;

    int num = 0;
    for (int i = 0; i < 27; i++) {
        a = i % 3;
        b = (i / 3) % 3;
        c = i / 9;
        if ((i == 4 || i == 10 || i == 12 || i == 14 || i == 16 || i == 22) &&
                !bmask[a][b][c] && !visited[a][b][c]) {
            visited[a][b][c] = true;
            num++;
            recursive_check(i, bmask, visited, CMINUS);
        }
    }
    return num;
}

void dilation(vector<node*>& leaves, vector<node*>& input, vector<node*>& dleaves, octree *oct) {

    vector<node*> neigh;
    dleaves.clear();

    // explode adjacent nodes to start the dilation
    for (vector<node*>::iterator i = leaves.begin(); i != leaves.end(); ++i) {
        node *n = *i;
        neighborhood(neigh, n);
        for (vector<node*>::iterator j = neigh.begin(); j != neigh.end(); ++j) {
            node *m = *j;
            if (m != NULL) {
                if (m->level() < n->level()) {
                    m->recAdd(n->level() - m->level());
                    oct->getLeaves(dleaves, m);
                }
            }
        }
    }

    // mark the voxel adjacent to the skel as candidates for skeleton selection
    for (vector<node*>::iterator i = input.begin(); i != input.end(); ++i) {
        node *n = *i;

        // we only add to the current skeleton those voxels that are
        // adjacent to the input (1 voxel dilation)
        fast_neighborhood(neigh, n);
        for (vector<node*>::iterator j = neigh.begin(); j != neigh.end(); ++j) {
            node *m = *j;
            if (m != NULL) {
                // if we didn't check this voxel
                if (!m->skel && m->dist == inf && !m->candidate) {
                    m->candidate = true;
                }
            }
        }
    }
}

/*
 * skeleton dilation: dilate the skeleton where no erosion was performed wrt
 * the original object, since this might suggest that we've been too
 *  aggressive in the past erosions. the new voxels are pushed into new_skel.
 */
bool skel_dilation(vector<node*>& skel, vector<node*>& new_skel, octree *oct) {

    vector<node*> new_leaves, neigh, neighm;
    bool dilation = false;

    for (vector<node*>::iterator i = skel.begin(); i != skel.end(); ++i) {
        node *n = *i;
        fast_neighborhood(neigh, n);
        for (vector<node*>::iterator j = neigh.begin(); j != neigh.end(); ++j) {
            node *m = *j;
            if (m != NULL) {
                // if this voxel belonged to the previous border dilation
                // then we might need to explode its neighbors
                if (m->dborder) {
                    neighborhood(neighm, m);
                    for (vector<node*>::iterator k = neighm.begin(); k != neighm.end(); ++k) {
                        node *l = *k;
                        if (l != NULL) {
                            if (l->level() < m->level() && !l->border()) {
                                l->recAdd(m->level() - l->level());
                                oct->getLeaves(new_leaves, l);
                            }
                        }
                    }
                }
            }
        }

        // we only add to the current skeleton those voxels that are
        // adjacent to the input (1 voxel dilation)
        // and only if the voxel was in the previous dilation border

        for (vector<node*>::iterator j = neigh.begin(); j != neigh.end(); ++j) {
            node *m = *j;
            if (m != NULL) {
                // if we didn't check this voxel and if it is a candidate
                if (!m->skel && m->dborder && !m->candidate) {
                    m->dborder = false;
                    m->candidate = true;
                    dilation = true;
                    new_leaves.push_back(m);
                }
            }
        }
    }

    for (vector<node*>::iterator i = new_leaves.begin(); i != new_leaves.end(); ++i) {
        node *n = *i;
        new_skel.push_back(n);
    }

    return dilation;
}

void set_dborder(vector<node*>& skel_input, vector<node*>& dborder) {

    vector<node*> neigh;
    dborder.clear();
    // mark the border of the dilation for further reference
    for (vector<node*>::iterator i = skel_input.begin(); i != skel_input.end(); ++i) {
        node *n = *i;
        if (n->skel) {
            neighborhood(neigh, n);
            for (vector<node*>::iterator j = neigh.begin(); j != neigh.end(); ++j) {
                node *m = *j;
                if (m != NULL) {
                    if (!m->skel && !m->dborder)
                        m->dborder = true;
                    dborder.push_back(m);
                }
            }
        }
    }
}

void reset_dborder(vector<node*>& dborder) {

    bool skip;
    vector<node*> neigh;
    for (vector<node*>::iterator it = dborder.begin(); it != dborder.end(); it++) {
        node *n = *it;
        skip = false;
        neighborhood(neigh, n);
        for (vector<node*>::iterator jt = neigh.begin(); jt != neigh.end(); jt++) {
            if (*jt != NULL)
                if ((*jt)->skel)
                    skip = true;
        }
        if (!skip)
            n->dborder = false;
    }
}

bool endpoint(node* n) {

    return
    (!n->neighbors[0]->skel && !n->neighbors[1]->skel && !n->neighbors[3]->skel) ||
    (!n->neighbors[0]->skel && !n->neighbors[4]->skel && !n->neighbors[5]->skel) ||

    (!n->neighbors[1]->skel && !n->neighbors[0]->skel && !n->neighbors[2]->skel) ||
    (!n->neighbors[1]->skel && !n->neighbors[4]->skel && !n->neighbors[5]->skel) ||

    (!n->neighbors[2]->skel && !n->neighbors[1]->skel && !n->neighbors[3]->skel) ||
    (!n->neighbors[2]->skel && !n->neighbors[4]->skel && !n->neighbors[5]->skel) ||

    (!n->neighbors[3]->skel && !n->neighbors[0]->skel && !n->neighbors[2]->skel) ||
    (!n->neighbors[3]->skel && !n->neighbors[4]->skel && !n->neighbors[5]->skel) ||

    (!n->neighbors[4]->skel && !n->neighbors[0]->skel && !n->neighbors[2]->skel) ||
    (!n->neighbors[4]->skel && !n->neighbors[1]->skel && !n->neighbors[3]->skel) ||

    (!n->neighbors[5]->skel && !n->neighbors[0]->skel && !n->neighbors[2]->skel) ||
    (!n->neighbors[5]->skel && !n->neighbors[1]->skel && !n->neighbors[3]->skel);
}

vector<node*> skeletonize(vector<node*>& input, double treshold, unsigned short int density_corr) {

    vector<node*> skel1, skel2;

    int cs, cm;
    double importance;

    for (vector<node*>::iterator i = input.begin(); i != input.end(); ++i) {
 
        cs = CStar(*i);
        cm = CMinus(*i);

        if (density_corr)
            importance = (*i)->div;
        else
            importance = (*i)->lap;

        if (cs == 0) {
            cerr << "\nWARNING #1: isolated voxel => ";
            cerr << **i << endl;
            (*i)->skel = false;
        } else if (cs == 1 && cm == 1 && importance > treshold) {
            (*i)->skel = false;
        } else {
            skel1.push_back(*i);
        }
    }

    if(density_corr)
        sort(skel1.begin(),skel1.end(),divOrd);
    else
        sort(skel1.begin(),skel1.end(),lapOrd);

    for (vector<node*>::iterator i = skel1.begin(); i != skel1.end(); ++i) {

        cs = CStar(*i);
        cm = CMinus(*i);

        if (density_corr)
            importance = (*i)->div;
        else
            importance = (*i)->lap;

        if (cs == 0) {
            cerr << "\nWARNING #1: isolated voxel => ";
            cerr << **i << endl;
            (*i)->skel = false;
        } else if (cs == 1 && cm == 1 && importance > treshold) {
            (*i)->skel = false;
        } else {
            skel2.push_back(*i);
        }
    }

    return skel2;
}

vector<node*> ensure_thinness(vector<node*>& input) {

    vector<node*> skel1;

    int cs, cm;

    for (vector<node*>::iterator i = input.begin(); i != input.end(); ++i) {

        cs = CStar(*i);
        cm = CMinus(*i);

        if (cs == 0) {
            cerr << "\nWARNING #2: isolated voxel => ";
            cerr << **i << endl;
            (*i)->skel = false;
        } else if (cs == 1 && cm == 1 && !endpoint(*i)) {
            (*i)->skel = false;
        } else {
            skel1.push_back(*i);
        }
    }

    return skel1;
}

void find_cc(node *n, vector<node*>& cc) {
    vector<node*> neigh;
    neighborhood(neigh, n);
    for (vector<node*>::iterator it = neigh.begin(); it != neigh.end(); ++it) {
        if (!(*it)->flag && (*it)->skel) {
            (*it)->flag = true;
            cc.push_back(*it);
            find_cc(*it, cc);
        }
    }
}

void keep_largest_cc(vector<node*>& interior) {

    vector<node*> lcc;
    for (vector<node*>::iterator it = interior.begin(); it != interior.end(); ++it) {
        node *n = *it;
        vector<node*> cc;
        if (!n->flag) {
            n->flag = true;
            cc.push_back(n);
            find_cc(n, cc);
            if (cc.size() > lcc.size())
                lcc = cc;
        }
    }
    for (vector<node*>::iterator it = interior.begin(); it != interior.end(); ++it) {
        node *n = *it;
        n->flag = false;
        n->skel = false;
    }

    interior.clear();
    for (vector<node*>::iterator it = lcc.begin(); it != lcc.end(); ++it) {
        node *n = *it;
        n->skel = true;
        interior.push_back(n);
    }
}