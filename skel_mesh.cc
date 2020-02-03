#include "skel_mesh.h"

vector<node*> refine_skeleton(vector<node*>& leaves, vector<node*>& skeleton, Tree& mesh){
    
    int s;
    vector<node*> neigh;
    double xoffset, yoffset, zoffset;
    
    // the gradient might have changed after the dilation, so we need
    // to recompute it here. since we have a complete neighborhood, we can
    // use common finite difference methods
    for(vector<node*>::iterator it=leaves.begin(); it!=leaves.end(); ++it){
        node *n=*it;             
		n->grad[0]=n->x-n->np[0];
        n->grad[1]=n->y-n->np[1];
        n->grad[2]=n->z-n->np[2];
        normalize(n->grad);
    }
    
    // re-align those voxels that were mistaken as being non-skeletal
    // warning: make sure not to change the topology of the skeleton!
   
    bool aligning=true;
    while(aligning){

        aligning=false;
        for(vector<node*>::iterator it=skeleton.begin(); it!=skeleton.end(); ++it){

            node *n=*it;       
            node *m;

            fast_neighborhood(neigh,n);

            for(vector<node*>::iterator it=neigh.begin(); it!=neigh.end(); ++it){
                m=*it;
                if(!m->skel && n->np==m->np && m->dist > n->dist && m->level()==n->level()){
                    m->skel=true;
                    if(CMinus(m)==1 && CStar(m)==1 && CMinus(n)==1 && CStar(n)==1){
                        n->skel=false;
                        aligning=true;
                    }
                    else
                        m->skel=false;                    
                }
            }
        }

        skeleton.clear();
        for(vector<node*>::iterator it=leaves.begin(); it!=leaves.end(); ++it){
            if((*it)->skel)
                skeleton.push_back(*it);
        }
    }

    // now move the voxels centers according to the local value of the gradient
    
    for(vector<node*>::iterator it=skeleton.begin(); it!=skeleton.end(); ++it){
        
        node *n=*it;
        node *m;

        fast_neighborhood(neigh,n);
        
        if(fabs(n->grad[0])<1e-5) xoffset = 0;
        else xoffset = (n->grad[0] < 0.0) ? -1 : 1;
        
        if(fabs(n->grad[1])<1e-5) yoffset = 0;
        else yoffset = (n->grad[1] < 0.0) ? -1 : 1;
        
        if(fabs(n->grad[2])<1e-5) zoffset = 0;
        else zoffset = (n->grad[2] < 0.0) ? -1 : 1;

        vector<node*> dir_neigh;
        
        if(xoffset!=0)
            dir_neigh.push_back(neigh[13 + xoffset * 3]);

        if(yoffset!=0)
            dir_neigh.push_back(neigh[13 + yoffset]);

        if(zoffset!=0)
            dir_neigh.push_back(neigh[13 + zoffset * 9]);
        
        if(xoffset!=0 && yoffset!=0)
            dir_neigh.push_back(neigh[13 + xoffset * 3 + yoffset]);
        
        if(xoffset!=0 && zoffset!=0)
            dir_neigh.push_back(neigh[13 + xoffset * 3 + zoffset * 9]);
        
        if(yoffset!=0 && zoffset!=0)
            dir_neigh.push_back(neigh[13 + yoffset + zoffset * 9]);
        
        if(xoffset!=0 && yoffset!=0 && zoffset!=0)
            dir_neigh.push_back(neigh[13 + xoffset * 3 + yoffset + zoffset * 9]);
        
        double cp[3];
        unsigned char counter=0; 
        vector<node*> endpoints;
        vector<vector<double> > points;
        vector<double> weights;
        
        //cerr << "===================\n";
        for(vector<node*>::iterator it=dir_neigh.begin(); it!=dir_neigh.end(); ++it){
            m=*it;
            if(dotProd(n->grad,m->grad)<=0){
                counter++;
                double nm[3]={n->x-m->x,n->y-m->y,n->z-m->z};
                
                normalize(nm);
                
                //crossProd(n->grad, nm, cp);
                endpoints.push_back(m);
                weights.push_back(1-dotProd(n->grad,nm));
                //weights.push_back(1-dotProd(cp,cp)/dotProd(n->grad,n->grad));                
            }
        }
        
        // if we couldn't improve the voxel location in any direction, just
        // skip the optimization of this voxel
        n->flag=false;
        if(counter<1){
            n->flag=true;
            continue;
        }
        
        for(vector<node*>::iterator it=endpoints.begin(); it!=endpoints.end(); ++it){
            
            node *nn=*it;
            gsl_matrix* A = gsl_matrix_calloc(5, 5);
            gsl_vector* b = gsl_vector_calloc(5);
            gsl_vector* x = gsl_vector_calloc(5);
            gsl_permutation* p = gsl_permutation_alloc(5);

            gsl_matrix_set(A, 0, 0, 0);
            gsl_matrix_set(A, 0, 1, 0);
            gsl_matrix_set(A, 0, 2, 0);
            gsl_matrix_set(A, 0, 3, 1);
            gsl_matrix_set(A, 0, 4, 1);

            gsl_matrix_set(A, 1, 0, 1);
            gsl_matrix_set(A, 1, 1, 0);
            gsl_matrix_set(A, 1, 2, 0);
            gsl_matrix_set(A, 1, 3, -n->x);
            gsl_matrix_set(A, 1, 4, -nn->x);

            gsl_matrix_set(A, 2, 0, 0);
            gsl_matrix_set(A, 2, 1, 1);
            gsl_matrix_set(A, 2, 2, 0);
            gsl_matrix_set(A, 2, 3, -n->y);
            gsl_matrix_set(A, 2, 4, -nn->y);

            gsl_matrix_set(A, 3, 0, 0);
            gsl_matrix_set(A, 3, 1, 0);
            gsl_matrix_set(A, 3, 2, 1);
            gsl_matrix_set(A, 3, 3, -n->z);
            gsl_matrix_set(A, 3, 4, -nn->z);

            gsl_matrix_set(A, 4, 0, 2*(nn->np[0]-n->np[0]));
            gsl_matrix_set(A, 4, 1, 2*(nn->np[1]-n->np[1]));
            gsl_matrix_set(A, 4, 2, 2*(nn->np[2]-n->np[2]));
            gsl_matrix_set(A, 4, 3, 0);
            gsl_matrix_set(A, 4, 4, 0);     

            gsl_vector_set(b, 0, 1);
            gsl_vector_set(b, 1, 0);
            gsl_vector_set(b, 2, 0);
            gsl_vector_set(b, 3, 0);
            gsl_vector_set(b, 4, nn->np[0]*nn->np[0]+nn->np[1]*nn->np[1]+
                    nn->np[2]*nn->np[2]-(n->np[0]*n->np[0]+n->np[1]*n->np[1]+
                    n->np[2]*n->np[2]));

            gsl_linalg_LU_decomp(A, p, &s);
            gsl_linalg_LU_solve(A, p, b, x);
            
            vector<double> pt;
            pt.push_back(gsl_vector_get(x,0));
            pt.push_back(gsl_vector_get(x,1));
            pt.push_back(gsl_vector_get(x,2));        
            points.push_back(pt);
            
            gsl_matrix_free(A);
            gsl_vector_free(b);
            gsl_vector_free(x);
            gsl_permutation_free(p);
        }
        
        double tot_weight=0.0;
        
        n->x=0;
        n->y=0;
        n->z=0;
                
        for(int i=0; i<points.size(); ++i){
            n->x+=1.0/weights[i]*points[i][0];
            n->y+=1.0/weights[i]*points[i][1];
            n->z+=1.0/weights[i]*points[i][2];          
            tot_weight+=1.0/weights[i];
        }
        
        n->x/=tot_weight;
        n->y/=tot_weight;
        n->z/=tot_weight;

        Point point_query(n->x,n->y,n->z);
        n->np = mesh.closest_point(point_query);
        n->dist = sqrt(CGAL::squared_distance(n->np,point_query));
        
        for(int i=0; i<3; ++i)
			n->grad[i]=point_query[i]-n->np[i];
        normalize(n->grad);
    }
    
    return skeleton;
}