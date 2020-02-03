#include<vector>
#include<queue>
#include<math.h>
#include<fstream>
#include<sstream>
#include"distTrans.h"
#include"misc.h"

using namespace std;

double OFFReader(char* inFile, Tree& mesh, double res){

	static ifstream in;
	list<Triangle> trs;

	// open model file
	in.open(inFile, ios::in);
	if (!in) {
		cerr << "Cannot open data file\n";
		exit(1);
	}

	string line;

	getline(in, line); // throw the 1st line away, it always says "OFF"

	// Read until we get to the first line that is not a comment or a blank line.
	// This line states the number of vertices, faces, and edges.
	while(getline(in, line)){

		if(line.size() == 0)
		  continue;
		if(line[0] == '#')
		  continue;

		// if we get to here, this is the info line
		break;
	}

	// At this point, the line we are interested in is stored in 'line'
	// We are only interested in vertices and faces.

	stringstream ss;
	ss << line;
	unsigned int nvert, nface;
	ss >> nvert >> nface;

	// read the vertices
	unsigned int vcount = 0;
	double pts[nvert][3];

	while(getline(in, line) && vcount < nvert){

		if(line.size() == 0)
			continue;

		stringstream ssvert;
		ssvert << line;
		double x,y,z;
		ssvert >> x >> y >> z;

		pts[vcount][0]=x;
		pts[vcount][1]=y;
		pts[vcount][2]=z;

		vcount++;
    }

	double min_x=inf, min_y=inf, min_z=inf, max_x=-inf, max_y=-inf, max_z=-inf, xres, yres, zres;
	for(unsigned int i=0; i<nvert; ++i){
		if(pts[i][0]<min_x)
			min_x = pts[i][0];
		if(pts[i][1]<min_y)
			min_y = pts[i][1];
		if(pts[i][2]<min_z)
			min_z = pts[i][2];
		if(pts[i][0]>max_x)
			max_x = pts[i][0];
		if(pts[i][1]>max_y)
			max_y = pts[i][1];
		if(pts[i][2]>max_z)
			max_z = pts[i][2];
	}
    
    // compute the length, width and depth of the mesh bounding box
	xres=(max_x-min_x);//*(1.0+0.1);//*rand()/(double)RAND_MAX);
	yres=(max_y-min_y);//*(1.0+0.1);//*rand()/(double)RAND_MAX);
	zres=(max_z-min_z);//*(1.0+0.1);//*rand()/(double)RAND_MAX);

    double max=res;
    res=res-3;
    
	if((xres>=yres) && (xres>=zres))
		res/=xres;
    
	else if((yres>=xres) && (yres>=zres))
		res/=yres;

	else
		res/=zres;

    // now, along the longest dimension, pts go from 0 to res
    // add a padding now..
	for(unsigned int i=0; i<vcount; ++i){
		pts[i][0]=(pts[i][0]-min_x)*res+3;
		pts[i][1]=(pts[i][1]-min_y)*res+3;
		pts[i][2]=(pts[i][2]-min_z)*res+3;
	}

	max_x=-inf;
	max_y=-inf;
	max_z=-inf;

	for(unsigned int i=0; i<nvert; ++i){
		if(pts[i][0]>max_x)
			max_x = pts[i][0];
		if(pts[i][1]>max_y)
			max_y = pts[i][1];
		if(pts[i][2]>max_z)
			max_z = pts[i][2];
	}

    // the discretized space has to be cubic
	double xoff, yoff, zoff;
	for(unsigned int i=0; i<nvert; ++i){
		if(max==max_x){
			yoff=(max-max_y)/2.0;
			zoff=(max-max_z)/2.0;
			pts[i][1]+=yoff;
			pts[i][2]+=zoff;
		}
		else if(max==max_y){
			xoff=(max-max_x)/2.0;
			zoff=(max-max_z)/2.0;
			pts[i][0]+=xoff;
			pts[i][2]+=zoff;
		}
		else{
			xoff=(max-max_x)/2.0;
			yoff=(max-max_y)/2.0;
			pts[i][0]+=xoff;
			pts[i][1]+=yoff;
		}
	}
    
	static ofstream out;
	stringstream ssname;
	ssname << inFile;

	string mesh_file=ssname.str()+"_mesh.ply";

	unsigned int fcount = 0;

	// read faces
	do{
		stringstream ssface;
		ssface << line;

		unsigned int nfv, v0, v1, v2;
		ssface >> nfv >> v0 >> v1 >> v2;
		out << nfv << " " << v0 << " " << v1 << " "<< v2 << "\n";

		if(nfv != 3){
			cerr << "File " << inFile << " contains a face with >3 (" << nfv << ") vertices.";
			exit(1);
		}

		Point pt1(pts[v0][0],pts[v0][1],pts[v0][2]);
		Point pt2(pts[v1][0],pts[v1][1],pts[v1][2]);
		Point pt3(pts[v2][0],pts[v2][1],pts[v2][2]);
		Triangle t(pt1,pt2,pt3);
		trs.push_back(t);

		fcount++;

	}while(getline(in, line) && fcount < nface);
    
    // constructs AABB tree
    mesh.rebuild(trs.begin(),trs.end());
    mesh.accelerate_distance_queries();

	in.close();
	return max;
}

void check_and_add(node* n, vector<node*>& skel_input){

	queue<node*> q;
	q.push(n);
	n->candidate=false;

	while(!q.empty()){
		// extract a node from the queue
		node *m=q.front();
		q.pop();

		// visit the node
		m->skel=true;
		if(CMinus(m)==1 && CStar(m)==1)
			skel_input.push_back(m);
		else
			m->skel=false;

		// mark the node as visited and add its neighbors to the queue
		for(int i=0; i<6; ++i){
			if(m->neighbors[i]->candidate){
				m->neighbors[i]->candidate=false;
				q.push(m->neighbors[i]);
			}
		}
	}
}

void dt(Tree& mesh, vector<node*>& leaves, vector<node*>& skel_input, int type){

	if(type==NORMAL){

		for(vector<node*>::iterator i=leaves.begin(); i!=leaves.end(); ++i){

			node *n=*i;

			Point point_query(n->x,n->y,n->z);
			Point endpt(0,0,0);

			n->np = mesh.closest_point(point_query);		
			n->dist=sqrt(CGAL::squared_distance(point_query,n->np));

			for(int i=0; i<3; ++i)
				n->grad[i]=point_query[i]-n->np[i];
			normalize(n->grad);

			Ray ray_query(point_query,endpt);

			// if the point is exterior, do not add it to the skeleton
			if(mesh.number_of_intersected_primitives(ray_query)%2==0){
				n->dist*=-1;
				n->grad[0]*=-1;
				n->grad[1]*=-1;
				n->grad[2]*=-1;
				n->dens=0.0;
				n->div=0.0;
				n->interior=false;
				n->skel=false;
			}

			// else add it only if it is simple
			else{
				n->interior=true;
				n->skel=true;
				skel_input.push_back(n);
			}
		}
	}

	else if(type==EXPLOSION){

		vector<node*> exterior;

		for(vector<node*>::iterator i=leaves.begin(); i!=leaves.end(); ++i){

			node *n=*i;

			Point point_query(n->x,n->y,n->z);
			Point endpt(0,0,0);
			/*if(n->parent->np!=NULL){
				closest_point = mesh.closest_point(point_query,*(n->parent->np));
			}
			else*/
				n->np = mesh.closest_point(point_query);		

			n->dist=sqrt(CGAL::squared_distance(point_query,n->np));
            
			for(int i=0; i<3; ++i)
				n->grad[i]=point_query[i]-n->np[i];
			normalize(n->grad);

		    Ray ray_query(point_query,endpt);

		    n->skel=false;
			if(mesh.number_of_intersected_primitives(ray_query)%2==0){
				n->dist*=-1;
				n->grad[0]*=-1;
				n->grad[1]*=-1;
				n->grad[2]*=-1;
				n->dens=0.0;
				n->div=0.0;
				n->interior=false;
				if(!n->border()){
					exterior.push_back(n);
				}
			}

		    else{
				n->interior=true;
				n->skel=true;
				skel_input.push_back(n);
			}
		}

		for(vector<node*>::iterator i=exterior.begin(); i!=exterior.end(); ++i){
			(*i)->skel=true;
		}

		sort(exterior.begin(), exterior.end(), distOrd);

		for(vector<node*>::iterator i=exterior.begin(); i!=exterior.end(); ++i){
			node *n=*i;
			if(CMinus(n)==1 && CStar(n)==1){
				n->skel=false;
			}
			else{
				skel_input.push_back(n);
			}
		}
	}

	else if(type==DILATION){

		vector<node*> candidates;

		for(vector<node*>::iterator i=leaves.begin(); i!=leaves.end(); ++i){

			node *n=*i;

			Point point_query(n->x,n->y,n->z);
			Point endpt(0,0,0);
			/*if(n->parent->np!=NULL){
				closest_point = mesh.closest_point(point_query,*(n->parent->np));
			}
			else*/
				n->np = mesh.closest_point(point_query);		

			n->dist=sqrt(CGAL::squared_distance(point_query,n->np));

			for(int i=0; i<3; ++i)
				n->grad[i]=point_query[i]-n->np[i];
			normalize(n->grad);

			Ray ray_query(point_query,endpt);

			// if the point is exterior, do not add it to the skeleton
			if(mesh.number_of_intersected_primitives(ray_query)%2==0){
				n->dist*=-1;
				n->grad[0]*=-1;
				n->grad[1]*=-1;
				n->grad[2]*=-1;
				n->dens=0.0;
				n->div=0.0;
				n->interior=false;
				n->skel=false;
				n->candidate=false;
			}

			// else add it only if it is simple
			else{
				n->interior=true;
				if(n->candidate){
					candidates.push_back(n);
				}
				else
					n->skel=false;
			}
		}

		for(vector<node*>::iterator i=candidates.begin(); i!=candidates.end(); ++i){
			node *n=*i;
			if(n->candidate && !n->skel){
				check_and_add(n,skel_input);
			}
		}
	}
}

void distmap2vtk(vector<node*> leaves, string filename){

	static ofstream dist_out;

    cout << "Saving "+filename+"_distmap.vtk..." << flush;

	string file=filename+"_distmap.vtk";
	dist_out.open(file.c_str(), ios::out);
	unsigned int num=0;
	double min=-inf;

	dist_out << "# vtk DataFile Version 2.0\ncvlab\nASCII\nDATASET UNSTRUCTURED_GRID\n";

	for(vector<node*>::iterator it=leaves.begin(); it!=leaves.end(); ++it){
		if((*it)->dist>0){
			num++;
			if(min<(*it)->dist)
				min=(*it)->dist;
		}
	}

	dist_out << "POINTS " << num << " double\n";

	for(vector<node*>::iterator it=leaves.begin(); it!=leaves.end(); ++it)
		if((*it)->dist>0)
			dist_out << std::fixed << std::setprecision(12) << (double)(*it)->x << " " << (double)(*it)->y << " " << (double)(*it)->z << endl;

	dist_out << "CELLS " << num << " " << 2*num << endl;

	for(unsigned int i=0; i<num; i++){
		dist_out << "1 " << i << endl;
	}

	dist_out << "CELL_TYPES " << num << endl;

	for(unsigned int i=0; i<num; i++){
		dist_out << "1" << endl;
	}

	dist_out << "POINT_DATA " << num << endl;
	dist_out << "VECTORS ScanColors unsigned_char\n";

	for(vector<node*>::iterator it=leaves.begin(); it!=leaves.end(); ++it){
		if((*it)->dist>0)
			dist_out << 255-(int)floor((*it)->dist/min*255) << " " << 255-(int)floor((*it)->dist/min*255) << " " << 255-(int)floor((*it)->dist/min*255)  << endl;
	}
    cout << " done." << endl;
    dist_out.close();
}
