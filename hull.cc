#include "hull.h"

using namespace std;

bool duplicated(vector<h_face>& faces, h_face* f){
	for(vector<h_face>::iterator it=faces.begin(); it!=faces.end(); ++it){
		if(*f==*it){
			return true;
		}
	}
	return false;
}

void overlapping_neigh(node *n, vector<node*>& neighbors, vector<char>& on, char id, char v_num){
	if(id==4 || id==10 || id==12 || id==14 || id==16 || id==22){

		char n1,n2,n3,n4,dir,c1,c2,c3,c4,c5,c6,c7,c8;
		bool b1=false,b2=false,b3=false,b4=false;

		if(id==4){n1=10;n2=12;n3=16;n4=14;dir=4;c1=1;c2=3;c3=7;c4=5;c5=0;c6=6;c7=8;c8=2;}

		if(id==10){n1=4;n2=12;n3=22;n4=14;dir=3;c1=1;c2=9;c3=19;c4=11;c5=0;c6=18;c7=20;c8=2;}

		if(id==12){n1=4;n2=10;n3=22;n4=16;dir=0;c1=3;c2=9;c3=21;c4=15;c5=0;c6=18;c7=24;c8=6;}

		if(id==16){n1=4;n2=12;n3=22;n4=14;dir=1;c1=7;c2=15;c3=25;c4=17;c5=6;c6=24;c7=26;c8=8;}

		if(id==14){n1=4;n2=10;n3=22;n4=16;dir=2;c1=5;c2=11;c3=23;c4=17;c5=2;c6=20;c7=26;c8=8;}

		if(id==22){n1=10;n2=12;n3=16;n4=14;dir=5;c1=19;c2=21;c3=25;c4=23;c5=18;c6=24;c7=26;c8=20;}

		if(neighbors[n1]!=NULL)
			if(neighbors[n1]->neighbors[dir]==neighbors[id]){b1=true;}
		if(neighbors[n2]!=NULL)
			if(neighbors[n2]->neighbors[dir]==neighbors[id]){b2=true;}
		if(neighbors[n3]!=NULL)
			if(neighbors[n3]->neighbors[dir]==neighbors[id]){b3=true;}
		if(neighbors[n4]!=NULL)
			if(neighbors[n4]->neighbors[dir]==neighbors[id]){b4=true;}

		if(b1){on[c1]=v_num;}
		if(b2){on[c2]=v_num;}
		if(b3){on[c3]=v_num;}
		if(b4){on[c4]=v_num;}
		if(b1&&b2){on[c5]=v_num;}
		if(b2&&b3){on[c6]=v_num;}
		if(b3&&b4){on[c7]=v_num;}
		if(b4&&b1){on[c8]=v_num;}
	}

	else if(id==1 || id==3 || id==5 || id==7 || id==9 || id==11 ||
			id==15 || id==17 || id==19 || id==21 || id==23 || id==25){

		char n1,n2,c1,c2;

		if(id==1){c1=2;c2=0;n1=2;n2=0;}
		if(id==3){c1=0;c2=6;n1=3;n2=1;}
		if(id==7){c1=8;c2=6;n1=2;n2=0;}
		if(id==5){c1=2;c2=8;n1=3;n2=1;}
		if(id==9){c1=0;c2=18;n1=4;n2=5;}
		if(id==15){c1=6;c2=24;n1=4;n2=5;}
		if(id==17){c1=8;c2=26;n1=4;n2=5;}
		if(id==11){c1=2;c2=20;n1=4;n2=5;}
		if(id==19){c1=20;c2=18;n1=2;n2=0;}
		if(id==21){c1=18;c2=24;n1=3;n2=1;}
		if(id==25){c1=26;c2=24;n1=2;n2=0;}
		if(id==23){c1=20;c2=26;n1=3;n2=1;}

		if(neighbors[c1]==NULL && isUberConsistent(neighbors[id]->neighbors[n1],n,neighbors)){on[c1]=v_num;}
		if(neighbors[c2]==NULL && isUberConsistent(neighbors[id]->neighbors[n2],n,neighbors)){on[c2]=v_num;}
	}
}

h_vertex::h_vertex(){
}

h_vertex::h_vertex(node *m, node *center){

	normal[0]=m->x-center->x;
	normal[1]=m->y-center->y;
	normal[2]=m->z-center->z;

	nl=length(normal);
	normalize(normal);

	w=0.0;
	n=m;
}

h_face::h_face(){}

h_face::h_face(h_vertex* v1, h_vertex* v2, h_vertex* v3){
	vts[0]=v1;
	vts[1]=v2;
	vts[2]=v3;
}

hull::hull(node *n,vector<node*> &neighbors){

	node *m;
	char v_num=0;
	vector<char> on(27,27);

	h_vertex v;

	// determine which neighbors are "covered" by other neighbors
	for(int i=0; i<27; ++i){
		if(i!=13){
			m=neighbors[i];
			if(m!=NULL){
				on[i]=v_num;
				v=h_vertex(m,n);
				vts.push_back(v);
				if(m->level()<n->level()){
					overlapping_neigh(n,neighbors,on,i,v_num);
				}
				v_num++;
			}
		}
	}

/*
	// just a check for error. TODO: remove it later
	for(int i=0; i<27; ++i){
		if(i!=13){
			if(on[i]==27)
				cerr << "There was an error during smoothing.\n";
		}
	}
*/

	// now for each face of the 27-neighborhood, create its triangulation
	h_face f;

	if(vts[on[0]].n!=vts[on[3]].n && vts[on[0]].n!=vts[on[4]].n && vts[on[3]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[0]],&vts[on[3]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[3]].n!=vts[on[6]].n && vts[on[4]].n!=vts[on[6]].n && vts[on[3]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[3]],&vts[on[6]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[6]].n!=vts[on[7]].n && vts[on[4]].n!=vts[on[7]].n && vts[on[6]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[6]],&vts[on[7]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[7]].n!=vts[on[8]].n && vts[on[4]].n!=vts[on[8]].n && vts[on[7]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[7]],&vts[on[8]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[8]].n!=vts[on[5]].n && vts[on[8]].n!=vts[on[4]].n && vts[on[4]].n!=vts[on[5]].n){
		f=h_face(&vts[on[4]],&vts[on[8]],&vts[on[5]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[5]].n!=vts[on[2]].n && vts[on[4]].n!=vts[on[2]].n && vts[on[5]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[5]],&vts[on[2]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[2]].n!=vts[on[1]].n && vts[on[4]].n!=vts[on[1]].n && vts[on[2]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[2]],&vts[on[1]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[1]].n!=vts[on[0]].n && vts[on[1]].n!=vts[on[4]].n && vts[on[4]].n!=vts[on[0]].n){
		f=h_face(&vts[on[4]],&vts[on[1]],&vts[on[0]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[18]].n!=vts[on[21]].n && vts[on[18]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[21]].n){
		f=h_face(&vts[on[22]],&vts[on[18]],&vts[on[21]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[21]].n!=vts[on[24]].n && vts[on[21]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[24]].n){
		f=h_face(&vts[on[22]],&vts[on[21]],&vts[on[24]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[24]].n!=vts[on[25]].n && vts[on[24]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[25]].n){
		f=h_face(&vts[on[22]],&vts[on[24]],&vts[on[25]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[25]].n!=vts[on[26]].n && vts[on[25]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[26]].n){
		f=h_face(&vts[on[22]],&vts[on[25]],&vts[on[26]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[26]].n!=vts[on[23]].n && vts[on[26]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[23]].n){
		f=h_face(&vts[on[22]],&vts[on[26]],&vts[on[23]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[23]].n!=vts[on[20]].n && vts[on[23]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[20]].n){
		f=h_face(&vts[on[22]],&vts[on[23]],&vts[on[20]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[20]].n!=vts[on[19]].n && vts[on[20]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[19]].n){
		f=h_face(&vts[on[22]],&vts[on[20]],&vts[on[19]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[19]].n!=vts[on[18]].n && vts[on[19]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[18]].n){
		f=h_face(&vts[on[22]],&vts[on[19]],&vts[on[18]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[0]].n!=vts[on[1]].n && vts[on[0]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[1]].n){
		f=h_face(&vts[on[10]],&vts[on[0]],&vts[on[1]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[1]].n!=vts[on[2]].n && vts[on[1]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[2]].n){
		f=h_face(&vts[on[10]],&vts[on[1]],&vts[on[2]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[2]].n!=vts[on[11]].n && vts[on[2]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[11]].n){
		f=h_face(&vts[on[10]],&vts[on[2]],&vts[on[11]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[11]].n!=vts[on[20]].n && vts[on[11]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[20]].n){
		f=h_face(&vts[on[10]],&vts[on[11]],&vts[on[20]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[20]].n!=vts[on[19]].n && vts[on[20]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[19]].n){
		f=h_face(&vts[on[10]],&vts[on[20]],&vts[on[19]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[19]].n!=vts[on[18]].n && vts[on[19]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[18]].n){
		f=h_face(&vts[on[10]],&vts[on[19]],&vts[on[18]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[18]].n!=vts[on[9]].n && vts[on[18]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[9]].n){
		f=h_face(&vts[on[10]],&vts[on[18]],&vts[on[9]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[9]].n!=vts[on[0]].n && vts[on[9]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[0]].n){
		f=h_face(&vts[on[10]],&vts[on[9]],&vts[on[0]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[6]].n!=vts[on[7]].n && vts[on[6]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[7]].n){
		f=h_face(&vts[on[16]],&vts[on[6]],&vts[on[7]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[7]].n!=vts[on[8]].n && vts[on[7]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[8]].n){
		f=h_face(&vts[on[16]],&vts[on[7]],&vts[on[8]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[8]].n!=vts[on[17]].n && vts[on[8]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[17]].n){
		f=h_face(&vts[on[16]],&vts[on[8]],&vts[on[17]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[17]].n!=vts[on[26]].n && vts[on[17]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[26]].n){
		f=h_face(&vts[on[16]],&vts[on[17]],&vts[on[26]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[26]].n!=vts[on[25]].n && vts[on[26]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[25]].n){
		f=h_face(&vts[on[16]],&vts[on[26]],&vts[on[25]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[25]].n!=vts[on[24]].n && vts[on[25]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[24]].n){
		f=h_face(&vts[on[16]],&vts[on[25]],&vts[on[24]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[24]].n!=vts[on[15]].n && vts[on[24]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[15]].n){
		f=h_face(&vts[on[16]],&vts[on[24]],&vts[on[15]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[15]].n!=vts[on[6]].n && vts[on[15]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[6]].n){
		f=h_face(&vts[on[16]],&vts[on[15]],&vts[on[6]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[0]].n!=vts[on[3]].n && vts[on[0]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[3]].n){
		f=h_face(&vts[on[12]],&vts[on[0]],&vts[on[3]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[3]].n!=vts[on[6]].n && vts[on[3]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[6]].n ){
		f=h_face(&vts[on[12]],&vts[on[3]],&vts[on[6]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[6]].n!=vts[on[15]].n && vts[on[6]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[15]].n){
		f=h_face(&vts[on[12]],&vts[on[6]],&vts[on[15]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[15]].n!=vts[on[24]].n && vts[on[15]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[24]].n){
		f=h_face(&vts[on[12]],&vts[on[15]],&vts[on[24]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[24]].n!=vts[on[21]].n && vts[on[24]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[21]].n){
		f=h_face(&vts[on[12]],&vts[on[24]],&vts[on[21]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[21]].n!=vts[on[18]].n && vts[on[21]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[18]].n){
		f=h_face(&vts[on[12]],&vts[on[21]],&vts[on[18]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[18]].n!=vts[on[9]].n && vts[on[18]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[9]].n){
		f=h_face(&vts[on[12]],&vts[on[18]],&vts[on[9]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[9]].n!=vts[on[0]].n && vts[on[9]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[0]].n){
		f=h_face(&vts[on[12]],&vts[on[9]],&vts[on[0]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[2]].n!=vts[on[5]].n && vts[on[2]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[5]].n){
		f=h_face(&vts[on[14]],&vts[on[2]],&vts[on[5]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[5]].n!=vts[on[8]].n && vts[on[5]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[8]].n){
		f=h_face(&vts[on[14]],&vts[on[5]],&vts[on[8]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[8]].n!=vts[on[17]].n && vts[on[8]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[17]].n){
		f=h_face(&vts[on[14]],&vts[on[8]],&vts[on[17]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[17]].n!=vts[on[26]].n && vts[on[17]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[26]].n){
		f=h_face(&vts[on[14]],&vts[on[17]],&vts[on[26]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[26]].n!=vts[on[23]].n && vts[on[26]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[23]].n){
		f=h_face(&vts[on[14]],&vts[on[26]],&vts[on[23]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[23]].n!=vts[on[20]].n && vts[on[23]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[20]].n){
		f=h_face(&vts[on[14]],&vts[on[23]],&vts[on[20]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[20]].n!=vts[on[11]].n && vts[on[20]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[11]].n){
		f=h_face(&vts[on[14]],&vts[on[20]],&vts[on[11]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[11]].n!=vts[on[2]].n && vts[on[11]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[2]].n){
		f=h_face(&vts[on[14]],&vts[on[11]],&vts[on[2]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
}

hull::hull(double x,double y,double z,vector<node*> &neighbors){

	node *m;
	node n;
	n.x=x;
	n.y=y;
	n.z=z;
	char v_num=0;
	char vn=0;
	vector<char> on(27,27);

	h_vertex v;

	// determine which neighbors are "covered" by other neighbors
	for(int i=0; i<27; ++i){
		m=neighbors[i];
		if(m!=NULL){
			on[i]=v_num;
			v=h_vertex(m,&n);
			vts.push_back(v);
			if(m->level()<neighbors[13]->level()){
				overlapping_neigh(neighbors[13],neighbors,on,i,v_num);
			}
			if(i==13)
				vn=v_num;
			v_num++;
		}
	}

	for(int i=0; i<27; ++i){
		if(on[i]==27)
			on[i]=vn;
	}

	// now for each face of the 27-neighborhood, create its triangulation
	h_face f;

	if(vts[on[0]].n!=vts[on[3]].n && vts[on[0]].n!=vts[on[4]].n && vts[on[3]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[0]],&vts[on[3]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[3]].n!=vts[on[6]].n && vts[on[4]].n!=vts[on[6]].n && vts[on[3]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[3]],&vts[on[6]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[6]].n!=vts[on[7]].n && vts[on[4]].n!=vts[on[7]].n && vts[on[6]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[6]],&vts[on[7]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[7]].n!=vts[on[8]].n && vts[on[4]].n!=vts[on[8]].n && vts[on[7]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[7]],&vts[on[8]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[8]].n!=vts[on[5]].n && vts[on[8]].n!=vts[on[4]].n && vts[on[4]].n!=vts[on[5]].n){
		f=h_face(&vts[on[4]],&vts[on[8]],&vts[on[5]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[5]].n!=vts[on[2]].n && vts[on[4]].n!=vts[on[2]].n && vts[on[5]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[5]],&vts[on[2]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[2]].n!=vts[on[1]].n && vts[on[4]].n!=vts[on[1]].n && vts[on[2]].n!=vts[on[4]].n){
		f=h_face(&vts[on[4]],&vts[on[2]],&vts[on[1]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[1]].n!=vts[on[0]].n && vts[on[1]].n!=vts[on[4]].n && vts[on[4]].n!=vts[on[0]].n){
		f=h_face(&vts[on[4]],&vts[on[1]],&vts[on[0]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[18]].n!=vts[on[21]].n && vts[on[18]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[21]].n){
		f=h_face(&vts[on[22]],&vts[on[18]],&vts[on[21]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[21]].n!=vts[on[24]].n && vts[on[21]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[24]].n){
		f=h_face(&vts[on[22]],&vts[on[21]],&vts[on[24]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[24]].n!=vts[on[25]].n && vts[on[24]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[25]].n){
		f=h_face(&vts[on[22]],&vts[on[24]],&vts[on[25]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[25]].n!=vts[on[26]].n && vts[on[25]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[26]].n){
		f=h_face(&vts[on[22]],&vts[on[25]],&vts[on[26]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[26]].n!=vts[on[23]].n && vts[on[26]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[23]].n){
		f=h_face(&vts[on[22]],&vts[on[26]],&vts[on[23]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[23]].n!=vts[on[20]].n && vts[on[23]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[20]].n){
		f=h_face(&vts[on[22]],&vts[on[23]],&vts[on[20]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[20]].n!=vts[on[19]].n && vts[on[20]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[19]].n){
		f=h_face(&vts[on[22]],&vts[on[20]],&vts[on[19]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[19]].n!=vts[on[18]].n && vts[on[19]].n!=vts[on[22]].n && vts[on[22]].n!=vts[on[18]].n){
		f=h_face(&vts[on[22]],&vts[on[19]],&vts[on[18]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[0]].n!=vts[on[1]].n && vts[on[0]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[1]].n){
		f=h_face(&vts[on[10]],&vts[on[0]],&vts[on[1]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[1]].n!=vts[on[2]].n && vts[on[1]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[2]].n){
		f=h_face(&vts[on[10]],&vts[on[1]],&vts[on[2]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[2]].n!=vts[on[11]].n && vts[on[2]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[11]].n){
		f=h_face(&vts[on[10]],&vts[on[2]],&vts[on[11]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[11]].n!=vts[on[20]].n && vts[on[11]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[20]].n){
		f=h_face(&vts[on[10]],&vts[on[11]],&vts[on[20]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[20]].n!=vts[on[19]].n && vts[on[20]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[19]].n){
		f=h_face(&vts[on[10]],&vts[on[20]],&vts[on[19]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[19]].n!=vts[on[18]].n && vts[on[19]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[18]].n){
		f=h_face(&vts[on[10]],&vts[on[19]],&vts[on[18]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[18]].n!=vts[on[9]].n && vts[on[18]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[9]].n){
		f=h_face(&vts[on[10]],&vts[on[18]],&vts[on[9]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[9]].n!=vts[on[0]].n && vts[on[9]].n!=vts[on[10]].n && vts[on[10]].n!=vts[on[0]].n){
		f=h_face(&vts[on[10]],&vts[on[9]],&vts[on[0]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[6]].n!=vts[on[7]].n && vts[on[6]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[7]].n){
		f=h_face(&vts[on[16]],&vts[on[6]],&vts[on[7]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[7]].n!=vts[on[8]].n && vts[on[7]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[8]].n){
		f=h_face(&vts[on[16]],&vts[on[7]],&vts[on[8]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[8]].n!=vts[on[17]].n && vts[on[8]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[17]].n){
		f=h_face(&vts[on[16]],&vts[on[8]],&vts[on[17]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[17]].n!=vts[on[26]].n && vts[on[17]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[26]].n){
		f=h_face(&vts[on[16]],&vts[on[17]],&vts[on[26]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[26]].n!=vts[on[25]].n && vts[on[26]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[25]].n){
		f=h_face(&vts[on[16]],&vts[on[26]],&vts[on[25]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[25]].n!=vts[on[24]].n && vts[on[25]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[24]].n){
		f=h_face(&vts[on[16]],&vts[on[25]],&vts[on[24]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[24]].n!=vts[on[15]].n && vts[on[24]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[15]].n){
		f=h_face(&vts[on[16]],&vts[on[24]],&vts[on[15]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[15]].n!=vts[on[6]].n && vts[on[15]].n!=vts[on[16]].n && vts[on[16]].n!=vts[on[6]].n){
		f=h_face(&vts[on[16]],&vts[on[15]],&vts[on[6]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[0]].n!=vts[on[3]].n && vts[on[0]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[3]].n){
		f=h_face(&vts[on[12]],&vts[on[0]],&vts[on[3]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[3]].n!=vts[on[6]].n && vts[on[3]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[6]].n ){
		f=h_face(&vts[on[12]],&vts[on[3]],&vts[on[6]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[6]].n!=vts[on[15]].n && vts[on[6]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[15]].n){
		f=h_face(&vts[on[12]],&vts[on[6]],&vts[on[15]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[15]].n!=vts[on[24]].n && vts[on[15]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[24]].n){
		f=h_face(&vts[on[12]],&vts[on[15]],&vts[on[24]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[24]].n!=vts[on[21]].n && vts[on[24]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[21]].n){
		f=h_face(&vts[on[12]],&vts[on[24]],&vts[on[21]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[21]].n!=vts[on[18]].n && vts[on[21]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[18]].n){
		f=h_face(&vts[on[12]],&vts[on[21]],&vts[on[18]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[18]].n!=vts[on[9]].n && vts[on[18]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[9]].n){
		f=h_face(&vts[on[12]],&vts[on[18]],&vts[on[9]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[9]].n!=vts[on[0]].n && vts[on[9]].n!=vts[on[12]].n && vts[on[12]].n!=vts[on[0]].n){
		f=h_face(&vts[on[12]],&vts[on[9]],&vts[on[0]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts[on[2]].n!=vts[on[5]].n && vts[on[2]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[5]].n){
		f=h_face(&vts[on[14]],&vts[on[2]],&vts[on[5]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[5]].n!=vts[on[8]].n && vts[on[5]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[8]].n){
		f=h_face(&vts[on[14]],&vts[on[5]],&vts[on[8]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[8]].n!=vts[on[17]].n && vts[on[8]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[17]].n){
		f=h_face(&vts[on[14]],&vts[on[8]],&vts[on[17]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[17]].n!=vts[on[26]].n && vts[on[17]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[26]].n){
		f=h_face(&vts[on[14]],&vts[on[17]],&vts[on[26]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[26]].n!=vts[on[23]].n && vts[on[26]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[23]].n){
		f=h_face(&vts[on[14]],&vts[on[26]],&vts[on[23]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[23]].n!=vts[on[20]].n && vts[on[23]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[20]].n){
		f=h_face(&vts[on[14]],&vts[on[23]],&vts[on[20]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[20]].n!=vts[on[11]].n && vts[on[20]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[11]].n){
		f=h_face(&vts[on[14]],&vts[on[20]],&vts[on[11]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}
	if(vts[on[11]].n!=vts[on[2]].n && vts[on[11]].n!=vts[on[14]].n && vts[on[14]].n!=vts[on[2]].n){
		f=h_face(&vts[on[14]],&vts[on[11]],&vts[on[2]]);
		if(!duplicated(faces,&f)) faces.push_back(f);
	}

	if(vts.size()<4)
		cerr << faces.size() << " " << vts.size() << endl;
}
