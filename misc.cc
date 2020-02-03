#include"misc.h"
#include <iostream>

using namespace std;

double dotProd(double v1[3], double v2[3]){
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

void crossProd(double v1[3], double v2[3], double cp[3]){
	cp[0] = v1[1]*v2[2]-v1[2]*v2[1];
	cp[1] = v1[2]*v2[0]-v1[0]*v2[2];
	cp[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

double length(double v1[3]){
	return sqrt(dotProd(v1,v1));
}

void normalize(double v[3]){
	double l=length(v);
	v[0]/=l;
	v[1]/=l;
	v[2]/=l;
}

double det(double b[3][3], int m){

	double sum = 0.0;
	double c[3][3];

	if(m==2){
		return b[0][0]*b[1][1] - b[0][1]*b[1][0];
	}

	for(int p=0; p<m; p++){
		int h = 0,k = 0;
		for(int i=1; i<m; i++){
			for(int j=0; j<m; j++){
				if(j==p)
					continue;
				c[h][k] = b[i][j];
				k++;
				if(k == m-1){
					h++;
					k = 0;
				}

			}
		}
		sum = sum + b[0][p]*pow(-1.0,p)*det(c,m-1);
	}
	return sum;
}
