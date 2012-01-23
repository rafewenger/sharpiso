/*
 *  sh_point_datastruct.cxx
 *  SHARPPOINT
 *
 *  Created by arindam bhattacharya on 11/6/11.
 *  Copyright 2011 Ohio State University. All rights reserved.
 *
 */

#include <iostream>
#include "sh_point_datastruct.h"
#include "sharpiso_types.h"
using namespace std;
using namespace SHARPISO;
using namespace sh_cube;

bool sh_cube::setup_shCube
(CUBE &cb,
 const GRADIENT_COORD_TYPE  gradients[],
 const SCALAR_TYPE isovalue,
 const SCALAR_TYPE scalar_vals[])
{
    
    float position[] = {0,0,0 ,1,0,0, 0,1,0, 1,1,0, 0,0,1, 1,0,1, 0,1,1, 1,1,1};
    cb.num_edges = 12;
    cb.num_pts = 8;
    cb.dim =3;
    cb.ne_intersect=0;
    int dim = cb.dim;
    
    //declare and initialize max and min scalar values
    SCALAR_TYPE max_sval, min_sval;
    
    max_sval = scalar_vals[0];
    min_sval = scalar_vals[1];
    
    for (int i=0; i<cb.num_pts; i++) {
        POINT p;
        //set scalar
        p.scalar = scalar_vals[i];
        
        if (scalar_vals[i] > max_sval)
            max_sval = scalar_vals[i];
        
        if(scalar_vals[i] < min_sval)
            min_sval = scalar_vals[i];
        
        for (int k=0; k<cb.dim; k++) {
            //set gradient
            p.grads[k]=gradients[i*dim+k];
            //set pos 
            p.pos[k]=position[i*dim+k];
        }
        cb.pts.push_back(p); 
    }
    
    if (!(isovalue >= min_sval && isovalue <= max_sval))
	{ 
        return false;
	}
    EDGE temp;
    temp.p1=cb.pts[0];
    temp.p2=cb.pts[1];
    cb.edges.push_back(temp); //1
    
    temp.p1=cb.pts[1];
    temp.p2=cb.pts[3];
    cb.edges.push_back(temp); //2
    
    temp.p1=cb.pts[3];
    temp.p2=cb.pts[2];
    cb.edges.push_back(temp); //3
    
    temp.p1=cb.pts[2];
    temp.p2=cb.pts[0];
    cb.edges.push_back(temp); //4
    
    temp.p1=cb.pts[4];
    temp.p2=cb.pts[5];
    cb.edges.push_back(temp); //5
    
    temp.p1=cb.pts[5];
    temp.p2=cb.pts[7];
    cb.edges.push_back(temp); //6
    
    temp.p1=cb.pts[7];
    temp.p2=cb.pts[6];
    cb.edges.push_back(temp); //7
    
    temp.p1=cb.pts[6];
    temp.p2=cb.pts[4];
    cb.edges.push_back(temp); //8
    
    temp.p1=cb.pts[0];
    temp.p2=cb.pts[4];
    cb.edges.push_back(temp); //9
    
    temp.p1=cb.pts[1];
    temp.p2=cb.pts[5];
    cb.edges.push_back(temp); //10
    
    temp.p1=cb.pts[3];
    temp.p2=cb.pts[7];
    cb.edges.push_back(temp); //11
    
    temp.p1=cb.pts[2];
    temp.p2=cb.pts[6];
    cb.edges.push_back(temp); //12
    
    return true;
};

/*
 PRINT POINT HELPER FUNCTION FOR DEBUGGING
 */
void POINT::print_point() const
{
	cout <<"\nscalar "<<scalar;
    cout <<" position ("<<pos[0]<<","<<pos[1]<<","<<pos[2]<<")";
    cout <<" gradient ("<<grads[0]<<","<<grads[1]<<","<<grads[2]<<")"<<endl;
};

void POINT::normalize_grads(){
    double sum=0;
    for (int i=0; i<3; i++) {
        sum = sum + grads[i]*grads[i];
    } 
    
    if (sum>0) {
        sum = sqrt(sum);
        for (int i=0; i<3; i++) {
            grads[i]=grads[i]/sum;
        }
    }
    
    
};
