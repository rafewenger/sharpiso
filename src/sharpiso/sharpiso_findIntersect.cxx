//FindIntersect.cxx
//

#include <iostream>
#include <cmath>
#include <limits>

#include "sharpiso_findIntersect.h"


using namespace std;
using namespace SHARPISO;


// HELPER FUNCTIONS:
// Findmax direction
int findmax(const float dir[]);



// FindIntersect  ,
// accepts as inputs  a point p[], a direction dir[].It returns a boolean which is TRUE
// if the Ray intersects the cube. It returns FALSE if the ray does not intersect the cube.
// If the bool is true , intersect[] returns the MID-point of intersection of the ray and the cueb.

bool SHARPISO::calculate_point_intersect
(const SCALAR_TYPE *p,
 const SCALAR_TYPE *dir,
 SCALAR_TYPE *intersect)
{
    
    int ind  = findmax(dir);
    //find t0 and t1
    SCALAR_TYPE t0 = ( -1.0*p[ind])/(dir[ind]);
    SCALAR_TYPE t1 = (1.0 - p[ind])/(dir[ind]);
    SCALAR_TYPE min_coord_j;
    SCALAR_TYPE max_coord_j;
    
    for (int j=0; j<3; j++)
    {
        if(j != ind)
        {
            min_coord_j = p[j] + t0*dir[j];
            max_coord_j = p[j] + t1*dir[j];
            
            //swap
            if (min_coord_j > max_coord_j)
            {
                swap(t0, t1);
                swap(min_coord_j, max_coord_j);
            }
            if (min_coord_j - 1.0 < THRESHOLD_CLAMP ) {
                min_coord_j = 1.0;
            }
            if (max_coord_j-0.0 > -THRESHOLD_CLAMP ) {
                max_coord_j=0.0;
            }
            
            if(min_coord_j > 1.0)
            {
                return false ;
            }
            else if(max_coord_j < 0.0)
            {
                return false;
            }
            else
            {
                if ((min_coord_j <= 0.0) && (max_coord_j > 0.0 ))
                {
                    t0  = (-p[j])/dir[j];
                }
                if ((min_coord_j < 1.0) &&(max_coord_j >= 1.0 ))
                {
                    t1  = (1.0 - p[j])/dir[j];
                }
            }
        }
    }
    
    SCALAR_TYPE endPt0[3];
    SCALAR_TYPE endPt1[3];
    for (int i=0; i<3; i++)
    {
        SCALAR_TYPE temp = p[i] + t0*dir[i];
        endPt0[i] = temp;
    }
    //t1
    for (int i=0; i<3; i++)
    {
        SCALAR_TYPE temp = p[i] + t1*dir[i];
        endPt1[i] = temp;
    }
    //Find the point of intersection.
    for (int i=0; i<3; i++) {
        intersect[i] = (endPt0[i] + endPt1[i])/2.0;
    }
    return true;
    
};


// Caculate point intesect between the ray and the cube.
// This version takes in the coordinate of the cube index and translates,
// to find the intersection with the unit cube and translate it back.

bool SHARPISO::calculate_point_intersect
(const COORD_TYPE cube_coord[],
 const SCALAR_TYPE *original_pt,
 const SCALAR_TYPE *dir,
 SCALAR_TYPE *intersect)
{
	SCALAR_TYPE p[DIM3] ={0.0};
	for (int i=0; i < DIM3; i++){
		p[i] = original_pt[i] - cube_coord[i];
	}
    
    int ind  = findmax(dir);
    //find t0 and t1
    SCALAR_TYPE t0 = ( -1.0*p[ind])/(dir[ind]);
    SCALAR_TYPE t1 = (1.0 - p[ind])/(dir[ind]);
    SCALAR_TYPE min_coord_j;
    SCALAR_TYPE max_coord_j;
    
    for (int j=0; j<3; j++)
    {
        if(j != ind)
        {
            min_coord_j = p[j] + t0*dir[j];
            max_coord_j = p[j] + t1*dir[j];
            //swap
            if (min_coord_j > max_coord_j)
            {
                swap(t0, t1);
                swap(min_coord_j, max_coord_j);
            }
            
            if (min_coord_j - 1.0 < THRESHOLD_CLAMP ) {
                min_coord_j = 1.0;
            }
            if (max_coord_j-0.0 > -THRESHOLD_CLAMP ) {
                max_coord_j=0.0;
            }
            
            if(min_coord_j > 1.0)
            {
                return false ;
            }
            else if(max_coord_j < 0.0)
            {
                return false;
            }
            else
            {
                
                if ((min_coord_j <= 0.0) && (max_coord_j > 0.0 ))
                {
                    t0  = (-p[j])/dir[j];
                }
                if ((min_coord_j < 1.0) &&(max_coord_j >= 1.0 ))
                {
                    t1  = (1.0 - p[j])/dir[j];
                }
            }
        }
    }
    
    SCALAR_TYPE endPt0[DIM3];
    SCALAR_TYPE endPt1[DIM3];
    for (int i=0; i<DIM3; i++)
    {
        SCALAR_TYPE temp = p[i] + t0*dir[i];
        endPt0[i] = temp;
    }
    //t1
    for (int i=0; i<DIM3; i++)
    {
        SCALAR_TYPE temp = p[i] + t1*dir[i];
        endPt1[i] = temp;
    }
    /* /// used for debug purpose
     cout <<" pt 1 x:["<< endPt0[0] <<"] y:["<< endPt0[1] <<"] z:["<< endPt0[2] <<"]"
     <<" t0: "<<t0<<endl;
     cout <<" pt 2 x:["<< endPt1[0] <<"] y:["<< endPt1[1] <<"] z:["<< endPt1[2] <<"]"
     <<" t1: "<<t1<<endl;
     */
    //Find the point of intersection.
    for (int i=0; i<DIM3; i++) {
        intersect[i] = (endPt0[i] + endPt1[i])/2.0;
        intersect[i] += cube_coord[i];
    }
    return true;
    
};




// Calculate the intersection with a larger cube.
// Takes an extra parameter of how big the cube is 

bool SHARPISO::calculate_point_intersect_complex
(const COORD_TYPE cube_coord[],
 const SCALAR_TYPE *original_pt,
 const SCALAR_TYPE *dir,
 const float th,
 SCALAR_TYPE *intersect)
{ 
	SCALAR_TYPE p[DIM3] ={0.0};
    
	for (int i = 0; i < DIM3; i++){
		p[i] = original_pt[i] - cube_coord[i];
	}
    
    int ind  = findmax(dir);    
    //find t0 and t1
    SCALAR_TYPE t0 = (-th -p[ind])/(dir[ind]);
    SCALAR_TYPE t1 = (1.0 + th - p[ind])/(dir[ind]);
    SCALAR_TYPE min_coord_j;
    SCALAR_TYPE max_coord_j;
    
    for (int j = 0; j < DIM3; j++)
    {
        if(j != ind)
        {
            min_coord_j = p[j] + t0*dir[j];
            max_coord_j = p[j] + t1*dir[j];
            //swap
            if (min_coord_j > max_coord_j)
            {
                swap(t0, t1);
                swap(min_coord_j, max_coord_j);
            }
            
            if (( min_coord_j - (1.0 + th)) > EPSILON)  
            {
                return false ;
            }
            else if (((0.0 - th) - max_coord_j) > EPSILON)
            {
                return false;
            }
            else
            {
                if ((min_coord_j <= (0.0 - th))
                    && (max_coord_j > (0.0 - th)))
                {  t0 = (-th - p[j])/dir[j]; }
                
                
                if ((min_coord_j < (1.0 + th)) &&
                    (max_coord_j >= (1.0 + th)))
                { t1  = (1.0 + th - p[j])/dir[j]; }
            }
        }
    }
    
    SCALAR_TYPE endPt0[3];
    SCALAR_TYPE endPt1[3];
    for (int i=0; i<3; i++)
    {
        endPt0[i] = p[i] + t0*dir[i];
        endPt1[i] = p[i] + t1*dir[i];
    }
    /*
     /// used for debug purpose
     cout <<" pt 1 ["<< endPt0[0] <<" "<< endPt0[1] <<" "<< endPt0[2] <<"]"
     <<" t0: "<<t0<<endl;
     cout <<" pt 2 ["<< endPt1[0] <<" "<< endPt1[1] <<"  "<< endPt1[2] <<"]"
     <<" t1: "<<t1<<endl;
     */
    
    //Find the point of intersection.
    for (int i = 0; i < DIM3; i++) {
        intersect[i] = (endPt0[i] + endPt1[i] ) / 2.0;
        intersect[i] = intersect[i] + cube_coord[i];
    }
    return true;    
};


// HELPER FUNCTIONS:
// find the max direction
int findmax(const SCALAR_TYPE dir[])
{
    int ind=0;
    SCALAR_TYPE max = abs(dir[ind]);
    for (int j=0;j<3;j++)
    {
        if (abs(dir[j]) >= max)
        {
            max = abs(dir[j]);
            ind = j;
        }
    }
    return ind;
}

