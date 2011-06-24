// This file implements the quaternion's.
// Quaternions have find the main application in object rotation as it has the ability to
// generate smoother animation/rotation compared to Euler angles and devoid of Gumbal lock.

/*
Design: I have not used the vector class. I have and I know I should do it, since else
code has lost the elgance :(. But the main reason for avoiding Vector class is to reduce
dependency. So in current sceneraio anyone can take this one file and plug into their code.
-- WHich is the main motivation for writing the Quaternion class. :)

Version Control: 
9/0  Fixed the bug in multiplication routine
8/26 Adding quad multiplication, matrix<->quat conversions and creating quat from axix and theta
8/25 Basic Quaternion structure is implemented and all operator overloading is done

*/

/*  ---- REFERENCE 
Virtual trackballs revisited - Henriksen, K.; Sporring, J.; Hornbaek, K.;
Visualization and Computer Graphics, IEEE Transactions on
Volume 10,  Issue 2,  Mar-Apr 2004 Page(s):206 - 216
Digital Object Identifier 10.1109/TVCG.2004.1260772 

Web Reference:
Exhaustive listing of quaternion reference:
http://home.att.net/~t.a.ell/QuatRef.htm

Basic intro article with animated images: www.darwin3d.com/gamedev/articles/col0398.pdf

All math operators : http://en.wikipedia.org/wiki/Quaternion
*/

// Author : Ketan Mehta / km223@msstate.edu / 8/25/05

// Quaternion is defined as Q = (u,v) -> ( r, (i,j,k) ) such that r is scalar part and i,j,k
// are complex part of the quaternion.

#ifndef _QUATERNION_H_
#define _QUATERNION_H_

#include <iostream>
#include <cmath>
#include "Quaternion.h"

#define DTOR            0.0174532925
#define RTOD            57.2957795
#define TWOPI           6.283185307179586476925287
#define PI              3.141592653589793238462643
#define PID2            1.570796326794896619231322
#define PID3            1.047197551196397746154214
#define PID4            0.785398163397448309615660
#define PIPI            9.869604401089358618834491
#define EPSILON         0.001

class Quaternion {
private:
    float r_;
    float i_, j_, k_;

public:
    // Assignment functions
    Quaternion();
    ~Quaternion() {} ;
    
    Quaternion(const Quaternion &);
    Quaternion(float, float, float, float);
    Quaternion(float);

    void assign(float, float, float, float);
    Quaternion & Quaternion::operator=(const Quaternion &);

    // access functions
    const float  r(void) const ;
    const float  i(void) const ;
    const float  j(void) const ;
    const float  k(void) const ;

    // operator functions
    Quaternion& operator +=(const Quaternion&);
    Quaternion& operator -=(const Quaternion&);
    Quaternion& operator *=(float);
    Quaternion& operator /=(float);

    Quaternion operator +(const Quaternion&) const;
    Quaternion operator -(const Quaternion&) const;
    Quaternion operator -() const;
    Quaternion operator *(const float);
    Quaternion operator /(float) const;
    Quaternion operator *(const Quaternion&);

    // basic math operations
    float length() const;
    Quaternion & normalize();
    float norm() const;
    Quaternion inv() const;
    Quaternion log(void) const;
    Quaternion exp(void) const;
    Quaternion conjugate();

    // mathematical operations
    float dot(const Quaternion& ) const;
    void ToRotationMatrix(float m[][4]);
    void FromRotationMatrix(const float m[][4]);
    
        // Create quaternion from the axis and the angle    
    friend Quaternion Axis2Quat(float, float v[]);
    friend Quaternion slerp(Quaternion , Quaternion , float);    
    friend Quaternion squad(Quaternion& , Quaternion& , Quaternion& , Quaternion& , float);

    // input , output functions
    friend std::ostream & operator<<(std::ostream &, const Quaternion &);
    friend std::istream & operator>>(std::istream &, Quaternion &);
};

        
// Assignment functions

inline Quaternion::Quaternion(void)
{
        r_ = 1;
        i_ = j_ = k_ = 0;
}

inline Quaternion::Quaternion(float complex)
{
        r_ = 1;
        i_ = j_ = k_ = complex; 
}

inline Quaternion::Quaternion(float r, float i, float j, float k)
{
        r_ = r;   i_ = i;    j_ = j;    k_ = k;
}

inline Quaternion::Quaternion(const Quaternion & q)
{
        r_ = q.r();
        i_ = q.i(); j_ = q.j() ; k_ = q.k();
}

inline void Quaternion::assign(float r, float x, float y, float z)
{
        r_ = r; i_ = x; j_ = y; k_ = z;
}
inline Quaternion& Quaternion::operator=(const Quaternion & q)
{
        r_ = q.r();  i_ = q.i();  j_ = q.j();  k_ = q.k();
        return *this;
}

// Access functions
inline const float  Quaternion::r(void) const {
        return r_;
}
inline const float  Quaternion::i(void) const {
        return i_;
}
inline const float  Quaternion::j(void) const {
        return j_;
}
inline const float  Quaternion::k(void) const {
        return k_;
}

// operator functions
inline Quaternion Quaternion::operator+(const Quaternion& q_rhs) const
{
        return Quaternion(r_ + q_rhs.r(), i_ + q_rhs.i(), j_ + q_rhs.j(), k_ + q_rhs.k());
}

inline Quaternion Quaternion::operator-(const Quaternion& q_rhs) const
{
        return Quaternion(r_ - q_rhs.r(), i_ - q_rhs.i(), j_ - q_rhs.j(), k_ - q_rhs.k());
}

// Create the new object and return it, instead of changing sign for given object. Flexibility !?
inline Quaternion Quaternion::operator-() const
{        
         return Quaternion(-r_, -i_, -j_, -k_);        
}
 
inline Quaternion& Quaternion::operator *=(float scale)
{
        r_ *= scale; 
        i_ *= scale;    j_ *= scale;    k_ *= scale;
        return *this;        
}
inline Quaternion Quaternion::operator *(const float scale)
{
        return Quaternion(r_*scale, i_*scale, j_*scale, k_*scale);
}

inline Quaternion& Quaternion::operator /=(float scale)
{
        r_ /= scale; 
        i_ /= scale;    j_ /= scale;    k_ /= scale;
        return *this;        
}
inline Quaternion Quaternion::operator /(float scale) const
{
        return Quaternion(r_/scale, i_/scale, j_/scale, k_/scale);
}

inline Quaternion& Quaternion::operator+=(const Quaternion& q_rhs)
{
        r_ += q_rhs.r(); 
        i_ += q_rhs.i(); j_ += q_rhs.j(); k_ += q_rhs.k();
        return *this;
}
inline Quaternion& Quaternion::operator-=(const Quaternion& q_rhs)
{
        r_ -= q_rhs.r(); 
        i_ -= q_rhs.i(); j_ -= q_rhs.j(); k_ -= q_rhs.k();
        return *this;
}
       
// mathematical operations
inline float Quaternion::length() const
{
        return sqrt(r_*r_ + i_*i_ + j_*j_ + k_*k_);
}

inline Quaternion& Quaternion::normalize()
{
       float len = length();
       if ( len > 0 )
              *this/=len; 
       return *this;
}

// Norm is defined the sum of squares in the Ken Shoemake's Quaternion tutorial
// which is different from that defined in the Realtime rendering from Moller.
// Here we return the sum of the squares.
float Quaternion::norm() const
{        
        return (r_*r_ + i_*i_ + j_*j_ + k_*k_);
//         return sqrt(r_*r_ + i_*i_ + j_*j_ + k_*k_);
}

Quaternion Quaternion::conjugate()
{
        return Quaternion(r_, -i_, -j_, -k_);
}

// Ref : Realtime rendering by Moller , pg 45, eq : 3.29
// Layout followed is different.
Quaternion Quaternion::operator*(const Quaternion& q_rhs)
{
        return Quaternion(
                r_ * q_rhs.r() - i_ * q_rhs.i() - j_ * q_rhs.j() - k_ * q_rhs.k() ,
                r_ * q_rhs.i() + i_ * q_rhs.r() + j_ * q_rhs.k() - k_ * q_rhs.j() ,
                r_ * q_rhs.j() + j_ * q_rhs.r() + k_ * q_rhs.i() - i_ * q_rhs.k() ,
                r_ * q_rhs.k() + k_ * q_rhs.r() + i_ * q_rhs.j() - j_ * q_rhs.i() 
        ).normalize();
}

float Quaternion::dot(const Quaternion& q_rhs) const
{
        return (  r_ * q_rhs.r()
                + i_ * q_rhs.i()
                + j_ * q_rhs.j()
                + k_ * q_rhs.k() );
        
}

// Quaternion log(q) returns theta*v, where v is the vector component of q
// and theta is the angle.. e^(theta * v ) - is unit quaternion
Quaternion Quaternion::log(void) const
{
        float len = length();
        if ( len < 1e-6 )
                return Quaternion(0.0, i_, j_, k_);
        else
        {                
                float coef = acos(r_)/len;                
                return Quaternion(0.0, i_*coef, j_*coef, k_*coef);
        }
}

// return exponential of the quaternion
inline Quaternion Quaternion::inv(void) const
{
        float coef = 1.0/norm();
        return Quaternion(r_*coef, -i_*coef, -i_*coef, -i_*coef);
}
 
// exp(q) = (cos(theta), sin(theta)*v), v is the vector component a
// ref: http://en.wikipedia.org/wiki/Quaternion - has the formula for general quaternion
// exp(q) = exp(r)(cos(theta), v*(sin(theta)/theta)) , where theta = sqrt(vx*vx + vy*vy + vz+vz)    // Made a friend function for the ease of invocation. Will see whether it makes sense NOT to make
// member function.   
Quaternion Quaternion::exp() const
{
        float theta = sqrt(i_*i_ + j_*j_ + k_*k_) ;
        
        if (theta < 1e-6 )
            return Quaternion(0.0, i_, j_, k_);
        else
        {
            float coef = sin(theta)/ theta;
            return Quaternion(cos(theta), coef*i_ ,coef*j_ ,coef*k_);
        }
}

// Given and axis vector and rotation angle compute the quaternion
Quaternion Axis2Quat(float theta, float v[3])
{
        float halfAng = 0.5 * theta ;
        
        // make sure that axix vector is normalized.        
        float len = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        v[0] /= len; v[1] /= len; v[2] /= len; 
        
        // scale the vector with sin(theta/2)
        float scale = sin(halfAng);
        v[0] *= scale;   v[1] *= scale;   v[2] *= scale;   
      
        return Quaternion(cos(halfAng), v[0], v[1], v[2]);
}

// Build the Rotation matrix given the quaternion 
// Ref : Realtime Rendering , page 48, eq : 3.44
// We are including the implementation of generic quaternion, instead of the unit quat.
// OpenGL is Column Major , here we have computed values as in equation. Going row wise.
//    and copied into the 2D array in column major format. Hence directly consumable by
//    OpenGL.

void Quaternion::ToRotationMatrix(float m[][4])
{

        float s = 2.0 / length();
 
        m[0][0] = 1 - s * ( j_*j_ + k_*k_ ) ;
        m[1][0] = s * ( i_*j_ - r_*k_ ) ;
        m[2][0] = s * ( i_*k_ + r_*j_ ) ;
        m[3][0] = 0.0 ;
        
        m[0][1] = s * ( i_*j_ + r_*k_ ) ;        
        m[1][1] = 1 - s * ( i_*i_ + k_*k_ ) ;
        m[2][1] = s * ( j_*k_ - r_*i_ ) ;
        m[3][1] = 0.0 ;
        
        m[0][2] = s * ( i_*k_ - r_*j_ ) ;
        m[1][2] = s * ( j_*k_ + r_*i_ ) ;
        m[2][2] = 1 - s * ( i_*i_ + j_*j_ ) ;
        m[3][2] = 0.0 ;
        
        m[0][3] = 0.0 ;        
        m[1][3] = 0.0 ;
        m[2][3] = 0.0 ;
        m[3][3] = 1.0 ;
      
//  std::cout << m[0][0]<< " " << m[0][1]<< " " << m[0][2]<< " " << m[0][3]<< std::endl;
//    std::cout << m[1][0]<< " " << m[1][1]<< " " << m[1][2]<< " " << m[1][3]<< std::endl;
//    std::cout << m[2][0]<< " " << m[2][1]<< " " << m[2][2]<< " " << m[2][3]<< std::endl;
//    std::cout << m[3][0]<< " " << m[3][1]<< " " << m[3][2]<< " " << m[3][3]<< std::endl 
//    << std::endl;
}

// This function is adopted from the Material from:
// Magic Software, Inc.
// http://www.magic-software.com
// Copyright (c) 2000, All Rights Reserved
//

// Code genelogy : First code base routines are given in following paper by Shoemake
// and then adopted by Watt and Watt and then by Dave Eberly (Wild Magic) and from there
// into lot of game engines :)

// Note :: Refer tutorial by Shoemake abt why returned quaternion can be Q or -Q.
// www.cs.wisc.edu/graphics/Courses/cs-838-2002/Papers/quatut.pdf , pg 8
// --- page extract
/* This last routine converts a rotation matrix to a unit quaternion, but it 
may not be the same as the one with which you created the matrix. Part 2 of 
Theorem 1 implies that quaternions are homogeneous coordinates for rotations. 
Thus q and  q give the same rotation matrix, and the extraction routine can 
return either one; sometimes the choice matters. Remember that SO(3) has the 
topology of 3-dimensional real projective space (RP3), but unit quaternions 
form a hypersphere (S3) in four dimensions, which is topologically different. 
The identification of opposite points is what makes up the difference.
*/

void Quaternion::FromRotationMatrix(const float kRot[][4] )
{
        // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
        // article "Quaternion Calculus and Fast Animation".
        
        float fTrace = kRot[0][0]+kRot[1][1]+kRot[2][2];
        float fRoot;
        
//         std::cout << kRot[0][0] << " " << kRot[0][1] << " " << kRot[0][2]  << std::endl;
//         std::cout << kRot[1][0] << " " << kRot[1][1] << " " << kRot[1][2]  << std::endl;
//         std::cout << kRot[2][0] << " " << kRot[2][1] << " " << kRot[2][2]  << std::endl;

        
        if ( fTrace > 0.0 )
        {
                // |w| > 1/2, may as well choose w > 1/2
                fRoot = sqrt(fTrace + 1.0);  // 2w
                r_ = 0.5*fRoot;
                
                fRoot = 0.5/fRoot;  // 1/(4w) == 1/ 4* sqrt(1/4(trace))
                i_ = (kRot[2][1]-kRot[1][2])*fRoot;
                j_ = (kRot[0][2]-kRot[2][0])*fRoot;
                k_ = (kRot[1][0]-kRot[0][1])*fRoot;
                
                return;
        }
        else
        {
                // |w| <= 1/2
                static size_t s_iNext[3] = { 1, 2, 0 };
                size_t i = 0;
                if ( kRot[1][1] > kRot[0][0] )   i = 1;
                if ( kRot[2][2] > kRot[i][i] )   i = 2;
                
                size_t j = s_iNext[i];
                size_t k = s_iNext[j];
        
                fRoot = sqrt(kRot[i][i] - kRot[j][j] - kRot[k][k] + 1.0);
                
                // Ease of addressing
                float* apkQuat[3] = { &i_, &j_, &k_ };
                
                *apkQuat[i] = 0.5*fRoot;
                fRoot = 0.5/fRoot;
                
                r_ = (kRot[j][k] - kRot[k][j]) * fRoot;
                *apkQuat[j] = (kRot[i][j] + kRot[j][i]) * fRoot;
                *apkQuat[k] = (kRot[i][k] + kRot[k][i]) * fRoot;
                
                return;
        }
}
    
// Spherical linear interpolation. Taking two unit quaternions on the unit sphere, slerp
//      interpolates between the two quats.
// Implemented as Eq 3.52 in Realtime rendering book.
Quaternion slerp(Quaternion q, Quaternion r, float t)
{
        float scaleQ, scaleR;
        
        // Ensure t is in the range [0,1]
        if ( t < 0 ) t = 0; 
        else if ( t > 1 ) t = 1;
        
        // Check for unit quaternion - remove it later as it should not have been here..
        if ( q.length() != 1 ) q.normalize();
        if ( r.length() != 1 ) r.normalize();
        
        float cos_theta = q.dot(r);  
        float theta = acos(cos_theta);      
        float invSin = 1.0 / sin(theta) ;
        
        // Check for inverting the rotation
        Quaternion val;
        
        // Travel along the shorter path. Ref : Adv Anim. by Watt & Watt
        if ( (1.0 + cos_theta) > EPSILON )
        {
                // If angle is not small use SLERP.
                if ( (1.0 - cos_theta) > EPSILON )
                {
                        scaleQ = sin( (1.0 -t)*theta ) * invSin ;
                        scaleR = sin(t*theta) * invSin;
                }
                else    // For small angles use LERP
                {       
//                         std::cout << " are we using LERP " << std::endl;
                        scaleQ = 1.0 - t;
                        scaleR = t;
                }
                val = q * scaleQ  + r * scaleR ;
        }
        else // This is a long way
        {
                // Clear the concept later...
                val.assign(r.k(), -r.j(), r.i(), -r.r());
                scaleQ = sin( (1.0 - t)*PID2 );
                scaleR = sin( t * PID2 );
                val*=scaleR;
                q*=scaleQ;
                val +=val;
        }
        val.normalize();
        return val;
}
     
// Perform symmetric quadric interpolation.
// Refer for graphcial explanation Watt: Adv Animation and Rendering techniques.
Quaternion squad(Quaternion& q, Quaternion& s_1, Quaternion& s_2, Quaternion& r, float t)
{
        Quaternion sl1 = slerp(q,r,t);
        Quaternion sl2 = slerp(s_1,s_2,t);        
        return Quaternion(slerp(sl1,sl2, 2*t*(1-t)) );
}
        
// input , output functions
inline std::ostream & operator<<(std::ostream & os, const Quaternion & q)
{
        os << "w[" << q.r() << "], v[" << q.i() << ", " << q.j() << ", " << q.k() << "]";
        return os;
}

inline std::istream & operator >>(std::istream & in, Quaternion & q)
{
        float r, x, y, z;
        in >> r >> x >> y >> z ;
        q.assign(r,x,y,z); 
        return in;
}

// Helper code : Code to compile the class and test each calls.

// int main(int argc, char ** argv)
// {
// 
//         Quaternion q(1, 0.4, 0.6, 0.8), q1(0.45), q2(0.3, 0.1, 0.5, 0.9), q3(0.0), q4(0.5);
//         
//         std::cout << " q1 " << q1 << " q2 " << q2 << std::endl;
//         std::cout << " q1 * q2 " << q1 * q2 << std::endl;
//         
//         float vec[3] = { 0, 1, 0 };
//         Quaternion qtmp = Axis2Quat(45.0,vec);
//         
//         std::cout << " Quaternion with y axis and 45 deg rotation " << qtmp << std::endl;
//         
//         float mat[4][4];
//         std::cout << " matrix for above quad is " << std::endl;
//         Quat2RotMatrix(q,mat);
//         std::cout << " quaternion from matrix is " << std::endl;
//         Quaternion qfm.FromRotationMatrix();
//         std::cout << qfm << std::endl;
//         
//         q3 = q1 + q2;
//         q4 = q3*0.75;
//         std::cout << "squad(..)  " << squad(q1,q2, q3,q4,0.8) << std::endl;
//         
//         std::cout << " q3 " << q3 << std::endl;
//         std::cout << " log(q3) " << q3.log() << std::endl;
//         std::cout << " exp(q3) " << q3.exp() << std::endl;
//         std::cout << " inv(q3) " << q3.inv() << std::endl;        
// 
// }

#endif // _QUATERNION_H_
