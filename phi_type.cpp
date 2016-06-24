#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phi_type.h"

using namespace Phi;

namespace Phi{

// class Vec2
Vec2 Vec2::gen(double _x, double _y){
	Vec2 v(_x, _y);
	return v;
}

Vec2::Vec2(){
	x = 0, y = 0;
}
Vec2::Vec2(double _x, double _y){
	x = _x, y = _y;
}
Vec2::~Vec2(){ ; }

Vec2 Vec2::operator - (){
	Vec2 v(-x, -y);
	return v;
}

double Vec2::operator + (){
	return fastSqrt(x*x + y*y);
}
double Vec2::length (){
	return fastSqrt(x*x + y*y);
}

double Vec2::invLength (){
	return fastInvSqrt(x*x + y*y);
}

Vec2& Vec2::normalize(){
	double n = fastInvSqrt(x*x + y*y);
	x *= n, y *= n;
	if(x == 0) y = 1;
	if(y == 0) x = 1;
	return *this;
}

Vec2& Vec2::operator = (Vec2 other){
	x = other.x, y = other.y;
	return *this;
}
Vec2& Vec2::operator () (Vec2 other){
	x = other.x, y = other.y;
	return *this;
}
Vec2& Vec2::operator () (double _x, double _y){
	x = _x, y = _y;
	return *this;
}

Vec2& Vec2::operator += (Vec2 other){
	x += other.x, y += other.y;
	return *this;
}
Vec2 Vec2::operator + (Vec2 other){
	Vec2 v;
	v.x = x + other.x;
	v.y = y + other.y;
	return v;
}

Vec2& Vec2::operator -= (Vec2 other){
	x -= other.x, y -= other.y;
	return *this;
}
Vec2 Vec2::operator - (Vec2 other){
	Vec2 v;
	v.x = x - other.x;
	v.y = y - other.y;
	return v;
}

Vec2& Vec2::operator *= (double value){
	x *= value, y *= value;
	return *this;
}
Vec2 Vec2::operator * (double value){
	Vec2 v;
	v.x = x * value;
	v.y = y * value;
	return v;
}

Vec2& Vec2::operator /= (double value){
	x /= value, y /= value;
	return *this;
}
Vec2 Vec2::operator / (double value){
	Vec2 v;
	v.x = x / value;
	v.y = y / value;
	return v;
}

double Vec2::X(double _v){ return x=_v; }
double Vec2::Y(double _v){ return y=_v; }
double Vec2::setX(double _v){ return x=_v; }
double Vec2::setY(double _v){ return y=_v; }
double Vec2::X(){ return x; }
double Vec2::Y(){ return y; }
double Vec2::getX(){ return x; }
double Vec2::getY(){ return y; }

double Vec2::operator * (Vec2 other){
	return x*other.x + y*other.y;
}
double Vec2::dot(Vec2 other){
	return x*other.x + y*other.y;
}
double Vec2::cross(Vec2 other){
	return x*other.y - y*other.x;
}

void Vec2::productEach(Vec2 v){
	x *= v.getX();
	y *= v.getY();
}

Vec2& Vec2::rotate(double angle){
	double tx=x, ty=y;
	double s=sin(angle*IRAD), c=cos(angle*IRAD);
	x = tx*c - ty*s;
	y = tx*s + ty*c;
	return *this;
}

Vec2& Vec2::rotate(double s, double c){
	double tx=x, ty=y;
	x = tx*c - ty*s;
	y = tx*s + ty*c;
	return *this;
}

Vec2& Vec2::rotateCW(){
	double t = -x;
	x = y, y = t;
	return *this;
}
Vec2& Vec2::rotateCCW(){
	double t = -y;
	y = x, x = t;
	return *this;
}

void Vec2::print(){
	printf("(%g, %g)", x, y);
}

// class Vec3
Vec3 Vec3::gen(double _x, double _y, double _z){
	Vec3 v(_x, _y, _z);
	return v;
}

Vec3::Vec3(){
	x = 0, y = 0, z = 0;
}
Vec3::Vec3(double _x, double _y, double _z){
	x = _x, y = _y, z = _z;
}
Vec3::~Vec3(){ ; }

Vec3 Vec3::operator - (){
	Vec3 v(-x, -y, -z);
	return v;
}

double Vec3::operator + (){
	return fastSqrt(x*x + y*y + z*z);
}
double Vec3::length (){
	return fastSqrt(x*x + y*y + z*z);
}

double Vec3::invLength (){
	return fastInvSqrt(x*x + y*y + z*z);
}

Vec3& Vec3::normalize(){
	double n = fastInvSqrt(x*x + y*y + z*z);
	x *= n, y *= n, z *= n;
	if(x == 0 && y == 0) z = 1;
	if(x == 0 && z == 0) y = 1;
	if(y == 0 && z == 0) x = 1;
	return *this;
}

Vec3& Vec3::operator = (Vec3 other){
	x = other.x, y = other.y, z = other.z;
	return *this;
}
Vec3& Vec3::operator () (Vec3 other){
	x = other.x, y = other.y, z = other.z;
	return *this;
}
Vec3& Vec3::operator () (double _x, double _y, double _z){
	x = _x, y = _y, z = _z;
	return *this;
}

Vec3& Vec3::operator += (Vec3 other){
	x += other.x, y += other.y, z += other.z;
	return *this;
}
Vec3 Vec3::operator + (Vec3 other){
	Vec3 v;
	v.x = x + other.x;
	v.y = y + other.y;
	v.z = z + other.z;
	return v;
}

Vec3& Vec3::operator -= (Vec3 other){
	x -= other.x, y -= other.y, z -= other.z;
	return *this;
}
Vec3 Vec3::operator - (Vec3 other){
	Vec3 v;
	v.x = x - other.x;
	v.y = y - other.y;
	v.z = z - other.z;
	return v;
}

Vec3& Vec3::operator *= (double value){
	x *= value, y *= value, z *= value;
	return *this;
}
Vec3 Vec3::operator * (double value){
	Vec3 v;
	v.x = x * value;
	v.y = y * value;
	v.z = z * value;
	return v;
}

Vec3& Vec3::operator /= (double value){
	x /= value, y /= value, z /= value;
	return *this;
}
Vec3 Vec3::operator / (double value){
	Vec3 v;
	v.x = x / value;
	v.y = y / value;
	v.z = z / value;
	return v;
}

double Vec3::X(double _v){ return x=_v; }
double Vec3::Y(double _v){ return y=_v; }
double Vec3::Z(double _v){ return z=_v; }
double Vec3::setX(double _v){ return x=_v; }
double Vec3::setY(double _v){ return y=_v; }
double Vec3::setZ(double _v){ return z=_v; }
double Vec3::X(){ return x; }
double Vec3::Y(){ return y; }
double Vec3::Z(){ return z; }
double Vec3::getX(){ return x; }
double Vec3::getY(){ return y; }
double Vec3::getZ(){ return z; }

double Vec3::operator * (Vec3 other){
	return x*other.x + y*other.y + z*other.z;
}
double Vec3::dot(Vec3 other){
	return x*other.x + y*other.y + z*other.z;
}
Vec3 Vec3::cross(Vec3 other){
	double _x=other.getX(),_y=other.getY(),_z=other.getZ();
	Vec3 v(y*_z-z*_y,z*_x-x*_z,x*_y-y*_x);
	return v;
}

void Vec3::productEach(Vec3 v){
	x *= v.getX();
	y *= v.getY();
	z *= v.getZ();
}

void Vec3::print(){
	printf("(%g, %g, %g)", x, y, z);
}

	Mat3::Mat3(){
		int i;
		for(i=0;i<9;i++) v[i] = 0;
	}
	Mat3::Mat3(double t[9]){
		memcpy(v, t, sizeof(double[9]));
	}
	Mat3::Mat3(double t[3][3]){
		int i, j;
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				v[i*3+j] = t[i][j];
	}
	Mat3::Mat3(const Mat3 &m){
		memcpy(this, &m, sizeof(Mat3));
	}

	Mat3::~Mat3(){ }

	Mat3& Mat3::operator = (Mat3 m){
		memcpy(this, &m, sizeof(Mat3));
		return *this;
	}
	Mat3& Mat3::operator = (double t[9]){
		memcpy(v, t, sizeof(double[9]));
		return *this;
	}
	Mat3& Mat3::operator = (double t[3][3]){
		int i, j;
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				v[i*3+j] = t[i][j];
		return *this;
	}

	Mat3& Mat3::operator () (Mat3 m){
		memcpy(this, &m, sizeof(Mat3));
		return *this;
	}
	Mat3& Mat3::operator () (double t[9]){
		memcpy(v, t, sizeof(double[9]));
		return *this;
	}
	Mat3& Mat3::operator () (double t[3][3]){
		int i, j;
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				v[i*3+j] = t[i][j];
		return *this;
	}

	Mat3& Mat3::operator += (Mat3 m){
		double *_v = m.get();
		int i;
		for(i=0;i<9;i++) v[i] += _v[i];
		return *this;
	}
	Mat3 Mat3::operator + (Mat3 m){
		double *_v = m.get(), r[9];
		int i;
		for(i=0;i<9;i++) r[i] = v[i] + _v[i];
		Mat3 mr(r);
		return mr;
	}

	Mat3& Mat3::operator -= (Mat3 m){
		double *_v = m.get();
		int i;
		for(i=0;i<9;i++) v[i] -= _v[i];
		return *this;
	}
	Mat3 Mat3::operator - (Mat3 m){
		double *_v = m.get(), r[9];
		int i;
		for(i=0;i<9;i++) r[i] = v[i] - _v[i];
		Mat3 mr(r);
		return mr;
	}

	Mat3& Mat3::operator *= (double d){
		int i;
		for(i=0;i<9;i++) v[i] *= d;
		return *this;
	}
	Mat3 Mat3::operator * (double d){
		double r[9];
		int i;
		for(i=0;i<9;i++) r[i] = v[i] * d;
		Mat3 mr(r);
		return mr;
	}

	Mat3& Mat3::operator /= (double d){
		int i;
		for(i=0;i<9;i++) v[i] /= d;
		return *this;
	}
	Mat3 Mat3::operator / (double d){
		double r[9];
		int i;
		for(i=0;i<9;i++) r[i] = v[i] / d;
		Mat3 mr(r);
		return mr;
	}

	Mat3 Mat3::operator * (Mat3 m){
		double *_v = m.get(), r[9];
		int i, j;
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				r[i+j*3] = v[0+j*3]*_v[i+0*3] + v[1+j*3]*_v[i+1*3] + v[2+j*3]*_v[i+2*3];
			}
		}
		Mat3 mr(r);
		return mr;
	}
	
	Mat3& Mat3::operator *= (Mat3 m){
		double *_v = m.get(), r[9];
		int i, j;
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				r[i+j*3] = v[0+j*3]*_v[i+0*3] + v[1+j*3]*_v[i+1*3] + v[2+j*3]*_v[i+2*3];
			}
		}
		(*this) = r;
		return *this;
	}

	Vec3 Mat3::operator * (Vec3 _v){
		double x=_v.getX(), y=_v.getY(), z=_v.getZ();
		Vec3 r;
		r.setX(x*v[0+0*3]+y*v[0+1*3]+z*v[0+2*3]);
		r.setY(x*v[1+0*3]+y*v[1+1*3]+z*v[1+2*3]);
		r.setZ(x*v[2+0*3]+y*v[2+1*3]+z*v[2+2*3]);
		return r;
	}

	Vec3 Mat3::product(Vec3 _v){
		double x=_v.getX(), y=_v.getY(), z=_v.getZ();
		Vec3 r;
		r.setX(x*v[0+0*3]+y*v[0+1*3]+z*v[0+2*3]);
		r.setY(x*v[1+0*3]+y*v[1+1*3]+z*v[1+2*3]);
		r.setZ(x*v[2+0*3]+y*v[2+1*3]+z*v[2+2*3]);
		return r;
	}

	double Mat3::operator + (){
		return v[0+0*3]*(v[1+1*3]*v[2+2*3]-v[1+2*3]*v[2+1*3])+v[0+1*3]*(v[1+2*3]*v[2+0*3]-v[1+0*3]*v[2+2*3])+v[0+2*3]*(v[1+0*3]*v[2+1*3]-v[1+1*3]*v[2+0*3]);
	}

	double Mat3::determinant(){
		return v[0+0*3]*(v[1+1*3]*v[2+2*3]-v[1+2*3]*v[2+1*3])+v[0+1*3]*(v[1+2*3]*v[2+0*3]-v[1+0*3]*v[2+2*3])+v[0+2*3]*(v[1+0*3]*v[2+1*3]-v[1+1*3]*v[2+0*3]);
	}

	Mat3& Mat3::transpose(){
		double t;
#define SWAP(A,B) t=(A),(A)=(B),(B)=t
		SWAP(v[1+0*3], v[0+1*3]);
		SWAP(v[2+0*3], v[0+2*3]);
		SWAP(v[2+1*3], v[1+2*3]);
#undef SWAP
		return *this;
	}

	Mat3& Mat3::invert(){
		double D = determinant();
		if(!D) return *this;
		double _v[9];
		int i, j;
#define I(a,b) (((a)%3)+((b)%3)*3)
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				_v[i+j*3] = v[I(i+1,j+1)]*v[I(i+2,j+2)]-v[I(i+2,j+1)]*v[I(i+1,j+2)];
			}
		}
		(*this) = _v;
		(*this) /= D;
		return transpose();
	}

	double* Mat3::get(){
		return v;
	}
	double Mat3::get(int i){
		return v[i];
	}
	double Mat3::get(int i, int j){
		return v[i+j*3];
	}

	void Mat3::print(){
		int i, j;
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				printf("%g ", v[i+j*3]);
			}
			if(i<2)printf("\n");
		}
	}

// class Momentum2
Momentum2::Momentum2(){
	density = volume = 0;

	linear.accel(0, 0);
	linear.velocity(0, 0);
	linear.position(0, 0);
	linear.momentum(0, 0);

	angular.accel = 0;
	angular.velocity = 0;
	angular.position = 0;
	angular.momentum = 0;

	temp_force(0, 0);
}
Momentum2::~Momentum2(){
}

void Momentum2::setDensity(double _value){
	density = _value;
}
double Momentum2::getDensity(){ return density; }
void Momentum2::setVolume(double _value){
	volume = _value;
}
double Momentum2::getVolume(){ return volume; }
double Momentum2::getMass(){ return density*volume; }

void Momentum2::setPosition(Vec2 v){
	linear.position = v;
}
Vec2 Momentum2::getPosition(){ return linear.position; }
void Momentum2::movePosition(Vec2 v){
	linear.position += v;
}

void Momentum2::setVelocity(Vec2 v){
	linear.velocity = v;
}
Vec2 Momentum2::getVelocity(){ return linear.velocity; }

void Momentum2::setAngularVelocity(Vec1 v){
	angular.velocity = v;
}
Vec1 Momentum2::getAngularVelocity(){ return angular.velocity; }

void Momentum2::setAngle(Vec1 v){
	angular.position = v;
	while(angular.position >= 360.0) angular.position -= 360.0;
	while(angular.position < 0.0) angular.position += 360.0;
}
Vec1 Momentum2::getAngle(){ return angular.position; }


void Momentum2::applyForce(Vec2 _vec){
	temp_force += _vec;
}
Vec2 Momentum2::getForce(){ return temp_force; }
void Momentum2::cleanForce(){
	temp_force(0, 0);
}

void Momentum2::update(double time_step){
	linear.accel = temp_force/getMass();
	linear.velocity += linear.accel*time_step;
	linear.position += linear.velocity*time_step;
	linear.momentum = linear.velocity*getMass();

	angular.accel = 0;
	angular.velocity += angular.accel*time_step;
	angular.position += angular.velocity*time_step;
	while(angular.position >= 360.0) angular.position -= 360.0;
	while(angular.position < 0.0) angular.position += 360.0;
	angular.momentum = angular.velocity*getMass();
	temp_force(0, 0);
}



void Momentum2::calculateVelocityAfterCollision(double e, double m1, Vec2 p1, Vec2 v1, Vec1 ap1, Vec1 av1, Vec2 pc, double m2, Vec2 p2, Vec2 v2, Vec1 ap2, Vec1 av2, Vec2 vn, Vec1 vd, Vec2 *vla, Vec1 *vaa){
	double angle = atan2(vn.getY(), vn.getX());
	double cosa = cos(angle), sina=sin(angle);
	/*double cosa1 = -sin(ap1*IRAD), sina1=cos(ap1*IRAD);
	  double cosa2 = -sin(ap2*IRAD), sina2=cos(ap2*IRAD);*/

	Vec2 vr1, vr2;

	vr1 = v1;//+Vec2::gen(cosa1,sina1)*av1*IRAD*(pc-p1).length();
	vr2 = v2;//+Vec2::gen(cosa2,sina2)*av2*IRAD*(pc-p2).length();

	double dv1 = v1.getX()*cosa + v1.getY()*sina;
	double dv2 = v2.getX()*cosa + v2.getY()*sina;

	double dva1x, dva1y;

	dva1x = (m1-m2*e)/(m1+m2)*dv1 + (m2+m2*e)/(m1+m2)*dv2;
	dva1y = v1.getY()*cosa - v1.getX()*sina;

	if(vla){
		(*vla)(
			dva1x*cosa - dva1y*sina,
			dva1x*sina + dva1y*cosa
			);
		vr1 = (vr1-*vla);
	}
	if(vaa){
		*vaa = 0;
		Vec1 delta = 0;
		Vec2 v_c = (pc-p1);
		Vec2 dv = vr1;
		Vec1 r = v_c.length();

		//axis;
		v_c.rotateCW();
		//Vec1 v = dv*v_c/(r*r*m1);
		Vec1 v = av1+(dv*v_c+vd)/(m1*r*r)*RAD;
		delta = v;

		*vaa = delta;
	}
}


Mat3 calculateInertiaTensor(double m, Vec3 r){
	// ry^2+rz^2, -rxry, -rxrz
	// -rxry, rz^2+rx^2, -ryrz
	// -rxrz, -ryrz, rx^2+ry^2
	double v[9];
	v[1] = v[3] = -r.X()*r.Y();
	v[2] = v[7] = -r.X()*r.Z();
	v[5] = v[8] = -r.Y()*r.Z();
	v[0] = r.Y()*r.Y()+r.Z()*r.Z();
	v[4] = r.Z()*r.Z()+r.X()*r.X();
	v[8] = r.X()*r.X()+r.Y()*r.Y();
	Mat3 mr(v);
	mr *= m;
	return mr;
}

Vec2 calculateLinearVelocity(Vec2 lv, Vec2 r, Vec1 av){
	Vec2 res = lv;
	Vec3 r3(r.X(), r.Y(), 0), av3(0, 0, av);
	r3 = av3.cross(r3);
	res += Vec2::gen(r3.X(), r3.Y());
	return res;
}

Vec2 calculateLinearVelocityByImpulse(Vec2 v, double m, Vec3 impulse, double i){
	Vec3 v3(v.X(), v.Y(), 0.0);
	v3 = (impulse / m)*i;
	v(v3.X(), v3.Y());
	return v;
}
Vec1 calculateAngularVelocityByImpulse(Vec1 v, Vec2 r, Mat3 invI, Vec3 impulse, double i){
	Vec3 v3(0.0, 0.0, v);
	Vec3 r3(r.X(), r.Y(), 0.0);
	v3 = (invI * (r3.cross(impulse)))*i;
	return v3.Z();
}

Vec3 calculateImpulse(Vec1 e, Vec1 c, Vec2 n, Mat3 invI1, Mat3 invI2, Vec2 r1, Vec2 r2, Vec1 m1, Vec1 m2){
    Vec3 j;
	Vec3 n3(n.X(), n.Y(), 0.0);
	Vec3 r13(r1.X(), r1.Y(), 0.0), r23(r2.X(), r2.Y(), 0.0);
	// j = (-(1+e)c)/(n((invI1)(r1Xn)Xr1) + n((invI2)(r2Xn)Xr2) + 1/m1 + 1/m2)
	j = n3 * ((c * (-1.0-e)) / (n3*(((invI1*(r13.cross(n3))).cross(r13))) + n3*(((invI2*(r23.cross(n3))).cross(r23))) + 1.0/m1 + 1.0/m2));
	return j;
}
	
}
