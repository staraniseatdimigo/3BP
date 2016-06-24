#ifndef __PHI_TYPE_H__
#define __PHI_TYPE_H__

#include "phi_math.h"

namespace Phi{

typedef double Vec1;
class Vec2;
typedef Vec2 Pos2;
class Vec3;
class Vec4;

class Mat2;
class Mat3;
class Mat4;

class LinearForce2;

class Vec2{
public:
	double x, y;
	
	static Vec2 gen(double, double);

	Vec2();
	Vec2(double, double);
	
	~Vec2();

	Vec2 operator - ();

	double operator + ();
	double length();

	double invLength();

	Vec2& normalize();

	Vec2& operator = (Vec2);
	Vec2& operator () (Vec2);
	Vec2& operator () (double, double);
	
	Vec2& operator += (Vec2);
	Vec2 operator + (Vec2);

	Vec2& operator -= (Vec2);
	Vec2 operator - (Vec2);

	Vec2& operator *= (double);
	Vec2 operator * (double);

	Vec2& operator /= (double);
	Vec2 operator / (double);

	double operator * (Vec2);
	double dot(Vec2);
	double cross(Vec2);

	void productEach(Vec2);

	Vec2& rotate(double);
	Vec2& rotate(double, double);
	Vec2& rotateCW();
	Vec2& rotateCCW();

	double X(double);
	double setX(double);
	double X();
	double getX();

	double Y(double);
	double setY(double);
	double Y();
	double getY();

	void print();
};

class Vec3{
public:
	double x, y, z;
	
	static Vec3 gen(double, double, double);

	Vec3();
	Vec3(double, double, double);
	
	~Vec3();

	Vec3 operator - ();

	double operator + ();
	double length();

	double invLength();

	Vec3& normalize();

	Vec3& operator = (Vec3);
	Vec3& operator () (Vec3);
	Vec3& operator () (double, double, double);
	
	Vec3& operator += (Vec3);
	Vec3 operator + (Vec3);

	Vec3& operator -= (Vec3);
	Vec3 operator - (Vec3);

	Vec3& operator *= (double);
	Vec3 operator * (double);

	Vec3& operator /= (double);
	Vec3 operator / (double);

	double operator * (Vec3);
	double dot(Vec3);
	Vec3 cross(Vec3);

	void productEach(Vec3);

	double X(double);
	double setX(double);
	double X();
	double getX();

	double Y(double);
	double setY(double);
	double Y();
	double getY();
	
	double Z(double);
	double setZ(double);
	double Z();
	double getZ();

	void print();
};

class Vec4{
	double x, y, z, w;
};

class Mat3{
private:
	double v[9];
public:
	Mat3();
	Mat3(double[9]);
	Mat3(double[3][3]);
	Mat3(const Mat3&);
	~Mat3();

	Mat3& operator = (Mat3);
	Mat3& operator = (double[9]);
	Mat3& operator = (double[3][3]);
	Mat3& operator () (Mat3);
	Mat3& operator () (double[9]);
	Mat3& operator () (double[3][3]);
	
	Mat3& operator += (Mat3);
	Mat3 operator + (Mat3);

	Mat3& operator -= (Mat3);
	Mat3 operator - (Mat3);

	Mat3& operator *= (double);
	Mat3 operator * (double);

	Mat3& operator /= (double);
	Mat3 operator / (double);

	Mat3 operator * (Mat3);
	Mat3& operator *= (Mat3);
	Mat3 product(Mat3);

	Vec3 operator * (Vec3);
	Vec3 product(Vec3);

	double operator + ();
	double determinant();

	Mat3& transpose();

	Mat3& invert();

	double* get();
	double get(int);
	double get(int,int);

	void print();
};

// 가속도, 속도, 위치, 밀도, 부피, 운동량, 

struct MomentumSet1{
	Vec1 accel;
	Vec1 velocity;
	Vec1 position;
	Vec1 momentum;
};

struct MomentumSet2{
	Vec2 accel;
	Vec2 velocity;
	Vec2 position;
	Vec2 momentum;
};

class Momentum2{
protected:
	double density, volume;
	MomentumSet2 linear;
	MomentumSet1 angular;
	
	Vec2 temp_force;

public:
	Momentum2();
	~Momentum2();

	void setDensity(double);
	double getDensity();
	void setVolume(double);
	double getVolume();
	double getMass();

	void setPosition(Vec2);
	Vec2 getPosition();
	void movePosition(Vec2);

	void setVelocity(Vec2);
	Vec2 getVelocity();

	void setAngularAccel(Vec1);
	Vec1 getAngularAccel();

	void setAngularVelocity(Vec1);
	Vec1 getAngularVelocity();

	void setAngle(Vec1);
	Vec1 getAngle();

	void applyForce(Vec2);
	Vec2 getForce();
	void cleanForce();

	void update(double);



	static void calculateVelocityAfterCollision(double e, double m1, Vec2 p1, Vec2 v1, Vec1 ap1, Vec1 av1, Vec2 pc, double m2, Vec2 p2, Vec2 v2, Vec1 ap2, Vec1 av2, Vec2 vn, Vec1 vd, Vec2 *vla, Vec1 *vaa);
};

Mat3 calculateInertiaTensor(double m, Vec3 r);
Vec2 calculateLinearVelocity(Vec2 lv, Vec2 r, Vec1 av);

Vec2 calculateLinearVelocityByImpulse(Vec2 v, double m, Vec3 impulse, double);
Vec1 calculateAngularVelocityByImpulse(Vec1 v, Vec2 r, Mat3 invI, Vec3 impulse, double);

Vec3 calculateImpulse(Vec1 e, Vec1 c, Vec2 n, Mat3 invI1, Mat3 invI2, Vec2 r1, Vec2 r2, Vec1 m1, Vec1 m2);
 

template <typename T>
struct List{
	T value;
	List<T> *next, *pre;

	List(T _value){
		value = _value;
		next = pre = this;
	}
	List(T _value, List<T> *p){
		value = _value;
		if(p){
			next = p->next;
			pre = p;
			p->next->pre = this;
			p->next = this;
		} else {
			next = pre = this;
		}
	}
	~List(){
		next->pre = pre;
		pre->next = next;
	}

	struct Iter{
		List<T> *p;

		T& operator *(){ return p->value; }
		void operator =(List<T> *o){ p=o; }
		void operator ++(){ p=p->next; }
		void operator ++(int){ p=p->next; }
		void operator --(){ p=p->pre; }
		void operator --(int){ p=p->pre; }
		bool operator ==(List<T> *o){ return o==p; }
		bool operator !=(List<T> *o){ return o!=p; }

		List<T>* get(){ return p; }
	};
};

}

#endif
