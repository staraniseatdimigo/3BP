/* IGRam-1Gst */
/* JinEunpa */
/* 2016-07-25 PM 21:02 */
#include "phi_math.h"
#include "phi_type.h"

#include <stdlib.h>
#include <stdio.h>

using namespace Phi;

#define G 1

class Planet {
public:
	double mass;
	double r;
	Vec3 p, v;

	Vec3 temp_force;

	Planet(double _mass, double _r, Vec3 _p, Vec3 _v) {
		mass = _mass;
		r = _r;
		p = _p;
		v = _v;
	}

	virtual ~Planet() {
	}

	void initializeTempForce() {
		temp_force(0, 0, 0);
	}

	void calculateTempForce(Planet *other, double timestamp) {
		Vec3 dir = ((*other).p - p);
		double dist = dir.length();
		dir.normalize();
		double m = other->mass * G;
		double fa = m / (dist * dist);
		dir *= fa;
		temp_force += dir;
	}

	void apply(double timestamp) {
		v += temp_force * timestamp;
		p += v * timestamp;
	}
};

#define PLANET_N 3
#define timestamp 1.0

Planet *planets[PLANET_N];

// Class added /////////////////////////////////////////////////////

class Experiment {
public:
    int[4] exGroup; // 실험그룹
    double dT;          // dT
    // 물리량 허용범위(질량, 초기위치, 반지름, 초기속도) + 추가바람
    // 가동시작시각(외부) 리눅스 코드 + 추가바람
    double T; // 소요시간
    double GT; // 그룹당 시간
    double unitT; // 1단위당 시간
    double maxT; // 안정판단한계시간
    char[15] drafter; // 작성자?
};

class Try {
public:
    Planet *planets[PLANET_N]; // 시행 초기정보
    /*
    double avgV;
    double avgP;
    Vec3 boonpoP;
    Vec3 boonpoV;
    */
};

class Result {
public:
    bool iscollided;
    int[2] cdPair;
    double cdTime; // 기준 - 시작 실험 시각(외부 시각이 아님)
    double[3] cdPoint; // 이거 좀 수정해줘라 충돌위치

    public Result(bool _iscollided, int[2] _cdPair, double _cdTime, double[3] _cdPoint) {
        iscollided = _iscollided;
        cdPair = _cdPair
        cdTime = _cdTime;
        cdPoint = _cdPoint;
    }

};

////////////////////////////////////////////////////////////////////////

void loop() {
	int running = 1;
	int i, j;
	while(running) {

            // collision detection
        for(i=1;i<PLANET_N;i++) {
			for(j=0;j<i;j++) {
				if(planets[i]->r + planets[j]->r <= Vec3::length(planets[i]->p - planets[j]->p))
				{
				    //collision detected

				}
			}
		}

        // Initialize
		for(i=0;i<PLANET_N;i++) {
			planets[i]->initializeTempForce();
		}
		// Calculation
		for(i=1;i<PLANET_N;i++) {
			for(j=0;j<i;j++) {
				planets[i]->calculateTempForce(planets[j], timestamp);
				planets[j]->calculateTempForce(planets[i], timestamp);
			}
		}
		// Render
		for(i=0;i<PLANET_N;i++) {
			planets[i]->apply(timestamp);
		}
	}
}

void input() {
	int i;
	double mass;
	Vec3 p, v;

	FILE *f = stdin;
	for(i=0;i<3;i++) {
		fscanf(stdin, "%lf", &mass);
		fscanf(stdin, "%lf %lf %lf", &p.x, &p.y, &p.z);
		fscanf(stdin, "%lf %lf %lf", &v.x, &v.y, &v.z);
		planets[i] = new Planet(mass, p, v);
	}
}

int main() {
	loop();
	return 0;
}

