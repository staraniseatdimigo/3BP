/* IGRam-1Gst */
/* JinEunpa */
/* 2016-07-25 PM 21:02 */
#include "phi_math.h"
#include "phi_type.h"

#include <stdlib.h>
#include <stdio.h>

using namespace Phi;

/* Gravity Constant */
#define G 1

class Planet {
public:
	double mass;
	double r;
	Vec3 p, v;

	Vec3 tempForce;

	Planet(Planet *other) {
		mass = other->mass;
		r = other->r;
		p = other->p;
		v = other->v;
	}

	Planet(double _mass, double _r, Vec3 _p, Vec3 _v) {
		mass = _mass;
		r = _r;
		p = _p;
		v = _v;
	}

	virtual ~Planet() {
	}

	void initializeTempForce() {
		tempForce(0, 0, 0);
	}

	void calculateTempForce(Planet *other, double timestamp) {
		Vec3 dir = ((*other).p - p);
		double dist = dir.length();
		dir.normalize();
		double m = other->mass * G;
		double fa = m / (dist * dist);
		dir *= fa;
		tempForce += dir;
	}

	void apply(double timestamp) {
		v += tempForce * timestamp;
		p += v * timestamp;
	}

	bool isCollided(Planet *other) {
		double distance = (this->p - other->p).length();
		double radiusSum = this->r + other->r;
		return distance <= radiusSum;
	}
};

#define PLANET_N 3
#define timestamp 1.0
#define timeout 10000000.0

Planet *planets[PLANET_N];


class Experiment {
public:
    int exGroup[4]; // 실험그룹
    double dT;          // dT
    // 물리량 허용범위(질량, 초기위치, 반지름, 초기속도) + 추가바람
    double T; // 소요시간
    double GT; // 그룹당 시간
    double unitT; // 1단위당 시간
    double maxT; // 안정판단한계시간
    char drafter[15]; // 작성자?
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

	void init(Planet **P) {
		for(int i=0;i<PLANET_N;i++) {
			this->planets[i] = new Planet(P[i]);
		}
	}
};

class CollisionResult {
public:
    int cdPair[2];
    double cdTime; // 기준 - 시작 실험 시각(외부 시각이 아님)
    Vec3 cdPoint; // 이건 고민좀

    CollisionResult(int _cdPair[], double _cdTime, Vec3 _cdPoint) {
        cdPair[0] = _cdPair[0];
		cdPair[1] = _cdPair[1];
        cdTime = _cdTime;
        cdPoint = _cdPoint;
    }
};

class Result {
public:
	Experiment E;
	Try T;

	/* If CollisionResult == NULL, No Collision */
	CollisionResult *C;

	Result()
		: E(), T() {
		C = NULL;
	}
};


/* Main Logic */
void loop(Result *R) {
	double currentTime = 0.0;
	int i, j;

	/* Backup Initial Status */
	R->T.init(planets);

	/* Main Loop */
	while(currentTime < timeout) {
		/* Check Collision Detection */
        for(i=1;i<PLANET_N;i++) {
			for(j=0;j<i;j++) {
				if(planets[i]->isCollided(planets[j])) {
					/* Is Collided */
					/* Write Collided Result */
					int cdPair[] = {i, j};
					Vec3 cp;
					R->C = new CollisionResult(cdPair, currentTime, cp);
					
					/* Forced End Loop */
					goto LOOP_END;
				}
			}
		}

        /* Initialize */
		for(i=0;i<PLANET_N;i++) {
			planets[i]->initializeTempForce();
		}
		
		/* Calculation */
		for(i=1;i<PLANET_N;i++) {
			for(j=0;j<i;j++) {
				planets[i]->calculateTempForce(planets[j], timestamp);
				planets[j]->calculateTempForce(planets[i], timestamp);
			}
		}
		
		/* Update Status */
		for(i=0;i<PLANET_N;i++) {
			planets[i]->apply(timestamp);
		}

		/* Time is gone */
		currentTime += timestamp;
	}

	/* Finalize */
LOOP_END:
}



/* Initialize Experiment */
void initExp(Result *R) {
	
}

/* Get Initial Value from STDIN */
void input(Result *R) {
	int i;
	double mass, radius;
	Vec3 p, v;

	FILE *f = stdin;
	for(i=0;i<3;i++) {
		fscanf(stdin, "%lf %lf", &mass, &radius);
		fscanf(stdin, "%lf %lf %lf", &p.x, &p.y, &p.z);
		fscanf(stdin, "%lf %lf %lf", &v.x, &v.y, &v.z);
		planets[i] = new Planet(mass, radius, p, v);
	}
}

int main() {
	Result r;
	initExp(&r);
	input(&r);
	loop(&r);
	return 0;
}

