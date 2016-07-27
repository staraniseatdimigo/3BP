/* IGRam-1Gst */
/* JinEunpa */
/* 2016-07-27 PM 20:48 */
#include "phi_math.h"
#include "phi_type.h"

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>

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

Planet *planets[PLANET_N];


class Experiment {
public:
    int exGroup[4]; // 실험그룹
    double dT;          // dT
    // 물리량 허용범위
    
    double mLimit;
    double pLimit;
    double vLimit;
    double rLimit;
	
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
	Vec3 c_mass;

	double timestamp = R->E.dT;
	double timeout = R->E.maxT;

	/* Backup Initial Status */
	// R->T.init(planets);

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
		
		/* adjusting to center of mass */
		for(i=0;i<PLANET_N;i++) {
			c_mass += planets[i]->p * planets[i]->mass;
		}
		
		for(i=0;i<PLANET_N;i++) {
			planets[i]->p = planets[i]->p - c_mass;
		}

		/* Check Velocity Limit */
		for(i=0;i<PLANET_N;i++) {
		}

		/* Time is gone */
		currentTime += timestamp;
	}

	/* Finalize */
LOOP_END:
	;
}



/* Initialize Experiment */
void initExp(Result *R, FILE *F) {
	fscanf(F, "%s", R->E.drafter);
	fscanf(F, "%lf", &(R->E.dT));
	fscanf(F, "%lf", &(R->E.maxT));
	fscanf(F, "%lf %lf %lf %lf", &(R->E.mLimit), &(R->E.pLimit), &(R->E.vLimit), &(R->E.rLimit));

	{
		int i;
		double mass, radius;
		Vec3 p, v;
		for(i=0;i<3;i++) {
			fscanf(F, "%lf %lf", &mass, &radius);
			fscanf(F, "%lf %lf %lf", &p.x, &p.y, &p.z);
			fscanf(F, "%lf %lf %lf", &v.x, &v.y, &v.z);
			planets[i] = new Planet(mass, radius, p, v);
		}
	}
	
	// I think it'd be better to make a file that can save E datas so that we can manage it easily.
	// Oh, I want to write in Korean. But it's github.한한국한국었한국어쓱한국어쓰곳한국어쓰고싶한국어쓰고싶다
	// 한국어 써 누가 쓰지 말래?
}


void writeResult(Result *R, FILE *f) {
	
	//static unsigned long long rnum = 0;
	//char filename[50];
	
	//sprintf(filename, "%x", rnum);
	
	//f = fopen(filename, "w");
	
	// output
	fprintf(f, "%d-%d-%d-%d", R->E.exGroup[0], R->E.exGroup[1], R->E.exGroup[2], R->E.exGroup[3]);
	fprintf(f, "%lf %lf %lf %lf %s", R->E.dT, R->E.T, R->E.GT, R->E.unitT, R->E.maxT, R->E.drafter);
	for(i=0;i<PLANET_N;i++) {
		fprintf(f, "%lf %lf", R->T.planets[i].mass, R->T.planets[i].r);
		fprintf(f, "%lf %lf %lf", R->T.planets[i].p.x,  R->T.planets[i].p.y,  R->T.planets[i].p.z);
		fprintf(f, "%lf %lf %lf", R->T.planets[i].v.x,  R->T.planets[i].v.y,  R->T.planets[i].v.z);
	}
	
	fprintf(f, "%d-%d", R->C->cdPair[0], R->C->cdPair[1]);
	fprintf(f, "%lf", R->C->cdTime);
	fprintf(f, "%lf %lf %lf", R->C->cdPoint.x, R->C->cdPoint.y, R->C->cdPoint.z);
	
	rnum++;
}

int main(int argc, char **argv) {
	if(argc < 3) {
		fprintf(stderr, "Usage: %s IN_FILE_NAME OUT_FILE_NAME\n", argv[0]);
		return 0;
	}
	
	Result r;
	FILE *f;

	fprintf(stderr, "%s: READ\n", argv[1]);

	/* Load Input Set */
	f = fopen(argv[1], "r");
	initExp(&r, f);
	fclose(f);

	fprintf(stderr, "%s: RUN\n", argv[1]);
	/* Run */
	loop(&r);

	fprintf(stderr, "%s: OUTPUT(%s)\n", argv[1], argv[2]);
	/* Output Result */
	f = fopen(argv[2], "w");
	writeResult(&r, f);
	fclose(f);
	return 0;
}

