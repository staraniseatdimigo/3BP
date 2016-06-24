/* IGRam-1Gst */

#include "phi_math.h"
#include "phi_type.h"

#include <stdlib.h>
#include <stdio.h>

using namespace Phi;

#define G 1

class Planet {
public:
	double mass;
	Vec3 p, v;
	
	Vec3 temp_force;

	Planet(double _mass, Vec3 _p, Vec3 _v) {
		mass = _mass;
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

void loop() {
	int running = 1;
	int i, j;
	while(running) {
		for(i=0;i<PLANET_N;i++) {
			planets[i]->initializeTempForce();
		}
		for(i=1;i<PLANET_N;i++) {
			for(j=0;j<i;j++) {
				planets[i]->calculateTempForce(planets[j], timestamp);
				planets[j]->calculateTempForce(planets[i], timestamp);
			}
		}
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

