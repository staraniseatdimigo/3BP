#include <math.h>
#include "phi_math.h"

using namespace Phi;

namespace Phi{

double fastInvSqrt(double var){
	double x2=var*0.5;
	const double tf = 1.5;
	union{long long i;double y;} t;
	t.y=var;
	t.i=0x5fe6ec85e7de30da-(t.i>>1);
	t.y*=tf-(x2*t.y*t.y);
	t.y*=tf-(x2*t.y*t.y);
	t.y*=tf-(x2*t.y*t.y);
	return t.y;
}

double fastSqrt(double var){
	double x2=var*0.5;
	const double tf = 1.5;
	union{long long i;double y;} t;
	t.y=var;
	t.i=0x5fe6ec85e7de30da-(t.i>>1);
	t.y*=tf-(x2*t.y*t.y);
	t.y*=tf-(x2*t.y*t.y);
	t.y*=tf-(x2*t.y*t.y);
	return t.y*var;
}

double abs(double v){
	if(v<0) return -v;
	return v;
}

}
