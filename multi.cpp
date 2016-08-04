/* IGRam-1Gst */
/* Multi Process Working */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>
#include <math.h>

#include "phi_math.h"
#include "phi_type.h"

using namespace Phi;


/* Limit and Step */
/* mass */
#define M_MIN 1.0
#define M_MAX 3.0
#define M_STEP 1.0
/* radius */
#define R_MIN 1.0
#define R_MAX 3.0
#define R_STEP 1.0
/* Velocity */
#define V_MIN 1.0
#define V_MAX 3.0
#define V_STEP 1.0
#define A1_MIN 0.0
#define A1_MAX 5.0
#define A1_STEP 1.0
#define A2_MIN -1.5
#define A2_MAX 1.5
#define A2_STEP 1.0

#define CAL_N(X) ((X##_MAX - X##_MIN) / X##_STEP + 1)
#define CAL_V(X, N) ((X##_MAX - X##_MIN) / X##_STEP * N + X##_MIN)

const long long int M_N = CAL_N(M);
const long long int R_N = CAL_N(R);
const long long int V_N = CAL_N(V);
const long long int A1_N = CAL_N(A1);
const long long int A2_N = CAL_N(A2);
const long long int CN = M_N * R_N * V_N * A1_N * A2_N;


char drafter[256];

/* Simultion Program Path */
char ex_path[256];

/* Current Experiment Number */
long long exp_num;

/* Input File Prefix */
char input_prefix[256];

/* Output File Prefix */
char output_prefix[256];

double timestamp, timeout;

int max_p;
struct Process {
	pid_t pid;
	long long exp_num;
} p[128];


double randDouble(double min, double max) {
	double d = max - min;
	int v = rand() % 32768;
	double dv = v / 32768.0;
	dv = min + dv * d;
	return dv;
}

struct Planet {
	double r, m;
	Vec3 p, v;
	
	Planet() {
		
	}

	Planet(int index, long long int num) {
		int mn = num % M_N; num /= M_N;
		int rn = num % R_N; num /= R_N;
		int v1n = num % V_N; num /= V_N;
		int v2n = num % A1_N; num /= A1_N;
		int v3n = num % A2_N; num /= A2_N;

		m = CAL_V(M, mn);
		r = CAL_V(R, rn);
		double d = CAL_V(V, v1n);
		double a1 = 2 * 3.141592 / 6.0 * CAL_V(A1, v2n);
		double a2 = 2 * 3.141592 / 6.0 * CAL_V(A2, v3n);
		
		v.x = cos(a1) * cos(a2);
		v.y = sin(a1) * cos(a2);
		v.z = sin(a2);
		v *= d;
		
		switch(index) {
		case 0: p.x = 0.0; p.y = 10.0; p.z = 0.0; break;
		case 1: p.x = -8.6; p.y = -5.0; p.z = 0.0; break;
		case 2: p.x = 8.6; p.y = -5.0; p.z = 0.0; break;
		}
	}

	void write(FILE *F) {
		/* radius */
		fprintf(F, "%lf\n", r);
		/* mass */
		fprintf(F, "%lf\n", m);
		/* Position */
		fprintf(F, "%lf %lf %lf\n", p.x, p.y, p.z);
		/* Velocity */
		fprintf(F, "%lf %lf %lf\n", v.x, v.y, v.z);
	}
};

void randomInput(FILE *F) {
	int i, j;
	/* drafter */
	
	fprintf(F, "%s\n", drafter);
	/* timestamp */
	fprintf(F, "%lf\n", timestamp);
	/* timeout */
	fprintf(F, "%lf\n", timeout);

	/* Limit */
	for(i=0;i<5;i++) {
		fprintf(F, "%lf ", 100000.0);
	}
	fprintf(F, "\n");

	

	Planet p[3];
	p[0] = Planet(0, exp_num % CN);
	p[1] = Planet(1, (exp_num / CN) % CN);
	p[2] = Planet(2, (exp_num / (CN * CN)) % CN);
	
	/* planets */
	for(i=0;i<3;i++) {
		p[i].write(F);
	}
}

void newProcess(int index) {
	pid_t pid;
	char buf_in[512];
	sprintf(buf_in, "%s%016llx", input_prefix, exp_num);

	/* Set Input File */
	FILE *input = fopen(buf_in, "w");
	randomInput(input);
	fclose(input);
	
	pid = fork();
	if(pid == -1) return;
	else if(pid == 0) { /* Child */
		char buf_out[512];
		sprintf(buf_out, "%s%016llx.out", output_prefix, exp_num);
		int result = execl(ex_path, ex_path, buf_in, buf_out, NULL);
		printf("%d\n", result);
	} else { /* Parent */
		p[index].pid = pid;
		p[index].exp_num = exp_num++;
	}
}

void deadProcess(int index) {
	char buf_in[512];
	sprintf(buf_in, "%s%016llx", input_prefix, p[index].exp_num);
	remove(buf_in);
	
	p[index].pid = -1;
	p[index].exp_num = 0;
}

void updateProcessState(int index) {
	pid_t pid = p[index].pid;

	/* Check Process */
	int r = kill(pid, 0);
	if(r == 0) {
		/* Process is Alive */
		return;
	} else {
		/* Process is Dead */
		deadProcess(index);
	}
}

/* Main Loop */
void loop() {
	int i;
	bool running = true;
	
	while(running) {
		for(i=0;i<max_p;i++) {
			if(p[i].pid > 0) { /* If Process is known as alive */
				/* Check Process State and Update State */
				updateProcessState(i);
			} else { /* No Current Process */
				/* New Process */
				newProcess(i);
			}
		}

		/* Sleep for 10 ms */
		usleep(10000);
		int status = 0;
		wait(&status);
	}
}

/* Initializer */
void init() {
	FILE *f = fopen("multi.init", "r");
	/* max process */
	fscanf(f, "%d ", &max_p);
	/* excutable path */
	fscanf(f, "%s ", ex_path);
	/* input prefix */
	fscanf(f, "%s ", input_prefix);
	/* output prefix */
	fscanf(f, "%s ", output_prefix);
	
	/* drafter */
	fscanf(f, "%s ", drafter);
	/* timestamp */
	fscanf(f, "%lf ", &timestamp);
	/* timeout */
	fscanf(f, "%lf ", &timeout);
	
	fclose(f);
}

void kill_all_process() {
	int i;
	for(i=0;i<max_p;i++) {
		if(p[i].pid > 0) {
			kill(p[i].pid, SIGKILL);
		}
	}
}

void sigint_handler() {
	kill_all_process();
	fprintf(stderr, "Interrupt");
	exit(0);
}

int main(int argc, char **argv) {
	signal(SIGINT, (__sighandler_t)sigint_handler);
	init();
	loop();
	return 0;
}

