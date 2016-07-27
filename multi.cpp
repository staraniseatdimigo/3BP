/* IGRam-1Gst */
/* Multi Process Working */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>

#include "phi_math.h"
#include "phi_type.h"

using namespace Phi;


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

void randomInput(FILE *F) {
	int i, j;
	/* drafter */
	fprintf(F, "%s\n", drafter);
	/* timestamp */
	fprintf(F, "%lf\n", timestamp);
	/* timeout */
	fprintf(F, "%lf\n", timeout);

	/* Limit */
	for(i=0;i<4;i++) {
		fprintf(F, "%lf ", 100000.0);
	}
	fprintf(F, "\n");
	
	/* planets */
	for(i=0;i<3;i++) {
		/* radius */
		fprintf(F, "%lf\n", randDouble(0.1, 20.0));
		/* mass */
		fprintf(F, "%lf\n", randDouble(1.0, 1000.0));
		/* Position */
		for(j=0;j<3;j++)
			fprintf(F, "%lf ", randDouble(-1000.0, 1000.0));
		fprintf(F, "\n");
		/* Velocity */
		for(j=0;j<3;j++)
			fprintf(F, "%lf ", randDouble(-100.0, 100.0));
		fprintf(F, "\n");
	}
}

void newProcess(int index) {
	pid_t pid;
	char buf_in[512];
	sprintf(buf_in, "%s%llx", input_prefix, exp_num);

	/* Set Input File */
	FILE *input = fopen(buf_in, "w");
	randomInput(input);
	fclose(input);
	
	pid = fork();
	if(pid == -1) return;
	else if(pid == 0) { /* Child */
		char buf_out[512];
		sprintf(buf_out, "%s%llx.out", output_prefix, exp_num);
		int result = execl(ex_path, ex_path, buf_in, buf_out, NULL);
		printf("%d\n", result);
	} else { /* Parent */
		p[index].pid = pid;
		p[index].exp_num = exp_num++;
	}
}

void deadProcess(int index) {
	char buf_in[512];
	sprintf(buf_in, "%s%lld", input_prefix, p[index].exp_num);
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

