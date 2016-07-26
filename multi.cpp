/* IGRam-1Gst */
/* Multi Process Working */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>

#include "phi_math.h"
#include "phi_type.h"

using namespace Phi;


/* Simultion Program Path */
char ex_path[256];

/* Current Experiment Number */
long long exp_num;

/* Input File Prefix */
char input_prefix[256];

/* Output File Prefix */
char output_prefix[256];

int max_p;
struct Process {
	pid_t pid;
	int exp_num;
} p[32];


void newProcess(int index) {
	pid_t pid;
	
	pid = fork();
	if(pid == -1) return;
	else if(pid == 0) { /* Child */
		char buf_in[512];
		char buf_out[512];
		sprintf(buf_in, "%s%lld", input_prefix, exp_num);
		sprintf(buf_out, "%s%lld.out", output_prefix, exp_num);
		execl(ex_path, ex_path, buf_in, buf_out, NULL);
	} else { /* Parent */
		p[index].pid = pid;
		p[index].exp_num = exp_num++;
	}
}

void updateProcessState(int index) {
	pid_t pid = p[index].pid;

	/* Check Process */
	if(kill(pid, 0) == 0) {
		/* Process is Alive */
		return;
	} else {
		/* Process is Dead */
		p[index].pid = -1;
		p[index].exp_num = 0;
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
	}
}

/* Initializer */
void init() {
	
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
}

