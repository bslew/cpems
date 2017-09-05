/*!
* \file taskControlServer.cpp
*
*  Created on: May 13, 2010
*      Author: blew
*/

/* A simple server in the internet domain using TCP
   The port number is passed as an argument */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include "Mscs-global-defs.h"

using namespace std;


/*!
	\brief reads the information from shared memory about the number of runs and names of the files that hold the run information
	\details
	@return number of runs also stored in the global variable MSCS_GLOBAL__RUN_SIMULATION_INFO
	Returns -1 if there was an error in reading the shared memory or if there are no running processes

	\date May 14, 2010, 1:56:46 PM
	\author Bartosz Lew
*/
int readRunsInfo();

/*!
	\brief reads the actual information about i'th run from run file
	\details
	@param i - run id number that should be from within range MSCS_GLOBAL__RUN_SIMULATION_INFO.st_thread to MSCS_GLOBAL__RUN_SIMULATION_INFO.en_thread
	@param simNum - holds the total number of number of simulations to be processed by i'th process
	@param stSim - holds the starting simulation id number
	@param enSim - holds the ending simulation id number
	@param currSim - holds the number of the currently processed simulation
	@return 0 if the information could be read from file and -1 if there was an error opening the run file

	This routine actually looks into files whose names are defined in the MSCS_GLOBAL__RUN_SIMULATION_INFO variable.
	The file should contain only four integer numbers separated by white space.
	numer_of_sumulations_on_current_node starting_sim ending_sim current_sim_id

	This information is next sent over tcp socket to the client as a text.

	\date May 14, 2010, 1:59:18 PM
	\author Bartosz Lew
*/
int getRunInfo(int i, int* simNum, int* stSim, int* enSim, int* currSim);

int check_command(char* buff) {
	if (strcmp(buff,"run status\r\n")==0) return 1;

	return 0;
}


void error(char *msg)
{
	perror(msg);
	exit(1);
}

int main(int argc, char *argv[]) {
	Mscs_initiate_global_variables();
	char s[200];
	string ss;
	int sockfd, newsockfd, portno, clilen;
	char buffer[256];
	struct sockaddr_in serv_addr, cli_addr;
	int n;

	int runsNo, currSim, simNum, stSim, enSim;

	if (argc < 2) {
		fprintf(stderr,"ERROR, no port provided\n");
		exit(1);
	}

	printf("* making a socket\n");
	sockfd = socket(AF_INET, SOCK_STREAM, 0); // get socket file descriptor
	if (sockfd < 0) error("ERROR opening socket");

	printf("* setting up the socket parameters\n");
	bzero((char *) &serv_addr, sizeof(serv_addr));
	portno = atoi(argv[1]);
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = INADDR_ANY;
	serv_addr.sin_port = htons(portno);

	// assign name (server address) to the scoket
	printf("* binding to a socket\n");
	if (bind(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0)  error("ERROR on binding");

	while (1) {
		runsNo=readRunsInfo();
		printf("runs number: %i\n",runsNo);
		printf("* waiting for connection\n");
		listen(sockfd,5); // enter a passiv state - waiting for connection, and set the maximal number of possible pending connections to 5
		// the program will stop here until the connection is made.

		clilen = sizeof(cli_addr);

		//
		// create a new connected socket from the list of pending sockets and assign the client addres and port to it
		//
		printf("* accepting incomming connection\n");
		newsockfd = accept(sockfd,
				(struct sockaddr *) &cli_addr,
				(socklen_t*)&clilen);
		// a new socket file descriptor is returned and it will be used for further communication full-duplex communication
		strcpy(s,inet_ntoa(cli_addr.sin_addr));
		printf(" -- an incomming connection from: %s\n", s);

		if (newsockfd < 0)  error("ERROR on accept");

		bzero(buffer,256);
		printf("* reading from the socket\n");
		n = read(newsockfd,buffer,255); // start reading from the socket
		if (n < 0) error("ERROR reading from socket");
		printf("Incoming command is: %s\n",buffer);

		if (check_command(buffer)) {
			printf("* sending answer to the socket\n");
			for (int i=0;i<runsNo;i++) {
				getRunInfo(i, &simNum, &stSim, &enSim, &currSim);
				bzero(s,sizeof(s));
				sprintf(s,"%i %i %i %i %i\n",i, simNum, stSim, enSim, currSim);
				ss=s;
				printf("sending: %i %i %i %i %i, total size: %i\n",i, simNum, stSim, enSim, currSim, ss.size());
				n = write(newsockfd,s,ss.size());
			}
		}
		else {
			printf("* unknown command\n");
		}
		if (n < 0) error("ERROR writing to socket");
		close(newsockfd);
	}
	return 0;
}



int readRunsInfo() {
	int fd;
	fd = shm_open(MSCS_GLOBAL__RUN_SIMULATION_SHM_KEY, O_RDONLY, 0444);
	if (fd == -1) {	 printf("error: shared memory not opened\n");   return -1; }
	//if (ftruncate(fd, sizeof(Mscs_run_control_t)) == -1) {	printf("error: not resizing\n");    return -1; }

	MSCS_GLOBAL__RUN_SIMULATION_INFO= (Mscs_run_control_t*)mmap(NULL, sizeof(Mscs_run_control_t), PROT_READ , MAP_SHARED, fd, 0);
	if (MSCS_GLOBAL__RUN_SIMULATION_INFO == MAP_FAILED) { printf("error: not mapping\n"); return -1; }
	return MSCS_GLOBAL__RUN_SIMULATION_INFO->procNum;
}

int getRunInfo(int i, int* simNum, int* stSim, int* enSim, int* currSim) {
	string tmps;
	if (MSCS_GLOBAL__RUN_SIMULATION_INFO->st_thread<0 && MSCS_GLOBAL__RUN_SIMULATION_INFO->en_thread<0) {
		tmps=MSCS_GLOBAL__RUN_SIMULATION_INFO->sim_run_file_pref;
		tmps+=MSCS_GLOBAL__RUN_SIMULATION_INFO->sim_run_file_suff;

	}
	else {
		char tmpch[4];
		sprintf(tmpch,"%i",i);
		tmps=MSCS_GLOBAL__RUN_SIMULATION_INFO->sim_run_file_pref;
		tmps+=tmpch;
		tmps+=MSCS_GLOBAL__RUN_SIMULATION_INFO->sim_run_file_suff;
	}

	FILE* f=fopen(tmps.c_str(),"r");
	if (f==NULL) { tmps="could not open the run file: "+tmps; perror(tmps.c_str()); return -1; }
	else printf("reading from file: %s\n",tmps.c_str());

	fscanf(f,"%i %i %i %i",simNum,stSim, enSim, currSim);
	return 0;
}
