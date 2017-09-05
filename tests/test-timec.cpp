/*!
 * \file test-timec.cpp
 *
 *  Project: cpeds
 *  Created on: Mar 8, 2013 11:25:47 AM
 *  Author: blew
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
//#include <sys/shm.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <string.h>

#define timestampSHMKEY "1979"
#define MBUFF_DEV_NAME "/dev/mbuff"
#define IOCTL_MBUFF_ALLOCATE 1
#define MBUFF_NAME_LEN 32

struct mbuff_request_struct{
	unsigned int flags;

	char name[MBUFF_NAME_LEN+1];

	size_t size;

	unsigned int reserved[4];
};

typedef struct {
		long sec;
		long nsec;
} systimestamp_t;

//typedef struct timeval systimestamp_t;
int SHmem;
void *SHmbuf;
pthread_mutex_t systimestampMutex = PTHREAD_MUTEX_INITIALIZER;
systimestamp_t *systimestamp;

void delay(int sec, int nsec){
	struct timespec time_delay;
	time_delay.tv_sec = sec;
	time_delay.tv_nsec = nsec;
	nanosleep(&time_delay,NULL);
}


inline void * mbuff_alloc_at(const char *name, int size, void * addr) {
	int fd;
	struct mbuff_request_struct req={0,"default",0,{0}};
	void * mbuf;

	if(name) strncpy(req.name,name,sizeof(req.name));
	req.name[sizeof(req.name)-1]='\0';
	req.size = size;
	if(( fd = open(MBUFF_DEV_NAME,O_RDWR) ) < 0 ){
		perror("open failed");
		return NULL;
	}
	size=ioctl(fd,IOCTL_MBUFF_ALLOCATE,&req);
	if(size<req.size)
		return NULL;
/* the type of first mmap's argument depends on libc version? This really
 * drives me crazy. Man mmap says "void * start" */
	mbuf=mmap(addr, size,PROT_WRITE|PROT_READ,MAP_SHARED|MAP_FILE,fd, 0);
	if( mbuf == (void *) -1) 
		mbuf=NULL;
	close(fd);
	return mbuf;
}

inline void * mbuff_alloc(const char *name, int size) {
	return mbuff_alloc_at(name, size, NULL);
}



int main() {
	
	systimestamp = (systimestamp_t*)mbuff_alloc(timestampSHMKEY,sizeof(systimestamp_t));
	if (systimestamp == NULL) return -1;
	
	//	SHmem = shm_open(timestampSHMKEY, O_WRONLY, 0660);
	//	SHmem=shmget(timestampSHMKEY,sizeof(systimestamp_t),0660 | IPC_CREAT);
	//	SHmbuf=shmat(SHmem,(char *)0,0);
	//	if ((long int)SHmbuf==-1)
	//	{
	//		perror("Pointer to the shared memory segment cannot be set.\n");
	//		exit(4);
	//	}
	//	else
	//		printf("Pointer to the shared memory set.\n");
	//	
	//	systimestamp=(systimestamp_t*)SHmbuf;

	
	
	
//	SHmem  = shm_open(timestampSHMKEY, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
//
//	if (SHmem==-1) {
//		perror("Shared memory segment cannot be opened.");
//		exit(-1);
//	}
//	else {
//		printf("Shared memory segment opened.\n");
//	}
//    if (ftruncate(SHmem, sizeof(systimestamp_t)) == -1) { 
//    	printf("shared memory: adjusting shared memory size failed\n"); 
//    	exit(-2);
//	}
//    systimestamp= (systimestamp_t*)mmap(NULL, sizeof(systimestamp_t), PROT_READ | PROT_WRITE, MAP_SHARED, SHmem, 0);
//    if (systimestamp == MAP_FAILED) { printf("shared memory: mapping failed\n"); exit(-3); }


	
	
	struct timeval  tv;
	struct timezone tz;
	while (1) {
		gettimeofday(&tv, &tz);
		printf("stamp: %.15lf\n",double(tv.tv_sec)+double(tv.tv_usec)*1e-6);
		pthread_mutex_lock(&systimestampMutex);
		systimestamp->sec=tv.tv_sec;
		systimestamp->nsec=tv.tv_usec*1000L;
		pthread_mutex_unlock(&systimestampMutex);
		delay(0,int(50e6));
	}
}
