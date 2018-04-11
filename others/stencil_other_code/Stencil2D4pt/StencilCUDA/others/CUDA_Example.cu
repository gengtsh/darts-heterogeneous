//Example CUDA code, written and commented by Jose Monsalve
//Taken from CUDA C/C++ Basics
//Supercomputing 2011 Tutorial
//by NVIDIA

/**
  This code executes c=a+b in a single thread in a GPU device.
It is a really simple code that is intended to show the memory 
movement between host and device, but not the division of the 
work between the differetn GPU threads.

First, we define the kernel (function that is executed in the 
GPU device), second, we initialize the values in the host (CPU)
and we allocate (reserve) some memory in the GPU). Then we move 
these values to the GPU device through explicitely memory movement
and we start the execution of the add kernel. Finally, after the 
computation is done, we copy the information back from the GPU and 
display it
**/

#include <stdio.h>
#include <stdlib.h>


//Simple add kernel, this function will be executed in the GPU device.
//C=A+B
__global__ void add(int *a, int *b, int *c)
{
	extern __shared__ int shared_mem[];
	int * shmem=shared_mem;
	shmem[threadIdx.x]=threadIdx.x;
	a[threadIdx.x]=shmem[threadIdx.x];
	b[threadIdx.x]=shmem[threadIdx.x];
	c[threadIdx.x]=a[threadIdx.x]+b[threadIdx.x];
}

//main func
int main(void)
{
	//Initializing the host variables
	//The *d_X are pointers that are not accessible in the host
	//directly but they represent a way to refer to the data in the 
	//GPU device
	int a[10], b[10], c[10];
	int *d_a, *d_b, *d_c;
	int size = sizeof(int)*10;

	//Initial values in the host 

	//Initializing memory in the GPU device, 
	//reserving the space, but there is not value yet (requires explicit movement)
	cudaMalloc ( (void **) &d_a, size );
	cudaMalloc ( (void **) &d_b, size );
	cudaMalloc ( (void **) &d_c, size );
	

	//moving the values that we just initialize in the host
	//to the GPU device (Explicit memory movement)
	cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);
	
	printf("Adding %d and %d in the device ... \n",a,b);
	
	//starting the kernel function
	add<<<1,10,10*sizeof(int)>>>(d_a,d_b,d_c);
	
	//Bringing back the result, which is stored in the GPU device 
	//And needs to be manually obtained.
	cudaMemcpy(c,d_c,size, cudaMemcpyDeviceToHost);

	printf("Result is \n");
	//printing the result
	for ( int i = 0; i < 10 ; i++)
		printf("%d\t",c[i]);
	printf("\n");
	
	//cleaning device.
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	return 0;
}
