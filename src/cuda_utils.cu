// Claro jijodebú, necesito un file .cu...

#include <iostream> // Required for standard input/output operations (e.g., printf, cout)
#include <cuda_runtime.h> // Required for CUDA runtime API functions (e.g., cudaMalloc, cudaMemcpy, cudaFree)
#include <cmath>

#include "device_launch_parameters.h"  // ??? -> ask for it.

// Def a checker:
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)
    
// Antes de todo, las __ctes__ & el kernel (__global__) de CUDA:
__constant__ int cant_particles;

//__constant__ float pos_centre[3];  // Camb por x, y , z...
__constant__ float x_centre;
__constant__ float y_centre;
__constant__ float z_centre;

__constant__ float scale;
__constant__ float softening;
__constant__ float grav_cte;
__constant__ float mass_centre;
__constant__ float delta_step;
// Estas variables las voy a tener que llamar del host (en cada step pls), así que son args del wrapper!!! (=/= kernel (!))

// Then, the function that is passed to the GPU a.k.a. "kernel" (will be executed on the 
// GPU by multiple threads in parallel; Each thread will process a single element of the array).
// Toca SOLO los arrays de pos y acc (!)
__global__ void integrate_CUDA(float* position, float* velocity, float* acceleration)
    {
    // Calculate the global index of the current element
    // threadIdx.x: Local thread ID within a block
    // blockIdx.x: Block ID within the grid
    // blockDim.x: Number of threads per block
    int tid = blockIdx.x * blockDim.x + threadIdx.x;  // "tid" == "thread ID"

    // "pos" such that pos[3*N] = x[N], pos[3*N + 1] = y[N], pos[3*N + 2] = z[N]  (!)
    float rMinusRjScaled[3];
    float dot, distance_ij3;

    // Check if the current index is within the bounds of the array, then do things:
    if (tid < cant_particles) {
        // copy & paste grav...
        rMinusRjScaled[0] = (position[3*tid + 0] - x_centre) * scale;
        rMinusRjScaled[1] = (position[3*tid + 1] - y_centre) * scale;
        rMinusRjScaled[2] = (position[3*tid + 2] - z_centre) * scale;

        dot = rMinusRjScaled[0] * rMinusRjScaled[0] + rMinusRjScaled[1] * rMinusRjScaled[1] +\
                rMinusRjScaled[2] * rMinusRjScaled[2];
        dot = sqrtf(dot);

        distance_ij3 = (dot + softening) * (dot + softening) * (dot + softening);

        // Updateo la gravedad:
        // acceleration += gravityTerm;
        acceleration[3*tid + 0] += -grav_cte * mass_centre * (rMinusRjScaled[0]/distance_ij3);
        acceleration[3*tid + 1] += -grav_cte * mass_centre * (rMinusRjScaled[1]/distance_ij3);
        acceleration[3*tid + 2] += -grav_cte * mass_centre * (rMinusRjScaled[2]/distance_ij3);

        // Acá le podría agregar que si la densidad de la part es muy alta,
        // que se vuele hacia afuera...
        
        // OJO! Me falto check for CFL condition:
		dot = (acceleration[3*tid +0] * acceleration[3*tid +0]) + (acceleration[3*tid +1] * acceleration[3*tid +1]) +\
			 (acceleration[3*tid +2] * acceleration[3*tid +2]);

		// Lo meto a dedo... (recall branchles...)
		bool limitExceeded = (dot > 1e+8f); //"mCflLimit2"
		if (limitExceeded)
		{
		  float cflScale = 1e+4f * rsqrtf(dot);
		  acceleration[3*tid +0] *= cflScale;
		  acceleration[3*tid +1] *= cflScale;
		  acceleration[3*tid +2] *= cflScale;
		}
        
        
        // Integro todo! LF-KDK: Only gravity (reset accel para el mid-step!)
		velocity[3*tid + 0] += acceleration[3*tid + 0] * delta_step * 0.5f;
		velocity[3*tid + 1] += acceleration[3*tid + 1] * delta_step * 0.5f;
		velocity[3*tid + 2] += acceleration[3*tid + 2] * delta_step * 0.5f;

		position[3*tid + 0] += velocity[3*tid + 0] * delta_step;
		position[3*tid + 1] += velocity[3*tid + 1] * delta_step;
		position[3*tid + 2] += velocity[3*tid + 2] * delta_step;

		// Nuevamente calc la accel...		  
		rMinusRjScaled[0] = (position[3*tid + 0] - x_centre) * scale;
        rMinusRjScaled[1] = (position[3*tid + 1] - y_centre) * scale;
        rMinusRjScaled[2] = (position[3*tid + 2] - z_centre) * scale;

		dot = rMinusRjScaled[0] * rMinusRjScaled[0] + rMinusRjScaled[1] * rMinusRjScaled[1] +\
					 rMinusRjScaled[2] * rMinusRjScaled[2];
		dot = sqrtf(dot);

		distance_ij3 = (dot + softening) * (dot + softening) * (dot + softening);

		// Updateo la gravedad:
		acceleration[3*tid + 0] = -grav_cte * mass_centre * (rMinusRjScaled[0]/distance_ij3);
		acceleration[3*tid + 1] = -grav_cte * mass_centre * (rMinusRjScaled[1]/distance_ij3);
		acceleration[3*tid + 2] = -grav_cte * mass_centre * (rMinusRjScaled[2]/distance_ij3);

		// Integro todo! LF-KDK: Only gravity
		velocity[3*tid + 0] += acceleration[3*tid + 0] * delta_step;
		velocity[3*tid + 1] += acceleration[3*tid + 1] * delta_step;
		velocity[3*tid + 2] += acceleration[3*tid + 2] * delta_step;

		// Y ya modifiqué todos los vectores!

    }
}

// Ese es el kernel principal, ahora falta el wrapper (que pide las ctes al host!):

// Host-side wrapper function
void launchMyKernel(float* h_position, float* h_velocity, float* h_acceleration,
				int h_cant_particles, float h_scale, float h_softening, float h_grav_cte,
				float h_mass_centre, float h_delta_step,
				float h_x_centre, float h_y_centre, float h_z_centre)
{
	const int N = h_cant_particles;
	// Pointer(s) to the array on the device (GPU)
	float* device_pos;
	float* device_vel;
	float* device_acc;

	// Allocate memory on the device
	// Crucial: Use CUDA_CHECK for all CUDA API calls for proper error handling.
	CUDA_CHECK(cudaMalloc(&device_pos, 3 * N * sizeof(float)));
	CUDA_CHECK(cudaMalloc(&device_vel, 3 * N * sizeof(float)));
	CUDA_CHECK(cudaMalloc(&device_acc, 3 * N * sizeof(float)));
	
	// Este COPY ocurre step by step...
	// Copy data from host to device
	// cudaMemcpy copies data between host and device.
	// Arguments: destination, source, size, direction (host to device)
	cudaMemcpy(device_pos, h_position, 3*N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(device_vel, h_velocity, 3*N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(device_acc, h_acceleration, 3*N * sizeof(float), cudaMemcpyHostToDevice);

	// Y esto ahora???
	
	// Copy the constant from host to device __constant__ memory
	// cudaMemcpyToSymbol is used to copy data to a __constant__ variable on the device.
	// Arguments: symbol name (address of the __constant__ variable), source, size
	cudaMemcpyToSymbol(cant_particles, &h_cant_particles, sizeof(int), 0, cudaMemcpyHostToDevice);
	
	//cudaMemcpyToSymbol(pos_centre, &h_pos_centre, sizeof(float) * 3, 0, cudaMemcpyHostToDevice);  // Cambio por x, y, z...
	cudaMemcpyToSymbol(x_centre, &h_x_centre, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(y_centre, &h_y_centre, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(z_centre, &h_z_centre, sizeof(float), 0, cudaMemcpyHostToDevice);
	
	cudaMemcpyToSymbol(scale, &h_scale, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(softening, &h_softening, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(grav_cte, &h_grav_cte, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(mass_centre, &h_mass_centre, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(delta_step, &h_delta_step, sizeof(float), 0, cudaMemcpyHostToDevice);	
	
	// Define grid and block dimensions for kernel execution
	// A block is a group of threads that can cooperate.
	// A grid is a collection of blocks.
	// Here, we use 256 threads per block.
	// The number of blocks is calculated to cover all elements.
	int threadsPerBlock = 256;
	int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock; // Ceiling division
	
	// Launch the kernel -> Ojo ctes...
	integrate_CUDA<<<blocksPerGrid, threadsPerBlock>>>(device_pos, device_vel, device_acc);
    
    // Synchronize the device to ensure all kernel operations are complete
    // cudaDeviceSynchronize blocks the CPU until all GPU tasks are finished.
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy data back from device to host
    // Arguments: destination, source, size, direction (device to host)
    CUDA_CHECK(cudaMemcpy(h_position, device_pos, 3 * N * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_velocity, device_vel, 3 * N * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_acceleration, device_acc, 3 * N * sizeof(float), cudaMemcpyDeviceToHost));

    // "h_position", etc. now contain the updated values from the GPU.

    // Free device memory
    // It's crucial to free allocated GPU memory to prevent memory leaks.
    CUDA_CHECK(cudaFree(device_pos));
    CUDA_CHECK(cudaFree(device_vel));
    CUDA_CHECK(cudaFree(device_acc));



}









