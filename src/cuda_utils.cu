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

// NEW: h & krnl-wise (todas con MemSymbolCpy()...)
__constant__ float h_krnl;
__constant__ float h_krnl2;
__constant__ float mHScaled9;
__constant__ float mKernel1Scaled;
__constant__ float mKernel2Scaled;
__constant__ float mKernel3Scaled;

// Also these:
__constant__ float mRho0;
__constant__ float mViscosityScalar;
__constant__ float mStiffness;

// Estas variables las voy a tener que llamar del host (en cada step pls), así que son args del wrapper!!! (=/= kernel (!))

/*

NEW: Como el findNeighbors() es muy hincha bolas, y quiero mergear todas las funcs para dárselas de comer a la GPU,
	¿Qué pasaría si transformo el problema a un N-body? Es decir, Todas las partículas le preguntan la distancia al resto,
	y a partir de los vecinos encontrados, computamos ON-THE-FLY las demás props hidrodinámicas (y de paso nos ahorramos
	guardar el indice de vecinos) ¿Mucho bardo? ¿Sería mejor hacer un "voxelize()" en la GPU? Probemos el approach O(n^2)...

*/

// Then, the function that is passed to the GPU a.k.a. "kernel" (will be executed on the 
// GPU by multiple threads in parallel; Each thread will process a single element of the array).
// Toca SOLO los arrays de pos y acc (!)

// NEW: separo la accel por vecinos de la grav central? Quiero aprovechar al máximo
// los blocks como working-units, y usar ~intrinsics de warps!!!
__global__ void neighbors_CUDA(float* position, float* velocity, float* acceleration,
							float* mass, float* density)
{
	// Thread ID - global!
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid >= cant_particles) return;

	// Within a block, min(tid) = 0; max(tid) = dim(block)
	int dim_block = blockDim.x;

	// 1st, necesito allocar la memoria que solo va a tocar este bloque (contigua!):
	// Each thread helps load a local reference.
    __shared__ float s_localPositions[3 * dim_block];  // (x,y,z)
	__shared__ float s_localVelocities[3 * dim_block];  // (vx,vy,vz)
	__shared__ float s_localAccels[3 * dim_block];  // (ax,ay,az)
	__shared__ float s_localMasses[dim_block];
	__shared__ float s_localDensities[dim_block];
	// If I want dynamic shared memory -> "extern" and pass it when calling the kernel

	// Thread index within the block
    int localThreadId = threadIdx.x;

    // Load local reference particles into shared memory
    // Fetch particles from global memory into shared memory.
    // This assumes local references are contiguous in global memory,
    // or that each block is given a pointer to its own local reference set.
    // For this example, let's assume each block is responsible for
    // a segment of the global input particles, and these are its "references".
    // Assume (localThreadId < numLocalReferencesPerBlock)...
	s_localPositions[3*localThreadId + 0] = position[3*tid + 0];
	s_localPositions[3*localThreadId + 1] = position[3*tid + 1];
	s_localPositions[3*localThreadId + 2] = position[3*tid + 2];

	s_localVelocities[3*localThreadId + 0] = velocity[3*tid + 0];
	s_localVelocities[3*localThreadId + 1] = velocity[3*tid + 1];
	s_localVelocities[3*localThreadId + 2] = velocity[3*tid + 2];

	s_localMasses[localThreadId] = mass[tid];
	s_localDensities[localThreadId] = density[tid];

	// I have to create an array of "neighbor_counts" to allow the warp to see the leader's neighbcount!!!
	// ...within the block, as always.
	__shared__ int s_neighbor_count[dim_block];
	s_neighbor_count[localThreadId] = 0;  // Now, all threads within a warp can access this.

    __syncthreads(); // Ensure all shared memory loads are complete before any thread uses them

	// vars def thread-wise.
	s_localDensities[localThreadId] = 0.f;
	float distance_ij3;
	float rMinusRjScaled[3];	

	// Main loop - block-wise (Ojo iters!)
	// Acá deberían ir de a 32 threads...
	for (int j = 0; j < dim_block; j++)
	{
		// Este break es thread-wise (privado), y tiene que venir al principio y NO al final (así
		// no queda ningún worker trabajando al pedo)
		if (s_neighbor_count[localThreadId] > 32) break;

		rMinusRjScaled[0] = (s_localPositions[3*localThreadId + 0] - s_localPositions[3*j + 0]) * scale;
        rMinusRjScaled[1] = (s_localPositions[3*localThreadId + 1] - s_localPositions[3*j + 1]) * scale;
        rMinusRjScaled[2] = (s_localPositions[3*localThreadId + 2] - s_localPositions[3*j + 2]) * scale;

		distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] +\
						rMinusRjScaled[1] * rMinusRjScaled[1] +\
						rMinusRjScaled[2] * rMinusRjScaled[2] + softening;  // This is SQUARED
		distance = sqrtf(distance_ij3);  // This is the true d_ij
		invDist = rsqrtf(distance_ij3);  // quick x^(-1/2)

		// From this, I can (branchless-ly) compute hydro properties (!)
		// Density:
		rightPart = h_krnl2 - distance_ij3;
		rightPart = rightPart * rightPart * rightPart;
		w = mKernel1Scaled * rightPart * (distance_ij3 < h_krnl2);  // true => 1, false => 0 (!)
		// apply weighted neighbor mass to our density
		s_localDensities[localThreadId] += (s_localMasses[localThreadId] * w);  // 0. if (d^2 >= h^2)

		// Pressure gradient:
		mj = s_localMasses[localThreadId + j];
		pj = (s_localDensities[localThreadId + j] - mRho0) * mStiffness;
		rhojInv = ((s_localDensities[localThreadId + j] > 0.0f) ? (1.0f / s_localDensities[localThreadId + j]) : 1.0f);
		rhojInv2 = rhojInv * rhojInv;

		centerPart = (h_krnl - distance) * (distance_ij3 < h_krnl2);  // true => 1, false => 0 (!)
		centerPart *= centerPart;  // 0 if d >= h
		centerPart *= mj * piDivRhoi2 * (pj * rhojInv2);  // 0 if d >= h

		pressureGradientContribution[0] = mKernel2Scaled * rMinusRjScaled[0] * invDist;
		pressureGradientContribution[1] = mKernel2Scaled * rMinusRjScaled[1] * invDist;
		pressureGradientContribution[2] = mKernel2Scaled * rMinusRjScaled[2] * invDist;

		// add pressure gradient contribution to pressure gradient
		pressureGradient[0] += pressureGradientContribution[0] * centerPart;  // 0 if d >= h
		pressureGradient[1] += pressureGradientContribution[1] * centerPart;
		pressureGradient[2] += pressureGradientContribution[2] * centerPart;

		// Viscosity (Reuso variables):
		centerPart = (h_krnl - distance) * (distance_ij3 < h_krnl2);  // true => 1, false => 0 (!)
		centerPart *= rhojInv * mj * mKernel3Scaled;  // 0 if d >= h

		// add contribution to viscous term (+0 if d >= h)
		viscousTerm[0] += (s_localVelocities[3*localThreadId + 0] - s_localVelocities[3*j + 0]) *\
						 centerPart * mViscosityScalar * rhoiInv;
		viscousTerm[1] += (s_localVelocities[3*localThreadId + 1] - s_localVelocities[3*j + 1]) *\
						 centerPart * mViscosityScalar * rhoiInv;
		viscousTerm[2] += (s_localVelocities[3*localThreadId + 2] - s_localVelocities[3*j + 2]) *\
						 centerPart * mViscosityScalar * rhoiInv;

		// Acá puedo hacer un mask (ya que en este loop vienen actuando warps de 32-threads)
		// y preg cuáles de ellos suman a valid neighbors:

		// Mask of active threads in the warp
    	unsigned int active_threads_mask = __activemask();

		// Perform a warp-wide ballot to see which threads meet the condition
		// (true (non-zero) or false (zero) for each thread).
		// The result 'condition_mask' is a 32-bit integer where bit 'i' is set
		// if thread 'i' in the warp met the condition.

		// Def the condition:
		bool valid_neighbor = (distance_ij3 < h_krnl2);
		// Ask them:
		unsigned int condition_mask = __ballot_sync(active_threads_mask, valid_neighbor);
		// __popc() to count how many threads in this warp met the condition:
    	int num_threads_meeting_condition = __popc(condition_mask);

		// Only the warp leader count its private counter!
		// =/= a LocalID == 0 !!!
		if (localThreadId % 32 == 0)
		{
			// Atomically add to a global or shared counter for this shared target particle
			// This ensures the sum for the warp's contributions is added only once
			atomicAdd(&s_neighbor_count[localThreadId], num_threads_meeting_condition);
		}

		// It's crucial to synchronize threads *after* the atomicAdd if *all* threads
        // need to see the *updated* shared count before the next iteration of the loop.
        // Without this, some threads might read an old value of s_shared_neighbor_count[LocalID]
        // at the start of the next iteration, potentially missing the break condition.
        __syncthreads(); // Synchronize all threads in the block to see the updated count

	}

	// Finish parte hydro, aparte haremos la accel central desp...
	// Assert que todos los del bloque terminaron de escribir previous to write globally!
	__syncthreads();

	// Every thread of the block writes in its own place the results from shared memory
    // to global memory.
	acceleration[3*tid + 0] = viscousTerm[0] - pressureGradient[0];
	acceleration[3*tid + 1] = viscousTerm[1] - pressureGradient[1];
	acceleration[3*tid + 2] = viscousTerm[2] - pressureGradient[2];

}


__global__ void integrate_CUDA(float* position, float* velocity, float* acceleration)
{
	// 1st, la 1ra parte grav que me falto desp de la hydro; desp integrate

	// Thread ID
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	// Init the vars needed...
	float distance_ij3, invDist, dot;
	float rMinusRjScaled[3];

	// Checker
	if (tid >= cant_particles) return;

	// Veamos la distancia al BH central (no dependo de las demás):
	rMinusRjScaled[0] = (position[3*tid + 0] - x_centre) * scale;
	rMinusRjScaled[1] = (position[3*tid + 1] - y_centre) * scale;
	rMinusRjScaled[2] = (position[3*tid + 2] - z_centre) * scale;

	distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] + rMinusRjScaled[1] * rMinusRjScaled[1] +\
				rMinusRjScaled[2] * rMinusRjScaled[2] + softening;

	invDist = rsqrtf(distance_ij3);  // quick x^(-1/2)
	invDist = invDist * invDist * invDist;  // dist^-3

	// Updateo la gravedad:
	// acceleration += gravityTerm;
	acceleration[3*tid + 0] += -grav_cte * mass_centre * (rMinusRjScaled[0] * invDist);
	acceleration[3*tid + 1] += -grav_cte * mass_centre * (rMinusRjScaled[1] * invDist);
	acceleration[3*tid + 2] += -grav_cte * mass_centre * (rMinusRjScaled[2] * invDist);
	
	// OJO! Me falto check for CFL condition:
	dot = (acceleration[3*tid +0] * acceleration[3*tid +0]) + (acceleration[3*tid +1] * acceleration[3*tid +1]) +\
		(acceleration[3*tid +2] * acceleration[3*tid +2]);

	// Lo meto a dedo... (recall branchles...)
	acceleration[3*tid +0] *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	acceleration[3*tid +1] *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	acceleration[3*tid +2] *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	// acc = acc if (|acc|^2 < CFL^2 => "x1"); acc = acc/|acc| * CFL if (|acc|^2 >= CFL^2)
	
	// ---------------------------------------------
	// Listo! Ahora, integro con LF-KDK -> Next kernel...	
	// ---------------------------------------------
		
	// Only gravity (reset accel para el mid-step!)
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

	distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] + rMinusRjScaled[1] * rMinusRjScaled[1] +\
				rMinusRjScaled[2] * rMinusRjScaled[2] + softening;

	invDist = rsqrtf(distance_ij3);  // quick x^(-1/2)
	invDist = invDist * invDist * invDist;  // dist^-3

	// Updateo la gravedad:
	// acceleration += gravityTerm;
	acceleration[3*tid + 0] += -grav_cte * mass_centre * (rMinusRjScaled[0] * invDist);
	acceleration[3*tid + 1] += -grav_cte * mass_centre * (rMinusRjScaled[1] * invDist);
	acceleration[3*tid + 2] += -grav_cte * mass_centre * (rMinusRjScaled[2] * invDist);

	// Integro todo! LF-KDK: Only gravity
	velocity[3*tid + 0] += acceleration[3*tid + 0] * delta_step;
	velocity[3*tid + 1] += acceleration[3*tid + 1] * delta_step;
	velocity[3*tid + 2] += acceleration[3*tid + 2] * delta_step;

	// Y ya modifiqué todos los vectores!
}


// Host-side wrapper function
void launchMyKernel(float* h_position, float* h_velocity, float* h_acceleration,
				float* h_mass, float* h_density,
				int h_cant_particles, float h_scale, float h_softening, float h_grav_cte,
				float h_mass_centre, float h_delta_step,
				float h_x_centre, float h_y_centre, float h_z_centre,
				float h_h_krnl, float h_h_krnl2, float h_mHScaled9, float h_mKernel1Scaled,
				float h_mKernel2Scaled, float h_mKernel3Scaled, float h_mRho0,
				float h_mViscosityScalar, float h_mStiffness)
{
	const int N = h_cant_particles;
	// Pointer(s) to the array on the device (GPU)
	float* device_pos;
	float* device_vel;
	float* device_acc;
	float* device_mass;
	float* device_dens;

	// Allocate memory on the device
	// Crucial: Use CUDA_CHECK for all CUDA API calls for proper error handling.
	CUDA_CHECK(cudaMalloc(&device_pos, 3 * N * sizeof(float)));
	CUDA_CHECK(cudaMalloc(&device_vel, 3 * N * sizeof(float)));
	CUDA_CHECK(cudaMalloc(&device_acc, 3 * N * sizeof(float)));
	CUDA_CHECK(cudaMalloc(&device_mass, N * sizeof(float)));
	CUDA_CHECK(cudaMalloc(&device_dens, N * sizeof(float)));
	
	// Este COPY ocurre step by step...
	// Copy data from host to device
	// cudaMemcpy copies data between host and device.
	// Arguments: destination, source, size, direction (host to device)
	cudaMemcpy(device_pos, h_position, 3*N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(device_vel, h_velocity, 3*N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(device_acc, h_acceleration, 3*N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(device_mass, h_mass, N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(device_dens, h_density, N * sizeof(float), cudaMemcpyHostToDevice);
	
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
	
	// Las nuevas!
	cudaMemcpyToSymbol(h_krnl2, &h_h_krnl2, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(h_krnl, &h_h_krnl, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(mHScaled9, &h_mHScaled9, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(mKernel1Scaled, &h_mKernel1Scaled, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(mKernel2Scaled, &h_mKernel2Scaled, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(mKernel3Scaled, &h_mKernel3Scaled, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(mRho0, &h_mRho0, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(mViscosityScalar, &h_mViscosityScalar, sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(mStiffness, &h_mStiffness, sizeof(float), 0, cudaMemcpyHostToDevice);


	// Define grid and block dimensions for kernel execution
	// A block is a group of threads that can cooperate.
	// A grid is a collection of blocks.
	// Here, we use 256 threads per block.
	// The number of blocks is calculated to cover all elements.
	int threadsPerBlock = 256;
	int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock; // Ceiling division
	
	// Launch the kernel(s)
	// 1st: neighbors!
	neighbors_CUDA<<<blocksPerGrid, threadsPerBlock>>>(device_pos, device_vel, device_acc,
														device_mass, device_dens);

	// Synchronize the device to ensure all kernel operations are complete
    // cudaDeviceSynchronize blocks the CPU until all GPU tasks are finished.
    CUDA_CHECK(cudaDeviceSynchronize());

	// 2nd: finish acc + integration
	integrate_CUDA<<<blocksPerGrid, threadsPerBlock>>>(device_pos, device_vel, device_acc);												
    
    // Synchronize the device to ensure all kernel operations are complete
    // cudaDeviceSynchronize blocks the CPU until all GPU tasks are finished.
    CUDA_CHECK(cudaDeviceSynchronize());


    // Copy data back from device to host
    // Arguments: destination, source, size, direction (device to host)
    CUDA_CHECK(cudaMemcpy(h_position, device_pos, 3 * N * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_velocity, device_vel, 3 * N * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_acceleration, device_acc, 3 * N * sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(h_mass, device_mass, N * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_density, device_dens, N * sizeof(float), cudaMemcpyDeviceToHost));

    // "h_position", etc. now contain the updated values from the GPU.

    // Free device memory
    // It's crucial to free allocated GPU memory to prevent memory leaks.
    CUDA_CHECK(cudaFree(device_pos));
    CUDA_CHECK(cudaFree(device_vel));
    CUDA_CHECK(cudaFree(device_acc));
	CUDA_CHECK(cudaFree(device_mass));
	CUDA_CHECK(cudaFree(device_dens));

}









