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

// New cte, size of the grid (!):
__constant__ int side_grid;

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

/*

// Re-do: launcheo muchos blocks buscando vecinos de UNA partícula a la vez,
// no estoy encontrando casi ningún vecino con el actual approach. Luego, el kernel
// se come la posición de la partícula i-ésima y barre de a muchas j-ésimas en simultáneo.
// Es muchísimo más parecido al cómputo de la dist al BH i.e., no shenanigans block-wise (!)

MUY lenta -> Volvamos a la voxelize() strategy but implemented in CUDA (homemade)

__global__ void neighbors_chunked_CUDA(
	float* position, float* mass, float* density,
	float x_i, float y_i, float z_i,
	float* density_i, int* cant_neighb)
{
	// Thread ID - global!
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= cant_particles) return;

	// Within a block, min(tid) = 0; max(tid) = dim(block)
	int dim_block = blockDim.x;
	// Thread index within the block
    int localThreadId = threadIdx.x;

	// Only 1 big extern shared:
    extern __shared__ char shared_buffer[];

	// Offsets: Calculate the size of each component in bytes:
    size_t vector_bytes = 3 * dim_block * sizeof(float);
    size_t scalar_bytes = dim_block * sizeof(float);

	// Calculate start addresses for each array within the buffer
    // The first array starts at the beginning
    char* current_ptr = shared_buffer;
    float* s_localPositions = (float*)current_ptr;
    current_ptr += vector_bytes; // Advance pointer past positions
    float* s_localMasses = (float*)current_ptr;
    current_ptr += scalar_bytes; // Advance pointer past masses
    float* s_localDensities = (float*)current_ptr;

	s_localPositions[3*localThreadId + 0] = position[3*tid + 0];
	s_localPositions[3*localThreadId + 1] = position[3*tid + 1];
	s_localPositions[3*localThreadId + 2] = position[3*tid + 2];

	s_localMasses[localThreadId] = mass[tid];
	s_localDensities[localThreadId] = 0.f;

    __syncthreads(); // Ensure all shared memory loads are complete before any thread uses them

	// TIENE que haber un if acá arriba de todo para que no se sigan ejecutando bloques de más:
	if (*cant_neighb > 32) return;

	// vars def thread-wise.
	float distance_ij3, rightPart, w;
	float rMinusRjScaled[3];	

	// Main loop - Cada hilo se encarga de 1 partícula (chunked):
	rMinusRjScaled[0] = (s_localPositions[3*localThreadId + 0] - x_i) * scale;
	rMinusRjScaled[1] = (s_localPositions[3*localThreadId + 1] - y_i) * scale;
	rMinusRjScaled[2] = (s_localPositions[3*localThreadId + 2] - z_i) * scale;

	distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] +\
					rMinusRjScaled[1] * rMinusRjScaled[1] +\
					rMinusRjScaled[2] * rMinusRjScaled[2] + softening;  // This is SQUARED

	// From this, I can (branchless-ly) compute hydro properties (!)
	// Density:
	rightPart = h_krnl2 - distance_ij3;
	rightPart = rightPart * rightPart * rightPart;
	w = mKernel1Scaled * rightPart * (distance_ij3 < h_krnl2);  // true => 1, false => 0 (!)
	// apply weighted neighbor mass to our density
	s_localDensities[localThreadId] += (s_localMasses[localThreadId] * w);  // 0. if (d^2 >= h^2)
	// Atomically add each density to the total sum to the global variable:
    atomicAdd(density_i, s_localDensities[localThreadId]);

	// Ahora si uso warp-wise functions!!!

	// Perform a warp-wide ballot to see which threads meet the condition
	// (true (non-zero) or false (zero) for each thread).
	// The result 'condition_mask' is a 32-bit integer where bit 'i' is set
	// if thread 'i' in the warp met the condition.

	// Mask of active threads in the warp
	unsigned int active_threads_mask = __activemask();
	// Def the condition:
	bool valid_neighbor = (distance_ij3 < h_krnl2);
	// Ask them:
	unsigned int condition_mask = __ballot_sync(active_threads_mask, valid_neighbor);
	// __popc() to count how many threads in this warp met the condition:
	int num_threads_meeting_condition = __popc(condition_mask);

	// Only the warp leader count its private counter!
	// =/= a LocalID == 0 !!!
	if (tid % 32 == 0)
	{
		// Atomically add to a global or shared counter for this shared target particle
		// This ensures the sum for the warp's contributions is added only once
		atomicAdd(cant_neighb, num_threads_meeting_condition);
	}

	// It's crucial to synchronize threads *after* the atomicAdd if *all* threads
	// need to see the *updated* shared count before the next iteration of the loop.
	// Without this, some threads might read an old value of s_shared_neighbor_count[LocalID]
	// at the start of the next iteration, potentially missing the break condition.
	__syncthreads(); // Synchronize all threads in the block to see the updated count

}

*/

// Voxelize, and then nighbors + dens & hydro...
__global__ void voxelize_CUDA(float* position, int* global_index)
{
	// Thread ID - global!
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid >= cant_particles) return;

	// Ojo, tengo que contar el tamaño de las celdas! (como en la simu):
	float cell_size = 2.0f * h_krnl;
	//float len_box = side_grid * (2.f * h_krnl);	

	// 1st, compute the grid_index of this cell:
	int idx_x = floor(position[3*tid + 0]/cell_size);
	int idx_y = floor(position[3*tid + 1]/cell_size);
	int idx_z = floor(position[3*tid + 2]/cell_size);

	// Ojo convencion de ejes, y es el vertical...
	int cell_index = idx_x + side_grid * idx_y + side_grid * side_grid * idx_z;

	global_index[tid] = cell_index;

	//printf("Thread %d: x_i = %.3f; idx_x = %d; lim_cell_0 = %.3f\n", tid, position[3*tid+0], idx_x, 2.f*h_krnl);
}


// Neigh + densities using the voxelization:
__global__ void neighbors_voxel_CUDA(float* position, float* mass, float* density,
								int* global_index)
{
	// Thread ID - global!
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid >= cant_particles) return;

	// vars def thread-wise.
	float distance_ij3, rightPart, w;
	float rMinusRjScaled[3];
	// Init
	int count_neighb = 0;
	int ieth_index, jeth_index;
	density[tid] = 0.f;
	ieth_index = global_index[tid];

	// OJO: Para la busqueda de vecinos, no quiero que toque potencialmente a todas las
	// particulas del sistema, así que 1ro cuento cuántas hay en su misma celda, y barro 
	// sólo esas!
	int shared_cell = 0;
	const int max_neighb = 512*4;  // Polemiquisimo...
	int list_neighb_idx[max_neighb];  // Polemiquisimo...
	for (int j=0; j < cant_particles; j++)
	{
		jeth_index = global_index[j];
		if (ieth_index == jeth_index)
		{
			list_neighb_idx[shared_cell] = j;  // Populo y sigo...
			shared_cell++;
		}
		//else {printf("Thread %d: dif cell = %d; current %d\n", tid, j, shared_cell);}

		if (shared_cell == max_neighb) break;
	}

	// Main loop - Cada loop toca solo la parte del array con los index prev calculated,
	// haciéndose cargo de 1 partícula...; Quiero barrer la lista de indexes creada (!)
	int j;
	for (int j_this = 0; j_this < shared_cell; j_this++)
	{
		j = list_neighb_idx[j_this];  // Por cosntruccion, (ieth_index == jeth_index)
		// j ahora es el indice GLOBAL
		// -> Que j =/= tid!!!
		if (j == tid) continue;
		
		rMinusRjScaled[0] = (position[3*tid + 0] - position[3*j + 0]) * scale;
		rMinusRjScaled[1] = (position[3*tid + 1] - position[3*j + 1]) * scale;
		rMinusRjScaled[2] = (position[3*tid + 2] - position[3*j + 2]) * scale;

		// -> Sin softening!
		distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] +\
						rMinusRjScaled[1] * rMinusRjScaled[1] +\
						rMinusRjScaled[2] * rMinusRjScaled[2];  // This is SQUARED

		// From this, I can (branchless-ly) compute hydro properties (!)
		// Density:
		rightPart = h_krnl2 - distance_ij3;
		rightPart = rightPart * rightPart * rightPart;
		w = mKernel1Scaled * rightPart * (distance_ij3 < h_krnl2);  // true => 1, false => 0 (!)
		// apply weighted neighbor mass to our density
		density[tid] += (mass[tid] * w);  // 0. if (d^2 >= h^2)

		count_neighb += 1 * (distance_ij3 < h_krnl2);  // 0. if (d^2 >= h^2)

		/*
		if (threadIdx.x == 0)
		{
			//printf("Thread %d: idx_x,i = %d; idx_x,j = %d\n", tid, ieth_index, global_index[j]);
			printf("Thread %d: distance_ij = %.3f; h^2 = %.3f\n", tid, distance_ij3, h_krnl2);
			printf("Thread %d: count_neighb = %d\n", tid, count_neighb);
		}
		*/

		if (count_neighb > 32) break;

		//if (distance_ij3 < h_krnl2) {printf("Thread %d: neighb found = %d\n", tid, count_neighb);}
	}

}


// Hydro using the voxelization (re-computo distancia, no hay con qué darle por ahora...)
__global__ void hydro_voxel_CUDA(float* position, float* velocity, float* acceleration,
							float* mass, float* density, int* global_index)
{
	// Thread ID - global!
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid >= cant_particles) return;

	// vars def thread-wise.
	float distance_ij3, distance, invDist;
	float pi, rhoiInv, rhoiInv2, piDivRhoi2;
	float centerPart, mj, pj, rhojInv, rhojInv2;
	float rMinusRjScaled[3];

	// Init
	int count_neighb = 0;
	int ieth_index, jeth_index;
	ieth_index = global_index[tid];

	// OJO: Para la busqueda de vecinos, no quiero que toque potencialmente a todas las
	// particulas del sistema, así que 1ro cuento cuántas hay en su misma celda, y barro 
	// sólo esas!
	int shared_cell = 0;
	const int max_neighb = 512*4;  // Polemiquisimo...
	int list_neighb_idx[max_neighb];  // Polemiquisimo...
	for (int j=0; j < cant_particles; j++)
	{
		jeth_index = global_index[j];
		if (ieth_index == jeth_index)
		{
			list_neighb_idx[shared_cell] = j;  // Populo y sigo...
			shared_cell++;
		}

		if (shared_cell == max_neighb) break;
	}

	pi = (density[tid] - mRho0) * mStiffness;  // One-liner...
	rhoiInv = ((pi > 0.0f) ? (1.0f / pi) : 1.0f);  // Better?
	rhoiInv2 = rhoiInv * rhoiInv;
	piDivRhoi2 = pi * rhoiInv2;

	// Pressure gradient
	float pressureGradient[3] = {0.0f, 0.0f, 0.0f};
	float pressureGradientContribution[3];
	// Viscous term
	float viscousTerm[3] = {0.0f, 0.0f, 0.0f};

	// A pseudo-softening:
	float pseudo_soft = 1e-3f;

	// Main loop - Cada loop toca solo la parte del array con los index prev calculated,
	// haciéndose cargo de 1 partícula...; Quiero barrer la lista de indexes creada (!)
	int j;
	for (int j_this = 0; j_this < shared_cell; j_this++)
	{
		j = list_neighb_idx[j_this];  // Por construccion, (ieth_index == jeth_index)
		// j ahora es el indice GLOBAL
		// j =/= tid!!!
		if (j == tid) continue;
		
		rMinusRjScaled[0] = (position[3*tid + 0] - position[3*j + 0]) * scale;
		rMinusRjScaled[1] = (position[3*tid + 1] - position[3*j + 1]) * scale;
		rMinusRjScaled[2] = (position[3*tid + 2] - position[3*j + 2]) * scale;

		// -> Sin softening!
		distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] +\
						rMinusRjScaled[1] * rMinusRjScaled[1] +\
						rMinusRjScaled[2] * rMinusRjScaled[2];  // This is SQUARED

		distance = sqrtf(distance_ij3);  // This is the true d_ij
		//invDist = rsqrtf(distance_ij3);  // quick x^(-1/2)
		// Cambio esto por un ~softening:
		invDist = 1.f/(distance + pseudo_soft);

		// Pressure gradient:
		mj = mass[j];
		pj = (density[j] - mRho0) * mStiffness;
		rhojInv = ((density[j] > 0.0f) ? (1.0f / density[j]) : 1.0f);
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
		// J - TID !!!!!!!!!!!
		viscousTerm[0] += (velocity[3*j + 0] - velocity[3*tid + 0]) *\
						centerPart * mViscosityScalar * rhoiInv;
		viscousTerm[1] += (velocity[3*j + 1] - velocity[3*tid + 1]) *\
						centerPart * mViscosityScalar * rhoiInv;
		viscousTerm[2] += (velocity[3*j + 2] - velocity[3*tid + 2]) *\
						centerPart * mViscosityScalar * rhoiInv;

		// Importante...
		count_neighb += 1 * (distance_ij3 < h_krnl2);  // 0. if (d^2 >= h^2)

		if (count_neighb > 32) break;
	}

	// Finish parte hydro, aparte haremos la accel central desp...
	// Every thread of the block writes in its own place the results from shared memory
    // to global memory. (no need for sync!)
	acceleration[3*tid + 0] = viscousTerm[0] - pressureGradient[0];
	acceleration[3*tid + 1] = viscousTerm[1] - pressureGradient[1];
	acceleration[3*tid + 2] = viscousTerm[2] - pressureGradient[2];

}



// So: 1st neigh + densities; 2nd hydro; 3rd grav + integ
__global__ void neighbors_CUDA(float* position, float* mass, float* density)
{
	// Thread ID - global!
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid >= cant_particles) return;

	// Within a block, min(tid) = 0; max(tid) = dim(block)
	int dim_block = blockDim.x;

	// 1st, necesito allocar la memoria que solo va a tocar este bloque (contigua!):
	// Each thread helps load a local reference.
	// __shared__ float s_localPositions[3 * dim_block];  // (x,y,z)
	// __shared__ float s_localMasses[dim_block];
	// __shared__ float s_localDensities[dim_block];

	// ...but only 1 big extern shared:
    extern __shared__ char shared_buffer[];

	// Offsets: Calculate the size of each component in bytes:
    size_t vector_bytes = 3 * dim_block * sizeof(float);
    size_t scalar_bytes = dim_block * sizeof(float);

	// Calculate start addresses for each array within the buffer
    // The first array starts at the beginning
    char* current_ptr = shared_buffer;
	// s_localPositions starts at the beginning
    float* s_localPositions = (float*)current_ptr;
    current_ptr += vector_bytes; // Advance pointer past positions
    float* s_localMasses = (float*)current_ptr;
    current_ptr += scalar_bytes; // Advance pointer past masses
    float* s_localDensities = (float*)current_ptr;

	// Now I can use s_localPositions, s_localMasses and s_localDensities
    // just like regular arrays within the kernel...

	// I have to create an array of "neighbor_counts" to allow the warp to see the leader's neighbcount!!!
	// ...within the block, as always.
	current_ptr += scalar_bytes; // Advance pointer past densities
    int* s_neighbor_count = (int*)current_ptr;

	// Thread index within the block
    int localThreadId = threadIdx.x;

	s_neighbor_count[localThreadId] = 0;  // Now, all threads within a warp can access this.

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

	s_localMasses[localThreadId] = mass[tid];
	s_localDensities[localThreadId] = 0.f;

    __syncthreads(); // Ensure all shared memory loads are complete before any thread uses them

	// vars def thread-wise.
	float distance_ij3, rightPart, w;
	float rMinusRjScaled[3];	

	// Main loop - block-wise (Ojo iters!)
	// Acá deberían ir de a 32 threads...  -> ¡Falso! Hay 32 threads que están haciendo esto
	// CON SUS PROPIAS PARTÍCULAS !!!
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

		// From this, I can (branchless-ly) compute hydro properties (!)
		// Density:
		rightPart = h_krnl2 - distance_ij3;
		rightPart = rightPart * rightPart * rightPart;
		w = mKernel1Scaled * rightPart * (distance_ij3 < h_krnl2);  // true => 1, false => 0 (!)
		// apply weighted neighbor mass to our density
		s_localDensities[localThreadId] += (s_localMasses[localThreadId] * w);  // 0. if (d^2 >= h^2)

		// Acá puedo hacer un mask (ya que en este loop vienen actuando warps de 32-threads)
		// y preg cuáles de ellos suman a valid neighbors:

		/*

		Como este loop lo está haciendo cada thread (y no de a warps), esto ya no sirve...

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

		*/

		// Then, back to the basics:
		s_neighbor_count[localThreadId] += 1 * (distance_ij3 < h_krnl2);

	}

	// Finish parte hydro, aparte haremos la accel central desp...
	// Assert que todos los del bloque terminaron de escribir previous to write globally!
	__syncthreads();

	// Every thread of the block writes in its own place the results from shared memory
    // to global memory.
	density[tid] = s_localDensities[localThreadId];

}

// Ahora si: 2nd hydro (re-computo distancia, no hay con qué darle por ahora...)
__global__ void hydro_CUDA(float* position, float* velocity, float* acceleration,
							float* mass, float* density)
{
	// Thread ID - global!
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid >= cant_particles) return;

	// Within a block, min(tid) = 0; max(tid) = dim(block)
	int dim_block = blockDim.x;

	// 1st, necesito allocar la memoria que solo va a tocar este bloque (contigua!):
	// Each thread helps load a local reference.
    /*
	__shared__ float s_localPositions[3 * dim_block];  // (x,y,z)
	__shared__ float s_localVelocities[3 * dim_block];  // (vx,vy,vz)
	__shared__ float s_localMasses[dim_block];
	__shared__ float s_localDensities[dim_block];
	*/
	
	// ...but only 1 big extern shared:
    extern __shared__ char shared_buffer[];

	// Offsets: Calculate the size of each component in bytes:
    size_t vector_bytes = 3 * dim_block * sizeof(float);
    size_t scalar_bytes = dim_block * sizeof(float);

	// Calculate start addresses for each array within the buffer
    // The first array starts at the beginning
    char* current_ptr = shared_buffer;
	// s_localPositions starts at the beginning
    float* s_localPositions = (float*)current_ptr;
    current_ptr += vector_bytes; // Advance pointer past positions
	float* s_localVelocities = (float*)current_ptr;
    current_ptr += vector_bytes; // Advance pointer past velocities
    float* s_localMasses = (float*)current_ptr;
    current_ptr += scalar_bytes; // Advance pointer past masses
    float* s_localDensities = (float*)current_ptr;
	current_ptr += scalar_bytes; // Advance pointer past densities
    int* s_neighbor_count = (int*)current_ptr;

	// I have to create an array of "neighbor_counts" to allow the warp to see the leader's neighbcount!!!
	// ...within the block, as always.
	/* __shared__ int s_neighbor_count[dim_block]; */

	// Thread index within the block
    int localThreadId = threadIdx.x;

	s_neighbor_count[localThreadId] = 0;  // Now, all threads within a warp can access this.

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

    __syncthreads(); // Ensure all shared memory loads are complete before any thread uses them

	// vars def thread-wise.
	float distance_ij3, distance, invDist;
	float pi, rhoiInv, rhoiInv2, piDivRhoi2;
	float centerPart, mj, pj, rhojInv, rhojInv2;
	float rMinusRjScaled[3];

	pi = (s_localDensities[localThreadId] - mRho0) * mStiffness;  // One-liner...
	rhoiInv = ((pi > 0.0f) ? (1.0f / pi) : 1.0f);  // Better?
	rhoiInv2 = rhoiInv * rhoiInv;
	piDivRhoi2 = pi * rhoiInv2;

	// Pressure gradient
	float pressureGradient[3] = {0.0f, 0.0f, 0.0f};
	float pressureGradientContribution[3];
	// Viscous term
	float viscousTerm[3] = {0.0f, 0.0f, 0.0f};

	// Main loop - block-wise (Ojo iters!)
	// Acá deberían ir de a 32 threads...  -> ¡Falso! Hay 32 threads que están haciendo esto
	// CON SUS PROPIAS PARTÍCULAS !!!
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

		/*

		Como este loop lo está haciendo cada thread (y no de a warps), esto ya no sirve...

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

		*/

		// Then, back to the basics:
		s_neighbor_count[localThreadId] += 1 * (distance_ij3 < h_krnl2);

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


// Lastly: grav + integration (!)
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

	// Softening va squared!
	distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] + rMinusRjScaled[1] * rMinusRjScaled[1] +\
				rMinusRjScaled[2] * rMinusRjScaled[2] + (softening * softening);

	invDist = rsqrtf(distance_ij3);  // quick x^(-1/2)
	invDist = invDist * invDist * invDist;  // dist^-3

	// Updateo la gravedad:
	// acceleration += gravityTerm;  -> Escribo mal a proposito, just to check (es +=)
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

	// Idem, soft^2...
	distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] + rMinusRjScaled[1] * rMinusRjScaled[1] +\
				rMinusRjScaled[2] * rMinusRjScaled[2] + (softening * softening);

	invDist = rsqrtf(distance_ij3);  // quick x^(-1/2)
	invDist = invDist * invDist * invDist;  // dist^-3

	// Updateo la gravedad (de 0!):
	acceleration[3*tid + 0] = -grav_cte * mass_centre * (rMinusRjScaled[0] * invDist);
	acceleration[3*tid + 1] = -grav_cte * mass_centre * (rMinusRjScaled[1] * invDist);
	acceleration[3*tid + 2] = -grav_cte * mass_centre * (rMinusRjScaled[2] * invDist);

	
	// OJO! Me falto check for CFL condition:
	dot = (acceleration[3*tid +0] * acceleration[3*tid +0]) + (acceleration[3*tid +1] * acceleration[3*tid +1]) +\
		(acceleration[3*tid +2] * acceleration[3*tid +2]);

	// Lo meto a dedo... (recall branchles...)
	acceleration[3*tid +0] *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	acceleration[3*tid +1] *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	acceleration[3*tid +2] *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));


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
				float h_mViscosityScalar, float h_mStiffness, int h_side_grid)
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

	// Lastly, the side of the grid! (I forgot it pls...)
	cudaMemcpyToSymbol(side_grid, &h_side_grid, sizeof(float), 0, cudaMemcpyHostToDevice);

	// Define grid and block dimensions for kernel execution
	// A block is a group of threads that can cooperate.
	// A grid is a collection of blocks.
	// Here, we use 256 threads per block.
	// The number of blocks is calculated to cover all elements.
	int threadsPerBlock = 256;
	int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock; // Ceiling division

	// Alloco la memoria del buffer (!)
	// Calculate component sizes in bytes
	//size_t vector_bytes = 3 * threadsPerBlock * sizeof(float);
	//size_t scalar_bytes = threadsPerBlock * sizeof(float);

	/*

	// Calculate total required shared memory bytes (including padding if applied)
	// pos, dens, mass, neighb
	size_t total_shared_mem_bytes_neighbors = vector_bytes + scalar_bytes + scalar_bytes + scalar_bytes;
	// pos, vel, dens, mass, neighb
	size_t total_shared_mem_bytes_hydro = vector_bytes + vector_bytes + scalar_bytes + scalar_bytes + scalar_bytes;
	// integrate no usa buffers (!)

	*/
	// --------------------------------------------------------
	/*

	Ahora si encuentra a los vecinos, pero se tarda x100 (!!!)

	// Try loop (for each particle) -> iterá sobre chunks de particles, pero sincronizá y cortá entre medio:
	// "chunk" == bloque!

	// OJO, necesito allocar las variables que voy a modificar dentro d eeste loop:
	// Device memory for accumulating results for the *current reference particle (i)*
    // These are single float/int values on the device, used for atomic adds
    float* d_current_density_contribution_ptr;
    int* d_current_neighb_count_ptr;
    CUDA_CHECK(cudaMalloc(&d_current_density_contribution_ptr, sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_current_neighb_count_ptr, sizeof(int)));

	for (int i = 0; i < N; i++)
	{
		float x_i = h_position[3*i + 0];  // Much slower?
		float y_i = h_position[3*i + 1];
		float z_i = h_position[3*i + 2];

		// Initialize the device-side accumulation variables for THIS particle 'i'
		float density_i = 0.f;
		int cant_neighb = 0;
		CUDA_CHECK(cudaMemcpy(d_current_density_contribution_ptr, &density_i, sizeof(float), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_current_neighb_count_ptr, &cant_neighb, sizeof(int), cudaMemcpyHostToDevice));

		// Quiero que cada bloque toque su propia memoria -> Vuelta al buffer
		neighbors_chunked_CUDA<<<blocksPerGrid, threadsPerBlock, total_shared_mem_bytes_neighbors>>>(
								device_pos, device_mass, device_dens,
								x_i, y_i, z_i,
								d_current_density_contribution_ptr,
								d_current_neighb_count_ptr);
								// Estos ults son atomicos
		
		// Hacé que todos se pongan al tanto de lo que pasa globalmente
		CUDA_CHECK(cudaDeviceSynchronize());

		// Copy the accumulated results for particle 'i' back to host
		float h_accumulated_density_i;
		int h_accumulated_neighb_count_i;
		CUDA_CHECK(cudaMemcpy(&h_accumulated_density_i, d_current_density_contribution_ptr, sizeof(float), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(&h_accumulated_neighb_count_i, d_current_neighb_count_ptr, sizeof(int), cudaMemcpyDeviceToHost));

		// Update the device_dens array for particle 'i' immediately.
        CUDA_CHECK(cudaMemcpy(&device_dens[i], &h_accumulated_density_i, sizeof(float), cudaMemcpyHostToDevice));

		// Check the exit condition based on accumulated neighbors
		if (h_accumulated_neighb_count_i > 32) continue;  // Next particle

	};

	*/

	// New, ft. voxelize() -> No necesito allcoar memoria (no uso buffers!). Pero si el global index
    int* d_global_index;
    CUDA_CHECK(cudaMalloc(&d_global_index, sizeof(int)));
	
	// Launch the kernel(s)
	// 1st: voxelize()
	voxelize_CUDA<<<blocksPerGrid, threadsPerBlock>>>(device_pos, d_global_index);
	// Synchronize the device to ensure all kernel operations are complete
    // cudaDeviceSynchronize blocks the CPU until all GPU tasks are finished.
    CUDA_CHECK(cudaDeviceSynchronize());

	// neighbors + density!
	neighbors_voxel_CUDA<<<blocksPerGrid, threadsPerBlock>>>(device_pos, device_mass, device_dens, d_global_index);
	// Synchronize the device to ensure all kernel operations are complete
    // cudaDeviceSynchronize blocks the CPU until all GPU tasks are finished.
    CUDA_CHECK(cudaDeviceSynchronize());

	// hydro!
	hydro_voxel_CUDA<<<blocksPerGrid, threadsPerBlock>>>(device_pos, device_vel, device_acc,
														device_mass, device_dens, d_global_index);
	// Synchronize the device to ensure all kernel operations are complete
    // cudaDeviceSynchronize blocks the CPU until all GPU tasks are finished.
    CUDA_CHECK(cudaDeviceSynchronize());

	/*
	// 2nd: hydro!
	hydro_CUDA<<<blocksPerGrid, threadsPerBlock, total_shared_mem_bytes_hydro>>>(device_pos, device_vel, device_acc,
														device_mass, device_dens);
	

	// Synchronize the device to ensure all kernel operations are complete
    // cudaDeviceSynchronize blocks the CPU until all GPU tasks are finished.
    CUDA_CHECK(cudaDeviceSynchronize());
	*/														

	// Finish acc + integration
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

    // Free device memory
    // It's crucial to free allocated GPU memory to prevent memory leaks.
    CUDA_CHECK(cudaFree(device_pos));
    CUDA_CHECK(cudaFree(device_vel));
    CUDA_CHECK(cudaFree(device_acc));
	CUDA_CHECK(cudaFree(device_mass));
	CUDA_CHECK(cudaFree(device_dens));

	// Free the newly created for the inner loop:
	//CUDA_CHECK(cudaFree(d_current_density_contribution_ptr));
	//CUDA_CHECK(cudaFree(d_current_neighb_count_ptr));

	CUDA_CHECK(cudaFree(d_global_index));

}









