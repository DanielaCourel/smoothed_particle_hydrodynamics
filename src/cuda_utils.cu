// Claro jijodebú, necesito un file .cu...

#include <unistd.h>
#include <fstream>
#include <iostream> // Required for standard input/output operations (e.g., printf, cout)
#include <cuda_runtime.h> // Required for CUDA runtime API functions (e.g., cudaMalloc, cudaMemcpy, cudaFree)
#include <cmath>

#include "device_launch_parameters.h"  // ??? -> ask for it.
#include "cuda_utils.h"

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/sequence.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/iterator/constant_iterator.h>


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

int countStemp7;

// ------------------------------------- FUNCIÓN INCIALIZADORA -------------------------------------
// Función que inicializa los valores para no tener que traerlos del host cada vez
DeviceData* initDeviceData(float* h_position, float* h_velocity, float* h_acceleration,
                           float* h_mass, float* h_density, int h_cant_particles,
                           float h_scale, float h_softening, float h_grav_cte,
                           float h_mass_centre, float h_delta_step,
                           float h_x_centre, float h_y_centre, float h_z_centre,
                           float h_h_krnl, float h_h_krnl2, float h_mHScaled9,
                           float h_mKernel1Scaled, float h_mKernel2Scaled,
                           float h_mKernel3Scaled, float h_mRho0,
                           float h_mViscosityScalar, float h_mStiffness, int h_side_grid) {

	countStemp7 = 0;

    DeviceData* devData = (DeviceData*)malloc(sizeof(DeviceData));

	devData->num_cells = h_side_grid * h_side_grid * h_side_grid;

    CUDA_CHECK(cudaMemcpyToSymbol(cant_particles, &h_cant_particles, sizeof(int), 0, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpyToSymbol(x_centre, &h_x_centre, sizeof(float), 0, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpyToSymbol(y_centre, &h_y_centre, sizeof(float), 0, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpyToSymbol(z_centre, &h_z_centre, sizeof(float), 0, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpyToSymbol(scale, &h_scale, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(softening, &h_softening, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(grav_cte, &h_grav_cte, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(mass_centre, &h_mass_centre, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(delta_step, &h_delta_step, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(h_krnl, &h_h_krnl, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(h_krnl2, &h_h_krnl2, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(mHScaled9, &h_mHScaled9, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(mKernel1Scaled, &h_mKernel1Scaled, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(mKernel2Scaled, &h_mKernel2Scaled, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(mKernel3Scaled, &h_mKernel3Scaled, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(mRho0, &h_mRho0, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(mViscosityScalar, &h_mViscosityScalar, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(mStiffness, &h_mStiffness, sizeof(float), 0, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(side_grid, &h_side_grid, sizeof(int), 0, cudaMemcpyHostToDevice));



    size_t f_size = sizeof(float) * h_cant_particles * 3;
    size_t s_size = sizeof(float) * h_cant_particles;
    size_t i_size = sizeof(int) * h_cant_particles;
	

    CUDA_CHECK(cudaMalloc(&(devData->d_position), f_size));
    CUDA_CHECK(cudaMalloc(&(devData->d_velocity), f_size));
    CUDA_CHECK(cudaMalloc(&(devData->d_acceleration), f_size));
    CUDA_CHECK(cudaMalloc(&(devData->d_mass), s_size));
    CUDA_CHECK(cudaMalloc(&(devData->d_density), s_size));
    CUDA_CHECK(cudaMalloc(&(devData->global_index), i_size));
	CUDA_CHECK(cudaMalloc(&(devData->d_particle_ids), i_size));

    CUDA_CHECK(cudaMemcpy(devData->d_position, h_position, f_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devData->d_velocity, h_velocity, f_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devData->d_acceleration, h_acceleration, f_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devData->d_mass, h_mass, s_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devData->d_density, h_density, s_size, cudaMemcpyHostToDevice));

    return devData;
}


// -------------------------- FUNCIONES AUXILIARES -------------------------------------

__global__ void set_counts_from_unique_cells(
    const int* d_temp_cells,
    const int* d_temp_counts,
    int* d_counts,
    int num_unique)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < num_unique) {
        int cell_id = d_temp_cells[i];
        d_counts[cell_id] = d_temp_counts[i];
    }
}

void prepare_sorted_particles(const DeviceData* devData, int h_cant_particles, int num_cells, SortedParticlesData& sortedData)
{
    // 1. Copiar global_index para no modificar el original
    thrust::device_vector<int> global_index_copy(devData->global_index, devData->global_index + h_cant_particles);

    // 2. Inicializar IDs de partículas
    sortedData.d_particle_ids.resize(h_cant_particles);
    thrust::sequence(sortedData.d_particle_ids.begin(), sortedData.d_particle_ids.end());

    // 3. Ordenar por celda (sin modificar devData->global_index)
    thrust::sort_by_key(
        global_index_copy.begin(),
        global_index_copy.end(),
        sortedData.d_particle_ids.begin()
    );

    // 4. Guardar los cell_ids ya ordenados
    sortedData.d_sorted_cell_ids = global_index_copy;

    // 5. Reducir: contar partículas por celda
    thrust::device_vector<int> d_temp_counts(h_cant_particles);
    thrust::device_vector<int> d_temp_cells(h_cant_particles);

	auto new_end = thrust::reduce_by_key(
	    sortedData.d_sorted_cell_ids.begin(),
	    sortedData.d_sorted_cell_ids.end(),
	    thrust::make_constant_iterator(1),
	    d_temp_cells.begin(),
	    d_temp_counts.begin()
	);
	int num_unique_cells = new_end.first - d_temp_cells.begin();

	// 5b. Zeros para todas las celdas
	thrust::device_vector<int> d_counts(num_cells, 0);

	// 5c. Lanzar kernel para asignar counts en sus posiciones
	set_counts_from_unique_cells<<<(num_unique_cells + 255)/256, 256>>>(
	    thrust::raw_pointer_cast(d_temp_cells.data()),
	    thrust::raw_pointer_cast(d_temp_counts.data()),
	    thrust::raw_pointer_cast(d_counts.data()),
	    num_unique_cells
	);

	// 6. Scan
	sortedData.d_cell_start.resize(num_cells + 1, 0);
	thrust::exclusive_scan(
	    d_counts.begin(),
	    d_counts.end(),
	    sortedData.d_cell_start.begin()
	);
}



// --------------------------------------------------- KERNELS --------------------------------------------------------

// Voxelize, and then nighbors + dens & hydro...
__global__ void voxelize_CUDA(float* position, int num_cells, int* global_index)
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

	if (idx_x < 0) idx_x= 0;
	if (idx_y < 0) idx_y= 0;
	if (idx_z < 0) idx_z= 0;
	if (idx_x >= side_grid) idx_x= side_grid-1;
	if (idx_y >= side_grid) idx_y= side_grid-1;
	if (idx_z >= side_grid) idx_z= side_grid-1;


	// Ojo convencion de ejes, y es el vertical...
	int cell_index = idx_x + side_grid * idx_y + side_grid * side_grid * idx_z;

	global_index[tid] = cell_index;

	//printf("Thread %d: x_i = %.3f; idx_x = %d; lim_cell_0 = %.3f\n", tid, position[3*tid+0], idx_x, 2.f*h_krnl);
}


// Neigh + densities using the voxelization:
__global__ void neighbors_voxel_CUDA(float* position, float* mass, float* density,
								int* global_index,
                                 int* d_particle_ids, int* d_sorted_cell_ids, int* d_cell_start)
{
	// Thread ID - global!
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid >= cant_particles) return;

	// vars def thread-wise.
	float distance_ij3, rightPart, w;
	float rMinusRjScaled[3];
	// Init
	int count_neighb = 0;
	density[tid] = 0.f;

    int cell_id = global_index[tid];
    int start = d_cell_start[cell_id];
    int end = d_cell_start[cell_id + 1];

	// Main loop - Cada loop toca solo la parte del array con los index prev calculated,
	// haciéndose cargo de 1 partícula...; Quiero barrer la lista de indexes creada (!)
	int j;
	for (int j_this = start; j_this < end; ++j_this)
    {
		j = d_particle_ids[j_this];  // Por cosntruccion, (ieth_index == jeth_index)
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
                                 float* mass, float* density, int* global_index,
                                 int* d_particle_ids, int* d_sorted_cell_ids, int* d_cell_start, int* d_flag)
{
	// Thread ID - global!
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= cant_particles) return;

	// vars def thread-wise.
    float distance_ij3, distance, invDist;
    float pi, rhoiInv, rhoiInv2, piDivRhoi2;
    float centerPart, mj, pj, rhojInv, rhojInv2;
    float rMinusRjScaled[3];

    // Presión y masa propia
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

    int count_neighb = 0;

    // Obtener celda y vecinos
    int cell_id = global_index[tid];
    int start = d_cell_start[cell_id];
    int end = d_cell_start[cell_id + 1];

    for (int idx = start; idx < end; ++idx)
    {
        int j = d_particle_ids[idx];  // índice real (global) de la partícula vecina
        if (j == tid) continue;

        int jeth_index = global_index[j];  // celda de j

        // Nos aseguramos que sea la misma celda
        if (jeth_index != cell_id) {printf("NUNCA DEBERÍA LLEGAR ACÁ cell_id es %d y jeth_index es %d\n", cell_id, jeth_index); atomicExch(d_flag, 1); continue;}

        rMinusRjScaled[0] = (position[3*tid + 0] - position[3*j + 0]) * scale;
        rMinusRjScaled[1] = (position[3*tid + 1] - position[3*j + 1]) * scale;
        rMinusRjScaled[2] = (position[3*tid + 2] - position[3*j + 2]) * scale;

        distance_ij3 = rMinusRjScaled[0]*rMinusRjScaled[0] +
                       rMinusRjScaled[1]*rMinusRjScaled[1] +
                       rMinusRjScaled[2]*rMinusRjScaled[2];

        distance = sqrtf(distance_ij3);
        invDist = 1.f / (distance + pseudo_soft);

        // Presión del vecino
        mj = mass[j];
        pj = (density[j] - mRho0) * mStiffness;
        rhojInv = (density[j] > 0.0f) ? (1.0f / density[j]) : 1.0f;
        rhojInv2 = rhojInv * rhojInv;

        // Gradiente de presión
        centerPart = (h_krnl - distance) * (distance_ij3 < h_krnl2);
        centerPart *= centerPart * mj * piDivRhoi2 * pj * rhojInv2;

        pressureGradientContribution[0] = mKernel2Scaled * rMinusRjScaled[0] * invDist;
        pressureGradientContribution[1] = mKernel2Scaled * rMinusRjScaled[1] * invDist;
        pressureGradientContribution[2] = mKernel2Scaled * rMinusRjScaled[2] * invDist;

        pressureGradient[0] += pressureGradientContribution[0] * centerPart;
        pressureGradient[1] += pressureGradientContribution[1] * centerPart;
        pressureGradient[2] += pressureGradientContribution[2] * centerPart;

        // Viscosidad
        centerPart = (h_krnl - distance) * (distance_ij3 < h_krnl2);
        centerPart *= rhojInv * mj * mKernel3Scaled;

        viscousTerm[0] += (velocity[3*j + 0] - velocity[3*tid + 0]) * centerPart * mViscosityScalar * rhoiInv;
        viscousTerm[1] += (velocity[3*j + 1] - velocity[3*tid + 1]) * centerPart * mViscosityScalar * rhoiInv;
        viscousTerm[2] += (velocity[3*j + 2] - velocity[3*tid + 2]) * centerPart * mViscosityScalar * rhoiInv;

        count_neighb += (distance_ij3 < h_krnl2);
        if (count_neighb > 32) break;
    }

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

//-------------------------------------- FUNCIÓN PRINCIPAL Y DESTRUCTURA -------------------------------------
// Host-side wrapper function
// (Hace 1000 steps, sin ir y volver...)
void launchMyKernel(DeviceData* devData, float* h_position, float* h_velocity, float* h_acceleration, float* h_mass, float* h_density, int h_cant_particles)
{
									int h_flag = 0;
									int* d_flag;
									cudaMalloc(&d_flag, sizeof(int));
									cudaMemcpy(d_flag, &h_flag, sizeof(int), cudaMemcpyHostToDevice);
									countStemp7++;
				
    int threadsPerBlock = 256;  // fiducial = 256
    int blocksPerGrid = (h_cant_particles + threadsPerBlock - 1) / threadsPerBlock;

    for (int step=0; step<1000; step++) {

    // Paso 1: voxelize
    voxelize_CUDA<<<blocksPerGrid, threadsPerBlock>>>(devData->d_position, devData->num_cells, devData->global_index);
									//cudaError_t err = cudaDeviceSynchronize();
									//if (err != cudaSuccess) {
									//    printf("0 CUDA error after set_counts_from_unique_cells: %s\n", cudaGetErrorString(err));
									//    exit(1); // o assert(false);
									//}
    CUDA_CHECK(cudaDeviceSynchronize());

    // Paso 2: ordenar y preparar estructura por celda
    SortedParticlesData sortedData;
    prepare_sorted_particles(devData, h_cant_particles, devData->num_cells, sortedData);

    // Paso 3: vecinos (usa ids ordenados y cell_start)
    neighbors_voxel_CUDA<<<blocksPerGrid, threadsPerBlock>>>(devData->d_position, devData->d_mass, devData->d_density, devData->global_index,
															thrust::raw_pointer_cast(sortedData.d_particle_ids.data()), thrust::raw_pointer_cast(sortedData.d_sorted_cell_ids.data()), thrust::raw_pointer_cast(sortedData.d_cell_start.data()));
									//cudaError_t err1 = cudaDeviceSynchronize();
									//if (err1 != cudaSuccess) {
									//    printf("1 CUDA error after set_counts_from_unique_cells: %s\n", cudaGetErrorString(err1));
									//    exit(1); // o assert(false);
									//}
    CUDA_CHECK(cudaDeviceSynchronize());/* 
									if (countStemp7 == 1793) {
									    cudaMemcpy(h_position, devData->d_position, 3 * h_cant_particles * sizeof(float), cudaMemcpyDeviceToHost);
										int* global_index = (int*)malloc(h_cant_particles * sizeof(int));
									    cudaMemcpy(global_index, devData->global_index, h_cant_particles * sizeof(int), cudaMemcpyDeviceToHost);
									    std::ofstream f("positions_step3705.dump");
									    for (int i = 0; i < h_cant_particles; ++i) {
											if (14030==global_index[i]){
									        	f << h_position[3*i] << " " << h_position[3*i+1] << " " << h_position[3*i+2] << "\n";
									        	f << global_index[i] << "\n";
											}
									    }										

									} */


    // Paso 4: hydro
    hydro_voxel_CUDA<<<blocksPerGrid, threadsPerBlock>>>(devData->d_position, devData->d_velocity, devData->d_acceleration, devData->d_mass, devData->d_density, devData->global_index,
														   thrust::raw_pointer_cast(sortedData.d_particle_ids.data()), thrust::raw_pointer_cast(sortedData.d_sorted_cell_ids.data()), thrust::raw_pointer_cast(sortedData.d_cell_start.data()), d_flag);
									//cudaError_t err2 = cudaDeviceSynchronize();
									//if (err2 != cudaSuccess) {
									//    printf("2 CUDA error after set_counts_from_unique_cells: %s, in the step: %d \n", cudaGetErrorString(err2), countStemp7);
									//    exit(1); // o assert(false);
									//}

    CUDA_CHECK(cudaDeviceSynchronize());

									//cudaError_t err3 = cudaDeviceSynchronize();
									//if (err3 != cudaSuccess) {
									//    printf("3 CUDA error after set_counts_from_unique_cells: %s\n", cudaGetErrorString(err3));
									//    exit(1); // o assert(false);
									//}
									cudaMemcpy(&h_flag, d_flag, sizeof(int), cudaMemcpyDeviceToHost);
									if (h_flag == 1) {
									    printf("Empieza a fallar en el step %d, el máxima número de celdas debería ser %d \n\n", countStemp7, devData->num_cells);
										exit(1	);
									}

									cudaFree(d_flag);

    // Paso 5: integración
    integrate_CUDA<<<blocksPerGrid, threadsPerBlock>>>(devData->d_position, devData->d_velocity, devData->d_acceleration);
    CUDA_CHECK(cudaDeviceSynchronize());

    }  // End itengration

    // Paso 6: copiar al host
    CUDA_CHECK(cudaMemcpy(h_position, devData->d_position, 3 * h_cant_particles * sizeof(float), cudaMemcpyDeviceToHost));
}


void cleanupDeviceData(DeviceData* devData) 
{
    cudaFree(devData->d_position);
    cudaFree(devData->d_velocity);
    cudaFree(devData->d_acceleration);
    cudaFree(devData->d_mass);
    cudaFree(devData->d_density);
	cudaFree(devData->global_index);

    free(devData);
}







