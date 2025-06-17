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

int countStemp;

// ------------------------------------- FUNCIÓN INCIALIZADORA -------------------------------------
// Función que inicializa los valores para no tener que traerlos del host cada vez
DeviceData* initDeviceData(float* h_position_x, float* h_position_y, float* h_position_z, 
                    	   float* h_velocity_x, float* h_velocity_y, float* h_velocity_z,
                    	   float* h_acceleration_x, float* h_acceleration_y, float* h_acceleration_z,
                           float* h_mass, float* h_density, int h_cant_particles,
                           float h_scale, float h_softening, float h_grav_cte,
                           float h_mass_centre, float h_delta_step,
                           float h_x_centre, float h_y_centre, float h_z_centre,
                           float h_h_krnl, float h_h_krnl2, float h_mHScaled9,
                           float h_mKernel1Scaled, float h_mKernel2Scaled,
                           float h_mKernel3Scaled, float h_mRho0,
                           float h_mViscosityScalar, float h_mStiffness, int h_side_grid) {

	countStemp = 0;

    DeviceData* devData = (DeviceData*)malloc(sizeof(DeviceData));

	devData->num_cells = h_side_grid * h_side_grid * h_side_grid;

	// Inicializar buffers de sorted_data
	devData->sorted_data.d_particle_ids.resize(h_cant_particles);
	devData->sorted_data.d_sorted_cell_ids.resize(h_cant_particles);
	devData->sorted_data.d_cell_start.resize(devData->num_cells + 1, 0);


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



    size_t s_size = sizeof(float) * h_cant_particles;
    size_t i_size = sizeof(int) * h_cant_particles;
	
	std::vector<float4> h_pos(h_cant_particles);
	std::vector<float4> h_vel(h_cant_particles);
	std::vector<float4> h_accel(h_cant_particles);
	
	for (int i = 0; i < h_cant_particles; ++i) {
	    h_pos[i] = make_float4(h_position_x[i], h_position_y[i], h_position_z[i], 0.0f);
	    h_vel[i] = make_float4(h_velocity_x[i], h_velocity_y[i], h_velocity_z[i], 0.0f);
	    h_accel[i] = make_float4(h_acceleration_x[i], h_acceleration_y[i], h_acceleration_z[i], 0.0f);
	}

	CUDA_CHECK(cudaMalloc(&(devData->d_pos), sizeof(float4) * h_cant_particles));
	CUDA_CHECK(cudaMalloc(&(devData->d_vel), sizeof(float4) * h_cant_particles));
	CUDA_CHECK(cudaMalloc(&(devData->d_accel), sizeof(float4) * h_cant_particles));

	CUDA_CHECK(cudaMemcpy(devData->d_pos, h_pos.data(), sizeof(float4) * h_cant_particles, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(devData->d_vel, h_vel.data(), sizeof(float4) * h_cant_particles, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(devData->d_accel, h_accel.data(), sizeof(float4) * h_cant_particles, cudaMemcpyHostToDevice));


    CUDA_CHECK(cudaMalloc(&(devData->d_mass), s_size));
    CUDA_CHECK(cudaMalloc(&(devData->d_density), s_size));
    CUDA_CHECK(cudaMalloc(&(devData->global_index), i_size));
	CUDA_CHECK(cudaMalloc(&(devData->d_particle_ids), i_size));
	

    CUDA_CHECK(cudaMemcpy(devData->d_mass, h_mass, s_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devData->d_density, h_density, s_size, cudaMemcpyHostToDevice));

    return devData;
}


// -------------------------- FUNCIONES AUXILIARES -------------------------------------

__device__ int hash(int x) 
{
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

__global__ void set_counts_from_unique_cells(const int* d_temp_cells, const int* d_temp_counts, int* d_counts, int num_unique)
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
	thrust::sequence(sortedData.d_particle_ids.begin(), sortedData.d_particle_ids.end());


    // 3. Ordenar por celda (sin modificar devData->global_index)
    thrust::sort_by_key(
        global_index_copy.begin(),
        global_index_copy.end(),
        sortedData.d_particle_ids.begin()
    );

    // 4. Guardar los cell_ids ya ordenados
    thrust::copy(global_index_copy.begin(), global_index_copy.end(), sortedData.d_sorted_cell_ids.begin());


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
	thrust::exclusive_scan(
	    d_counts.begin(),
	    d_counts.end(),
	    sortedData.d_cell_start.begin()
	);

}



// --------------------------------------------------- KERNELS --------------------------------------------------------

// Voxelize, and then nighbors + dens & hydro...
__global__ void voxelize_CUDA(float4* pos, int num_cells, int* global_index)
{
	// Thread ID - global!
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid >= cant_particles) return;

	// Ojo, tengo que contar el tamaño de las celdas! (como en la simu):
	float cell_size = 2.0f * h_krnl;
	//float len_box = side_grid * (2.f * h_krnl);	

	// 1st, compute the grid_index of this cell:
	int idx_x = floor(pos[tid].x/cell_size);
	int idx_y = floor(pos[tid].y/cell_size);
	int idx_z = floor(pos[tid].z/cell_size);

	if (idx_x < 0) idx_x= 0;
	if (idx_y < 0) idx_y= 0;
	if (idx_z < 0) idx_z= 0;
	if (idx_x >= side_grid) idx_x= side_grid-1;
	if (idx_y >= side_grid) idx_y= side_grid-1;
	if (idx_z >= side_grid) idx_z= side_grid-1;


	// Ojo convencion de ejes, y es el vertical...
	int cell_index = idx_x + side_grid * idx_y + side_grid * side_grid * idx_z;

	global_index[tid] = cell_index;

}


// Neigh + densities using the voxelization:
__global__ void neighbors_voxel_CUDA(float4* pos, float* mass, float* density,
                                     int* global_index, int* d_particle_ids,
                                     int* d_sorted_cell_ids, int* d_cell_start, int step)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= cant_particles) return;

	float4 pos_i = __ldg(&pos[tid]);
	int cell_id = global_index[tid];

	int start = d_cell_start[cell_id];
	int end   = d_cell_start[cell_id + 1];
	int len   = end - start;

	if (len <= 1) return;

	int offset = hash(step) % len;

	float dens = 0.0f;
	int count_neighb = 0;

	for (int i = 0; i < len; ++i) {
		int j_this = start + ((i + offset) % len);
		int j = d_particle_ids[j_this];
		if (j == tid) continue;

		float4 pos_j = __ldg(&pos[j]);

		float dx = (pos_i.x - pos_j.x) * scale;
		float dy = (pos_i.y - pos_j.y) * scale;
		float dz = (pos_i.z - pos_j.z) * scale;

		float dist2 = dx*dx + dy*dy + dz*dz;

		if (dist2 >= h_krnl2) continue;

		float diff = h_krnl2 - dist2;
		float w = mKernel1Scaled * diff * diff * diff;

		dens += __ldg(&mass[j]) * w;

		if (++count_neighb > 32) break;
	}

	density[tid] = dens;
}


// Hydro using the voxelization (re-computo distancia, no hay con qué darle por ahora...)
__global__ void hydro_voxel_CUDA(float4* pos, float4* vel, float4* accel,
                                 float* mass, float* density, int* global_index,
                                 int* d_particle_ids, int* d_sorted_cell_ids,
                                 int* d_cell_start, int step)
{
    int tid_global = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid_global >= cant_particles) return;

    int tid = d_particle_ids[tid_global];

    // Datos de la partícula actual
    float4 pos_i = __ldg(&pos[tid]);
    float4 vel_i = __ldg(&vel[tid]);
    float rhoi = __ldg(&density[tid]);
    float pi = (rhoi - mRho0) * mStiffness;
    float rhoiInv = (pi > 0.0f) ? (1.0f / pi) : 1.0f;
    float rhoiInv2 = rhoiInv * rhoiInv;
    float piDivRhoi2 = pi * rhoiInv2;

    // Acumuladores
    float pressureGradX = 0.0f, pressureGradY = 0.0f, pressureGradZ = 0.0f;
    float viscousX = 0.0f, viscousY = 0.0f, viscousZ = 0.0f;
    float pseudo_soft = 1e-3f;

    int count_neighb = 0;

    // Acceso a vecinos
    int cell_id = global_index[tid];
    int start = d_cell_start[cell_id];
    int end = d_cell_start[cell_id + 1];
    int len = end - start;
    if (len <= 1) return;

    int offset = hash(step) % len;

    for (int i = 0; i < len; ++i) {
        int j_this = start + ((i + offset) % len);
        int j = d_particle_ids[j_this];
        if (j == tid) continue;

        float4 pos_j = __ldg(&pos[j]);
        float dx = (pos_i.x - pos_j.x) * scale;
        float dy = (pos_i.y - pos_j.y) * scale;
        float dz = (pos_i.z - pos_j.z) * scale;

        float dist2 = dx*dx + dy*dy + dz*dz;
        if (dist2 >= h_krnl2) continue;

        float invDist = rsqrtf(dist2 + pseudo_soft);
        float dist = 1.0f / invDist;

        float mj = __ldg(&mass[j]);
        float rhoj = __ldg(&density[j]);
        float pj = (rhoj - mRho0) * mStiffness;
        float rhojInv = (rhoj > 0.0f) ? (1.0f / rhoj) : 1.0f;
        float rhojInv2 = rhojInv * rhojInv;

        // Gradiente de presión
        float w = (h_krnl - dist);
        float factor = w * w * mj * piDivRhoi2 * pj * rhojInv2 * mKernel2Scaled * invDist;

        pressureGradX += dx * factor;
        pressureGradY += dy * factor;
        pressureGradZ += dz * factor;

        // Viscosidad
        float visc = w * rhojInv * mj * mKernel3Scaled * mViscosityScalar * rhoiInv;
        float4 vel_j = __ldg(&vel[j]);

        viscousX += (vel_j.x - vel_i.x) * visc;
        viscousY += (vel_j.y - vel_i.y) * visc;
        viscousZ += (vel_j.z - vel_i.z) * visc;

		if (++count_neighb > 32) break;
    }

    accel[tid].x = viscousX - pressureGradX;
    accel[tid].y = viscousY - pressureGradY;
    accel[tid].z = viscousZ - pressureGradZ;
}



// Lastly: grav + integration (!)
__global__ void integrate_CUDA(float4* pos, float4* vel, float4* accel)
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
	rMinusRjScaled[0] = (pos[tid].x - x_centre) * scale;
	rMinusRjScaled[1] = (pos[tid].y - y_centre) * scale;
	rMinusRjScaled[2] = (pos[tid].z - z_centre) * scale;

	// Softening va squared!
	distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] + rMinusRjScaled[1] * rMinusRjScaled[1] +\
				rMinusRjScaled[2] * rMinusRjScaled[2] + (softening * softening);

	invDist = rsqrtf(distance_ij3);  // quick x^(-1/2)
	invDist = invDist * invDist * invDist;  // dist^-3

	// Updateo la gravedad:
	// acceleration += gravityTerm;  -> Escribo mal a proposito, just to check (es +=)
	accel[tid].x += -grav_cte * mass_centre * (rMinusRjScaled[0] * invDist);
	accel[tid].y += -grav_cte * mass_centre * (rMinusRjScaled[1] * invDist);
	accel[tid].z += -grav_cte * mass_centre * (rMinusRjScaled[2] * invDist);
	
	// OJO! Me falto check for CFL condition:
	dot = (accel[tid].x * accel[tid].x) + (accel[tid].y * accel[tid].y) +\
		(accel[tid].z * accel[tid].z);

	// Lo meto a dedo... (recall branchles...)
	accel[tid].x *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	accel[tid].y *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	accel[tid].z *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	// acc = acc if (|acc|^2 < CFL^2 => "x1"); acc = acc/|acc| * CFL if (|acc|^2 >= CFL^2)
	
	// ---------------------------------------------
	// Listo! Ahora, integro con LF-KDK -> Next kernel...	
	// ---------------------------------------------
		
	// Only gravity (reset accel para el mid-step!)
	vel[tid].x += accel[tid].x * delta_step * 0.5f;
	vel[tid].y += accel[tid].y * delta_step * 0.5f;
	vel[tid].z += accel[tid].z * delta_step * 0.5f;

	pos[tid].x += vel[tid].x * delta_step;
	pos[tid].y += vel[tid].y * delta_step;
	pos[tid].z += vel[tid].z * delta_step;

	// Nuevamente calc la accel...		  
	rMinusRjScaled[0] = (pos[tid].x - x_centre) * scale;
	rMinusRjScaled[1] = (pos[tid].y - y_centre) * scale;
	rMinusRjScaled[2] = (pos[tid].z - z_centre) * scale;

	// Idem, soft^2...
	distance_ij3 = rMinusRjScaled[0] * rMinusRjScaled[0] + rMinusRjScaled[1] * rMinusRjScaled[1] +\
				rMinusRjScaled[2] * rMinusRjScaled[2] + (softening * softening);

	invDist = rsqrtf(distance_ij3);  // quick x^(-1/2)
	invDist = invDist * invDist * invDist;  // dist^-3

	// Updateo la gravedad (de 0!):
	accel[tid].x = -grav_cte * mass_centre * (rMinusRjScaled[0] * invDist);
	accel[tid].y = -grav_cte * mass_centre * (rMinusRjScaled[1] * invDist);
	accel[tid].z = -grav_cte * mass_centre * (rMinusRjScaled[2] * invDist);

	
	// OJO! Me falto check for CFL condition:
	dot = (accel[tid].x * accel[tid].x) + (accel[tid].y * accel[tid].y) +\
		(accel[tid].z * accel[tid].z);

	// Lo meto a dedo... (recall branchles...)
	accel[tid].x *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	accel[tid].y *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));
	accel[tid].z *= (1.f * (dot < 1e+8f) + (1e+4f * rsqrtf(dot) * (dot >= 1e+8f)));


	// Integro todo! LF-KDK: Only gravity
	vel[tid].x += accel[tid].x * delta_step;
	vel[tid].y += accel[tid].y * delta_step;
	vel[tid].z += accel[tid].z * delta_step;

	// Y ya modifiqué todos los vectores!
}

//-------------------------------------- FUNCIÓN PRINCIPAL Y DESTRUCTURA -------------------------------------
// Host-side wrapper function
void launchMyKernel(DeviceData* devData, float* h_position_x, float* h_position_y, float* h_position_z, 
				    float* h_mass, float* h_density, int h_cant_particles)
{
	countStemp++;
    int threadsPerBlock = 512;
    int blocksPerGrid = (h_cant_particles + threadsPerBlock - 1) / threadsPerBlock;

    // Paso 1: voxelize
    voxelize_CUDA<<<blocksPerGrid, threadsPerBlock>>>(devData->d_pos, devData->num_cells, devData->global_index);
    CUDA_CHECK(cudaDeviceSynchronize());

    // Paso 2: ordenar y preparar estructura por celda
    //SortedParticlesData sortedData;
    prepare_sorted_particles(devData, h_cant_particles, devData->num_cells, devData->sorted_data);


    // Paso 3: vecinos (usa ids ordenados y cell_start)
    neighbors_voxel_CUDA<<<blocksPerGrid, threadsPerBlock>>>(devData->d_pos, devData->d_mass, devData->d_density, devData->global_index,devData->sorted_data.raw_particle_ids(), thrust::raw_pointer_cast(devData->sorted_data.d_sorted_cell_ids.data()), thrust::raw_pointer_cast(devData->sorted_data.d_cell_start.data()), countStemp);
	CUDA_CHECK(cudaDeviceSynchronize());


    // Paso 4: hydro
    hydro_voxel_CUDA<<<blocksPerGrid, threadsPerBlock>>>( devData->d_pos, devData->d_vel, devData->d_accel, 
												 devData->d_mass, devData->d_density, devData->global_index, devData->sorted_data.raw_particle_ids(), thrust::raw_pointer_cast(devData->sorted_data.d_sorted_cell_ids.data()), thrust::raw_pointer_cast(devData->sorted_data.d_cell_start.data()), countStemp);
    CUDA_CHECK(cudaDeviceSynchronize());

    // Paso 5: integración
    integrate_CUDA<<<blocksPerGrid, threadsPerBlock>>>(devData->d_pos, devData->d_vel, devData->d_accel);
    CUDA_CHECK(cudaDeviceSynchronize());

    // Paso 6: copiar al host

	if (countStemp % 2 == 0) {
		float4* h_pos = (float4*)malloc(sizeof(float4) * h_cant_particles);
		cudaMemcpy(h_pos, devData->d_pos, sizeof(float4) * h_cant_particles, cudaMemcpyDeviceToHost);

		// Separar en componentes si necesitás
		for (int i = 0; i < h_cant_particles; ++i) {
			h_position_x[i] = h_pos[i].x;
			h_position_y[i] = h_pos[i].y;
			h_position_z[i] = h_pos[i].z;
		}

		free(h_pos);
	}
}


void cleanupDeviceData(DeviceData* devData) 
{
    cudaFree(devData->d_pos);
    cudaFree(devData->d_vel);
    cudaFree(devData->d_accel);
    cudaFree(devData->d_mass);
    cudaFree(devData->d_density);
	cudaFree(devData->global_index);

    free(devData);
}







