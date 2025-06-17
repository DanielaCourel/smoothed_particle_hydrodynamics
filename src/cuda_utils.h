#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/sequence.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>

// Struct para almacenar arreglos uxiliares de vecinos
struct SortedParticlesData{
    thrust::device_vector<int> d_particle_ids;
    thrust::device_vector<int> d_sorted_cell_ids;
    thrust::device_vector<int> d_cell_start;

    int* raw_particle_ids() { return thrust::raw_pointer_cast(d_particle_ids.data()); }
    int* raw_cell_start()   { return thrust::raw_pointer_cast(d_cell_start.data()); }
};


// Struct para encapsular datos del dispositivo
typedef struct {
    __align__(16) float4* d_pos;
    __align__(16) float4* d_vel;
    __align__(16) float4* d_accel;
    float* d_mass;
    float* d_density;
    int* global_index;
    int* d_particle_ids;
    int num_cells;

    SortedParticlesData sorted_data;

} DeviceData;

// Declare host-side wrapper function(s)
void launchMyKernel(DeviceData* devData, float* h_position_x, float* h_position_y, float* h_position_z, 
				    float* h_mass, float* h_density, int h_cant_particles);

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
                           float h_mViscosityScalar, float h_mStiffness, int h_side_grid);

void cleanupDeviceData(DeviceData* devData);

#endif // CUDA_UTILS_H
