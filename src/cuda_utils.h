#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

// Declare host-side wrapper function(s)
void launchMyKernel(float* h_position, float* h_velocity, float* h_acceleration,
				int h_cant_particles, float h_scale, float h_softening, float h_grav_cte,
				float h_mass_centre, float h_delta_step,
				float h_x_centre, float h_y_centre, float h_z_centre);

#endif // CUDA_UTILS_H
