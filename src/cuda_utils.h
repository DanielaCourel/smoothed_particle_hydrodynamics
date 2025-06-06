#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

// Declare host-side wrapper function(s)
void launchMyKernel(float* h_position, float* h_velocity, float* h_acceleration,
				float* h_mass, float* h_density,
				int h_cant_particles, float h_scale, float h_softening, float h_grav_cte,
				float h_mass_centre, float h_delta_step,
				float h_x_centre, float h_y_centre, float h_z_centre,
				float h_h_krnl, float h_h_krnl2, float h_mHScaled9, float h_mKernel1Scaled,
				float h_mKernel2Scaled, float h_mKernel3Scaled, float h_mRho0,
				float h_mViscosityScalar, float h_mStiffness, int h_side_grid);

#endif // CUDA_UTILS_H
