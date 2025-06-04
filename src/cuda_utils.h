#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <GL/glew.h> // Required for GLuint
#include <cuda_runtime.h> // Required for cudaGraphicsResource

// Forward declaration of the CUDA graphics resource
struct cudaGraphicsResource;

#ifdef __cplusplus
extern "C" {
#endif

/*
Initializes CUDA-OpenGL interoperability.
Creates an OpenGL Vertex Buffer Object (VBO) and registers it with CUDA.
numParticles The total number of particles in the simulation.
cuda_vbo_resource_out Pointer to a cudaGraphicsResource* where the registered CUDA resource will be stored.
GLuint The ID of the created OpenGL VBO. Returns 0 on error.
*/

GLuint initCudaGLInterop(int numParticles, struct cudaGraphicsResource** cuda_vbo_resource_out);

/*
Updates particle positions on the GPU using a CUDA kernel.
Maps the OpenGL VBO, launches the kernel to write into it, and then unmaps it.

cuda_vbo_resource The registered CUDA graphics resource for the VBO.
numParticles The total number of particles.
deltaTime The time step for the simulation.
*/

void updateParticlePositionsOnGPU(struct cudaGraphicsResource* cuda_vbo_resource, int numParticles, float deltaTime);
// CAMBIAR UNA POR LA OTRA, OJO CON ESA STRUCT

// Declare host-side wrapper function(s)
void launchMyKernel(float* h_position, float* h_velocity, float* h_acceleration,
				float* h_mass, float* h_density,
				int h_cant_particles, float h_scale, float h_softening, float h_grav_cte,
				float h_mass_centre, float h_delta_step,
				float h_x_centre, float h_y_centre, float h_z_centre,
				float h_h_krnl, float h_h_krnl2, float h_mHScaled9, float h_mKernel1Scaled,
				float h_mKernel2Scaled, float h_mKernel3Scaled, float h_mRho0,
				float h_mViscosityScalar, float h_mStiffness);

/*
Cleans up CUDA-OpenGL interoperability resources.
Unregisters the CUDA resource and deletes the OpenGL VBO.

vboID The ID of the OpenGL VBO to delete.
cuda_vbo_resource The registered CUDA graphics resource to unregister.
*/

void cleanupCudaGLInterop(GLuint vboID, struct cudaGraphicsResource* cuda_vbo_resource);

#ifdef __cplusplus
}
#endif

#endif // CUDA_UTILS_H
