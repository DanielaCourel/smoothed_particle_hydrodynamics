// base
#include "particle.h"
//#include <vector>  // The vector header is already included in particle.h ...?

Particle::Particle(size_t numParticles)
 : mMass(numParticles, 0.0f),
   mDensity(numParticles, 0.0f),

   mPositionX(numParticles, 0.0f),
   mPositionY(numParticles, 0.0f),
   mPositionZ(numParticles, 0.0f),
   mVelocityX(numParticles, 0.0f),
   mVelocityY(numParticles, 0.0f),
   mVelocityZ(numParticles, 0.0f),
   mAccelerationX(numParticles, 0.0f),
   mAccelerationY(numParticles, 0.0f),
   mAccelerationZ(numParticles, 0.0f),

   mNeighborCount(numParticles, 0)
   //Maybe internalEnergy?
   // Accessing the position of the i-th particle:
   // float x = mPosition[i * 3];
   // float y = mPosition[i * 3 + 1];
   // float z = mPosition[i * 3 + 2];
{
   
}