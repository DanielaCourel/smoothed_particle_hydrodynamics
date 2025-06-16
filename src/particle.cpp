// base
#include "particle.h"
//#include <vector>  // The vector header is already included in particle.h

Particle::Particle(size_t numParticles)
 : mMass(numParticles, 0.0f),
   mDensity(numParticles, 0.0f),
   mPosition_x(numParticles, 0.0f),
   mPosition_y(numParticles, 0.0f),
   mPosition_z(numParticles, 0.0f),
   mVelocity_x(numParticles, 0.0f),
   mVelocity_y(numParticles, 0.0f),
   mVelocity_z(numParticles, 0.0f),
   mAcceleration_x(numParticles, 0.0f),
   mAcceleration_y(numParticles, 0.0f),
   mAcceleration_z(numParticles, 0.0f),
   mNeighborCount(numParticles, 0)
{
   
}