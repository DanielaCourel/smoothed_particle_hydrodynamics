// base
#include "particle.h"

/* Es demasiado grande al pedo. No hace falta guardar las inversas
de las densidades, más vale las calculemos on-the-fly; ¿Ídem Energía/AngMom? */
Particle::Particle(size_t numParticles)
 : mMass(numParticles, 0.0f),
   mDensity(numParticles, 0.0f),
   mPosition(numParticles, vec3(0.0f, 0.0f, 0.0f)),
   mVelocity(numParticles, vec3(0.0f, 0.0f, 0.0f)),
   mAcceleration(numParticles, vec3(0.0f, 0.0f, 0.0f)),
   mNeighborCount(numParticles, 0)
   //,
   //mDensityInverse(0.0f),
   //mDensityInverse2(0.0f),
   //mPressure(0.0f),
   //mPotentialEnergy(0.0f),
   //mKineticEnergy(0.0f),
   //mAngularMomentum(0.0f, 0.0f, 0.0f)
{
   
}