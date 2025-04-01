// base
#include "particle.h"

/* Es demasiado grande al pedo. No hace falta guardar las inversas
de las densidades, más vale las calculemos on-the-fly; ¿Ídem Energía/AngMom? */
Particle::Particle()
 : mMass(1.0f),
   mDensity(0.0f)//,
   //mDensityInverse(0.0f),
   //mDensityInverse2(0.0f),
   //mPressure(0.0f),
   //mPotentialEnergy(0.0f),
   //mKineticEnergy(0.0f),
   //mAngularMomentum(0.0f, 0.0f, 0.0f)
{
   mVelocity.set(0.0f, 0.0f, 0.0f);
   mAcceleration.set(0.0f, 0.0f, 0.0f);
}


