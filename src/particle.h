#ifndef PARTICLE_H
#define PARTICLE_H

// #include "vec3.h"
#include <math.h>  // Lo unico que necesitaba de vec3.h ...
#include <vector>

// Quiero en vez de Pos[3N] -> X[N], Y[N], Z[N]
class Particle
{
public:

   Particle(size_t numParticles);
   ~Particle() = default;
   std::vector<float> mMass;
   std::vector<float> mDensity;

   std::vector<float> mPositionX;
   std::vector<float> mPositionY;
   std::vector<float> mPositionZ;
   std::vector<float> mVelocityX;
   std::vector<float> mVelocityY;
   std::vector<float> mVelocityZ;
   std::vector<float> mAccelerationX;
   std::vector<float> mAccelerationY;
   std::vector<float> mAccelerationZ;
   
   std::vector<int> mNeighborCount;

};

#endif // PARTICLE_H
