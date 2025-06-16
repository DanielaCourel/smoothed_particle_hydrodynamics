#ifndef PARTICLE_H
#define PARTICLE_H

#include "vec3.h"
#include <vector>

class Particle
{
public:

   Particle(size_t numParticles);
   ~Particle() = default;
   std::vector<float> mMass;
   std::vector<float> mDensity;
   std::vector<float> mPosition_x;
   std::vector<float> mPosition_y;
   std::vector<float> mPosition_z;
   std::vector<float> mVelocity_x;
   std::vector<float> mVelocity_y;
   std::vector<float> mVelocity_z;
   std::vector<float> mAcceleration_x;
   std::vector<float> mAcceleration_y;
   std::vector<float> mAcceleration_z;
   std::vector<int> mNeighborCount;

};

#endif // PARTICLE_H
