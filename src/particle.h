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
   std::vector<vec3> mPosition;
   std::vector<vec3> mVelocity;
   std::vector<vec3> mAcceleration;
   std::vector<int> mNeighborCount;

};

#endif // PARTICLE_H
