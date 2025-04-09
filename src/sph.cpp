// base
#include "sph.h"

// sph
#include "particle.h"

// Qt
#include <QElapsedTimer>

// cmath
#include <math.h>

// openmp
#include <omp.h>

#include <QDateTime>

// write
#include <iostream>
#include <fstream>
#include <sys/stat.h> 
#include <sys/types.h> // write
#include <iostream>
#include <fstream>
#include <sys/stat.h>

#ifndef M
#define M 16
#endif

// Unidades: [km/s pc M_sun Myr]...
// ¿Define a cada step estos valores?
SPH::SPH()
 : mParticleCount(0),
   mGridCellCount(0),
   mRho0(0.0f),
   mStopped(false),
   mPaused(false),
   mKineticEnergyTotal(0.0f),
   mPotentialEnergyTotal(0.0f),
   mAngularMomentumTotal(vec3(0.0f, 0.0f, 0.0f))
{
   // grid
   float h = 0.1f;  // OG = 3.34
   mSimulationScale = 1.0f;  // OG = 4e-3
   mSimulationScaleInverse = 1.0f / mSimulationScale;
   mH = h;
   mH2 = pow(h, 2);
   mHTimes2 = h * 2.0f;
   mHTimes2Inv = 1.0f / mHTimes2;
   mHScaled  = h * mSimulationScale;
   mHScaled2 = pow(h * mSimulationScale, 2);
   mHScaled6 = pow(h * mSimulationScale, 6);
   mHScaled9 = pow(h * mSimulationScale, 9);  // Creo que está bien mantener estos floats así
   // ¿Pero no conviene definirlos como constantes?
   mParticleCount = M * 1024;
   mGridCellsX = 32;  // OG 32...
   mGridCellsY = 32;
   mGridCellsZ = 32;
   mGridCellCount = mGridCellsX * mGridCellsY * mGridCellsZ;
   mCellSize = 2.0f * h;
   mMaxX = mCellSize * mGridCellsX;
   mMaxY = mCellSize * mGridCellsY;
   mMaxZ = mCellSize * mGridCellsZ;

   float time_simu = 1.0f;  // [Myr]
   mTimeStep = 0.001f;
   totalSteps = (int)round(time_simu/mTimeStep);

   // physics
   mRho0 = 0.1f;  // Check qué debería ser para el H_1 + He_2 pristino...
   mStiffness = 0.1f;  // idk
   mGravity = vec3(0.0f, 0.0f, 0.0f);
   mViscosityScalar = 1.0f;  // 1e+1~2 == nice disk formation (!!!)
   mDamping = 0.001f;  // Deberíamos "tirar" las que se escapen (En vez de checkear boundaries...)
   // Deberíamos definir acá la const de grav, el softening, la masa central y su pos?
   mGravConstant = 4.3009e-3f;  // En pc (km/s)^2 / M_sun
   mCentralMass = 1e+5f;  // As we wish
   mCentralPos = vec3(mMaxX * 0.5f, mMaxY * 0.5f, mMaxZ * 0.5f);
   mSoftening = mHScaled * 0.5;  // Ojo, se recomienda que sea = h...

   float mass = 1.0f;  // Cada star = 1 M_sun (Pero orig son "part de gas"...)
   mCflLimit = 10000.0f;  // Esto no debería existir...
   mCflLimit2 = mCflLimit * mCflLimit;

   // smoothing kernels
   mKernel1Scaled = 315.0f / (64.0f * (float)(M_PI) * mHScaled9);
   mKernel2Scaled = -45.0f / ((float)(M_PI) * mHScaled6);  // M_PI as float?
   mKernel3Scaled = -mKernel2Scaled;

   // Valor fiducial = 32
   mExamineCount = 32;

   mSrcParticles = new Particle[mParticleCount];
   mVoxelIds= new int[mParticleCount];
   mVoxelCoords= new vec3i[mParticleCount];

   // Para difs masas (estaria bueno ver que onda...)
   for (int i = 0; i < mParticleCount; i++)
   {
      mSrcParticles[i].mMass = mass;
   }

   mGrid = new QList<uint16_t>[mGridCellCount];

   mNeighbors = new uint16_t[mParticleCount*mExamineCount];
   mNeighborDistancesScaled = new float[mParticleCount*mExamineCount];

   // randomize particle start positions -> Cambiar por orta config...
   // initParticlePositionsRandom();
   initParticlePolitionsSphere();
}

SPH::~SPH()
{
   stopSimulation();
   delete[] mSrcParticles;
   delete[] mVoxelIds;
   delete[] mVoxelCoords;
   delete[] mGrid;
   delete[] mNeighbors;
   delete[] mNeighborDistancesScaled;
   quit();
   wait();
}


bool SPH::isStopped() const
{
   mMutex.lock();
   bool stopped = mStopped;
   mMutex.unlock();

   return stopped;
}


bool SPH::isPaused() const
{
   mMutex.lock();
   bool paused = mPaused;
   mMutex.unlock();

   return paused;
}



void SPH::run()
{
   int stepCount = 0;

   // Create directory ./out
   const char *path = "out";
   int result = mkdir(path, 0777);
   if (result == 0)
      std::cout << "Directory created" << std::endl;
   else
      std::cout << "Directory already exists" << std::endl;

   // Create files
   std::ofstream outfile1("out/energy.txt");
   outfile1 << "Step, Kinetic Energy, Potential Energy, Total Energy" << std::endl;
   std::ofstream outfile2("out/angularmomentum.txt");
   outfile2 << "Step, Angular Momentum" << std::endl;
   std::ofstream outfile3("out/timing.txt");
   outfile3 << "Step, Voxelize, Find Neighbors, Compute Density, Compute Pressure, Compute Acceleration, Integrate" << std::endl;


   while(!isStopped() && stepCount <= totalSteps)
   {
      if (!isPaused())
      {
         step();
         outfile1 << stepCount << ", " << mKineticEnergyTotal << ", " << mPotentialEnergyTotal << ", " << mKineticEnergyTotal + mPotentialEnergyTotal << std::endl;
         outfile2 << stepCount << ", " << mAngularMomentumTotal.length() << std::endl;
         outfile3 << stepCount << ", " << timeVoxelize << ", " << timeFindNeighbors << ", " << timeComputeDensity << ", " << timeComputePressure << ", " << timeComputeAcceleration << ", " << timeIntegrate << std::endl;
         stepCount++;
      }
   }

   outfile1.close();
   outfile2.close();
   outfile3.close();
}


void SPH::step()
{
   timeVoxelize = 0;
   timeFindNeighbors = 0;
   timeComputeDensity = 0;
   timeComputePressure = 0;
   timeComputeAcceleration = 0;
   timeIntegrate = 0;
   QElapsedTimer t;
   mKineticEnergyTotal = 0.0f;
   mPotentialEnergyTotal = 0.0f;
   mAngularMomentumTotal = vec3(0.0f, 0.0f, 0.0f);

   // put particles into voxel grid
   t.start();
   voxelizeParticles();
   timeVoxelize = t.nsecsElapsed() / 1000000;

   // find neighboring particles
   t.start();
   // #pragma omp parallel for
   for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];
      const vec3i& voxel= mVoxelCoords[particleIndex];

      // neighbors for this particle
      uint16_t* neighbors= &mNeighbors[particleIndex*mExamineCount];

      // a) examine a local region of 2x2x2 voxels using the spatial index to
      //    look up the particles that we are going to examine
      //    -> create a neighbor map based on the analaysis
      findNeighbors(particle, particleIndex, neighbors, voxel.x, voxel.y, voxel.z);
   }
   timeFindNeighbors = t.nsecsElapsed() / 1000000;

   // compute density
   //    we only compute interactions with 32 particles.
   //    -> compute the interaction and the physics (with these 32 particles)
   t.start();
   // #pragma omp parallel for

   // Maybe se puede hacer durante el anterior loop?
   for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];

      // neighbors for this particle
      uint16_t* neighbors= &mNeighbors[particleIndex*mExamineCount];
      float* neighborDistances= &mNeighborDistancesScaled[particleIndex*mExamineCount];

      computeDensity(particle, neighbors, neighborDistances);
   }
   timeComputeDensity = t.nsecsElapsed() / 1000000;

   // compute pressure
   t.start();
   // #pragma omp parallel for

   // Maybe se puede hacer durante el anterior loop?
   /* for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];

      computePressure(particle);
   } */ // Skip because now are one-liners
   timeComputePressure = t.nsecsElapsed() / 1000000;

   // compute acceleration
   t.start();
   // #pragma omp parallel for

   // Maybe se puede hacer durante el anterior loop?
   for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];

      // neighbors for this particle
      uint16_t* neighbors= &mNeighbors[particleIndex*mExamineCount];
      float* neighborDistances= &mNeighborDistancesScaled[particleIndex*mExamineCount];

      computeAcceleration(particle, neighbors, neighborDistances);
   }
   timeComputeAcceleration = t.nsecsElapsed() / 1000000;

   // integrate
   t.start();
   // #pragma omp parallel for

   // Maybe se puede hacer durante el anterior loop? Ojo con el orden...
   for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];

      integrate(particle);
      // Calc acá T y W del sistema (no guardar las energías en las Particles)

      mKineticEnergyTotal += 0.5f * particle->mMass * particle->mVelocity.length2();
      //mKineticEnergyTotal += particle->mKineticEnergy;
      
      mAngularMomentumTotal += (particle->mMass * \
         (particle->mPosition - mCentralPos).cross(particle->mVelocity));
      //mAngularMomentumTotal += particle->mAngularMomentum;
   }
   timeIntegrate = t.nsecsElapsed() / 1000000;

   emit updateElapsed(
      timeVoxelize,
      timeFindNeighbors,
      timeComputeDensity,
      timeComputePressure,
      timeComputeAcceleration,
      timeIntegrate
   );

   emit stepFinished();
}


void SPH::pauseResume()
{
   mMutex.lock();
   mPaused = !mPaused;
   mMutex.unlock();
}


void SPH::stopSimulation()
{
   mMutex.lock();
   mStopped = true;
   mMutex.unlock();
}


// Sobra... -> Cambiar por otras cond inic
void SPH::initParticlePositionsRandom()
{
   srand(QDateTime::currentMSecsSinceEpoch() % 1000);

   for (int i = 0; i < mParticleCount; i++)
   {
      float x = rand() / (float)RAND_MAX;
      float y = rand() / (float)RAND_MAX;
      float z = rand() / (float)RAND_MAX;

      x *= mGridCellsX * mHTimes2 * 0.1f;
      y *= mGridCellsY * mHTimes2 * 0.75f;
      z *= mGridCellsZ * mHTimes2;

      if (x == (float)mGridCellsX)
         x -= 0.00001f;
      if (y == (float)mGridCellsY)
         y -= 0.00001f;
      if (z == (float)mGridCellsZ)
         z -= 0.00001f;

      mSrcParticles[i].mPosition.set(x, y, z);
   }

   // just set up random directions
   for (int i = 0; i < mParticleCount; i++)
   {
      
      // have a range from -1 to 1
      float x = ((rand() / (float)RAND_MAX) * 2.0f) - 1.0f;
      float y = ((rand() / (float)RAND_MAX) * 2.0f) - 1.0f;
      float z = ((rand() / (float)RAND_MAX) * 2.0f) - 1.0f;

      mSrcParticles[i].mVelocity.set(x, y, z);
   }
}


void SPH::initParticlePolitionsSphere()
{
   // Fix seed:
   //srand(QDateTime::currentMSecsSinceEpoch() % 1000);
   srand(42);

   float dist = 0.0f;

   float x = 0.0f;
   float y = 0.0f;
   float z = 0.0f;

   vec3 sphereCenter;
   sphereCenter.set(
      mMaxX * 0.5f,
      mMaxY * 0.5f,
      mMaxZ * 0.5f
   );

   float radius = 2.0f;
   float phi;  // El ang acimutal para la v_tangencial. (atan2(y,x))
   float v_x_inic, v_y_inic, v_z_inic;  // El hdp puso a y como la comp vertical...
                              // (no quiero v_inic en "z" (que aca es "y"))

   int a = mParticleCount;

   for (int i = 0; i < a; i++)
   {
      do
      {
         x = rand() / (float)RAND_MAX;
         y = rand() / (float)RAND_MAX;
         z = rand() / (float)RAND_MAX;

         x *= mGridCellsX * mHTimes2;
         y *= mGridCellsY * mHTimes2;
         z *= mGridCellsZ * mHTimes2;

         if (x == (float)mGridCellsX)
            x -= 0.00001f;
         if (y == (float)mGridCellsY)
            y -= 0.00001f;
         if (z == (float)mGridCellsZ)
            z -= 0.00001f;

         dist = (vec3(x,y,z) - sphereCenter).length();
      }
      while (dist > radius);

      mSrcParticles[i].mPosition.set(x, y, z);
      phi = atan2(z - mMaxZ * 0.5f, x - mMaxX * 0.5f);  // Acomodar por el centro de la esfera!
      v_x_inic = 20.0f * pow(dist + mHScaled*0.5, -0.5) * -sin(phi);  // L ~ r (!) -> Check
      v_z_inic = 20.0f * pow(dist + mHScaled*0.5, -0.5) * cos(phi);
      // Some random movements on "y" (z)
      v_y_inic = ((rand() / (float)RAND_MAX) * 0.5f) - 0.25f;
      mSrcParticles[i].mVelocity.set(v_x_inic, v_y_inic, v_z_inic);
   }

}



void SPH::clearGrid()
{
   for (int i = 0; i < mGridCellCount; i++)
   {
      mGrid[i].clear();
   }
}


void SPH::voxelizeParticles()
{
   // Esto sobre si no usamos OMP, no?
   //omp_lock_t writelock;
   //omp_init_lock(&writelock);

   clearGrid();

   // #pragma omp parallel for
   for (int i = 0; i < mParticleCount; i++)
   {
      Particle* particle = &mSrcParticles[i];

      // compute a scalar voxel id from a position
      vec3 pos = particle->mPosition;

      // the voxel size is 2h * 2h * 2h
      // to get the voxel just divide the coordinate by 2h
      // then convert this x, y, z index to a serial number
      // between 0 - maximum number of voxels
      int voxelX = (int)floor(pos.x * mHTimes2Inv);
      int voxelY = (int)floor(pos.y * mHTimes2Inv);
      int voxelZ = (int)floor(pos.z * mHTimes2Inv);

      // it has been seen the positions can run slightly out of bounds for
      // one solver step. so the positions are temporarily fixed here.
      if (voxelX < 0) voxelX= 0;
      if (voxelY < 0) voxelY= 0;
      if (voxelZ < 0) voxelZ= 0;
      if (voxelX >= mGridCellsX) voxelX= mGridCellsX-1;
      if (voxelY >= mGridCellsY) voxelY= mGridCellsY-1;
      if (voxelZ >= mGridCellsZ) voxelZ= mGridCellsZ-1;

      // don't write into particle but into separate memory
      mVoxelCoords[i].x= voxelX;
      mVoxelCoords[i].y= voxelY;
      mVoxelCoords[i].z= voxelZ;

      int voxelId = computeVoxelId(voxelX, voxelY, voxelZ);

      // the lock performs terribly! (~17ms instead of 3ms)
      //      omp_set_lock(&writelock);
      //      mGrid[voxelId].push_back(particle);
      //      omp_unset_lock(&writelock);

      // instead, remember the voxelId and put the particle in later
      mVoxelIds[i]= voxelId;
   }

   // put each particle into according voxel (sequential)
   for (int i = 0; i < mParticleCount; i++)
   {
       //holaxd
       mGrid[ mVoxelIds[i] ].push_back(i);
   }
}


void SPH::findNeighbors(Particle* particle, int particleIndex, uint16_t* neighbors, int voxelX, int voxelY, int voxelZ)
{
   float xOrientation = 0.0f;
   float yOrientation = 0.0f;
   float zOrientation = 0.0f;

   int x = 0;
   int y = 0;
   int z = 0;

   int particleOffset = 0;
   int particleIterateDirection = 0;
   int neighborIndex = 0;
   bool enoughNeighborsFound = false;

   vec3 pos = particle->mPosition;

   // get voxel orientation within voxel
   //
   //
   // 3d illustration
   //
   //    1 voxel:
   //      ____________
   //     /     /     /|
   //    /     /     / |
   //   +-----+-----+  |  the interaction radius of particle x may
   //   |     |  x  |  |  reach into 7 other voxels
   //   |     |     | /|  => 8 voxels to check including 'self'
   //   +-----+-----+/ |
   //   |     |     |  |
   //   |     |     | /
   //   +-----+-----+/
   //
   //   <---< h >--->
   //
   //
   // 2d illustration of 1 slice
   //
   //    +-----+-----+-----+-----
   //    |     |/////|/////|
   //    |     |/////|/////|
   //    +-----+-----+-----+-----
   //    |     |////x|/////|
   //    |     |/////|/////|
   //    +-----+-----+-----+-----
   //    |     |     |     |
   //    |     |     |     |
   //    +-----+-----+-----+-----
   //    |     |     |     |
   //    |     |     |     |

   // this gives us the relative position; i.e the orientation within a voxel
   xOrientation = pos.x - (voxelX * mHTimes2);
   yOrientation = pos.y - (voxelY * mHTimes2);
   zOrientation = pos.z - (voxelZ * mHTimes2);

   // get neighbour voxels
   x = 0;
   y = 0;
   z = 0;

   // retrieve location within voxel
   // 1 voxel side ^= 2h
   // => check in which half of the voxel a particle is located
   //    btw, we just ignore the fact a particle can be positioned exactly on a
   //    center axis of a voxel
   (xOrientation > mH) ? x++ : x--;
   (yOrientation > mH) ? y++ : y--;
   (zOrientation > mH) ? z++ : z--;

   // neighbour voxels
   int vx[8];
   int vy[8];
   int vz[8];

   // same slice
   vx[0] = voxelX;
   vy[0] = voxelY;
   vz[0] = voxelZ;

   // distance 1
   vx[1] = voxelX + x;
   vy[1] = voxelY;
   vz[1] = voxelZ;

   vx[2] = voxelX;
   vy[2] = voxelY + y;
   vz[2] = voxelZ;

   vx[3] = voxelX;
   vy[3] = voxelY;
   vz[3] = voxelZ + z;

   // distance 2
   vx[3] = voxelX + x;
   vy[3] = voxelY + y;
   vz[3] = voxelZ;

   vx[5] = voxelX + x;
   vy[5] = voxelY;
   vz[5] = voxelZ + z;

   vx[6] = voxelX;
   vy[6] = voxelY + y;
   vz[6] = voxelZ + z;

   // distance 3
   vx[7] = voxelX + x;
   vy[7] = voxelY + y;
   vz[7] = voxelZ + z;

   int vxi;
   int vyi;
   int vzi;

   for (int voxelIndex = 0; voxelIndex < 8; voxelIndex++)
   {
      vxi = vx[voxelIndex];
      vyi = vy[voxelIndex];
      vzi = vz[voxelIndex];

      // check if voxels can be processed
      if (
            vxi > 0 && vxi < mGridCellsX
         && vyi > 0 && vyi < mGridCellsY
         && vzi > 0 && vzi < mGridCellsZ
      )
      {
         // +-----+-----+-----+-----+-----+-----+-----+-----
         // |     |     |     |     |     |     |     |
         // |     |     |    0|   *1|     |     |     |
         // +-----+-----+-----+-----+-----+-----+-----+-----
         // |     |     | *  x|  * *|     |     |     |
         // |     |     |  * 2| * *4|     |     |     |
         // +-----+-----+-----+-----+-----+-----+-----+-----
         // |     |     |     |     |     |     |     |
         // |     |     |     |     |     |     |     |
         // +-----+-----+-----+-----+-----+-----+-----+-----

         const QList<uint16_t>& voxel = mGrid[computeVoxelId(vxi, vyi, vzi)];

         if (!voxel.isEmpty())
         {
            // randomize idea:
            //
            // for every voxel n
            //
            //    pick a random start particle within voxel n
            //    randomly go up or down from this start particle
            //    check if that particle is within the interaction radius
            //    increment 'found counter'
            //
            //    if (32 particles found),
            //       break
            //
            //    optimization: store dot product in order to avoid
            //                  computing it twice

            // TODO: that's neither fast no a good idea
            //       if there's only 1 particle nearby, the code below is pretty pointless

            // B: OJO CON ESTE RAND, QUE SE LLAMA MUCHÍSIMO
            particleOffset = rand() % voxel.length();
            particleIterateDirection = (particleIndex % 2) ? -1 : 1;

            int i = 0;
            while (true)
            {
               int nextIndex = particleOffset + i * particleIterateDirection;

               // leave if we're out out the voxel's bounds
               if (nextIndex < 0 || nextIndex > voxel.length() - 1)
                  break;

               uint16_t realIndex = voxel[nextIndex];
               Particle* neighbor = &mSrcParticles[realIndex];
               i++;

               Particle* validNeighbor = evaluateNeighbor(particle, neighbor);

               if (validNeighbor)
               {
                  neighbors[neighborIndex] = realIndex;
                  neighborIndex++;
               }

               // leave if we have sufficient neighbor particles
               enoughNeighborsFound = (neighborIndex > mExamineCount - 1);
               if (enoughNeighborsFound)
                  break;
            }
         }
      }

      // no need to process any other voxels
      if (enoughNeighborsFound)
         break;
   }

   particle->mNeighborCount = neighborIndex;
}


// ??
Particle* SPH::evaluateNeighbor(
   Particle* current,
   Particle* neighbor
)
{
   Particle* validNeighbor = 0;

   if (current != neighbor)
   {
      // we save the sqrt here and use mInteractionRadius2 (h^2) instead
      // of mInteractionRadius
      //
      // float distance = sqrt(dot);
      // if (distance < mInteractionRadius)

      vec3 dist = current->mPosition - neighbor->mPosition;
      float dot = dist * dist;

      // the dot product is unscaled and so is MH2;
      // so there's no need to add any simulation scale here
      if (dot < mH2)
      {
         validNeighbor = neighbor;
      }
   }

   return validNeighbor;
}



void SPH::computeDensity(Particle* particle, uint16_t* neighbors, float* neighborDistances)
{
   float density = 0.0f;
   float mass = 0.0f;
   vec3 pos = particle->mPosition;
   float w = 0.0f;
   float rightPart = 0.0f;

   for (int neighborIndex = 0; neighborIndex < particle->mNeighborCount; neighborIndex++)
   {
      uint16_t realIndex = neighbors[neighborIndex];
      Particle* neighbor = &mSrcParticles[realIndex];

      if (!neighbor)
         break;

      if (neighbor != particle)
      {
         // add mass of neighbor
         mass = neighbor->mMass;

         // apply smoothing kernel to mass
         
         vec3 dist = pos - neighbor->mPosition;
         float dot = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
         float distance = sqrt(dot);
         float distanceScaled = distance * mSimulationScale;
         neighborDistances[neighborIndex] = distanceScaled;

         if (distanceScaled > mHScaled)
         {
            w = 0.0f;
         }
         else
         {
            // the dot product is not used here since it's not properly scaled
            rightPart = (mHScaled2 - (distanceScaled * distanceScaled));
            rightPart = (rightPart * rightPart * rightPart);

            //           315
            // w =  -------------  * rightPart
            //      64 * PI * h^9
            w = mKernel1Scaled * rightPart;

            // apply weighted neighbor mass to our density
            density += (mass * w);
         }
      }
   }

   particle->mDensity = density;

   // we'll need those for the acceleration evaluation later
   // that's why we precompute them now

   // B: Quizás es mejor calcularlo en donde lo necesitemos (achicar el tamaño del struct)

   //float inverseDensity = ((density > 0.0f) ? (1.0f / density) : 1.0f);
   //particle->mDensityInverse = inverseDensity;
   //particle->mDensityInverse2 = inverseDensity * inverseDensity;
}

// Skip this... Calc the pressure on-the-fly, maybe change later (diff EoS)
void SPH::computePressure(Particle* particle)
{
   // rho0: resting density
   // float deltaRho = particle->mDensity - mRho0;
   // float p = mStiffness * deltaRho;
   // particle->mPressure = p;
}


void SPH::computeAcceleration(Particle* p, uint16_t* neighbors, float* neighborDistances)
{

   // viscosity term
   //
   //     u = a scalar coefficient that says how viscous is the fluid
   //         small value => water
   //         large value => sirup
   //
   //     the viscosity term is approximately
   //
   //         u
   //        ----  *  the sum of the masses at various points j
   //        rho      multiplied by the velocity between two points i and j
   //           i     divided by the pressure of j                     __
   //                 multiplied by the laplacian of a smoothin kernel \/^2W

   Particle* neighbor = 0;
   float distanceToNeighborScaled = 0.0f;

   // OJO, hay muchas cosas que no guardamos dentro de partículas:
   //float pi = p->mPressure;
   float pi = (p->mDensity - mRho0) * mStiffness;  // One-liner...
   //inverseDensity = ((density > 0.0f) ? (1.0f / density) : 1.0f);
   //particle->mDensityInverse = inverseDensity;
   //particle->mDensityInverse2 = inverseDensity * inverseDensity;
   //float rhoiInv = 1.0f / p->mDensity;  // One-liner (if rho > 0, which should be the case...; add a "softening"?)
   float rhoiInv = ((pi > 0.0f) ? (1.0f / pi) : 1.0f);  // Better?
   float rhoiInv2 = rhoiInv * rhoiInv;
   float piDivRhoi2 = pi * rhoiInv2;
   vec3 r = p->mPosition;
   vec3 vi = p->mVelocity;

   float pj = 0.0f;
   float rhoj = 0.0f;
   float rhojInv = 0.0f;
   float rhojInv2 = 0.0f;
   float mj = 0.0f;
   vec3 rj;
   vec3 vj;
   vec3 rMinusRj;
   vec3 rMinusRjScaled;

   // pressure gradient...
   vec3 pressureGradient(0.0f, 0.0f, 0.0f);
   vec3 pressureGradientContribution;

   // ...and viscous term
   vec3 viscousTerm(0.0f, 0.0f, 0.0f);
   vec3 viscousTermContribution;

   // are added to the final acceleration
   vec3 acceleration(0.0f, 0.0f, 0.0f);

   // Acaso llama a cada rato al "neighbor count" o se entiende que es un numero fijo?
   for (int neighborIndex = 0; neighborIndex < p->mNeighborCount; neighborIndex++)
   {
      uint16_t realIndex = neighbors[neighborIndex];
      neighbor = &mSrcParticles[realIndex];

      pj = (neighbor->mDensity - mRho0) * mStiffness;  // One-liner...
      rhoj = neighbor->mDensity;
      //rhojInv = 1.0f / rhoj;  // One-liner ¿Why not this way @ 1st time?
      rhojInv = ((rhoj > 0.0f) ? (1.0f / rhoj) : 1.0f);
      rhojInv2 = rhojInv * rhojInv;
      rj = neighbor->mPosition;
      vj = neighbor->mVelocity;
      mj = neighbor->mMass;

      // pressure gradient
      rMinusRj = (r - rj);
      rMinusRjScaled = rMinusRj * mSimulationScale;
      distanceToNeighborScaled = neighborDistances[neighborIndex];
      pressureGradientContribution = rMinusRjScaled / distanceToNeighborScaled;

      //   -45
      // --------
      // PI * h^6
      pressureGradientContribution *= mKernel2Scaled;

      // (h - ||r - r  ||)^2
      //             j
      float centerPart = (mHScaled - distanceToNeighborScaled);
      centerPart = centerPart * centerPart;
      pressureGradientContribution *= centerPart;

      float factor = mj * piDivRhoi2 * (pj * rhojInv2);
      pressureGradientContribution *= factor;

      // add pressure gradient contribution to pressure gradient
      pressureGradient += pressureGradientContribution;


      // viscosity
      viscousTermContribution = vj - vi;
      viscousTermContribution *= rhojInv;
      viscousTermContribution *= mj;

      // (4)
      //      45
      //   -------- * (h - ||r - r  ||)
      //   PI * h^6               b
      viscousTermContribution *= mKernel3Scaled;
      viscousTermContribution *= (mHScaled - distanceToNeighborScaled);

      // add contribution to viscous term
      viscousTerm += viscousTermContribution;
   }

   viscousTerm *= (mViscosityScalar * rhoiInv);

   // potential optimization:
   // multiplication with rhoiInv could be done here in just one go
   acceleration = viscousTerm - pressureGradient;

   // New vecs (grav):
   vec3 gravityTerm(0.0f, 0.0f, 0.0f);
   vec3 gravityTermContribution;
   float distance_ij3;  // Needed...
   // LASTLY: add a point-mass accel @ the center (e.g. central Black-Hole/NSC)
   // *once per particle. acc = -G M /r^3 (r^, porque apunta al centro)
   rMinusRj = (r - mCentralPos);
   rMinusRjScaled = rMinusRj * mSimulationScale;
   distance_ij3 = pow((rMinusRjScaled.length() + mSoftening), 3);
   gravityTermContribution = rMinusRjScaled/distance_ij3;
   gravityTermContribution *= -mGravConstant * mCentralMass;
   gravityTerm += gravityTermContribution;  // Try it...

   // Updateo la gravedad:
   acceleration += gravityTerm;

   // check CFL condition
   float dot =
        acceleration.x * acceleration.x
      + acceleration.y * acceleration.y
      + acceleration.z * acceleration.z;

   bool limitExceeded = (dot > mCflLimit2);
   if (limitExceeded)
   {
      float length = sqrt(dot);
      float cflScale = mCflLimit / length;
      acceleration *= cflScale;
   }

   // Energía potencial sería G * Mcentral * m_i/r_i
   mPotentialEnergyTotal += -mGravConstant * mCentralMass * p->mMass / (rMinusRjScaled.length());
   // + softening);  // B: Without soft (i.e. without a Plummer equivalent)

   p->mAcceleration = acceleration;
}


void SPH::integrate(Particle* p)
{
   vec3 acceleration = p->mAcceleration;
   vec3 position = p->mPosition;
   vec3 velocity = p->mVelocity;

   // apply external forces
   //acceleration += mGravity;  // Esto ya lo hacemos en el otro loop ¿Debería hacerlo acá por LF-KDK?

   float posTimeStep = mTimeStep * mSimulationScaleInverse;
   // LF-KDK:
   // Claro, pero necesito calcular la aceleración de nuevo... ¿Only gravity? 
   vec3 velocity_halfstep;
   velocity_halfstep = velocity + (acceleration * mTimeStep*0.5);
   vec3 newPosition = position + (velocity_halfstep * posTimeStep);
   // new half-accel grav (half-step)
   vec3 newAcceleration;
   // Chanchada, re def estas cosas acá...
   // Lo mejor sería handle de otra forma el loop accel + integration en el "step()".
   float distance_ij3;
   vec3 rMinusRj  = (newPosition - mCentralPos);
   vec3 rMinusRjScaled = rMinusRj * mSimulationScale;
   distance_ij3 = pow((rMinusRjScaled.length() + mSoftening), 3);
   newAcceleration = rMinusRjScaled/distance_ij3;
   newAcceleration *= -mGravConstant * mCentralMass;

   // NO repito CFL condition...

   vec3 newVelocity = velocity_halfstep + (newAcceleration * mTimeStep);

   // let particles bounce back if they collide with the grid boundaries
   // handleBoundaryConditions(
   //    position,
   //    &newVelocity,
   //    posTimeStep,
   //    &newPosition
   // );  // B: Actually, no... Let them escape...

   p->mVelocity = newVelocity;
   p->mPosition = newPosition;  // Order matters?

   // OJO: queremos achicar el struct. Calc la energy when needed...

   // update kinetic energy
   //p->mKineticEnergy = 0.5f * p->mMass * (newVelocity.x * newVelocity.x + newVelocity.y * newVelocity.y + newVelocity.z * newVelocity.z);

   // update angular momentum m * (r x v)
   // B: OJO, tiene que ser respecto a la masa central, no al origen...
   //p->mAngularMomentum = p->mMass * (rMinusRj.cross(newVelocity));

   // B: It's negative because we are using "y" as the vertical. Como las v_inic (rot) las def según
   // (x, z) => ^x X ^z = -^y...
}


void SPH::handleBoundaryConditions(
   vec3 position,
   vec3* newVelocity,
   float timeStep,
   vec3* newPosition
)
{
   // x coord
   if (newPosition->x < 0.0f)
   {
      vec3 normal(1, 0, 0);
      float intersectionDistance = -position.x / newVelocity->x;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }
   else if (newPosition->x > mMaxX)
   {
      vec3 normal(-1, 0, 0);
      float intersectionDistance = (mMaxX - position.x) / newVelocity->x;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }

   // y coord
   if (newPosition->y < 0.0f)
   {
      vec3 normal(0, 1, 0);
      float intersectionDistance = -position.y / newVelocity->y;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }
   else if (newPosition->y > mMaxY)
   {
      vec3 normal(0, -1, 0);
      float intersectionDistance = (mMaxY - position.y) / newVelocity->y;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }

   // z coord
   if (newPosition->z < 0.0f)
   {
      vec3 normal(0, 0, 1);
      float intersectionDistance = -position.z / newVelocity->z;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }
   else if (newPosition->z > mMaxZ)
   {
      vec3 normal(0, 0, -1);
      float intersectionDistance = (mMaxZ - position.z) / newVelocity->z;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }
}


void SPH::applyBoundary(
      vec3 position,
      float timeStep,
      vec3* newPosition,
      float intersectionDistance,
   vec3 normal,
   vec3* newVelocity
)
{
   
   vec3 intersection = position + (*newVelocity * intersectionDistance);

   float dotProduct =
        newVelocity->x * normal.x
      + newVelocity->y * normal.y
      + newVelocity->z * normal.z;

   vec3 reflection = *newVelocity - (normal * dotProduct * 2.0f);

   float remaining = timeStep - intersectionDistance;

   // apply boundaries
   *newVelocity = reflection;
   *newPosition = intersection + reflection * (remaining * mDamping);
}


int SPH::computeVoxelId(int voxelX, int voxelY, int voxelZ)
{
   return (voxelZ * mGridCellsY + voxelY) * mGridCellsX + voxelX;
}


void SPH::clearNeighbors()
{
   memClear32(mNeighbors, mParticleCount * mExamineCount * sizeof(Particle*));
}


void SPH::memClear32(void* dst, int len)
{
   unsigned int* dst32= (unsigned int*)dst;
   len>>=2;
   while (len--)
      *dst32++= 0;
}


float SPH::getCellSize() const
{
   return mCellSize;
}


Particle* SPH::getParticles()
{
   return mSrcParticles;
}


int SPH::getParticleCount() const
{
   return mParticleCount;
}


void SPH::getGridCellCounts(int &x, int &y, int &z)
{
   x = mGridCellsX;
   y = mGridCellsY;
   z = mGridCellsZ;
}


void SPH::getParticleBounds(float &x, float &y, float &z)
{
   x = mMaxX;
   y = mMaxY;
   z = mMaxZ;
}


float SPH::getInteractionRadius2() const
{
   return mHScaled2;
}


QList<uint16_t>* SPH::getGrid()
{
   return mGrid;
}



vec3 SPH::getGravity() const
{
   return mGravity;
}


void SPH::setGravity(const vec3 &gravity)
{
   mGravity = gravity;
}


float SPH::getCflLimit() const
{
   return mCflLimit;
}


void SPH::setCflLimit(float cflLimit)
{
   mCflLimit = cflLimit;
   mCflLimit2 = mCflLimit * mCflLimit;
}


float SPH::getDamping() const
{
   return mDamping;
}


void SPH::setDamping(float damping)
{
   mDamping = damping;
}


float SPH::getTimeStep() const
{
   return mTimeStep;
}


void SPH::setTimeStep(float timeStep)
{
   mTimeStep = timeStep;
}


float SPH::getViscosityScalar() const
{
   return mViscosityScalar;
}


void SPH::setViscosityScalar(float viscosityScalar)
{
   mViscosityScalar = viscosityScalar;
}


float SPH::getStiffness() const
{
   return mStiffness;
}


void SPH::setStiffness(float stiffness)
{
   mStiffness = stiffness;
}
