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

#include <immintrin.h>

// Include CUDA_UTILS...
#include "cuda_utils.h" // Include the header for CUDA wrapper functions

#ifndef M
#define M 16
#endif
#define K 8

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
   mH2 = h * h;
   mHTimes2 = h * 2.0f;
   mHTimes2Inv = 1.0f / mHTimes2;
   mHScaled  = h * mSimulationScale;
   mHScaled2 = pow(h * mSimulationScale, 2);
   mHScaled6 = pow(h * mSimulationScale, 6);
   mHScaled9 = pow(h * mSimulationScale, 9);  // Creo que está bien mantener estos floats así
   // ¿Pero no conviene definirlos como constantes?
   mParticleCount = M * 1024;
   mGridCellsX = 24;  // OG 32...; Try 24 or whatever
   mGridCellsY = 24;
   mGridCellsZ = 24;
   mGridCellCount = mGridCellsX * mGridCellsY * mGridCellsZ;
   mCellSize = 2.0f * h;
   mMaxX = mCellSize * mGridCellsX;
   mMaxY = mCellSize * mGridCellsY;
   mMaxZ = mCellSize * mGridCellsZ;

   float time_simu = 1.0f;  // [Myr]
   mTimeStep = 1e-4f;
   totalSteps = (int)round(time_simu/mTimeStep);

   // physics
   mRho0 = 1e+2f;  // Check qué debería ser para el H_1 + He_2 pristino...
   mStiffness = 1e-1f;  // idk
   mGravity = vec3(0.0f, 0.0f, 0.0f);
   mViscosityScalar = 1e+1f;  // 1e+1~2 == nice disk formation (!!!)
   mDamping = 0.001f;  // Deberíamos "tirar" las que se escapen (En vez de checkear boundaries...)
   // Deberíamos definir acá la const de grav, el softening, la masa central y su pos?
   mGravConstant = 4.3009e-3f;  // En pc (km/s)^2 / M_sun
   mCentralMass = 1e+5f;  // As we wish
   //mCentralPos = vec3(mMaxX * 0.5f, mMaxY * 0.5f, mMaxZ * 0.5f);
   mCentralPos[0] = mMaxX * 0.5f;
   mCentralPos[1] = mMaxY * 0.5f;
   mCentralPos[2] = mMaxZ * 0.5f;
   mSoftening = mHScaled; // * 0.5;  // Ojo, se recomienda que sea = h...

   // Deberia depender de la cant de parts...
   float mass = 0.25f;  // Cada star = 1 M_sun (Pero orig son "part de gas"...)
   mCflLimit = 1e+4f;  // Esto no debería existir...
   mCflLimit2 = mCflLimit * mCflLimit;

   // smoothing kernels
   mKernel1Scaled = 315.0f / (64.0f * (float)(M_PI) * mHScaled9);
   mKernel2Scaled = -45.0f / ((float)(M_PI) * mHScaled6);  // M_PI as float?
   mKernel3Scaled = -mKernel2Scaled;

   // Valor fiducial = 32
   mExamineCount = 32;

   mSrcParticles = new Particle(mParticleCount);
   mVoxelIds= new int[mParticleCount];
   mVoxelCoords= new vec3i[mParticleCount];

   // Para difs masas (estaria bueno ver que onda...)
   for (int i = 0; i < mParticleCount; i++)
   {
      mSrcParticles->mMass[i] = mass;
   }

   mGrid = new QList<uint32_t>[mGridCellCount];

   mNeighbors = new uint32_t[mParticleCount*mExamineCount];
   //mNeighborDistancesScaled = new float[mParticleCount*mExamineCount];  // De más...
   // Memory-bounded => Increase complexity (recalc dist to neighbors)

   // randomize particle start positions -> Cambiar por otra config...
   // initParticlePositionsRandom();
   initParticlePolitionsSphere();

}

SPH::~SPH()
{
   stopSimulation();
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
   
   /*
   
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
   std::ofstream outfile4("out/neighbors.txt");
   
   */

   DeviceData* devData = initDeviceData(mSrcParticles->mPosition_x.data(), mSrcParticles->mPosition_y.data(), mSrcParticles->mPosition_z.data(), 
                                        mSrcParticles->mVelocity_x.data(), mSrcParticles->mVelocity_y.data(), mSrcParticles->mVelocity_z.data(),
                                        mSrcParticles->mAcceleration_x.data(), mSrcParticles->mAcceleration_y.data(), mSrcParticles->mAcceleration_z.data(),
                                        mSrcParticles->mMass.data(), mSrcParticles->mDensity.data(),
                                        mParticleCount, mSimulationScale, mSoftening, mGravConstant,
                                        mCentralMass, mTimeStep, mCentralPos[0], mCentralPos[1], mCentralPos[2],
                                        mH, mH2, mHScaled9, mKernel1Scaled, mKernel2Scaled, mKernel3Scaled,
                                        mRho0, mViscosityScalar, mStiffness, mGridCellsX);

   while(!isStopped() && stepCount <= totalSteps)
   {
      if (!isPaused())
      {
         launchMyKernel(devData, mSrcParticles->mPosition_x.data(), mSrcParticles->mPosition_y.data(), mSrcParticles->mPosition_z.data(), 
                        mSrcParticles->mMass.data(), mSrcParticles->mDensity.data(), mParticleCount);
         
         /*
         outfile1 << stepCount << ", " << mKineticEnergyTotal << ", " << mPotentialEnergyTotal << ", " << mKineticEnergyTotal + mPotentialEnergyTotal << std::endl;
         outfile2 << stepCount << ", " << mAngularMomentumTotal.length() << std::endl;
         outfile3 << stepCount << ", " << timeVoxelize << ", " << timeFindNeighbors << ", " << timeComputeDensity << ", " << timeComputePressure << ", " << timeComputeAcceleration << ", " << timeIntegrate << std::endl;
         */
         stepCount++;
         //if (stepCount > 1) break;
      }

   }
   cleanupDeviceData(devData);

   /*
   outfile1.close();
   outfile2.close();
   outfile3.close();
   outfile4.close();
   */
   
}

// El main loop va a la GPU!
void SPH::step()
{
   timeVoxelize = 0;
   timeFindNeighbors = 0;
   timeComputeDensity = 0;
   timeComputePressure = 0;
   timeComputeAcceleration = 0;
   timeIntegrate = 0;
   QElapsedTimer t;

   /* Ahora no!

   // put particles into voxel grid
   t.start();
   voxelizeParticles();
   timeVoxelize = t.nsecsElapsed() / 1000000;

   // time all the //-zone as a whole:
   t.start();

   #pragma omp parallel 
   {
      // find neighboring particles
      #pragma omp for schedule(guided)
      for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
      {
         const vec3i& voxel= mVoxelCoords[particleIndex];

         // neighbors for this particle
         uint32_t* neighbors= &mNeighbors[particleIndex*mExamineCount];
         // Quiero calc 2 veces la dist a neighbors (higher complexity (!))
         // float* neighborDistances= &mNeighborDistancesScaled[particleIndex*mExamineCount];

         //findNeighbors(particleIndex, neighbors, voxel.x, voxel.y, voxel.z, neighborDistances);
         findNeighbors(particleIndex, neighbors, voxel.x, voxel.y, voxel.z);

         //computeDensity(particleIndex, neighbors, neighborDistances);
         computeDensity(particleIndex, neighbors);
      }

      // compute acceleration
      #pragma omp for schedule(guided)
      for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
      {
         // neighbors for this particle
         uint32_t* neighbors= &mNeighbors[particleIndex*mExamineCount];
         //float* neighborDistances= &mNeighborDistancesScaled[particleIndex*mExamineCount];

         //computeAcceleration(particleIndex, neighbors, neighborDistances);
         computeAcceleration(particleIndex, neighbors);
      }

   } // End parallel -> gravity + integrate va a la GPU!
   
   timeFindNeighbors = t.nsecsElapsed() / 1000000;

   */
   
	// time all the GPU-zone as a whole:
	t.start();
	
   // Kernel! (se encarga de todo...)
   // Son "std::vector<float>", asi que así pedimos los punteros a la 1ra direc de memoria
   // (También podría hacer "&vector[0]"...)
/*    launchMyKernel(mSrcParticles->mPosition.data(), mSrcParticles->mVelocity.data(), mSrcParticles->mAcceleration.data(),
         mSrcParticles->mMass.data(), mSrcParticles->mDensity.data()); */
   // Ojo con la escala del grideado... (cada celda mide 2*h (!))

				
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


// Sobra... -> CAMBIAR POR CUALQUIER OTRA IC.
void SPH::initParticlePositionsRandom()
{
   // srand(QDateTime::currentMSecsSinceEpoch() % 1000);

   // for (int i = 0; i < mParticleCount; i++)
   // {
   //    float x = rand() / (float)RAND_MAX;
   //    float y = rand() / (float)RAND_MAX;
   //    float z = rand() / (float)RAND_MAX;

   //    x *= mGridCellsX * mHTimes2 * 0.1f;
   //    y *= mGridCellsY * mHTimes2 * 0.75f;
   //    z *= mGridCellsZ * mHTimes2;

   //    if (x == (float)mGridCellsX)
   //       x -= 0.00001f;
   //    if (y == (float)mGridCellsY)
   //       y -= 0.00001f;
   //    if (z == (float)mGridCellsZ)
   //       z -= 0.00001f;
   //    mSrcParticles->mPosition[i].set(x, y, z);
   // }

   // // just set up random directions
   // for (int i = 0; i < mParticleCount; i++)
   // {
      
   //    // have a range from -1 to 1
   //    float x = ((rand() / (float)RAND_MAX) * 2.0f) - 1.0f;
   //    float y = ((rand() / (float)RAND_MAX) * 2.0f) - 1.0f;
   //    float z = ((rand() / (float)RAND_MAX) * 2.0f) - 1.0f;

   //    mSrcParticles->mVelocity[i].set(x, y, z);
   // }
}


// ~disc
void SPH::initParticlePolitionsSphere()
{
   // Fix seed:
   //srand(QDateTime::currentMSecsSinceEpoch() % 1000);
   srand(42);

   float dist = 0.0f;

   float x = 0.0f;
   float y = 0.0f;
   float z = 0.0f;

   float sphereCenter_x = mMaxX * 0.5f;
   float sphereCenter_y = mMaxY * 0.5f;
   float sphereCenter_z = mMaxZ * 0.5f;

   float radius = 2.0f;
   float phi;  // El ang acimutal para la v_tangencial. (atan2(y,x))
   float v_x_inic, v_y_inic, v_z_inic;  // El hdp puso a y como la comp vertical...
                              // (no quiero v_inic en "z" (que aca es "y"))

   for (int i = 0; i < mParticleCount; i++)
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

         //dist = (vec3(x,y,z) - sphereCenter).length();
         dist = (x - sphereCenter_x) * (x - sphereCenter_x) +\
                (y - sphereCenter_y) * (y - sphereCenter_y) +\
                (z - sphereCenter_z) * (z - sphereCenter_z);
         dist = sqrt(dist);
      }
      while (dist > radius);

      mSrcParticles->mPosition_x[i] = x;
      mSrcParticles->mPosition_y[i] = y;
      mSrcParticles->mPosition_z[i] = z;

      phi = atan2(z - mMaxZ * 0.5f, x - mMaxX * 0.5f);  // Acomodar por el centro de la esfera!
      v_x_inic = 20.0f * pow(dist + mHScaled*0.5, -0.5) * -sin(phi);  // a = 20.0
      v_z_inic = 20.0f * pow(dist + mHScaled*0.5, -0.5) * cos(phi);  // a = 20.0
      // Some random movements on "y" (z)
      // Change it a little:
      if (phi > 0 && phi < 3.14) {
         v_y_inic = ((rand() / (float)RAND_MAX) * 6.5f);
      }
      else {
         v_y_inic = -((rand() / (float)RAND_MAX) * 6.5f);
      };
      

      //mSrcParticles->mVelocity[i].set(v_x_inic, v_y_inic, v_z_inic);
      mSrcParticles->mVelocity_x[i] = v_x_inic;
      mSrcParticles->mVelocity_y[i] = v_y_inic;
      mSrcParticles->mVelocity_z[i] = v_z_inic;
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
   clearGrid();

   #pragma omp parallel for
   for (int i = 0; i < mParticleCount; i++)
   {
      // compute a scalar voxel id from a position
      //vec3 pos = mSrcParticles->mPosition[i];
      float pos[3];
      pos[0] = mSrcParticles->mPosition_x[i];
      pos[1] = mSrcParticles->mPosition_y[i];
      pos[2] = mSrcParticles->mPosition_z[i];

      int voxelX = (int)floor(pos[0] * mHTimes2Inv);
      int voxelY = (int)floor(pos[1] * mHTimes2Inv);
      int voxelZ = (int)floor(pos[2] * mHTimes2Inv);
      
      /*
      
      B: Me quise hacer el canchero y comentar esto, pero es FUNDAMENTAL para evitar un segfault...
      
      Go branchless:
      
      // it has been seen the positions can run slightly out of bounds for
      // one solver step. so the positions are temporarily fixed here.
      
      if (voxelX < 0) voxelX= 0;
      if (voxelY < 0) voxelY= 0;
      if (voxelZ < 0) voxelZ= 0;
      if (voxelX >= mGridCellsX) voxelX= mGridCellsX-1;
      if (voxelY >= mGridCellsY) voxelY= mGridCellsY-1;
      if (voxelZ >= mGridCellsZ) voxelZ= mGridCellsZ-1;
      
      */
      
      // Go branchless:
      voxelX = 0 * (voxelX < 0) + 0 * (voxelX >= mGridCellsX) + voxelX * ((voxelX > 0) && (voxelX < mGridCellsX));
      voxelY = 0 * (voxelY < 0) + 0 * (voxelY >= mGridCellsY) + voxelY * ((voxelY > 0) && (voxelY < mGridCellsY));
      voxelZ = 0 * (voxelZ < 0) + 0 * (voxelZ >= mGridCellsZ) + voxelZ * ((voxelZ > 0) && (voxelZ < mGridCellsZ));

      // don't write into particle but into separate memory
      mVoxelCoords[i].x= voxelX;
      mVoxelCoords[i].y= voxelY;
      mVoxelCoords[i].z= voxelZ;

      int voxelId = computeVoxelId(voxelX, voxelY, voxelZ);

      mVoxelIds[i]= voxelId;
   }

   // put each particle into according voxel (sequential)
   for (int i = 0; i < mParticleCount; i++)
   {
       //holaxd;    // B: Not parallel?
       mGrid[ mVoxelIds[i] ].push_back(i);
   }
}


void SPH::findNeighbors(int particleIndex, uint32_t* neighbors, int voxelX, int voxelY, int voxelZ)
{
   float xOrientation = 0.0f;
   float yOrientation = 0.0f;
   float zOrientation = 0.0f;

   int x = 0;
   int y = 0;
   int z = 0;

   int neighborIndex = 0;
   bool enoughNeighborsFound = false;

   float pos[3];
   pos[0] = mSrcParticles->mPosition_x[particleIndex];
   pos[1] = mSrcParticles->mPosition_y[particleIndex];
   pos[2] = mSrcParticles->mPosition_z[particleIndex];

   // this gives us the relative position; i.e the orientation within a voxel
   xOrientation = pos[0] - (voxelX * mHTimes2);
   yOrientation = pos[1] - (voxelY * mHTimes2);
   zOrientation = pos[2] - (voxelZ * mHTimes2);

   // get neighbour voxels
   x = 0;
   y = 0;
   z = 0;

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

   // Para ahorrarse tirar randoms, entre 0 y 4~5 -> LCG?
   int linear_cong_gen;
   int almost_a_random = 0;
   
   float pos_neighbor[3];
   float dot;

   // Variables para usar dentro del loop
   __m256i zeros = _mm256_setzero_si256();
   __m256i ka = _mm256_set_epi32(8,8,8,8,8,8,8,8);
   __m256i jota = _mm256_set_epi32(7,6,5,4,3,2,1,0);  //valor para los 8 jota's
   __m256i ones = _mm256_set1_epi32(1);
   __m256 mH2vec = _mm256_set1_ps(mH2);

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

         const QList<uint32_t>& voxel = mGrid[computeVoxelId(vxi, vyi, vzi)];
         __m256i voxelLength = _mm256_set1_epi32(voxel.length());

         if (!voxel.isEmpty())
         {
            linear_cong_gen = (1664525*(particleIndex + almost_a_random) + 1013904223) % 4294967296;
            const int particleOffset = linear_cong_gen % voxel.length();  //rand() % voxel.length();
            almost_a_random++;
            const int particleIterateDirection = (particleIndex % 2) ? -1 : 1;
            __m256i particleOffsetVec = _mm256_set1_epi32(particleOffset);
            __m256i particleIterateDirectionVec = _mm256_set1_epi32(particleIterateDirection);
            __m256i particleId = _mm256_set1_epi32(particleIndex);

            __m256i iii = _mm256_setzero_si256();
            int maxSteps = (voxel.length() + K - 1) / K;

            for (int step = 0; step < maxSteps; ++step)
            {
               //(off + jota) + i * pID
               __m256i nextIndexs = _mm256_add_epi32(particleOffsetVec, jota);
               __m256i nextIndexsaux = _mm256_mullo_epi32(iii, particleIterateDirectionVec);
               nextIndexs = _mm256_add_epi32(nextIndexs, nextIndexsaux);
                              
               __m256i voxelLength = _mm256_set1_epi32(voxel.length());

               __m256i tooSmall = _mm256_cmpgt_epi32(zeros, nextIndexs);           // nextIndexs < 0
               __m256i tooBig = _mm256_cmpgt_epi32(nextIndexs, _mm256_sub_epi32(voxelLength, _mm256_set1_epi32(1)));  // nextIndexs >= voxel.length()

               __m256i invalid = _mm256_or_si256(tooSmall, tooBig);  // invalid = (nextIndexs < 0 || nextIndexs >= voxel.length())

               __m256i valid = _mm256_cmpeq_epi32(invalid, ones);
               int cmp = _mm256_movemask_ps(_mm256_castsi256_ps(invalid));
               __m256i same = _mm256_cmpeq_epi32(nextIndexs, particleId);

               if (cmp != 0)
                  break;

               uint32_t realIndex[K];
               

               alignas(32) int invalidarray[8];
               _mm256_store_si256((__m256i*)invalidarray, invalid);
               alignas(32) int samearray[8];
               _mm256_store_si256((__m256i*)samearray, same);
               alignas(32) int nextIndexsarray[8];
               _mm256_store_si256((__m256i*)nextIndexsarray, nextIndexs);

               #pragma omp simd  // este no se vectoriza :'v
               for (int j = 0; j < K; j++) {
                     int value = invalidarray[j] ? 1 : 0;
                     int same = samearray[j];
                     int idx = nextIndexsarray[j];
                     realIndex[j] = value ? -2 : (same ? -1 : voxel[idx]);
   
               }

               iii = _mm256_add_epi32(iii,ka);

               int validMask[K];
               float dotVals[K];
               int realNeighbors[K];

               #pragma omp simd // este loop se vectoriza bien
               for (int j = 0; j < K; j++) {
                  int idx = realIndex[j];
                  int isValid = (idx >= 0);

                  float dx = pos[0] - mSrcParticles->mPosition_x[idx];
                  float dy = pos[1] - mSrcParticles->mPosition_y[idx];
                  float dz = pos[2] - mSrcParticles->mPosition_z[idx];

                  dx *= isValid;
                  dy *= isValid;
                  dz *= isValid;

                  dotVals[j] = dx * dx + dy * dy + dz * dz;
                  realNeighbors[j] = idx;
                  validMask[j] = isValid;
               }

               __m256 dotValsV = _mm256_loadu_ps(dotVals);
               __m256 cmp1 = _mm256_cmp_ps(dotValsV, mH2vec, _CMP_LT_OQ); // dotVals < mH2

               int simdValid[K];
               #pragma omp simd // este loop se vectoriza bien
               for (int j = 0; j < K; j++) {
                   simdValid[j] = validMask[j] ? 0xFFFFFFFF : 0x00000000;
               }
               __m256 validMaskV = _mm256_castsi256_ps(_mm256_loadu_si256((__m256i*)simdValid));

               __m256 mask = _mm256_and_ps(validMaskV, cmp1);
               int bitmask = _mm256_movemask_ps(mask);

               #pragma omp simd  // este no se vectoriza :'v
               for (int j = 0; j < K; j++) {
                  if (bitmask & (1 << j)) {
                     neighbors[neighborIndex] = realNeighbors[j];
                     //neighborDistances[neighborIndex] = sqrtf(dotVals[j]) * mSimulationScale;
                     neighborIndex++;
                  }
               }
            
               enoughNeighborsFound = (neighborIndex > mExamineCount - K);
               if (enoughNeighborsFound){
                  break;
               }
            }
         }
      }
   
         if (enoughNeighborsFound)
            break;
      }

   mSrcParticles->mNeighborCount[particleIndex] = neighborIndex;
}


// ??
int SPH::evaluateNeighbor(
   int current,
   int neighbor
)
{
   int validNeighbor = 0;

   // if (current != neighbor)
   // {
   //    vec3 dist = mSrcParticles->mPosition[current] - mSrcParticles->mPosition[neighbor];
   //    float dot = dist * dist;

   //    // the dot product is unscaled and so is MH2;
   //    // so there's no need to add any simulation scale here
   //    if (dot < mH2)
   //    {
   //       validNeighbor = neighbor;
   //    }
   // }

   return validNeighbor;
}



// Try less branches? Recall video...
void SPH::computeDensity(int particleIndex, uint32_t* neighbors)
{
   float density = 0.0f;
   float mass = 0.0f;
   float w = 0.0f;
   float rightPart = 0.0f;
   float distanceScaled;
   
   float ri[3], rj[3];
   ri[0] = mSrcParticles->mPosition_x[particleIndex];
   ri[1] = mSrcParticles->mPosition_y[particleIndex];
   ri[2] = mSrcParticles->mPosition_z[particleIndex];

   for (int neighborIndex = 0; neighborIndex < mSrcParticles->mNeighborCount[particleIndex]; neighborIndex++)
   {
      uint32_t realIndex = neighbors[neighborIndex];

	  // No debería ser necesario esto... (checked before)
      if(realIndex >= mParticleCount)
         break;
      
      if (realIndex != particleIndex)
      {
		   rj[0] = mSrcParticles->mPosition_x[realIndex];
	      rj[1] = mSrcParticles->mPosition_y[realIndex];
	      rj[2] = mSrcParticles->mPosition_z[realIndex];
         // add mass of neighbor
         mass = mSrcParticles->mMass[realIndex];
         
         //distanceScaled = neighborDistances[neighborIndex];  // Again -> Higher complexity
         distanceScaled = (ri[0] - rj[0]) * (ri[0] - rj[0]) +\
     					(ri[1] - rj[1]) * (ri[1] - rj[1]) +\
     					(ri[2] - rj[2]) * (ri[2] - rj[2]);
     	 // Esta ^2!
         
         // apply smoothing kernel to mass; // Branchless?
         rightPart = mHScaled2 - distanceScaled;
         rightPart = rightPart * rightPart * rightPart;
         w = mKernel1Scaled * rightPart * (distanceScaled < mHScaled2);  // true => 1, false => 0 (!)
         // apply weighted neighbor mass to our density
         density += (mass * w);  // 0. if (d^2 > h^2)
         
         //           315
         // w =  -------------  * rightPart
         //      64 * PI * h^9
         
      }
   }

   mSrcParticles->mDensity[particleIndex] = density;
   
}

// Skip this... Calc the pressure on-the-fly, maybe change later (diff EoS)
void SPH::computePressure(int particle)
{
   // rho0: resting density
   // float deltaRho = particle->mDensity - mRho0;
   // float p = mStiffness * deltaRho;
   // particle->mPressure = p;
}


void SPH::computeAcceleration(int particleIndex, uint32_t* neighbors)
{
   Particle* neighbor = 0;
   float distanceToNeighborScaled = 0.0f;

   // OJO, hay muchas cosas que no guardamos dentro de partículas:
   //float pi = p->mPressure;
   float pi = (mSrcParticles->mDensity[particleIndex] - mRho0) * mStiffness;  // One-liner...
   float rhoiInv = ((pi > 0.0f) ? (1.0f / pi) : 1.0f);  // Better?
   float rhoiInv2 = rhoiInv * rhoiInv;
   float piDivRhoi2 = pi * rhoiInv2;
   
   float r[3];
   r[0] = mSrcParticles->mPosition_x[particleIndex];
   r[1] = mSrcParticles->mPosition_y[particleIndex];
   r[2] = mSrcParticles->mPosition_z[particleIndex];
   
   float vi[3];
   vi[0] = mSrcParticles->mVelocity_x[particleIndex];
   vi[1] = mSrcParticles->mVelocity_y[particleIndex];
   vi[2] = mSrcParticles->mVelocity_z[particleIndex];

   float pj = 0.0f;
   float rhoj = 0.0f;
   float rhojInv = 0.0f;
   float rhojInv2 = 0.0f;
   float mj = 0.0f;

   float rj[3];
   float vj[3];

   float rMinusRjScaled[3];

   // pressure gradient...
   float pressureGradient[3] = {0.0f, 0.0f, 0.0f};
   float pressureGradientContribution[3];

   // ...and viscous term
   float viscousTerm[3] = {0.0f, 0.0f, 0.0f};

   // are added to the final acceleration
   float acceleration[3] = {0.0f, 0.0f, 0.0f};

   float centerPart;
   int upper_bound_loop = mSrcParticles->mNeighborCount[particleIndex];
   for (int neighborIndex = 0; neighborIndex < upper_bound_loop; neighborIndex++)
   {
      uint32_t realIndex = neighbors[neighborIndex];

	  rhoj = mSrcParticles->mDensity[realIndex];
      pj = (rhoj - mRho0) * mStiffness;  // One-liner...
      
      rhojInv = 1.0f / rhoj;  // One-liner
      rhojInv = ((rhoj > 0.0f) ? (1.0f / rhoj) : 1.0f);
      rhojInv2 = rhojInv * rhojInv;
      
      rj[0] = mSrcParticles->mPosition_x[realIndex];
      rj[1] = mSrcParticles->mPosition_y[realIndex];
      rj[2] = mSrcParticles->mPosition_z[realIndex];
      
      vj[0] = mSrcParticles->mVelocity_x[realIndex];
      vj[1] = mSrcParticles->mVelocity_y[realIndex];
      vj[2] = mSrcParticles->mVelocity_z[realIndex];

      mj = mSrcParticles->mMass[realIndex];

      // pressure gradient
      rMinusRjScaled[0] = (r[0] - rj[0]) * mSimulationScale;
      rMinusRjScaled[1] = (r[1] - rj[1]) * mSimulationScale;
      rMinusRjScaled[2] = (r[2] - rj[2]) * mSimulationScale;
      
      distanceToNeighborScaled = (rMinusRjScaled[0] * rMinusRjScaled[0]) +\
						      (rMinusRjScaled[1] * rMinusRjScaled[1]) +\
						      (rMinusRjScaled[2] * rMinusRjScaled[2]);
						      
	  distanceToNeighborScaled = sqrtf(distanceToNeighborScaled);

      // Ya sabemos que la distancie > 0 (cuando definimos vecinos validos). However,
      // let's add a ~softening; Use rsqrtf()?
      pressureGradientContribution[0] = mKernel2Scaled * rMinusRjScaled[0] / (distanceToNeighborScaled + 0.01);
      pressureGradientContribution[1] = mKernel2Scaled * rMinusRjScaled[1] / (distanceToNeighborScaled + 0.01);
      pressureGradientContribution[2] = mKernel2Scaled * rMinusRjScaled[2] / (distanceToNeighborScaled + 0.01);

      centerPart = (mHScaled - distanceToNeighborScaled);
      centerPart *= centerPart;
      centerPart *= mj * piDivRhoi2 * (pj * rhojInv2);

      // add pressure gradient contribution to pressure gradient
      pressureGradient[0] += pressureGradientContribution[0] * centerPart;
      pressureGradient[1] += pressureGradientContribution[1] * centerPart;
      pressureGradient[2] += pressureGradientContribution[2] * centerPart;

      // viscosity
      //viscousTermContribution *= (mHScaled - distanceToNeighborScaled); -> "centerpart" (!!!)
      // Reuso variables:
      centerPart = (mHScaled - distanceToNeighborScaled);
      centerPart *= rhojInv * mj * mKernel3Scaled;

      // add contribution to viscous term
      //viscousTerm += viscousTermContribution;
      viscousTerm[0] += (vj[0] - vi[0]) * centerPart;
      viscousTerm[1] += (vj[1] - vi[1]) * centerPart;
      viscousTerm[2] += (vj[2] - vi[2]) * centerPart;

      // Hago acá el viscosityscalar & rho_i^-1
      viscousTerm[0] *= mViscosityScalar * rhoiInv;
      viscousTerm[1] *= mViscosityScalar * rhoiInv;
      viscousTerm[2] *= mViscosityScalar * rhoiInv;

   }

   acceleration[0] = viscousTerm[0] - pressureGradient[0];
   acceleration[1] = viscousTerm[1] - pressureGradient[1];
   acceleration[2] = viscousTerm[2] - pressureGradient[2];


   // WIP: IF rho_i > rho_threshold => kick cinetico!


   /* 
   
   Esto ahora se lo doy a integrate!
	
   // New vecs (grav):
   float gravityTerm[3] = {0.0f, 0.0f, 0.0f};
   float distance_ij3;
   // LASTLY: add a point-mass accel @ the center (e.g. central Black-Hole/NSC)
   // *once per particle. acc = -G M /r^3 (r^, porque apunta al centro)
   rMinusRjScaled[0] = (r[0] - mCentralPos[0]) * mSimulationScale;
   rMinusRjScaled[1] = (r[1] - mCentralPos[1]) * mSimulationScale;
   rMinusRjScaled[2] = (r[2] - mCentralPos[2]) * mSimulationScale;

   float dot = (rMinusRjScaled[0] * rMinusRjScaled[0]) + (rMinusRjScaled[1] * rMinusRjScaled[1]) +\
         (rMinusRjScaled[2] * rMinusRjScaled[2]);
   dot = sqrt(dot);
   distance_ij3 = (dot + mSoftening) * (dot + mSoftening) * (dot + mSoftening);

   // gravityTermContribution = rMinusRjScaled/distance_ij3;
   // gravityTermContribution *= -mGravConstant * mCentralMass;
   gravityTerm[0] = rMinusRjScaled[0]/distance_ij3;
   gravityTerm[1] = rMinusRjScaled[1]/distance_ij3;
   gravityTerm[2] = rMinusRjScaled[2]/distance_ij3;
   // Updateo la gravedad:
   // acceleration += gravityTerm;
   acceleration[0] += -mGravConstant * mCentralMass * gravityTerm[0];
   acceleration[1] += -mGravConstant * mCentralMass * gravityTerm[1];
   acceleration[2] += -mGravConstant * mCentralMass * gravityTerm[2];

   // check CFL condition
   dot = (acceleration[0] * acceleration[0]) + (acceleration[1] * acceleration[1]) +\
         (acceleration[2] * acceleration[2]);

   bool limitExceeded = (dot > mCflLimit2);
   if (limitExceeded)
   {
      float length = sqrt(dot);
      float cflScale = mCflLimit / length;
      acceleration[0] *= cflScale;
      acceleration[1] *= cflScale;
      acceleration[2] *= cflScale;
   }
   
   */

   // Updateo ESTO (la parte hidro), para que desp GPU se encargue de grav + solver...
   mSrcParticles->mAcceleration_x[particleIndex] = acceleration[0];
   mSrcParticles->mAcceleration_y[particleIndex] = acceleration[1];
   mSrcParticles->mAcceleration_z[particleIndex] = acceleration[2];
}


// Deprecated (por el kernel CUDA)
void SPH::integrate(int particleIndex)
{   
   // vec3 position = mSrcParticles->mPosition[particleIndex];
   float position[3];
   position[0] = mSrcParticles->mPosition_x[particleIndex];
   position[1] = mSrcParticles->mPosition_y[particleIndex];
   position[2] = mSrcParticles->mPosition_z[particleIndex];
   // vec3 velocity = mSrcParticles->mVelocity[particleIndex];
   float velocity[3];
   velocity[0] = mSrcParticles->mVelocity_x[particleIndex];
   velocity[1] = mSrcParticles->mVelocity_y[particleIndex];
   velocity[2] = mSrcParticles->mVelocity_z[particleIndex];
   // vec3 acceleration = mSrcParticles->mAcceleration[particleIndex];
   float acceleration[3];
   acceleration[0] = mSrcParticles->mAcceleration_x[particleIndex];
   acceleration[1] = mSrcParticles->mAcceleration_y[particleIndex];
   acceleration[2] = mSrcParticles->mAcceleration_z[particleIndex];

   float mass_here = mSrcParticles->mMass[particleIndex];
   float posTimeStep = mTimeStep * mSimulationScaleInverse;  // ??

   // LF-KDK: Only gravity

   //vec3 velocity_halfstep;
   float velocity_halfstep[3];
   velocity_halfstep[0] = velocity[0] + (acceleration[0] * mTimeStep * 0.5f);
   velocity_halfstep[1] = velocity[1] + (acceleration[1] * mTimeStep * 0.5f);
   velocity_halfstep[2] = velocity[2] + (acceleration[2] * mTimeStep * 0.5f);

   //vec3 newPosition = position + (velocity_halfstep * posTimeStep);
   float newPosition[3];
   newPosition[0] = position[0] + (velocity_halfstep[0] * posTimeStep);
   newPosition[1] = position[1] + (velocity_halfstep[1] * posTimeStep);
   newPosition[2] = position[2] + (velocity_halfstep[2] * posTimeStep);

   // copy & paste grav...
   float rMinusRjScaled[3];
   rMinusRjScaled[0] = (newPosition[0] - mCentralPos[0]) * mSimulationScale;
   rMinusRjScaled[1] = (newPosition[1] - mCentralPos[1]) * mSimulationScale;
   rMinusRjScaled[2] = (newPosition[2] - mCentralPos[2]) * mSimulationScale;

   float dot = rMinusRjScaled[0] * rMinusRjScaled[0] + rMinusRjScaled[1] * rMinusRjScaled[1] +\
         rMinusRjScaled[2] * rMinusRjScaled[2];
   dot = sqrtf(dot);

   float distance_ij3;
   distance_ij3 = (dot + mSoftening) * (dot + mSoftening) * (dot + mSoftening);

   // Updateo la gravedad:
   // acceleration += gravityTerm;
   acceleration[0] = -mGravConstant * mCentralMass * (rMinusRjScaled[0]/distance_ij3);
   acceleration[1] = -mGravConstant * mCentralMass * (rMinusRjScaled[1]/distance_ij3);
   acceleration[2] = -mGravConstant * mCentralMass * (rMinusRjScaled[2]/distance_ij3);

   //vec3 newVelocity = velocity_halfstep + (acceleration * mTimeStep);
   float newVelocity[3];
   newVelocity[0] = velocity_halfstep[0] + (acceleration[0] * mTimeStep);
   newVelocity[1] = velocity_halfstep[1] + (acceleration[1] * mTimeStep);
   newVelocity[2] = velocity_halfstep[2] + (acceleration[2] * mTimeStep);

   dot = newVelocity[0] * newVelocity[0] + newVelocity[1] * newVelocity[1] +\
         newVelocity[2] * newVelocity[2];

   // Muchos NaNs... Skip them:
   /* if (dot > 0)
   {
      // Calc acá T, W y L del sistema (no guardar las energías en las Particles)
      #pragma omp atomic
      mKineticEnergyTotal += 0.5f * mass_here * dot;

      // Energía potencial sería G * Mcentral * m_i/r_i
      #pragma omp atomic
      mPotentialEnergyTotal -= mGravConstant * mCentralMass * mass_here / distance_ij3;
      // + softening);  // B: Without soft (i.e. without a Plummer equivalent)

      // WIP
      //mAngularMomentumTotal += (mass_here * (newPosition - mCentralPos).cross(newVelocity));

   } */

   mSrcParticles->mPosition_x[particleIndex] = newPosition[0];
   mSrcParticles->mPosition_y[particleIndex] = newPosition[1];
   mSrcParticles->mPosition_z[particleIndex] = newPosition[2];

   mSrcParticles->mVelocity_x[particleIndex] = newVelocity[0];
   mSrcParticles->mVelocity_y[particleIndex] = newVelocity[1];
   mSrcParticles->mVelocity_z[particleIndex] = newVelocity[2];
}


// Maybe en algun momento lo podriamos usar (otras IC)
void SPH::handleBoundaryConditions(
   vec3 position,
   vec3* newVelocity,
   float timeStep,
   vec3* newPosition
)
{

/*

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
   
*/

}


// Idem prev...
void SPH::applyBoundary(
      vec3 position,
      float timeStep,
      vec3* newPosition,
      float intersectionDistance,
   vec3 normal,
   vec3* newVelocity
)
{

/*

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
   
*/

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


QList<uint32_t>* SPH::getGrid()
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
