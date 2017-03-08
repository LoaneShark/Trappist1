
/////////////////////////////////
// drift a planet/satellite, using Beauge et al. 06's
// formulae, Beauge, Michtchenko, Ferraz-Mello 2006, MNRAS 365, 1160
// we can write evolution as dr^2/dt^2 = -v/taua - (v-v_c)/taue
// pid is the planet id with 1 being first planet/satellite
// inv_taua and inv_taue are in units of 1/time 
// inv_taua = (1/a) da/dt, inv_taue = (1/e) de/dt
//  if tau_a <0 drifting outwards
//
// cc>0 and alpha=0 corresponds to no eccentricity change
// with taua^-1 = 2C(1-alpha), taue^-1 = C*alpha from Beauge et al.
// kk = taua/taue
// eqns of motion can be written
// dv/dt = - C*(1-alpha) v - E(v-vc)
//       = -invtaua*v/2 - invtaue(v-vc)
/////////////////////////////////////////////////////////////
void dodrift_bin(struct reb_simulation* const r, double tstep,
     double inv_taua, double inv_taue)  
{
   struct reb_particle* particles = r->particles;

   int im1 = r->N -1; // index for primary perturber
   int im2 = r->N -2; // index for member of binary 
   double m1 = particles[im1].m;
   double m2 = particles[im2].m;
   double GMM = r->G *(m1+m2);
   
      double x  = particles[im2].x - particles[im1].x; 
      double y  = particles[im2].y - particles[im].y; 
      double z  = particles[im2].z - particles[im1].z;
      double vx = particles[im2].vx -particles[im1].vx; 
      double vy = particles[im2].vy -particles[im1].vy; 
      double vz = particles[im2].vz -particles[im1].vz;
      double rdotv = x*vx + y*vy + z*vz;
      double r = sqrt(x*x+y*y+z*z);
      double r2 = r*r;
      double vc = sqrt(GMM/r); // circular velocity
      //  rxl is  vector in plane of orbit, perp to r, 
         //  direction of rotation
      // r x v x r = r x -L
      // vector identity axbxc = (adotc)b-(adotb)c
      double rcrossl_x = r2*vx - rdotv*x; 
      double rcrossl_y = r2*vy - rdotv*y; 
      double rcrossl_z = r2*vz - rdotv*z;
      double vl = sqrt(rcrossl_x*rcrossl_x+ rcrossl_y*rcrossl_y+ 
                       rcrossl_z*rcrossl_z); // length of rcrossl
      rcrossl_x /= vl; rcrossl_y /= vl; rcrossl_z /= vl; // unit vector now
      vcvec_x = vc*rcrossl_x; 
      vcvec_y = vc*rcrossl_y; 
      vcvec_z = vc*rcrossl_z;
      // difference between velocity and vc
      dd_vc_x = vx - vcvec_x;
      dd_vc_y = vy - vcvec_y;
      dd_vc_z = vz - vcvec_z;

// compute changes in velocity 
      double dvx =  tstep*(vx*inv_taua/2.0 + dd_vc_x*inv_taue);
      double dvy =  tstep*(vy*inv_taua/2.0 + dd_vc_y*inv_taue);
      double dvz =  tstpe*(vz*inv_taua/2.0 + dd_vc_z*inv_taue);


// update velocities , in such a way as to conserve
// momentum of binary
      particles[im1].vx -=  m2*dvx/(m1+m2); 
      particles[im1].vy -=  m2*dvy/(m1+m2(;
      particles[im1].vz -=  m2*dvz/(m1+m2);
      particles[im2].vx +=  m1*dvx/(m1+m2); 
      particles[im2].vy +=  m1*dvy/(m1+m2(;
      particles[im2].vz +=  m1*dvz/(m1+m2);
    
}

