/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// epsilon1: attractive region (sigma < 2^1/6); epsilon2: repulsive region.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_lj_twopiece.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLJ2P::PairLJ2P(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJ2P::~PairLJ2P()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon1);
    memory->destroy(epsilon2);
    memory->destroy(sigma);
    memory->destroy(sigmacut2);
    memory->destroy(eshift);
    memory->destroy(lj11);
    memory->destroy(lj12);
    memory->destroy(lj13);
    memory->destroy(lj14);
    memory->destroy(lj21);
    memory->destroy(lj22);
    memory->destroy(lj23);
    memory->destroy(lj24);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJ2P::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        if (rsq >= sigmacut2[itype][jtype]) {
          forcelj = r6inv * (lj11[itype][jtype]*r6inv - lj12[itype][jtype]);
        }
        else {
          forcelj = r6inv * (lj21[itype][jtype]*r6inv - lj22[itype][jtype]);
        }
        //forcelj = r6inv * (lj11[itype][jtype]*r6inv - lj12[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;
        /* if (fpair > 1000) {
          printf("f %g rsq %g i %d j %d itype %d jtype %d\n", fpair, rsq, i, j, itype, jtype);
        } */

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq >= sigmacut2[itype][jtype]) {
          evdwl = r6inv*(lj13[itype][jtype]*r6inv-lj14[itype][jtype]) -
            offset[itype][jtype];
          }
          else {
            evdwl = r6inv*(lj23[itype][jtype]*r6inv-lj24[itype][jtype]) -
            offset[itype][jtype] + eshift[itype][jtype];
          }
          //evdwl = r6inv*(lj13[itype][jtype]*r6inv-lj14[itype][jtype]) -
          //  offset[itype][jtype];
          evdwl *= factor_lj;
          /* if (evdwl > 1000) {
            printf("e %g rsq %g i %d j %d itype %d jtype %d\n", evdwl, rsq, i, j, itype, jtype);
          } */
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJ2P::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon1,n+1,n+1,"pair:epsilon1");
  memory->create(epsilon2,n+1,n+1,"pair:epsilon2");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(sigmacut2,n+1,n+1,"pair:sigmacut2");
  memory->create(eshift,n+1,n+1,"pair:eshift");
  memory->create(lj11,n+1,n+1,"pair:lj11");
  memory->create(lj12,n+1,n+1,"pair:lj12");
  memory->create(lj13,n+1,n+1,"pair:lj13");
  memory->create(lj14,n+1,n+1,"pair:lj14");
  memory->create(lj21,n+1,n+1,"pair:lj21");
  memory->create(lj22,n+1,n+1,"pair:lj22");
  memory->create(lj23,n+1,n+1,"pair:lj23");
  memory->create(lj24,n+1,n+1,"pair:lj24");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJ2P::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJ2P::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon1_one = force->numeric(FLERR,arg[2]);
  double epsilon2_one = force->numeric(FLERR,arg[3]);
  double sigma_one = force->numeric(FLERR,arg[4]);

  double cut_one = cut_global;
  if (narg == 6) cut_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon1[i][j] = epsilon1_one;
      epsilon2[i][j] = epsilon2_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJ2P::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    printf("setflag %d %d", i, j);
    epsilon1[i][j] = mix_energy(epsilon1[i][i],epsilon1[j][j],
                               sigma[i][i],sigma[j][j]);
    epsilon2[i][j] = mix_energy(epsilon2[i][i],epsilon2[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }
  //printf("sigma %d %d %g\n", i, j, sigma[i][j]);
  lj11[i][j] = 48.0 * epsilon1[i][j] * pow(sigma[i][j],12.0);
  lj12[i][j] = 24.0 * epsilon1[i][j] * pow(sigma[i][j],6.0);
  lj13[i][j] = 4.0 * epsilon1[i][j] * pow(sigma[i][j],12.0);
  lj14[i][j] = 4.0 * epsilon1[i][j] * pow(sigma[i][j],6.0);

  lj21[i][j] = 48.0 * epsilon2[i][j] * pow(sigma[i][j],12.0);
  lj22[i][j] = 24.0 * epsilon2[i][j] * pow(sigma[i][j],6.0);
  lj23[i][j] = 4.0 * epsilon2[i][j] * pow(sigma[i][j],12.0);
  lj24[i][j] = 4.0 * epsilon2[i][j] * pow(sigma[i][j],6.0);

  eshift[i][j] = epsilon2[i][j] - epsilon1[i][j];

  if (offset_flag && (cut[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon1[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  lj11[j][i] = lj11[i][j];
  lj12[j][i] = lj12[i][j];
  lj13[j][i] = lj13[i][j];
  lj14[j][i] = lj14[i][j];
  lj21[j][i] = lj21[i][j];
  lj22[j][i] = lj22[i][j];
  lj23[j][i] = lj23[i][j];
  lj24[j][i] = lj24[i][j];
  eshift[j][i] = eshift[i][j];
  offset[j][i] = offset[i][j];
  sigmacut2[i][j] = sigma[i][j]*sigma[i][j]*cbrt(2.0);
  sigmacut2[j][i] = sigmacut2[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon1[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon1[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJ2P::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon1[i][j],sizeof(double),1,fp);
        fwrite(&epsilon2[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJ2P::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon1[i][j],sizeof(double),1,fp);
          fread(&epsilon2[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJ2P::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJ2P::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJ2P::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,epsilon1[i][i],epsilon2[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJ2P::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,epsilon1[i][j],epsilon2[i][j],sigma[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJ2P::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r2inv,r6inv,forcelj,philj;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  if (rsq >= sigmacut2[itype][jtype]) {
    forcelj = r6inv * (lj11[itype][jtype]*r6inv - lj12[itype][jtype]);
    philj = r6inv*(lj13[itype][jtype]*r6inv-lj14[itype][jtype]) -
    offset[itype][jtype];
  }
  else {
    forcelj = r6inv * (lj21[itype][jtype]*r6inv - lj22[itype][jtype]);
    philj = r6inv*(lj23[itype][jtype]*r6inv-lj24[itype][jtype]) -
    offset[itype][jtype] + eshift[itype][jtype];
  }
  fforce = factor_lj*forcelj*r2inv;
  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

void *PairLJ2P::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon1") == 0) return (void *) epsilon1;
  if (strcmp(str,"epsilon2") == 0) return (void *) epsilon2;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return NULL;
}
