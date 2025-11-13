/*! \file Y3Dfd.c
 *  \brief Y nodal forces
 *
 *  \todo
 *  finish documentation
 *
 *  Copyright (C) 2008, Queen Mary University of London (QMUL) &
 *  Imperial College of Science, Technology and Medicine (ICSTM).
 *  All rights reserved. Implemented by Prof Antonio Munjiza &
 *  Dr Jiansheng Xiang.
 *
 *  This code is part of the Virtual Geoscience Workbench (VGW)
 *  developed jointly by ICSTM and QMUL through two related parallel
 *  projects at ICSTM and QMUL respectively funded by EPSRC.
 *
 *  This code is provided by copyright holders under the GNU Lesser
 *  General Public License (LGPL). It is open source code; you can
 *  redistribute it and/or modify it under the terms of the GNU Lesser
 *  General Public License version 3.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
 *  the GNU Lesser General Public License for more details,
 *  http://www.gnu.org/licenses/lgpl-3.0.txt.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this code; if not, write to:
 *
 *  Dr Jiansheng Xiang   <j.xiang@imperial.ac.uk>    \n
 *  Prof Antonio Munjiza <a.munjiza@qmul.ac.uk>      \n
 *  Dr John-Paul Latham  <j.p.latham@imperial.ac.uk> \n
 *
 */


/**********************************************************************/
/**********************************************************************/


#include "Yproto.h"


/* check for negative values - Z - what significance?
*/
#if (defined(DEBUG) || defined(_DEBUG))
# define DBG_NEG(d,i) if((d)<=R0){ CHRw(stdout,"Negative Jacobian for element "); INTw(stdout,i,5); CHRwcr(stdout); }
#else
# define DBG_NEG(d,i)
#endif


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


/*calculate fragment size distribution*/

static INLINE void calFragVol(INT nelem,INT *i1elfrg,INT *NElemSet, DBL *VolElemSet,
		DBL *d1elvol, INT **i2eljp,INT *i1elpr,INT *i1pejp,INT **i2ejjp,
        INT *i1elchprm, INT *i1nochprm, INT **i2elto, DBL *d1nccx,DBL *d1nccy,DBL *d1nccz,
	DBL *d1nvcx,DBL *d1nvcy,DBL *d1nvcz, INT nelno,INT ichprm,INT iremoval,
	INT indexv,DBL dcrmvx,DBL dcrmvy,DBL dcrmvz,
       	DBL **d2nfrm,INT *i1elfrm, INT *i1frgrm,DBL *d1ndisp,DBL dcfrmca,
	INT *i1nopr, INT *i1bnvx,INT *i1bnvy,INT *i1bnvz,INT *i1rmgrav, DBL *d1mart)
{
    INT i,j,k,istart,isurf,jsurf,flag;
    int NoIncrease;
    INT ielem,jelem,ijoint;
    DBL *minix, *miniy,*miniz, dis_cen,*minal, temp;
    INT *i1remov,*i1removal,*i1remov_bc;

//CHRw(stdout,"enter calFragVol");
    istart=0;

      for(i=istart;i<nelem;i++)
    {
        i1elfrg[i]=-1;
 //      INTw(stdout,i,10);
    }

    *NElemSet=istart;
      for(i=istart;i<nelem;i++)
    {
	    if(i1pejp[i1elpr[i]]>0)
	    {
	    NoIncrease = 1;
        if(i1elfrg[i]==-1)
        {    i1elfrg[i]=*NElemSet;
          do
          {
              NoIncrease = 0;
              for(j=istart;j<nelem;j++)
              {
		      if(i1pejp[i1elpr[j]]>0)
            {
                  if((i1elfrg[j]==*NElemSet))
                  {
                      for(isurf=0;isurf<4;isurf++)
                        {
                            ijoint=i2eljp[j][isurf];
                            if(ijoint>=0)
                            {
                                ielem=i2ejjp[ijoint][0];
                                jelem=i2ejjp[ijoint][1];
                                if((ielem>=0)&&(jelem>=0)&&(ielem!=jelem))
                                {
                                if((ielem==j)&&(i1elfrg[jelem]==-1))
                                {
                                    i1elfrg[jelem]=*NElemSet;
                                    NoIncrease=1;
                                }
                                else if((jelem==j)&&(i1elfrg[ielem]==-1))
                                {
                                    i1elfrg[ielem]=*NElemSet;
                                    NoIncrease=1;
                                }
                                }
                            }
                            }

                        }

            }
         }
       }
	  while(NoIncrease==1);
         *NElemSet=*NElemSet+1;
        }
    }
    }
//INTw(stdout, *NElemSet,10);
    minix=TalDBL1(*NElemSet+1);
    miniy=TalDBL1(*NElemSet+1);
    miniz=TalDBL1(*NElemSet+1);
    minal=TalDBL1(*NElemSet+1);
    i1remov=TalINT1(*NElemSet+1);
    i1removal=TalINT1(*NElemSet+1);
    i1remov_bc=TalINT1(*NElemSet+1);
    for(i=0;i<*NElemSet;i++)
        {
            VolElemSet[i]=0.0;
            minix[i]=BEPSILON;
            miniy[i]=BEPSILON;
            miniz[i]=BEPSILON;
	    minal[i]=EPSILON;
            i1remov[i]=-1;
	    i1removal[i]=-1;
	    i1remov_bc[i]=-1;
        }


    for(i=istart;i<nelem;i++)
        {

            if(i1pejp[i1elpr[i]]>0)
            {
                j=i1elfrg[i];
                VolElemSet[j] = VolElemSet[j]+d1elvol[i]/d1mart[i] ;
		if(i1elchprm[i]!=1)
		{
                for(k=0;k<nelno-1;k++)
                {
                    minix[j]=MINIM(minix[j],d1nccx[i2elto[i][k]]);
                    miniy[j]=MINIM(miniy[j],d1nccy[i2elto[i][k]]);
                    miniz[j]=MINIM(miniz[j],d1nccz[i2elto[i][k]]);
		    if(iremoval)
		    {
			    i1rmgrav[i2elto[i][k]]=1;
		    temp=SQR(d1nccx[i2elto[i][k]]-d2nfrm[0][i2elto[i][k]])+
			    SQR(d1nccy[i2elto[i][k]]-d2nfrm[1][i2elto[i][k]])+
			    SQR(d1nccz[i2elto[i][k]]-d2nfrm[2][i2elto[i][k]]);
		    temp=SQRT(temp);
		    minal[j]=MAXIM(minal[j],temp);
		    d1ndisp[i2elto[i][k]]=temp;
		    if((i1bnvx[i1nopr[i2elto[i][k]]]>0)||
	              (i1bnvy[i1nopr[i2elto[i][k]]]>0)||
	    	      (i1bnvz[i1nopr[i2elto[i][k]]]>0))
		    {
			    i1remov_bc[j]=1;
		    }
		    }
		    dis_cen=SQRT(d1nccx[i2elto[i][k]]*d1nccx[i2elto[i][k]]+
				    d1nccz[i2elto[i][k]]*d1nccz[i2elto[i][k]]);
		    if(dis_cen>0.03)i1remov[j]=1;
                }
		}
	    }
        }
    NoIncrease=0;
    for(i=0;i<*NElemSet;i++)
        {
		i1frgrm[i]=0;
            if((miniz[i]>0.025)&&(ichprm==1))
	    {
		    i1remov[i]=1;
		    NoIncrease=1;
	    }
	    if((minal[i]>dcfrmca)&&(iremoval==1)&&(i1remov_bc[i]!=1))
	    {
		    i1removal[i]=1;
		    i1frgrm[i]=1;
	    }
        }
				for(i=istart;i<nelem;i++)
		        {

		                j=i1elfrg[i];
				i1elfrm[i]=0;
				if(i1removal[j]>0)
				{
			             i1elfrm[i]=1;
				}
				if((i1remov_bc[j]==1)&&(iremoval==1))
				{
				for(k=0;k<nelno-1;k++)
		                {
	                        i1rmgrav[i2elto[i][k]]=0;
				}
				}
				if((i1remov_bc[j]!=1)&&(indexv==1))
			 {
			 for(k=0;k<nelno-1;k++)
			 {
			 d1nvcx[i2elto[i][k]]=d1nvcx[i2elto[i][k]]+dcrmvx;
			 d1nvcy[i2elto[i][k]]=d1nvcy[i2elto[i][k]]+dcrmvy;
			 d1nvcz[i2elto[i][k]]=d1nvcz[i2elto[i][k]]+dcrmvz;
			 }
			 }

			}
    if(NoIncrease)
    {
    for(i=istart;i<nelem;i++)
        {
            if((i1pejp[i1elpr[i]]>0)&&(i1elchprm[i]!=1))
            {
                j=i1elfrg[i];
		if(i1removal[j]>0)
		{
	             i1elfrm[i]=1;
		}

                if(i1remov[j]>0)
                {
                    i1elchprm[i]=1;
                for(k=0;k<nelno-1;k++)
                {
                    i1nochprm[i2elto[i][k]]=1;
                    d1nvcx[i2elto[i][k]]=R0;
                    d1nvcy[i2elto[i][k]]=R0;
                    d1nvcz[i2elto[i][k]]=R0;

                }
                }
	    }
        }
    }
    FREE(i1remov_bc);
    FREE(i1removal);
    FREE(i1remov);
    FREE(miniz);
    FREE(miniy);
    FREE(minix);
}

/* Calculate the Jacobian matrix xj
*/
static INLINE void jacobi(DBL d2xj[NDIME][NDIME], DBL **d2ncc, DBL **d2dsh, INT *i1elto)
{
  INT i, j, k, ip;
  DBL *d1dsh;
  DBL *d1ncc, *d1xj;

  for(i=0; i<NDIME; i++)
  {
    d1ncc = d2ncc[i];
    d1xj  = d2xj[i];

    for(j=0; j<NDIME; j++)
    {
      d1dsh   = d2dsh[j];
      d1xj[j] = R0;

      for(k=0; k<NNODE; k++)
      {
        ip = i1elto[k];
        d1xj[j] += d1ncc[ip]*d1dsh[k];
      }
    }
  }
}


/**********************************************************************/
/**********************************************************************/


/* Calculate the Jacobian matrices xj, xjci
*/
static INLINE void jacobi2(DBL d2xj[NDIME][NDIME], DBL **d2ncc, DBL d2xjci[NDIME][NDIME], DBL **d2nci, DBL **d2dsh, INT *i1elto)
{
  INT i, j, k, ip;
  DBL *d1dsh;
  DBL *d1ncc, *d1xj;
  DBL *d1nci, *d1xjci;

  for(i=0; i<NDIME; i++)
  {
    d1ncc  = d2ncc[i];
    d1nci  = d2nci[i];
    d1xj   = d2xj[i];
    d1xjci = d2xjci[i];

    for(j=0; j<NDIME; j++)
    {
      d1dsh     = d2dsh[j];
      d1xj[j]   = R0;
      d1xjci[j] = R0;

      for(k=0; k<NNODE; k++)
      {
        ip = i1elto[k];
        d1xj[j]   += d1ncc[ip]*d1dsh[k];
        d1xjci[j] += d1nci[ip]*d1dsh[k];
      }
    }
  }
}


/**********************************************************************/
/**********************************************************************/


/* setup for gradient calculations
*/
static INLINE void grad(DBL d2nd[NDIME][NNODE], DBL d2xjinv[NDIME][NDIME], DBL **d2dsh)
{
  INT i;

  for(i=0; i<NNODE; i++)
  {
    d2nd[0][i] = (d2xjinv[0][0]*d2dsh[0][i]) + (d2xjinv[1][0]*d2dsh[1][i]) + (d2xjinv[2][0]*d2dsh[2][i]);
    d2nd[1][i] = (d2xjinv[0][1]*d2dsh[0][i]) + (d2xjinv[1][1]*d2dsh[1][i]) + (d2xjinv[2][1]*d2dsh[2][i]);
    d2nd[2][i] = (d2xjinv[0][2]*d2dsh[0][i]) + (d2xjinv[1][2]*d2dsh[1][i]) + (d2xjinv[2][2]*d2dsh[2][i]);
  }
}


/**********************************************************************/
/**********************************************************************/

/* calculate the Cartesian gradient d2finv
*/
static INLINE void ctgrad(DBL d2nd[NDIME][NNODE], DBL d2finv[NDIME][NDIME], DBL **d2nci, INT *i1elto)
{
  INT i, j, k, ip;
  DBL *d1nd;
  DBL *d1finv, *d1nci;

  /* d2nd must be complete
  */
  for(i=0; i<NDIME; i++)
  {
    d1finv = d2finv[i];
    d1nci  = d2nci[i];

    for(j=0; j<NDIME; j++)
    {
      d1nd      = d2nd[j];
      d1finv[j] = R0;

      for(k=0; k<NNODE; k++)
      {
        ip = i1elto[k];
        d1finv[j] += d1nci[ip]*d1nd[k];
      }
    }
  }
}


/**********************************************************************/
/**********************************************************************/


/* calculate the Cartesian gradient d2finv, and velocity gradient L
*/
static INLINE void ctvgrad(DBL d2nd[NDIME][NNODE], DBL d2finv[NDIME][NDIME], DBL **d2nci, DBL d2vgrad[NDIME][NDIME], DBL **d2nvc, INT *i1elto)
{
  INT i, j, k, ip;
  DBL *d1nd;
  DBL *d1finv,  *d1nci;
  DBL *d1vgrad, *d1nvc;

  /* d2nd must be complete DBL d2lvg[NDIME][NDIME],
  */
  for(i=0; i<NDIME; i++)
  {
    d1finv = d2finv[i];
    d1nci  = d2nci[i];

    d1vgrad = d2vgrad[i];
    d1nvc   = d2nvc[i];

    for(j=0; j<NDIME; j++)
    {
      d1nd       = d2nd[j];
      d1finv[j]  = R0;
      d1vgrad[j] = R0;

      for(k=0; k<NNODE; k++)
      {
        ip = i1elto[k];
        d1finv[j]  += d1nci[ip]*d1nd[k];
        d1vgrad[j] += d1nvc[ip]*d1nd[k];
      }
    }
  }
}


/**********************************************************************/
/**********************************************************************/


/* calculate calculate tensor B, Cauchy stress T, tcs
*/
static INLINE void cauchy(DBL T[NDIME][NDIME], DBL L[NDIME][NDIME], DBL d2ftens[NDIME][NDIME],
                 //         DBL ***d3tcs,
                          DBL dpemu, DBL dpela, DBL detf, DBL detf0, DBL dpeks, INT ielem)
{
  INT i, j;
  DBL scale1, scale2, tweak, D, B;

  scale1 = pow(detf0/(R1/detf),R2/R3);
  scale2 = dpemu/detf0;
  tweak = (dpela*log(detf0)-dpemu)/detf0;

  for(i=0; i<NDIME; i++)
  {
    /* hopefully the compiler will unroll and expand this entirely */
    for(j=0; j<NDIME; j++)
    {
      /* left Cauchy-Green strain tensor B */
      B = ((d2ftens[i][0]*d2ftens[j][0]) + (d2ftens[i][1]*d2ftens[j][1]) + (d2ftens[i][2]*d2ftens[j][2])) * scale1;

      /* deformation rate D (NOTE: [j][i] access forces L to be complete prior) */
      D = (L[i][j] + L[j][i])*RP5;

      /* Cauchy stress T */
      T[i][j] = B*scale2 + D*dpeks;
    }

    /* Cauchy stress along diagonal of T (Neo-Hookean's law) */
    T[i][i] += tweak;
  }

}


/**********************************************************************/
/**********************************************************************/


/* calculate surface traction d3st[ig] (tmp1)
*/
static INLINE void stract(DBL d2st[NDIME][NNODE], DBL T[NDIME][NDIME], DBL d2nd[NDIME][NNODE], DBL detjb)
{
  INT i;

  for(i=0; i<NNODE; i++)
  {
    d2st[0][i] = ((d2nd[0][i]*T[0][0]) + (d2nd[1][i]*T[0][1]) + (d2nd[2][i]*T[0][2]))*detjb;
    d2st[1][i] = ((d2nd[0][i]*T[1][0]) + (d2nd[1][i]*T[1][1]) + (d2nd[2][i]*T[1][2]))*detjb;
    d2st[2][i] = ((d2nd[0][i]*T[2][0]) + (d2nd[1][i]*T[2][1]) + (d2nd[2][i]*T[2][2]))*detjb;
  }
}


/**********************************************************************/
/**********************************************************************/


/* calculate nodal force
*/
static INLINE void nforce(DBL *d1nftx, DBL *d1nfty, DBL *d1nftz, DBL d3st[NGRSH][NDIME][NNODE], INT *i1elto)
{
  INT i, ip;

  for(i=0; i<NNODE; i++)
  {
    ip = i1elto[i];

    /* nodal total current force */
    d1nftx[ip] -= WEIGHT*(d3st[1][0][i] + d3st[2][0][i] + d3st[3][0][i] + d3st[4][0][i]);
    d1nfty[ip] -= WEIGHT*(d3st[1][1][i] + d3st[2][1][i] + d3st[3][1][i] + d3st[4][1][i]);
    d1nftz[ip] -= WEIGHT*(d3st[1][2][i] + d3st[2][2][i] + d3st[3][2][i] + d3st[4][2][i]);

    /* nodal mass current translation */
  }
}


/**********************************************************************/
/**********************************************************************/


/* small strain elastic tetrahetra
*/
static void Yfd3TET10ELS(
            INT nelem,
            DBL *dpeks,    DBL *dpela,    DBL *dpemu,
            DBL **d2ncc,  DBL **d2nci,  DBL **d2nft,    DBL **d2nvc,
            INT *i1elpr,  INT **i2elto,
            DBL ***d3dsh, INT nthreads, INT nblock, INT nnopo, INT nmode
            )
{

#pragma omp parallel default(none) shared(nelem,dpeks,dpela,dpemu,d2ncc,d2nci,\
d2nft,d2nvc,i1elpr,i2elto,d3dsh,nnopo,nthreads,nblock,nmode)
    {
        DBL d2nd[NDIME][NNODE], d2finv[NDIME][NDIME], d2ftens[NDIME][NDIME];
        DBL d2xj[NDIME][NDIME], d2xjci[NDIME][NDIME], d2xjinv[NDIME][NDIME];
        DBL L[NDIME][NDIME], T[NDIME][NDIME];
        DBL d3st[NGRSH][NDIME][NNODE];
        DBL voli, detci, detj, detf, detf0;
        INT ielem, ig;
        INT *i1elto;
        INT iprop;
        INT i,ifrom,j,ito;
        int tid;
        DBL **d2nfs;

        d2nfs=TalDBL2(4,nnopo+1);
        TzDBL2(d2nfs,4,nnopo+1);

#pragma omp for schedule(guided) //nowait
  for(ielem=0; ielem<nelem; ielem++)
  {
    iprop=i1elpr[ielem];
    i1elto = i2elto[ielem];

    /* simpler calculations for ig=0 */
    ig = 0;
    /* calculate Jacobian matrix xj */
    jacobi(d2xj, d2ncc, d3dsh[ig], i1elto);
    voli = R0;
    /* calculate inverse of xj (debug: warn if determinant negative) */
    YMATINV3(d2xj, d2xjinv, detj); DBG_NEG(detj, ielem)
    /* calculate Cartesian gradients nd, finv */
    grad(d2nd, d2xjinv, d3dsh[ig]);
    ctgrad(d2nd, d2finv, d2nci, i1elto);
    /* calculate determinant of finv */
    YMATDET3(d2finv, detf);

    /* store determinant */
    detf0=R1/detf;

    /* iterate overy Gauss points using deformation (shape) gradient */
    for(ig=1;ig<NGRSH;ig++)
    {
      /* full calculations for ig>0 */
      /* calculate the Jacobian matrixes xj, xjci */
      jacobi2(d2xj, d2ncc, d2xjci, d2nci, d3dsh[ig], i1elto);
      /* calculate determinant of xjci */
      YMATDET3(d2xjci,detci);
	    voli += ABS(detci)*WEIGHT/R6;
      /* calculate inverse and determinant of xj (debug: warn if negative) */
      YMATINV3(d2xj,d2xjinv,detj); DBG_NEG(detj,ielem)
      /* calculate Cartesian gradients nd, finv, and velocity gradient L */
      grad(d2nd, d2xjinv, d3dsh[ig]);
      ctvgrad(d2nd, d2finv, d2nci, L, d2nvc, i1elto);
      /* calculate inverse and determinant of finv */
      YMATINV3(d2finv, d2ftens, detf);

      /* calculate left Cauchy-Green strain tensor B, deformation rate D, Cauchy stress T, net stress tcs */
      cauchy(T, L, d2ftens, dpemu[iprop], dpela[iprop], detf, detf0, dpeks[iprop], ielem);

      /* end result for ig=1,2,3,4 */
      stract(d3st[ig], T, d2nd, detj/R6);
    }

    /* calculate nodal force */
    nforce(d2nfs[0], d2nfs[1], d2nfs[2], d3st, i1elto);
  }
        tid = omp_get_thread_num();
        for(i=0;i<nthreads;i++)
        {
            ifrom=tid+i;
              if(ifrom>=nthreads)ifrom=ifrom-nthreads;
              if((nmode!=0)&&((ifrom+nmode)>=nthreads))
              {
                  ito=ifrom+1;
                  ifrom=ifrom*nblock-nthreads+ifrom+nmode;
                  ito=ito*nblock-nthreads+ito+nmode;
              }
              else
              {
                  ifrom=ifrom*nblock;
                  ito=ifrom+nblock;
              }
            for(j=ifrom;j<ito;j++)
            {
                d2nft[0][j]+=d2nfs[0][j];
                d2nft[1][j]+=d2nfs[1][j];
                d2nft[2][j]+=d2nfs[2][j];
            }

#pragma omp barrier
        }
        FREE(d2nfs);

    }
}

static INLINE DBL FdftensTet4(DBL F[3][3],DBL L[3][3],DBL *d1nccx,  DBL *d1nccy,  DBL *d1nccz,
                    DBL *d1ncix,  DBL  *d1nciy,
                    DBL *d1nciz,DBL *d1nvcx,  DBL *d1nvcy,  DBL *d1nvcz, INT *i1elto, DBL *d1elvol)
{
    INT j,k; int i,ii;
    DBL FX[3][3]; /* current local base  also delta ux/delta X */
    DBL F0inv[3][3]; /* global base in initial local base */
    DBL FXinv[3][3]; /* global base in current local base */
    DBL LX[3][3]; /* velocity gradient in current local base = delta x/delta X */
    DBL voli,volc,detf;
k = i1elto[0];
for(i=1;i<4;i++)
{   j=i1elto[i];
    ii=i-1;
    F[0][ii]=d1ncix[j]-d1ncix[k]; /* init. base */
    F[1][ii]=d1nciy[j]-d1nciy[k];
    F[2][ii]=d1nciz[j]-d1nciz[k];
    FX[0][ii]=d1nccx[j]-d1nccx[k]; /* curr. base */
    FX[1][ii]=d1nccy[j]-d1nccy[k];
    FX[2][ii]=d1nccz[j]-d1nccz[k];
    LX[0][ii]=d1nvcx[j]-d1nvcx[k]; /* vel. grad. */
    LX[1][ii]=d1nvcy[j]-d1nvcy[k];
    LX[2][ii]=d1nvcz[j]-d1nvcz[k];
}

YMATINV3(F,F0inv,voli);         /* global base in initial local coordinates    */

YMATINV3(FX,FXinv,volc);         /* global base in current local coordinates    */

    *d1elvol=volc/6.0;
detf=volc/voli;

    __builtin_prefetch(&L[0][0],1,1);
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            F[i][j]=R0;
            L[i][j]=R0;
        }
    }
for(i=0;i<3;i++)
{
    for(k=0;k<3;k++)
    {
        for(j=0;j<3;j++)
       {

			F[i][j]=F[i][j]+FX[i][k]*F0inv[k][j]; /* deformation gradient  */
            L[i][j]=L[i][j]+LX[i][k]*FXinv[k][j]; /* velocity gradient     */
        }
    }
}
    return detf;
}


static INLINE void CauchyTet4_pore_pressure(DBL F[3][3],DBL L[3][3],DBL dpemu,DBL dpeks,DBL dpela,DBL detf,DBL T[3][3], DBL D[3][3])
{
    DBL  B[3][3]; /* left Cauchy-Green strain tensor */
//    DBL  D[3][3]; /* rate of deformation (stretching) tensor */
    int i,j,k;
    DBL dups_u; /* [ielem] poisson ratio - undrained  */
    DBL dups; /* [ielem] poisson ratio - drained  */
    DBL ALPHA;
  DBL dredf,SK,Tm;

          dups = dpela / (2 * (dpela + dpemu));
          dups_u = dups/0.8;  /* using 20% or lower approaximate at the moment - Devonian Shale; Hydrostone -   */
    dups_u=0.3839;

  SK=0.644559576;
//  SK=R1;
  /* Biot coefficient -- input from .y in the future */
          ALPHA=3.0*(dups_u - dups)/ ((1.0 - 2.0 * dups) * (1.0 + dups_u));

for(i=0;i<3;i++)
{
    for(j=0;j<3;j++)
    {
        B[i][j]=R0;
        for(k=0;k<3;k++)
        {
            B[i][j]=B[i][j]+F[i][k]*F[j][k];  /* left Cauchy-Green strain  */
        }
        D[i][j]=RP5*(L[i][j]+L[j][i]);      /* rate of deformation       */
    }
}
/*Calculate Cauchy stress*/

for(i=0;i<3;i++)
{
    for(j=0;j<3;j++)
    {
        T[i][j]=(dpemu/detf)*B[i][j]+dpeks*D[i][j];
    }
    T[i][i]=T[i][i]+(dpela*log(detf)-dpemu)/detf;
}
    Tm=(T[0][0]+T[1][1]+T[2][2])/3.0;
      Tm=-(detf-1.0)*2.0*SK*dpemu*(1.0 + dups_u)/(3.0*(1.0-2.0*dups_u));
      Tm=R0;
for(i=0;i<3;i++)
  {
        T[i][i]=T[i][i]-ALPHA/SK*Tm;

  }

Tm=-ALPHA/SK*Tm;
}
static INLINE void CauchyTet4(DBL F[3][3],DBL L[3][3],DBL dpemu,DBL dpeks,DBL dpela,DBL detf,
		DBL T[3][3],DBL D[3][3], DBL *dem0n,DBL dem0,DBL demf,DBL dpedm, DBL *deldam)
{
    DBL  B[3][3]; /* left Cauchy-Green strain tensor */
    int i,j,k;
    DBL df;
for(i=0;i<3;i++)
{
    for(j=0;j<3;j++)
    {
        B[i][j]=R0;
        for(k=0;k<3;k++)
        {
			B[i][j]=B[i][j]+F[i][k]*F[j][k];  /* left Cauchy-Green strain  */
        }
        D[i][j]=RP5*(L[i][j]+L[j][i]);      /* rate of deformation       */
    }
}
/*Calculate Cauchy stress*/

/*calculate damage factor */
df=R0;
*deldam=df;

/*Calculate Cauchy stress*/
for(i=0;i<3;i++)
{
    for(j=0;j<3;j++)
    {
        T[i][j]=(R1-df)*(dpemu/detf)*B[i][j]+dpeks*D[i][j];
    }
    T[i][i]=T[i][i]+(R1-df)*(dpela*log(detf)-dpemu)/detf;
}
}
static INLINE void nforceTet4(DBL T[3][3],DBL *d1nccx,DBL *d1nccy,DBL *d1nccz,
                              DBL *d1nfcx,DBL *d1nfcy,DBL *d1nfcz,INT *i1elto,
			      DBL *d1elfsx,  DBL *d1elfsy,  DBL *d1elfsz)
{
    DBL nx,ny,nz;
    int i;
    INT ii,jj,kk,ll;
    DBL tx,ty,tz;
    int p1[4]={1,2,3,0};
    int p2[4]={2,3,0,1};
    int p3[4]={3,0,1,2};

for(i=0;i<4;i++)      /* Nodal Forces */
{
    ii=i1elto[i];
    jj=i1elto[p1[i]];
    kk=i1elto[p2[i]];
    ll=i1elto[p3[i]];
                nx=((d1nccy[kk]-d1nccy[jj])*
                    (d1nccz[ll]-d1nccz[jj])-
                    (d1nccy[ll]-d1nccy[jj])*
                    (d1nccz[kk]-d1nccz[jj]))/R6;
                ny=((d1nccz[kk]-d1nccz[jj])*
                    (d1nccx[ll]-d1nccx[jj])-
                    (d1nccx[kk]-d1nccx[jj])*
                    (d1nccz[ll]-d1nccz[jj]))/R6;
                nz=((d1nccx[kk]-d1nccx[jj])*
                    (d1nccy[ll]-d1nccy[jj])-
                    (d1nccy[kk]-d1nccy[jj])*
                    (d1nccx[ll]-d1nccx[jj]))/R6;
  //   __builtin_prefetch(&T[0][0],0,2);
                tx=(T[0][0]*nx+T[0][1]*ny+T[0][2]*nz);
                ty=(T[1][0]*nx+T[1][1]*ny+T[1][2]*nz);
                tz=(T[2][0]*nx+T[2][1]*ny+T[2][2]*nz);
                //          tm=dpero*voli/24;

                //        d1nmct[ii]=d1nmct[ii]+tm;
		if((i==1)||(i==3))
		{
			tx=-tx;ty=-ty;tz=-tz;
		}
			d1elfsx[i]=tx;
			d1elfsy[i]=ty;
			d1elfsz[i]=tz;

}
}
static void Yfd3TET4ELS(  /* small strain elastic 4-noded tetrahedra  */
			  INT    nelem, INT *i1elpr,
			  DBL *d1peks,    DBL *d1pela,    DBL *d1pemu,
			  DBL *d1nccx,  DBL *d1nccy,  DBL *d1nccz,  DBL *d1ncix,  DBL  *d1nciy,
			  DBL *d1nciz,  DBL **d2nft,DBL *d1nvcx,  DBL *d1nvcy,  DBL *d1nvcz,
                          INT **i2elto, DBL *d1elvol, INT *i1remov,
                           DBL *d1m0n,DBL *d1m0,DBL *d1mf,  DBL *d1pedm, DBL *d1eldam, DBL ***d3tst

			  )
{INT j,k,ii;
DBL **d2elfsx,  **d2elfsy, **d2elfsz;
d2elfsx=TalDBL2(nelem+1,4);
d2elfsy=TalDBL2(nelem+1,4);
d2elfsz=TalDBL2(nelem+1,4);
#pragma omp parallel  default(none) shared(nelem,i1elpr,d1peks,d1pela,d1pemu,d1nccx,d1nccy,d1nccz,d1ncix,\
d1nciy,d1nciz,d2nft,d1nvcx,d1nvcy,d1nvcz,i2elto,d1elvol,i1remov,d1m0n,d1m0,d1mf,d1eldam,d1pedm,d3tst,\
d2elfsx, d2elfsy, d2elfsz)
    {
        DBL  F[3][3]; /* deformation gradient in global base delta ux/delta x */
        DBL  L[3][3]; /* velocity gradient in global base  delta vx/delta x    */
        DBL  T[3][3]; /* Cauchy stress */
        DBL  D[3][3]; /* rate of deformation*/
        DBL  detf;
        INT *i1elto;
        INT iprop,i;
        INT m,n;


#pragma omp for schedule(guided) nowait
//#endif
  for(i=0;i<nelem;i++)
  {
      if(i1remov[i]!=1)
      {
      iprop=i1elpr[i];
      i1elto=i2elto[i];

      detf=FdftensTet4(F,L,d1nccx,d1nccy,d1nccz,d1ncix,d1nciy,d1nciz,d1nvcx,d1nvcy,d1nvcz,i1elto,&d1elvol[i]);
      CauchyTet4(F,L,d1pemu[iprop],d1peks[iprop],d1pela[iprop],detf,T,D,&d1m0n[i],d1m0[iprop],d1mf[iprop],d1pedm[iprop],&d1eldam[i]);
//      CauchyTet4_pore_pressure(F,L,d1pemu[iprop],d1peks[iprop],d1pela[iprop],detf,T,D);
for(m=0;m<3;m++)
      {
        for(n=0;n<3;n++)
        {
            d3tst[i][m][n]=D[m][n];
        }
      }
      nforceTet4(T,d1nccx,d1nccy,d1nccz,d2nft[0],d2nft[1],d2nft[2],i1elto,d2elfsx[i], d2elfsy[i], d2elfsz[i]);
      }
    }
        
 }

for(j=0;j<nelem;j++)
{
for(k=0;k<4;k++)      /* Nodal Forces */
{
 ii=i2elto[j][k];
 d2nft[0][ii]+=d2elfsx[j][k];
 d2nft[1][ii]+=d2elfsy[j][k];
 d2nft[2][ii]+=d2elfsz[j][k];
}
}
FREE(d2elfsz);
FREE(d2elfsy);
FREE(d2elfsx);
}

static INLINE void S_N_direction(DBL o[3],DBL s[3],DBL normal[3], DBL sx[3], DBL sy[3],DBL sz[3],DBL *area,DBL *el,
                                 INT in[6],
                                 DBL *d1nccx,DBL *d1nccy, DBL *d1nccz,INT ieljo)
{

    DBL prjx[3],prjy[3],prjz[3];   /*projection of first 3 joint element nodes on middle plane*/
    DBL mpx[6],mpy[6],mpz[6];   /*middle points of the 6 edges: 0-01, 1-12, 2-20, 3-54, 4-43, 5-35*/
    DBL midx[3],midy[3],midz[3];   /*coordinates of middle points; 0-05, 1-14, 2-23*/
    DBL midy10,midz20,midz10,midy20,midx20,midx10,midx12,midy12,midz12;
    DBL tmp,tmp1,tmp2,tmp3,mpx03,mpy03,mpz03,mpx14,mpy14,mpz14,mpx25,mpy25,mpz25;
    DBL small,h;
     small=EPSILON;
midx[0]=RP5*(d1nccx[in[0]]+d1nccx[in[5]]);
midy[0]=RP5*(d1nccy[in[0]]+d1nccy[in[5]]);
midz[0]=RP5*(d1nccz[in[0]]+d1nccz[in[5]]);
midx[1]=RP5*(d1nccx[in[1]]+d1nccx[in[4]]);
midy[1]=RP5*(d1nccy[in[1]]+d1nccy[in[4]]);
midz[1]=RP5*(d1nccz[in[1]]+d1nccz[in[4]]);
midx[2]=RP5*(d1nccx[in[2]]+d1nccx[in[3]]);
midy[2]=RP5*(d1nccy[in[2]]+d1nccy[in[3]]);
midz[2]=RP5*(d1nccz[in[2]]+d1nccz[in[3]]);

mpx[0]=RP5*(d1nccx[in[0]]+d1nccx[in[1]]);
mpy[0]=RP5*(d1nccy[in[0]]+d1nccy[in[1]]);
mpz[0]=RP5*(d1nccz[in[0]]+d1nccz[in[1]]);
mpx[1]=RP5*(d1nccx[in[1]]+d1nccx[in[2]]);
mpy[1]=RP5*(d1nccy[in[1]]+d1nccy[in[2]]);
mpz[1]=RP5*(d1nccz[in[1]]+d1nccz[in[2]]);
mpx[2]=RP5*(d1nccx[in[2]]+d1nccx[in[0]]);
mpy[2]=RP5*(d1nccy[in[2]]+d1nccy[in[0]]);
mpz[2]=RP5*(d1nccz[in[2]]+d1nccz[in[0]]);
mpx[3]=RP5*(d1nccx[in[5]]+d1nccx[in[4]]);
mpy[3]=RP5*(d1nccy[in[5]]+d1nccy[in[4]]);
mpz[3]=RP5*(d1nccz[in[5]]+d1nccz[in[4]]);
mpx[4]=RP5*(d1nccx[in[4]]+d1nccx[in[3]]);
mpy[4]=RP5*(d1nccy[in[4]]+d1nccy[in[3]]);
mpz[4]=RP5*(d1nccz[in[4]]+d1nccz[in[3]]);
mpx[5]=RP5*(d1nccx[in[3]]+d1nccx[in[5]]);
mpy[5]=RP5*(d1nccy[in[3]]+d1nccy[in[5]]);
mpz[5]=RP5*(d1nccz[in[3]]+d1nccz[in[5]]);

/* 2*area of this triangle */
midy10=midy[1]-midy[0];
midz20=midz[2]-midz[0];
midz10=midz[1]-midz[0];
midy20=midy[2]-midy[0];
midx20=midx[2]-midx[0];
midx10=midx[1]-midx[0];
midx12=midx[1]-midx[2];
midy12=midy[1]-midy[2];
midz12=midz[1]-midz[2];

tmp=midy10*midz20-midz10*midy20;
tmp1=tmp*tmp;
tmp=midz10*midx20-midx10*midz20;
tmp2=tmp*tmp;
tmp=midx10*midy20-midy10*midx20;
tmp3=tmp*tmp;

h=SQRT(tmp1+tmp2+tmp3);

/* in numerical integration, the area used in the calculation of node force is 1/6 of the triangle area */
*area=h/12.0;

/* nx,ny,nz can be positive or negative */
if((ieljo==0)||(ieljo==2))
{
    
    normal[0]=(midy10*midz20-midz10*midy20)/(h+small);
    normal[1]=(midz10*midx20-midx10*midz20)/(h+small);
    normal[2]=(midx10*midy20-midy10*midx20)/(h+small);
}
else
{
    
    normal[0]=(midz10*midy20-midy10*midz20)/(h+small);
    normal[1]=(midx10*midz20-midz10*midx20)/(h+small);
    normal[2]=(midy10*midx20-midx10*midy20)/(h+small);
}


*el=(SQRT(midx10*midx10+midy10*midy10
         +midz10*midz10)
    +SQRT(midx12*midx12+midy12*midy12
          +midz12*midz12)
    +SQRT(midx20*midx20+midy20*midy20
          +midz20*midz20))/R3;

mpx03=mpx[0]-mpx[3];
mpy03=mpy[0]-mpy[3];
mpz03=mpz[0]-mpz[3];
mpx14=mpx[1]-mpx[4];
mpy14=mpy[1]-mpy[4];
mpz14=mpz[1]-mpz[4];
mpx25=mpx[2]-mpx[5];
mpy25=mpy[2]-mpy[5];
mpz25=mpz[2]-mpz[5];

o[0]=normal[0]*mpx03+normal[1]*mpy03+normal[2]*mpz03;
o[1]=normal[0]*mpx14+normal[1]*mpy14+normal[2]*mpz14;
o[2]=normal[0]*mpx25+normal[1]*mpy25+normal[2]*mpz25;

prjx[0]=mpx03-o[0]*normal[0];
prjy[0]=mpy03-o[0]*normal[1];
prjz[0]=mpz03-o[0]*normal[2];
prjx[1]=mpx14-o[1]*normal[0];
prjy[1]=mpy14-o[1]*normal[1];
prjz[1]=mpz14-o[1]*normal[2];
prjx[2]=mpx25-o[2]*normal[0];
prjy[2]=mpy25-o[2]*normal[1];
prjz[2]=mpz25-o[2]*normal[2];

/* sx[i],sy[i],sz[i] can be positive or negative */
tmp=SQRT(prjx[0]*prjx[0]+prjy[0]*prjy[0]+prjz[0]*prjz[0]);
if(DABS(prjx[0])<EPSILON)
{
    sx[0]=R0;
}
else
{
    sx[0]=prjx[0]/tmp;
}
if(DABS(prjy[0])<EPSILON)
{
    sy[0]=R0;
}
else
{
    sy[0]=prjy[0]/tmp;
}
if(DABS(prjz[0])<EPSILON)
{
    sz[0]=R0;
}
else
{
    sz[0]=prjz[0]/tmp;
}
tmp=SQRT(prjx[1]*prjx[1]+prjy[1]*prjy[1]+prjz[1]*prjz[1]);
if(DABS(prjx[1])<EPSILON)
{
    sx[1]=R0;
}
else
{
    sx[1]=prjx[1]/tmp;
}
if(DABS(prjy[1])<EPSILON)
{
    sy[1]=R0;
}
else
{
    sy[1]=prjy[1]/tmp;
}
if(DABS(prjz[1])<EPSILON)
{
    sz[1]=R0;
}
else
{
    sz[1]=prjz[1]/tmp;
}
tmp=SQRT(prjx[2]*prjx[2]+prjy[2]*prjy[2]+prjz[2]*prjz[2]);

if(DABS(prjx[2])<EPSILON)
{
    sx[2]=R0;
}
else
{
    sx[2]=prjx[2]/tmp;
}
if(DABS(prjy[2])<EPSILON)
{
    sy[2]=R0;
}
else
{
    sy[2]=prjy[2]/tmp;
}
if(DABS(prjz[2])<EPSILON)
{
    sz[2]=R0;
}
else
{
    sz[2]=prjz[2]/tmp;
}


s[0]=sx[0]*mpx03+sy[0]*mpy03+sz[0]*mpz03;
s[1]=sx[1]*mpx14+sy[1]*mpy14+sz[1]*mpz14;
s[2]=sx[2]*mpx25+sy[2]*mpy25+sz[2]*mpz25;

}




static INLINE void Sigma_tau(DBL o[3],DBL s[3], DBL sx[3], DBL sy[3], DBL sz[3], DBL normal[3],
                             DBL area, DBL el, DBL *d1nftx,DBL *d1nfty,DBL *d1nftz,
                                   DBL dpeft, DBL  dpegfn, DBL dpegfs,DBL dpepe,
                             DBL dpcoh, DBL dpicf, INT **i2eljp, INT **i2ejjp, INT *i1elpr, INT *i1elty,
                             INT ielem, INT *i1elbe, int tid, INT *index, INT in[6],INT iprop,
                             INT *nelem_t, INT **i2elto_t, INT *i1bd,DBL *d1nbp,DBL dbkp,
                             DBL *d1ejfcx,DBL *d1ejfcy,DBL *d1ejfcz,INT *i1eljt, INT *i1elncstep, INT ncstep,
			     DBL *d1ejfjx,DBL *d1ejfjy,DBL *d1ejfjz)
{
    DBL odis,sdis,sabs,op,sp,ot,st,z,sigma,tau,sigma_tmp,tmp1,tmp2,tmp,tmp3;
    DBL dpefa=0.63;
    DBL dpefb=1.8;
    DBL dpefc=6.0;
    DBL dpefm=0.0;

    DBL dpefs;
    int integ,nfail,nfail_n,nfail_s,nnfail,iz,i;
DBL d2nforce[3][6]={{0.0}};
INT nelem_t_c;

op=R2*el*dpeft/dpepe;
ot=MAXIM((R2*op),(R3*dpegfn/dpeft));   /*need further investigation*/

nfail=0;
nfail_n=0;
nfail_s=0;
nnfail=0;
//numerical integration: three integration points
for(integ=0;integ<3;integ++)
{
    
    /* Mohr-Coulomb failure criterion with the tension cut-off */
    odis=o[integ];
    sigma_tmp=EPSILON;
    if(odis<R0)
    {
        sigma_tmp=R2*odis*dpeft/op;
    }
    if(sigma_tmp>R0)
    {
        dpefs=dpcoh;
    }
    else
    {
        dpefs=dpcoh-dpicf*sigma_tmp;
    }

    sp=R2*el*dpefs/dpepe;
    st=MAXIM((R2*sp),(R3*dpegfs/dpefs));   /*need further investigation*/

    odis=o[integ];
    sdis=s[integ];
    sabs=DABS(sdis);
    tmp1=(odis-op)/ot;
    tmp2=(sabs-sp)/st;
    iz=1;
    if((odis>op)&&(sabs>sp))
    {
        z=SQRT(tmp1*tmp1+tmp2*tmp2);
        if(z>=R1)
        {
            nfail_n=nfail_n+1;
            nfail_s=nfail_s+1;
            nnfail=nnfail+1;

        }

    }
    else if(odis>op)
    {
        z=tmp1;
        if(z>=R1)
        {
            nfail_n=nfail_n+1;
            //                         nfail_s=nfail_s+1;
        }

    }
    else if(sabs>sp)
    {
        z=tmp2;
        if(z>=R1)
        {
            //                         nfail_n=nfail_n+1;
            nfail_s=nfail_s+1;
        }

    }
    else
    {
        z=R0;
        iz=0;
    }

    if(z>=R1)
    {
        nfail=nfail+1;
        if((nfail>1)&&(i1elpr[ielem]>=0))
        {
            i1elpr[ielem]=iprop-YIPROPMAX;
            if(nnfail>1)
            {
                i1elty[ielem]=5;
            }
            else if ((nfail_s>1)&&(nfail_n<1))
            {
                i1elty[ielem]=4; //shear
            }
            else if ((nfail_s<1)&&(nfail_n>1))
            {
                i1elty[ielem]=3; //normal
            }
            else
            {
                i1elty[ielem]=5; //shear-rotation
            }

            //-PY changed it for 3D_fracture_coupling_with_multiphase -ao commented this again
            

            for(i=0;i<4;i++)
            {
                if(i2eljp[i2ejjp[ielem][0]][i]==ielem)
                {
                    i2eljp[i2ejjp[ielem][0]][i]=-(i2eljp[i2ejjp[ielem][0]][i]+YIPROPMAX);
                }
                if(i2eljp[i2ejjp[ielem][1]][i]==ielem)
                {
                    i2eljp[i2ejjp[ielem][1]][i]=-(i2eljp[i2ejjp[ielem][1]][i]+YIPROPMAX);
                }
            }
#pragma omp critical
{
    i1eljt[i2ejjp[ielem][0]]=1;
    i1eljt[i2ejjp[ielem][1]]=1;
    i1elncstep[i2ejjp[ielem][0]]=ncstep;
    i1elncstep[i2ejjp[ielem][1]]=ncstep;

    nelem_t_c=*nelem_t;
    *nelem_t=*nelem_t+2;
}
    i1bd[nelem_t_c]=0;
    d1nbp[nelem_t_c]=dbkp;
    i2elto_t[nelem_t_c][0]=in[0];
    i2elto_t[nelem_t_c][1]=in[1];
    i2elto_t[nelem_t_c][2]=in[2];
    nelem_t_c=nelem_t_c+1;

    i1bd[nelem_t_c]=0;
    d1nbp[nelem_t_c]=dbkp;
    i2elto_t[nelem_t_c][0]=in[3];
    i2elto_t[nelem_t_c][1]=in[4];
    i2elto_t[nelem_t_c][2]=in[5];

//#pragma omp critical
{
//            if(i1elbe[i2ejjp[ielem][0]]<YIPROPMAX)
            {
            #pragma omp atomic
            i1elbe[i2ejjp[ielem][0]]+=YIPROPMAX;
            }
//            if(i1elbe[i2ejjp[ielem][1]]<YIPROPMAX)
            {
            #pragma omp atomic
            i1elbe[i2ejjp[ielem][1]]+=YIPROPMAX;
            }
}
            
            index[tid]=1;



        }
        z=R1;
        iz=2;
    }

    if(iz==0)z=R1;
        else if(iz==2)z=R0;
            else
                z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*
                   exp(z*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))
                *(dpefa*(R1-z)+dpefb*pow((R1-z),dpefc));

                tmp=odis/op;
                if(odis<R0)
                {
                    sigma=R2*tmp*dpeft;
                }
                else if(odis>op)
                {
                    sigma=dpeft*z;
                    //             nsoft=nsoft+1;
                }
                else
                {
                    sigma=(R2*tmp-tmp*tmp)*z*dpeft;
                }


    tmp=sabs/sp;
    if((sigma>R0)&&(sabs>sp))
    {
        tau=z*dpefs;
    }
    else if(sigma>R0)
    {
        tau=(R2*tmp-tmp*tmp)*z*dpefs;
    }
    else if(sabs>sp)
    {
        tau=z*dpefs-dpefm*sigma;
    }
    else
    {
        tau=(R2*tmp-tmp*tmp)*(z*dpefs-dpefm*sigma); //-ao getting NaNs sometimes
    }

    if(sdis<R0)
    {
        tau=-tau;
    }


    if(integ==0)
    {
        tmp1=area*(tau*sx[0]+sigma*normal[0]);
        tmp2=area*(tau*sy[0]+sigma*normal[1]);
        tmp3=area*(tau*sz[0]+sigma*normal[2]);

        d2nforce[0][0]=d2nforce[0][0]-tmp1;
        d2nforce[0][5]=d2nforce[0][5]+tmp1;
        d2nforce[1][0]=d2nforce[1][0]-tmp2;
        d2nforce[1][5]=d2nforce[1][5]+tmp2;
        d2nforce[2][0]=d2nforce[2][0]-tmp3;
        d2nforce[2][5]=d2nforce[2][5]+tmp3;
        d2nforce[0][1]=d2nforce[0][1]-tmp1;
        d2nforce[0][4]=d2nforce[0][4]+tmp1;
        d2nforce[1][1]=d2nforce[1][1]-tmp2;
        d2nforce[1][4]=d2nforce[1][4]+tmp2;
        d2nforce[2][1]=d2nforce[2][1]-tmp3;
        d2nforce[2][4]=d2nforce[2][4]+tmp3;
    
    }
    else if(integ==1)
    {
        tmp1=area*(tau*sx[1]+sigma*normal[0]);
        tmp2=area*(tau*sy[1]+sigma*normal[1]);
        tmp3=area*(tau*sz[1]+sigma*normal[2]);
        d2nforce[0][1]=d2nforce[0][1]-tmp1;
        d2nforce[0][4]=d2nforce[0][4]+tmp1;
        d2nforce[1][1]=d2nforce[1][1]-tmp2;
        d2nforce[1][4]=d2nforce[1][4]+tmp2;
        d2nforce[2][1]=d2nforce[2][1]-tmp3;
        d2nforce[2][4]=d2nforce[2][4]+tmp3;
        d2nforce[0][2]=d2nforce[0][2]-tmp1;
        d2nforce[0][3]=d2nforce[0][3]+tmp1;
        d2nforce[1][2]=d2nforce[1][2]-tmp2;
        d2nforce[1][3]=d2nforce[1][3]+tmp2;
        d2nforce[2][2]=d2nforce[2][2]-tmp3;
        d2nforce[2][3]=d2nforce[2][3]+tmp3;

    }
    else
    {
        tmp1=area*(tau*sx[2]+sigma*normal[0]);
        tmp2=area*(tau*sy[2]+sigma*normal[1]);
        tmp3=area*(tau*sz[2]+sigma*normal[2]);
        d2nforce[0][2]=d2nforce[0][2]-tmp1;
        d2nforce[0][3]=d2nforce[0][3]+tmp1;
        d2nforce[1][2]=d2nforce[1][2]-tmp2;
        d2nforce[1][3]=d2nforce[1][3]+tmp2;
        d2nforce[2][2]=d2nforce[2][2]-tmp3;
        d2nforce[2][3]=d2nforce[2][3]+tmp3;
        d2nforce[0][0]=d2nforce[0][0]-tmp1;
        d2nforce[0][5]=d2nforce[0][5]+tmp1;
        d2nforce[1][0]=d2nforce[1][0]-tmp2;
        d2nforce[1][5]=d2nforce[1][5]+tmp2;
        d2nforce[2][0]=d2nforce[2][0]-tmp3;
        d2nforce[2][5]=d2nforce[2][5]+tmp3;

    }

}
    if(i1elpr[ielem]<0)
    {
          d1ejfcx[ielem]=d2nforce[0][3]+d2nforce[0][4]+d2nforce[0][5];
          d1ejfcy[ielem]=d2nforce[1][3]+d2nforce[1][4]+d2nforce[1][5];
          d1ejfcz[ielem]=d2nforce[2][3]+d2nforce[2][4]+d2nforce[2][5];

    }


    for(i=0;i<6;i++)
       {
	       d1ejfjx[i]=d2nforce[0][i];
	       d1ejfjy[i]=d2nforce[1][i];
	       d1ejfjz[i]=d2nforce[2][i];
	   
         }


}

static void Yfd3TET4JOINT(		/* joint element for 4-node tetrahedra */
     INT nelem,
     INT *i1ptyp,
      DBL *dpefti, DBL *dpegfni, DBL *dpegfsi,DBL *dpepei,
     DBL *d1nccx,DBL *d1nccy,DBL *d1nccz,
     DBL *d1nftx,DBL *d1nfty,DBL *d1nftz,
     DBL *d1nvcx,DBL *d1nvcy,DBL *d1nvcz,
     INT *i1elpr,INT **i2elto,INT *i1eljo,
     INT **i2eljp, INT **i2ejjp,
     DBL *d0iedi,
     INT **i2elfr, INT *i1elbe,INT *i1joint,INT njoint,INT *ifrac,
     DBL *dpcohi, DBL *dpicf, INT *i1elty, INT *i1elst,
    INT *nelem_t, INT **i2elto_t, INT *i1bd, DBL *d1nbp, DBL dbkp,
    DBL *d1ejfcx,DBL *d1ejfcy,DBL *d1ejfcz,INT *i1defect,int nthreads, INT *i1eljt,DBL *d1ejmr,
                          INT *i1elncstep, INT ncstep,DBL ***d3tst, DBL **d2jst
     )

{
    INT *index,i,*i1elpr_old,ii;
    DBL **d2ejfjx, **d2ejfjy, **d2ejfjz;
    index=TalINT1(nthreads+1);
    TzINT1(index,nthreads+1);
    i1elpr_old=TalINT1(nelem+1);
    d2ejfjx=TalDBL2(nelem+1,6);
    d2ejfjy=TalDBL2(nelem+1,6);
    d2ejfjz=TalDBL2(nelem+1,6);

#pragma omp parallel  default(none) shared(nelem,dpefti,dpegfni,dpegfsi,dpepei,d1nccx,\
d1nccy,d1nccz,d1nftx,d1nfty,d1nftz,d1nvcx,d1nvcy,d1nvcz,i1elpr,i2elto,i1eljo,i2eljp,i2ejjp,\
d0iedi,i2elfr,i1elbe,i1ptyp,i1elst,d1ejfcx,d1ejfcy,d1ejfcz,i1elncstep,ncstep,d1ejmr,\
dpcohi,dpicf,i1elty,index,i1joint,njoint,nelem_t,i2elto_t,i1bd,d1nbp,dbkp,i1defect,i1eljt,d3tst,d2jst,stdout,\
d2ejfjx,d2ejfjy,d2ejfjz,i1elpr_old)
    {
     DBL area,el,normal[3];
     INT ielem,in[6];
     INT *i1elto;
    DBL o[3],s[3], sx[3], sy[3], sz[3];
     INT j,jj,iprop;
     DBL dpegfn,dpegfs,dpepe,dpeft,dpcoh;
     DBL fact;
     int tid;
  DBL dpeftdif, dpcohdif, dpegfndif, dpegfsdif;

	 tid = omp_get_thread_num();
      index[tid]=0;
 //   #pragma omp for schedule(dynamic,64) nowait
    #pragma omp for schedule(guided) nowait

        for(ielem=0;ielem<nelem;ielem++)
	  {

          iprop=i1elpr[ielem];
	  i1elpr_old[ielem]=iprop;

if(iprop>=0)
{
for(jj=0;jj<6;jj++)
           {
	       in[jj]=i2elto[ielem][jj];
           }
}
          if((iprop>=0)&&(i2elto[ielem][0]!=i2elto[ielem][5]))
          {
       S_N_direction(o,s,normal,sx,sy,sz,&area,&el,in,
		       d1nccx,d1nccy, d1nccz,i1eljo[ielem]);
       
          //Strain rate effect of Tensile strength
	        dpeftdif=R1;
            dpeft=dpeftdif*dpefti[iprop];
            dpcohdif=R1;
//Strain rate effect of cohesive

            dpcoh=dpcohdif*dpcohi[iprop];
            dpegfn=dpegfni[iprop]*dpeftdif;
            dpegfs=dpegfsi[iprop]*dpcohdif;

	  /* fact=0.1 for Sadjad's BD */
              if(i1defect[ielem]==1)
              {
		      fact=0.1;
                  dpeft=dpeft*fact;
                  dpcoh=dpcoh*fact;
                  dpepe=dpepei[iprop];///d1ejmr[ielem];
                  dpegfn=dpegfn*fact;
                  dpegfs=dpegfs*fact;
              }
              else
              {
           dpeft=dpeft;
           dpcoh=dpcoh;
           dpepe=dpepei[iprop];// /d1ejmr[ielem];
           dpegfn=dpegfn;
           dpegfs=dpegfs;
              }

      //    }
         for(jj=0;jj<6;jj++)
           {
	       in[jj]=i2elto[ielem][jj];
           }
                 Sigma_tau(o,s,sx,sy,sz, normal,
                           area, el, d1nftx,d1nfty,d1nftz,
                           dpeft, dpegfn, dpegfs,dpepe,
                           dpcoh, dpicf[iprop],i2eljp,i2ejjp, i1elpr,i1elty,
                           ielem, i1elbe, tid, index,in,iprop,
                           nelem_t,i2elto_t,i1bd,d1nbp,dbkp,
                           d1ejfcx,d1ejfcy,d1ejfcz,i1eljt,i1elncstep,ncstep,
			   d2ejfjx[ielem],d2ejfjy[ielem],d2ejfjz[ielem]);

              if(i1elpr[ielem]<0)i1elst[ielem]=0;
	  }
     }
    }

for(i=0;i<nelem;i++)
{
          if((i1elpr_old[i]>=0)&&(i2elto[i][0]!=i2elto[i][5]))
         {

for(ii=0;ii<6;ii++)
{
	d1nftx[i2elto[i][ii]]+=d2ejfjx[i][ii];
	d1nfty[i2elto[i][ii]]+=d2ejfjy[i][ii];
	d1nftz[i2elto[i][ii]]+=d2ejfjz[i][ii];
}
	 }
}
    for(i=0;i<nthreads;i++)
    {
        if(index[i])
        {
            (*d0iedi)=2000.0;
            *ifrac=1;
            break;
        }
    }
    FREE(d2ejfjz);
    FREE(d2ejfjy);
    FREE(d2ejfjx);
    FREE(i1elpr_old);
    FREE(index);

}

/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief calculate nodal forces
 *  \param[in]     yde element database
 *  \param[in,out] ydn nodal database
 *  \param[in]     ydp property database
 *  \param[in]     ydk constants database
 *  \par Details:
 *  Yfd() ... concepts described in chapter 4, I think
 *
    \verbatim
    for each tetraheron
      for eash Gauss point
        calculate Jacobian matrix
        calculate Cartesian gradients
        calculate tensor
        calculate velocity gradients
        calculate Cauchy stress
        calculate some unnamed value for this iteration
      calculate the nodal forces
    \endverbatim
 *
 *
 */
void Yfd(YDE yde, YDN ydn, YDP ydp, YDK ydk,YDI ydi,YDC ydc,YDX ydx,YDJ ydj, YDJP ydjp,
		YDB ydb,int nthreads)
{
    INT i,j,ielem,nblock,iprop,nmode,inopo,index,indexv;

    nblock=(INT)(ydn->nnopo/nthreads);
    nmode=ydn->nnopo-nblock*nthreads;



    for(i=0;i<ydn->nnopo;i++)
    {
        for(j=0;j<3;j++)
        {
            ydn->d2nfc[j][i]=R0;
            ydn->d2nft[j][i]=R0;

        }
    }
    for(i=0;i<yde->nelem;i++)
    {yde->i1eljt[i]=-1;
    }

    if(ydp->i1ptyp[0] == YTE3TET10ELS) /* small strain elastic tetrahetra  */
    {
      Yfd3TET10ELS(
        yde->nelemi,
        ydp->d1peks, ydp->d1pela, ydp->d1pemu,
        ydn->d2ncc, ydn->d2nci, ydn->d2nft, ydn->d2nvc,
        yde->i1elpr, yde->i2elto,
        ydk->d3dsh,
        nthreads,nblock,ydn->nnopo,nmode
        );

    }
  else if((ydp->i1ptyp[0])==(YTE3TET4ELS))
	  {
		Yfd3TET4ELS(  // small strain elastic tetrahetra
		yde->nelemi, yde->i1elpr,
		ydp->d1peks,ydp->d1pela,ydp->d1pemu,
		ydn->d2ncc[0],ydn->d2ncc[1],
		ydn->d2ncc[2],ydn->d2nci[0],ydn->d2nci[1],
		ydn->d2nci[2],ydn->d2nft,ydn->d2nvc[0],
		ydn->d2nvc[1],ydn->d2nvc[2],
		yde->i2elto,yde->d1elvol,yde->i1elchprm,
                yde->d1pem0n,ydp->d1pem0,ydp->d1pemf,ydp->d1pedm,
		yde->d1eldam,yde->d3tst);
	}

	       Yfd3TET4JOINT(	// joint element for 4-node tetrahedra
		    ydj->nelem,
                    ydp->i1ptyp,
		    ydjp->d1peft,
		    ydjp->d1pegfn,ydjp->d1pegfs,ydjp->d1pepe,
		    ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
		    ydn->d2nft[0],ydn->d2nft[1],ydn->d2nft[2],
		    ydn->d2nvc[0],ydn->d2nvc[1],ydn->d2nvc[2],
		    ydj->i1elpr,ydj->i2elto,ydj->i1eljo,
		    yde->i2eljp, ydj->i2ejjp,
            &(ydi->diedi),
            yde->i2elfr,yde->i1elbe,yde->i1joint,yde->njoint_intact,&(yde->ifrac),
		    ydjp->d1pcoh, ydjp->d1picf,ydj->i1elty,ydj->i1elst,
            &(ydx->nelem_t),ydx->i2elto_t,ydx->i1jtp,ydx->d1nbp,ydx->dbkp,
            ydj->d2ejfc[0],ydj->d2ejfc[1],ydj->d2ejfc[2],
			 ydj->i1defect, nthreads,yde->i1eljt,ydj->d1ejmr,
                         yde->i1elncstep,ydc->ncstep,yde->d3tst,ydj->d2jst);

    if(yde->ifrac)
    {
    for(ielem=0;ielem<ydj->nelem;ielem++)
         {
             if(ydj->i1elst[ielem]==0)
   {
       ydj->i1elst[ielem]=1;
    for(i=0;i<6;i++)
    {
        inopo=ydj->i2elto[ielem][i];
        if(ydn->i1nobf[inopo]<=0)
	{
        ydn->i1nobfse[ydn->nnobf]=inopo;
        ydn->nnobf++;
	}
        ydn->i1nobf[inopo]+=YIPROPMAX;

    }
   }
    }
        
    yde->ifrac=0;
    ydc->ifrag=1;
    }

       if(((ydc->ncstep%ydc->icoutf)==0&&(ydc->ifrag==1))||
	    ((ydc->ncstep%ydc->icoutf)==0&&(ydc->ncstep>=ydc->mcstep))||
		(ydc->ncstep==ydc->mcstep))
       {
	       index=0;
	       if(ydc->ncstep>=ydc->mcstep)index=1;
               indexv=0;
               if(ydc->ncstep==ydc->mcstep)indexv=1;

	       calFragVol(yde->nelemi,yde->i1elfrg,&(yde->nfrag), yde->d1volfrg,
	       yde->d1emct,yde->i2eljp,yde->i1elpr,ydp->i1pejp,ydj->i2ejjp,
                yde->i1elchprm, ydn->i1nochprm, yde->i2elto, ydn->d2ncc[0],
                ydn->d2ncc[1],ydn->d2ncc[2],ydn->d2nvc[0],
                      ydn->d2nvc[1],ydn->d2nvc[2],yde->nelno,ydc->ichprm,
		      index,indexv,ydc->dcrmvx,ydc->dcrmvy,ydc->dcrmvz,
		      ydn->d2nfrm,yde->i1elfrm,yde->i1frgrm,ydn->d1ndisp,
		      ydc->dcfrmca,ydn->i1nopr,ydb->i1bnvx,ydb->i1bnvy,ydb->i1bnvz,ydn->i1rmgrav,yde->d1mart);
	       ydc->ifrag=0;
       }
       
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/
