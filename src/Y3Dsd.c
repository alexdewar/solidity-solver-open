/*! \file Y3Dsd.c
 *  \brief solver database for Y3D model
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


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


/* calculate current nodal force and change in velocity
*/
static INLINE void nfvel(
  DBL *d1bnfx, DBL *d1bnfy, DBL *d1bnfz,
  DBL *d1nfcx, DBL *d1nfcy, DBL *d1nfcz,
  DBL *d1nftx, DBL *d1nfty, DBL *d1nftz,
  DBL *d1kax,  DBL *d1kay,  DBL *d1kaz,
  INT *i1bnvx, INT *i1bnvy, INT *i1bnvz,
  DPT *dp1dv,
  INT nnopo,
  DBL *d1nmct,
  INT *i1nopr,
  int izero,INT *i1norem)
{
    INT inopo, iprop;
    DBL nmct;
#pragma omp parallel for default(none) shared(d1bnfx,d1bnfy,d1bnfz,d1nfcx,d1nfcy,\
d1nfcz,d1nftx, d1nfty,d1nftz,d1kax,d1kay,d1kaz,i1bnvx,i1bnvy,i1bnvz,\
dp1dv,nnopo,d1nmct,i1nopr,izero,i1norem) private(inopo,iprop,nmct) schedule(guided)
  for(inopo=0; inopo<nnopo; inopo++)
  {  if(i1norem[inopo]==1)continue;
    iprop = i1nopr[inopo];
    nmct  = d1nmct[inopo];
 
    if(i1bnvx[iprop]||DABS(nmct)<EPSILON)
    {
     d1nftx[inopo]=R0;
     dp1dv[inopo][0] =R0;
    }
    else
    {
     d1nftx[inopo] += (d1bnfx[iprop] + d1nfcx[inopo]+d1kax[iprop]*nmct);
     dp1dv[inopo][0] = d1nftx[inopo]/nmct;
    }

    if(i1bnvy[iprop]||DABS(nmct)<EPSILON)
    {
     d1nfty[inopo]=R0;
     dp1dv[inopo][1] =R0;
    }
    else
    {
     d1nfty[inopo] += (d1bnfy[iprop] + d1nfcy[inopo]+d1kay[iprop]*nmct);
     dp1dv[inopo][1] = d1nfty[inopo]/nmct;
    }

    if(i1bnvz[iprop]||DABS(nmct)<EPSILON)
    {
     d1nftz[inopo]=R0;
     dp1dv[inopo][2] =R0;
    }
    else
    {
     d1nftz[inopo] += (d1bnfz[iprop] + d1nfcz[inopo]+d1kaz[iprop]*nmct);
     dp1dv[inopo][2] = d1nftz[inopo]/nmct;
 
    }

      if(izero)
     {
      d1nftx[inopo]=R0;
      d1nfty[inopo]=R0;
      d1nftz[inopo]=R0;
      d1nfcx[inopo]=R0;
      d1nfcy[inopo]=R0;
      d1nfcz[inopo]=R0;
     }
  }
}


/**********************************************************************/
/**********************************************************************/


/* one-pass solver
*/
static INLINE void solve(
  DBL *d1nftx0, DBL *d1nfty0, DBL *d1nftz0,
  DBL *d1nftx1, DBL *d1nfty1, DBL *d1nftz1,
  DBL *d1nftx2, DBL *d1nfty2, DBL *d1nftz2,
  INT *i1bnvx, INT *i1bnvy, INT *i1bnvz,
  DPT *dp1dv,
  DPT *dp1am,
  INT nnopo,
  DBL *d1nmct,
  INT nelem,
  DBL *d1emct,
  INT **i2elto,
  DBL **d2csmm,
  DBL dcurelx,
  INT *i1nopr, INT *i1elrem, INT *i1norem)
{
    /* zero acceleration work array (old method was less cohesive)
     */
    TzDPT1(dp1am, nnopo);
#pragma omp parallel default(none) shared(d1nftx0,d1nfty0,d1nftz0,d1nftx1,d1nfty1,\
d1nftz1,d1nftx2,d1nfty2,d1nftz2,i1bnvx,i1bnvy,i1bnvz,dp1dv,dp1am,nnopo,d1nmct,nelem,\
d1emct,i2elto,d2csmm,dcurelx,i1nopr,i1elrem,i1norem)
    {
  INT ielem, inopo, i, j, ip, jp, *i1elto,iprop;
  DBL emct, mass, *d1csmm;

  /* calculate force*/
#pragma omp for schedule(static,64) //nowait
  for(ielem=0; ielem<nelem; ielem++)
  {
      if(i1elrem[ielem]==1)continue;
    i1elto = i2elto[ielem];
    emct   = d1emct[ielem];

	  for(i=0; i<NELNO; i++)
	  {
      ip     = i1elto[i];
      d1csmm = d2csmm[i];

		  for(j=0; j<NELNO; j++)
		  {
        jp   = i1elto[j];
        mass = d1csmm[j]*emct;

			  dp1am[ip][0] += mass*dp1dv[jp][0];
			  dp1am[ip][1] += mass*dp1dv[jp][1];
			  dp1am[ip][2] += mass*dp1dv[jp][2];
		  }
	  }
  }

  /* calculate velocity*/
    #pragma omp for schedule(guided) nowait
  for(inopo=0; inopo<nnopo; inopo++)
  {
      if(i1norem[inopo]==1)continue;
if(d1nmct[inopo]>EPSILON)
{
 iprop=i1nopr[inopo];
 if(i1bnvx[iprop]) /* supplied velocity     */
 dp1am[inopo][0]=R0;

 if(i1bnvy[iprop]) /* supplied velocity     */
 dp1am[inopo][1]=R0;

 if(i1bnvz[iprop]) /* supplied velocity     */
 dp1am[inopo][2]=R0;

    d1nftx2[inopo] = d1nftx0[inopo] + (d1nftx1[inopo] - dp1am[inopo][0])*dcurelx;
    d1nfty2[inopo] = d1nfty0[inopo] + (d1nfty1[inopo] - dp1am[inopo][1])*dcurelx;
    d1nftz2[inopo] = d1nftz0[inopo] + (d1nftz1[inopo] - dp1am[inopo][2])*dcurelx;

    dp1dv[inopo][0] = d1nftx2[inopo]/d1nmct[inopo];
    dp1dv[inopo][1] = d1nfty2[inopo]/d1nmct[inopo];
    dp1dv[inopo][2] = d1nftz2[inopo]/d1nmct[inopo];
  }
    }
}
}


/**********************************************************************/
/**********************************************************************/


/* calculate updated nodal velocity and coordinates*/
static INLINE void nupdate(
  DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
  DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
  DBL *d1bnvx, DBL *d1bnvy, DBL *d1bnvz,
  INT *i1bnvx, INT *i1bnvy, INT *i1bnvz,
//  DBL *d1kax,  DBL *d1kay,  DBL *d1kaz,
  DPT *dp1dv,
  INT nnopo,
  INT *i1nopr,
  DBL dcstec,INT *i1norem)
{
      INT inopo, iprop;
#pragma omp parallel for default(none) shared(d1nccx,d1nccy,d1nccz,\
d1nvcx,d1nvcy,d1nvcz,d1bnvx,d1bnvy,d1bnvz,i1bnvx,i1bnvy,i1bnvz,\
dp1dv,nnopo,i1nopr,dcstec,i1norem) private(inopo,iprop) schedule(guided)
  for(inopo=0; inopo<nnopo; inopo++)
  {
    if(i1norem[inopo]==1)continue;
    iprop = i1nopr[inopo];

    /* current velocity */
    d1nvcx[inopo] += dp1dv[inopo][0]*dcstec;
    d1nvcy[inopo] += dp1dv[inopo][1]*dcstec;
    d1nvcz[inopo] += dp1dv[inopo][2]*dcstec;

    /* override using supplied velocity (turns out this is quicker than if/else) */
    if(i1bnvx[iprop]) d1nvcx[inopo] = d1bnvx[iprop];
    if(i1bnvy[iprop]) d1nvcy[inopo] = d1bnvy[iprop];
    if(i1bnvz[iprop]) d1nvcz[inopo] = d1bnvz[iprop];

    /* new coordinates */
    d1nccx[inopo] += d1nvcx[inopo]*dcstec;
    d1nccy[inopo] += d1nvcy[inopo]*dcstec;
    d1nccz[inopo] += d1nvcz[inopo]*dcstec;
  }
}


/**********************************************************************/
/**********************************************************************/


/* mechanical solver for 3D nodes with x,y,z d.o.f.*/
static void Ysd3MEC_multi_pass(
					  DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
            DBL *d1nfcx, DBL *d1nfcy, DBL *d1nfcz,
					  DBL *d1nftx, DBL *d1nfty, DBL *d1nftz,
            DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
            INT *i1bnvx, INT *i1bnvy, INT *i1bnvz,
					  DBL *d1bnvx, DBL *d1bnvy, DBL *d1bnvz,
            DBL *d1bnfx, DBL *d1bnfy, DBL *d1bnfz,
					  DBL *d1kax,  DBL *d1kay,  DBL *d1kaz,
            DBL **d2csmm,
            INT nnopo,
            DBL *d1nmct,
            INT *i1nopr,
            INT *i1nobf,
            DBL *d1emct,
            INT nelem,
            INT **i2elto,
            INT initer,
            DBL dcstec,
            DBL dcurelx,int izero,
            INT *i1elrem, INT *i1norem
					  )
{
  DBL *d1nftx1, *d1nfty1, *d1nftz1;
  DPT *dp1dv, *dp1am;

  /* allocate temporary working arrays (could be static in YDK) */
  d1nftx1 = TalDBL1(nnopo);
  d1nfty1 = TalDBL1(nnopo);
  d1nftz1 = TalDBL1(nnopo);
  dp1dv   = TalDPT1(nnopo);
  dp1am   = TalDPT1(nnopo);

  /* calculate current nodal force nft, and change in velocity dv */
  nfvel(d1bnfx, d1bnfy, d1bnfz,
        d1nfcx, d1nfcy, d1nfcz,
        d1nftx, d1nfty, d1nftz,
        d1kax,  d1kay,  d1kaz,
        i1bnvx, i1bnvy, i1bnvz,
        dp1dv,
        nnopo,
        d1nmct,
        i1nopr,izero, i1norem);

  /* solver pass 1 (or could have simpler solve1() and solve2()) */
if(initer>0)
{
  solve(d1nftx,  d1nfty,  d1nftz,
        d1nftx,  d1nfty,  d1nftz,
        d1nftx1, d1nfty1, d1nftz1,
        i1bnvx, i1bnvy, i1bnvz,
        dp1dv,
        dp1am,
        nnopo,
        d1nmct,
        nelem,
        d1emct,
        i2elto,
        d2csmm,
        dcurelx,
        i1nopr, i1elrem,i1norem);

  /* solver pass 2 (this could be compile-time constant: #if(INITER==2) */
  if(initer == 2)
  {
    solve(d1nftx,  d1nfty,  d1nftz,
          d1nftx1, d1nfty1, d1nftz1,
          d1nftx1, d1nfty1, d1nftz1,
          i1bnvx, i1bnvy, i1bnvz,
          dp1dv,
          dp1am,
          nnopo,
          d1nmct,
          nelem,
          d1emct,
          i2elto,
          d2csmm,
          dcurelx,
          i1nopr, i1elrem,i1norem);
  }
}
  /* update nodal velocity nvc and coordinates ncc */
  nupdate(d1nccx, d1nccy, d1nccz,
          d1nvcx, d1nvcy, d1nvcz,
          d1bnvx, d1bnvy, d1bnvz,
          i1bnvx, i1bnvy, i1bnvz,
          dp1dv,
          nnopo,
          i1nopr,
          dcstec,i1norem);

  /* release temporary working arrays (could be static in YDK) */
  FREE(dp1am);
  FREE(dp1dv);
  FREE(d1nftz1);
  FREE(d1nfty1);
  FREE(d1nftx1);
}

/**********************************************************************/
/**********************************************************************/


/* mechanical solver for 3D nodes with x,y,z d.o.f.*/
static void Ysd3MEC(
            DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
            DBL *d1nfcx, DBL *d1nfcy, DBL *d1nfcz,
            DBL *d1nftx, DBL *d1nfty, DBL *d1nftz,
            DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
            INT *i1bnvx, INT *i1bnvy, INT *i1bnvz,
            DBL *d1bnvx, DBL *d1bnvy, DBL *d1bnvz,
            DBL *d1bnfx, DBL *d1bnfy, DBL *d1bnfz,
            DBL *d1kax,  DBL *d1kay,  DBL *d1kaz,
            INT nnopo,
            DBL *d1nmct,
            INT *i1nopr,
            INT *i1nobf,
            INT nelem,
            INT **i2elto,
            DBL dcstec,
            int izero, INT *i1norem,
	    INT *i1rmgrav,DBL dcrmgravx,DBL dcrmgravy,DBL dcrmgravz,
	    int index
                                          )
{
  /* allocate temporary working arrays (could be static in YDK) */

/* calculate current nodal force and change in velocity*/
    INT inopo, iprop;
    DBL nmct,dp1dv[3],dbeta;
    dbeta=0.0005;

#pragma omp parallel for default(none) shared(d1bnfx,d1bnfy,d1bnfz,d1nfcx,d1nfcy,\
d1nfcz,d1nftx, d1nfty,d1nftz,d1kax,d1kay,d1kaz,i1bnvx,i1bnvy,i1bnvz,\
nnopo,d1nmct,i1nopr,izero,dcstec,d1nccx,d1nccy,d1nccz,\
d1nvcx,d1nvcy,d1nvcz,d1bnvx,d1bnvy,d1bnvz, i1norem,i1rmgrav, dcrmgravx, dcrmgravy,\
dcrmgravz,index) private(inopo,iprop,nmct,dp1dv,dbeta)  schedule(guided)
  for(inopo=0; inopo<nnopo; inopo++)
  {
          nmct  = d1nmct[inopo];
      if(nmct<EPSILON||i1norem[inopo]==1)continue;
    iprop = i1nopr[inopo];


    if(i1bnvx[iprop]||DABS(nmct)<EPSILON)
    {
     d1nftx[inopo]=R0;
     dp1dv[0] =R0;
    }
    else
    {
     d1nftx[inopo] += (d1bnfx[iprop] + d1nfcx[inopo]+d1kax[iprop]*nmct
		     -4795.0*dbeta*nmct*d1nvcx[inopo])+(double)(i1rmgrav[inopo])*dcrmgravx*nmct;
     dp1dv[0] = d1nftx[inopo]/nmct;
    }

    if(i1bnvy[iprop]||DABS(nmct)<EPSILON)
    {
     d1nfty[inopo]=R0;
     dp1dv[1] =R0;
    }
    else
    {
     d1nfty[inopo] += (d1bnfy[iprop] + d1nfcy[inopo]+d1kay[iprop]*nmct
		     -4795.0*dbeta*nmct*d1nvcy[inopo])+(double)(i1rmgrav[inopo])*dcrmgravy*nmct;
     dp1dv[1] = d1nfty[inopo]/nmct;
    }

    if(i1bnvz[iprop]||DABS(nmct)<EPSILON)
    {
     d1nftz[inopo]=R0;
     dp1dv[2] =R0;
    }
    else
    {
     d1nftz[inopo] += (d1bnfz[iprop] + d1nfcz[inopo]+d1kaz[iprop]*nmct
		     -4795.0*dbeta*nmct*d1nvcz[inopo])+(double)(i1rmgrav[inopo])*dcrmgravz*nmct;
     dp1dv[2] = d1nftz[inopo]/nmct;
    }

      if(izero)
     {
      d1nftx[inopo]=R0;
      d1nfty[inopo]=R0;
      d1nftz[inopo]=R0;
      d1nfcx[inopo]=R0;
      d1nfcy[inopo]=R0;
      d1nfcz[inopo]=R0;
     }


  /* update nodal velocity nvc and coordinates ncc */
    /* current velocity */
    d1nvcx[inopo] += dp1dv[0]*dcstec;
    d1nvcy[inopo] += dp1dv[1]*dcstec;
    d1nvcz[inopo] += dp1dv[2]*dcstec;

    /* override using supplied velocity (turns out this is quicker than if/else) */
    if(i1bnvx[iprop]) d1nvcx[inopo] = d1bnvx[iprop];
    if(i1bnvy[iprop]) d1nvcy[inopo] = d1bnvy[iprop];
    if(i1bnvz[iprop]) d1nvcz[inopo] = d1bnvz[iprop];
if(index==1)
{
    d1nvcx[inopo] = d1nvcx[inopo]*(double)(i1rmgrav[inopo]);
    d1nvcy[inopo] = d1nvcy[inopo]*(double)(i1rmgrav[inopo]);
    d1nvcz[inopo] = d1nvcz[inopo]*(double)(i1rmgrav[inopo]);
}
    /* new coordinates */
    d1nccx[inopo] += d1nvcx[inopo]*dcstec;
    d1nccy[inopo] += d1nvcy[inopo]*dcstec;
    d1nccz[inopo] += d1nvcz[inopo]*dcstec;

}

}


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief explicit solver of equations
 *  \param[in]     ydc control database
 *  \param[in]     yde element database
 *  \param[in,out] ydn nodal database
 *  \param[out]    ydo output database
 *  \param[in]     ydb boundary condition database
 *  \param[in]     ydk constants (for d2csmm consistent mass matrix)
 *  \par Details:
 *  Ysd() calls Ysd3MEC() which is the mechanical solver for 3D nodes
 *  with x,y,z d.o.f..  Adjustments are made afterward for the history
 *  variable of total kinetic energy.
 *
 */
void Ysd(YDC ydc, YDE yde, YDN ydn, YDO ydo, YDB ydb, YDK ydk)
{
  INT i, ihys,j;
  DBL Ek, ohys, stprev, ohyp,avx,avy,avz,av,bind;
    int izero;
    int index;
    for(i=0;i<yde->nelem;i++)
    {
	    if(yde->d1eldam[i]>=0.959)
	    {       bind=1.0;
		    avx=(ydn->d2nvc[0][yde->i2elto[i][0]]+ydn->d2nvc[0][yde->i2elto[i][1]]+
		    ydn->d2nvc[0][yde->i2elto[i][2]]+ydn->d2nvc[0][yde->i2elto[i][3]])/4.0;
		    avy=(ydn->d2nvc[1][yde->i2elto[i][0]]+ydn->d2nvc[1][yde->i2elto[i][1]]+                   
			ydn->d2nvc[1][yde->i2elto[i][2]]+ydn->d2nvc[1][yde->i2elto[i][3]])/4.0;
		    avz=(ydn->d2nvc[2][yde->i2elto[i][0]]+ydn->d2nvc[2][yde->i2elto[i][1]]+
                     ydn->d2nvc[2][yde->i2elto[i][2]]+ydn->d2nvc[2][yde->i2elto[i][3]])/4.0;
		    av=SQRT(avx*avx+avy*avy+avz*avz);
		    if(av>bind)avx=avx*bind/av;
		    if(av>bind)avy=avy*bind/av;
		    if(av>bind)avz=avz*bind/av;

		    for(j=0;j<4;j++)
		    {
//			    ydn->i1ndam[yde->i2elto[i][j]]=1;
                ydn->d2nft[0][yde->i2elto[i][j]]=ydn->d2nft[0][yde->i2elto[i][j]]/1000.0;
			    ydn->d2nft[1][yde->i2elto[i][j]]=ydn->d2nft[1][yde->i2elto[i][j]]/1000.0;
			    ydn->d2nft[2][yde->i2elto[i][j]]=ydn->d2nft[2][yde->i2elto[i][j]]/1000.0;
			    ydn->d2nvc[0][yde->i2elto[i][j]]=avx;
			    ydn->d2nvc[1][yde->i2elto[i][j]]=avy;
			    ydn->d2nvc[2][yde->i2elto[i][j]]=avz;
		    }
	    }
    }
    izero=1;
    if((ydc->ncstep%ydc->icoutf)==0)izero=0;
    index=0;
    if(ydc->ncstep>=ydc->mcstep)
    {
	    index=1;

    }
  /* mechanical solver for 3D nodes with x,y,z d.o.f.  */
  if(ydc->initer>0)
  Ysd3MEC_multi_pass(
    ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
    ydn->d2nfc[0], ydn->d2nfc[1], ydn->d2nfc[2],
    ydn->d2nft[0], ydn->d2nft[1], ydn->d2nft[2],
    ydn->d2nvc[0], ydn->d2nvc[1], ydn->d2nvc[2],
    ydb->i1bnvx,   ydb->i1bnvy,   ydb->i1bnvz,
    ydb->d1bnvx,   ydb->d1bnvy,   ydb->d1bnvz,
    ydb->d1bnfx,   ydb->d1bnfy,   ydb->d1bnfz,
    ydk->d2ka[0],  ydk->d2ka[1],  ydk->d2ka[2],
    ydk->d2csmm,
    ydn->nnopo,
    ydn->d1nmct,
    ydn->i1nopr,
    ydn->i1nobf,
    yde->d1emct,
    yde->nelemi,
    yde->i2elto,
    ydc->initer,
    ydc->dcstec,
    ydc->dcurelx,
    izero,
    yde->i1elchprm,ydn->i1nochprm
    );
   else
   Ysd3MEC(
    ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
    ydn->d2nfc[0], ydn->d2nfc[1], ydn->d2nfc[2],
    ydn->d2nft[0], ydn->d2nft[1], ydn->d2nft[2],
    ydn->d2nvc[0], ydn->d2nvc[1], ydn->d2nvc[2],
    ydb->i1bnvx,   ydb->i1bnvy,   ydb->i1bnvz,
    ydb->d1bnvx,   ydb->d1bnvy,   ydb->d1bnvz,
    ydb->d1bnfx,   ydb->d1bnfy,   ydb->d1bnfz,
    ydk->d2ka[0],  ydk->d2ka[1],  ydk->d2ka[2],
    ydn->nnopo,
    ydn->d1nmct,
    ydn->i1nopr,
    ydn->i1nobf,
    yde->nelemi,
    yde->i2elto,
    ydc->dcstec,
    izero,
    ydn->i1nochprm,ydn->i1rmgrav,
    ydc->dcrmgrax,ydc->dcrmgray,ydc->dcrmgraz,index

    );

    if((ydc->ncstep%1000)==0)
    {
  for(ihys=0; ihys<ydo->nohys; ihys++) /* get history variables */
  {
    if(ydo->i1ohyt[ihys] == YFLEK) /* total kinetic energy */
    {
      Ek=R0;
      for(i=0; i<ydn->nnopo; i++)
      {
          if(ydn->i1noel[i]==1)
          {
        Ek += RP5*(SQR(ydn->d2nvc[0][i]) + SQR(ydn->d2nvc[1][i]) + SQR(ydn->d2nvc[2][i]))*ydn->d1nmct[i];
          }
      }

      /* isolated macros to avoid redundant evaluation */
      ohys = ABS(ydo->d1ohys[ihys]);
      stprev = MAXIM(EPSILON, ohys);
      ohyp = R1 - Ek/stprev;
      ohyp = ABS(ohyp);

      if(ohyp > ydo->dohyp)
      {
        ydo->d1ohys[ihys] = Ek;
        ydo->d1ohyt[ihys] = ydc->dctime;
      }
    }
    }
  }
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/
