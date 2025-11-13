/*! \file Y3Dinit.c
 *  \brief Y constants initialization
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


static INLINE DBL Ycd3TETminilength( double *d1nccx, double *d1nccy, double *d1nccz,
                                      INT iele, INT **i2elto)
{
    DBL xorig,yorig,zorig,dct[4];
    DBL xc[4],yc[4],zc[4],xnc,ync,znc,xtmp,ytmp,ztmp,xct[4],yct[4],zct[4];
    int ipc[4],i,j,icons;
    INT ii;
    int p1[4] = { 1, 3, 3, 1 };
    int p2[4] = { 2, 2, 0, 0 };
    int p3[4] = { 3, 0, 1, 2 };
    DBL le_mini,dnor;
    for(i=0;i<4;i++)
    {ii=i2elto[iele][i];
        xct[i]=d1nccx[ii];
        yct[i]=d1nccy[ii];
        zct[i]=d1nccz[ii];
 
    }
    for(icons=0;icons<4;icons++)
    {
        ipc[0]=icons;
        ipc[1]=p1[icons];
        ipc[2]=p2[icons];
        ipc[3]=p3[icons];
        
        for(i=0;i<4;i++)
        {
            xc[i]=xct[ipc[i]];
            yc[i]=yct[ipc[i]];
            zc[i]=zct[ipc[i]];
        }
        
        xorig=xc[0]; yorig=yc[0]; zorig=zc[0];
        for(i=0;i<4;i++)
        { xc[i]=xc[i]-xorig; yc[i]=yc[i]-yorig; zc[i]=zc[i]-zorig;
        }
        /* contactor normal, e-base and target points in e-base */
        V3DCro(xnc,ync,znc,xc[1],yc[1],zc[1],xc[2],yc[2],zc[2]);
                V3DNor(dnor,xnc,ync,znc);
            xtmp=xc[3]; ytmp=yc[3]; ztmp=zc[3];
            V3DDot(dct[icons],xnc,ync,znc,xtmp,ytmp,ztmp);
        }
    le_mini=MINIM(ABS(dct[0]),ABS(dct[1]));
                  le_mini=MINIM(le_mini,ABS(dct[2]));
                  le_mini=MINIM(le_mini,ABS(dct[3]));
 
    return le_mini;
}


/* calculate consistent mass matrix
*/
static void CalConMassM(DBL **d2csmm)
{
  DBL d1shapf[10][4];
  DBL d1xi[4]   = { alpha, beta,  beta,  beta };
  DBL d1eta[4]  = { beta,  alpha, beta,  beta };
  DBL d1zeta[4] = { beta,  beta,  alpha, beta };
  INT ig,ng,i,j;

  ng=4; /* NOTE:  ng != NGRSH */

  for(i=0;i<NNODE;i++)
  {
    for(j=0;j<NNODE;j++)
    {
      d2csmm[i][j]=R0;
    }
  }
 
  for(ig=0;ig<ng;ig++)
  {
    d1shapf[0][ig]=(R1-R2*d1xi[ig]-R2*d1eta[ig]-R2*d1zeta[ig])*(R1-d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d1shapf[1][ig]=d1xi[ig]*(R2*d1xi[ig]-R1);
    d1shapf[2][ig]=d1eta[ig]*(R2*d1eta[ig]-R1);
    d1shapf[3][ig]=d1zeta[ig]*(R2*d1zeta[ig]-R1);
    d1shapf[4][ig]=R4*d1xi[ig]*(R1-d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d1shapf[5][ig]=R4*d1xi[ig]*d1eta[ig];
    d1shapf[6][ig]=R4*d1eta[ig]*(R1-d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d1shapf[7][ig]=R4*d1zeta[ig]*(R1-d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d1shapf[8][ig]=R4*d1xi[ig]*d1zeta[ig];
    d1shapf[9][ig]=R4*d1zeta[ig]*d1eta[ig];

    for(i=0;i<NNODE;i++)
    {
      for(j=0;j<NNODE;j++)
      {
        d2csmm[i][j] += WEIGHT*d1shapf[i][ig]*d1shapf[j][ig];
      }
    }
  }
}


/**********************************************************************/
/**********************************************************************/


/* calculate lumped mass matrix
*/
static void CalLumMassM(DBL **d2csmm, DBL *d1lpmm)
{
  DBL sum;
  INT i;

  for(sum=R0,i=0;i<NNODE;i++) sum += d2csmm[i][i];   /* sum of diagonal */

  for(i=0;i<NNODE;i++) d1lpmm[i] = d2csmm[i][i]/sum; /* normalized diagonal */
}


/**********************************************************************/
/**********************************************************************/


/* quadratic shape function (deformation gradient?)
*/
static void SHAPEGRAD(DBL ***d3dsh)
{
  DBL d1xi[NGRSH]   = { RP25, alpha, beta,  beta,  beta };
  DBL d1eta[NGRSH]  = { RP25, beta,  alpha, beta,  beta };
  DBL d1zeta[NGRSH] = { RP25, beta,  beta,  alpha, beta };
  INT ig;

  for(ig=0;ig<NGRSH;ig++)
  {    	
    d3dsh[ig][0][0]=R4*(d1xi[ig]+d1eta[ig]+d1zeta[ig])-R3;
    d3dsh[ig][0][1]=R4*d1xi[ig]-R1;
    d3dsh[ig][0][2]=R0;
    d3dsh[ig][0][3]=R0;
    d3dsh[ig][0][4]=R4*(R1-R2*d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d3dsh[ig][0][5]=R4*d1eta[ig];
    d3dsh[ig][0][6]=-R4*d1eta[ig];
    d3dsh[ig][0][7]=-R4*d1zeta[ig];
    d3dsh[ig][0][8]=R4*d1zeta[ig];
    d3dsh[ig][0][9]=R0;

    d3dsh[ig][1][0]=R4*(d1xi[ig]+d1eta[ig]+d1zeta[ig])-R3;
    d3dsh[ig][1][1]=R0;
    d3dsh[ig][1][2]=R4*d1eta[ig]-R1;
    d3dsh[ig][1][3]=R0;
    d3dsh[ig][1][4]=-R4*d1xi[ig];
    d3dsh[ig][1][5]=R4*d1xi[ig];
    d3dsh[ig][1][6]=R4*(R1-d1xi[ig]-R2*d1eta[ig]-d1zeta[ig]);
    d3dsh[ig][1][7]=-R4*d1zeta[ig];
    d3dsh[ig][1][8]=R0;
    d3dsh[ig][1][9]=R4*d1zeta[ig];

    d3dsh[ig][2][0]=R4*(d1xi[ig]+d1eta[ig]+d1zeta[ig])-R3;
    d3dsh[ig][2][1]=R0;
    d3dsh[ig][2][2]=R0;
    d3dsh[ig][2][3]=R4*d1zeta[ig]-R1;
    d3dsh[ig][2][4]=-R4*d1xi[ig];
    d3dsh[ig][2][5]=R0;
    d3dsh[ig][2][6]=-R4*d1eta[ig];
    d3dsh[ig][2][7]=R4*(R1-d1xi[ig]-d1eta[ig]-R2*d1zeta[ig]);
    d3dsh[ig][2][8]=R4*d1xi[ig];
    d3dsh[ig][2][9]=R4*d1eta[ig];
  }
}


/**********************************************************************/
/**********************************************************************/


static int testIEEE(void)
{
  double zero_flt = 0.0;
  double zero_bits = 42.0;
  memset(&zero_bits, 0, sizeof(double));
  return (0 == memcmp(&zero_flt, &zero_bits, sizeof(double))); /* (zero_flt == zero_bits) */
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


/* small strain elastic tetrahetra
 */
static void Ymass3TET10ELS(
                         INT nelem, INT iprop,DBL dpero,
                         DBL **d2ncc,  DBL **d2nci,
                         INT *i1elpr,  INT **i2elto,
                           DBL *d1emct, DBL  *d1nmct, DBL *d1lpmm,
                           DBL ***d3dsh
                         )
{
        DBL d2xj[NDIME][NDIME], d2xjci[NDIME][NDIME];
        DBL voli, detci;
        INT ielem, ig;
        INT *i1elto;
        INT i, ip;

        for(ielem=0; ielem<nelem; ielem++)
        {
            if(i1elpr[ielem] != iprop) continue;
            i1elto = i2elto[ielem];
            
            /* simpler calculations for ig=0 */
            ig = 0;
            /* calculate Jacobian matrix xj */
//            jacobi(d2xj, d2ncc, d3dsh[ig], i1elto);
            voli = R0;
            /* calculate inverse of xj (debug: warn if determinant negative) */
            /* iterate overy Gauss points using deformation (shape) gradient */
            for(ig=1;ig<NGRSH;ig++)
            {
                /* full calculations for ig>0 */
                /* calculate the Jacobian matrixes xj, xjci */
                jacobi2(d2xj, d2ncc, d2xjci, d2nci, d3dsh[ig], i1elto);
                /* calculate determinant of xjci */
                YMATDET3(d2xjci,detci);
                voli += ABS(detci)*WEIGHT/R6;

            }
            
            /* calculate nodal force */
                d1emct[ielem] = dpero*voli; //Z property ro - density
            
            for(i=0; i<NNODE; i++)
            {
                ip = i1elto[i];
                /* nodal mass current translation */
                d1nmct[ip] += d1emct[ielem]*d1lpmm[i];
            }
        }

}


static void Ymass3TET4el(  /* small strain elastic 4-noded tetrahedra  */
                        INT    nelem,  INT    iprop,
                        DBL dpero,DBL *d1ncix,  DBL  *d1nciy,
                        DBL *d1nciz,DBL *d1emct, DBL  *d1nmct,INT *i1elpr,  
			INT **i2elto, DBL *d1mart, DBL dcstec, DBL lamda, DBL mu, DBL pepe, DBL dcmart
                    )
{ DBL voli;

    DBL F0[3][3]; /* initial local base */
    DBL F0inv[3][3]; /* global base in initial local base */
    INT ielem,i,*i1elto,itmp;
    DBL tm,le,dcstec_el,ifa,penalty;
    ifa=0.2;
        for(ielem=0;ielem<nelem;ielem++)
        {
            if(i1elpr[ielem]==iprop)
            {
                i1elto = i2elto[ielem];
                for(i=1;i<4;i++)
                {
                    F0[0][i-1]=d1ncix[i1elto[i]]-d1ncix[i1elto[0]]; /* init. base */
                    F0[1][i-1]=d1nciy[i1elto[i]]-d1nciy[i1elto[0]];
                    F0[2][i-1]=d1nciz[i1elto[i]]-d1nciz[i1elto[0]];
                }
                YMATINV3(F0,F0inv,voli);         /* global base in initial local coordinates    */
                
                /*Calculate Cauchy stress*/

                if(voli<R0)
                {
        //        CHRw(stdout,"direction change  ");
                itmp=i2elto[ielem][2];
                i2elto[ielem][2]=i2elto[ielem][3];
                i2elto[ielem][3]=itmp;
                voli=-voli;
                }
                le=Ycd3TETminilength( d1ncix, d1nciy, d1nciz,ielem, i2elto);
		penalty=MAXIM(pepe,lamda+2.0*mu);
		dcstec_el=dcmart*le*SQRT(dpero/penalty);
		if(dcstec>dcstec_el)d1mart[ielem]=dcstec/dcstec_el;
		d1mart[ielem]=d1mart[ielem]*d1mart[ielem];
                d1emct[ielem] = d1mart[ielem]*dpero*voli/6.0; //Z property ro - density

                for(i=0;i<4;i++)  
                {

                    tm=dpero*voli/24.0*d1mart[ielem];
                        d1nmct[i1elto[i]]=d1nmct[i1elto[i]]+tm;
                }
		
            }
        }
    }


static void Ymass3JOINT(  /* small strain elastic 6-noded joint  */
                                   INT nelem, DBL *d1pepe,INT **i2ejjp,
                                   DBL *d1emct,DBL *d1nmct,INT *i1elpr,
                                   DBL *d1mart,INT **i2elto,DBL *d1picf,
				   DBL *d1jpmr,DBL dcstec
                         )
{
	DBL dcstec_jp,m1,m2, factor;
	INT ielem,n1,n2;
	for(ielem=0;ielem<nelem;ielem++)
	  {
	if(i2elto[ielem][0]!=i2elto[ielem][5])
	{
	        n1=i2ejjp[ielem][0];
	        n2=i2ejjp[ielem][1];
		m1=d1emct[n1];
		m2=d1emct[n2];
		factor=R1;
		if(d1picf[i1elpr[ielem]]>R1)factor=d1picf[i1elpr[ielem]];
          dcstec_jp=SQRT(m1*m2/(factor*(m1+m2)*d1pepe[i1elpr[ielem]]))*R2;
	  if(dcstec>dcstec_jp)
	  {
		  factor=dcstec/dcstec_jp;
		  factor=factor*factor;
//		  d1jpmr[ielem]=factor;
                  if(factor>d1mart[n1])d1mart[n1]=factor;
                  if(factor>d1mart[n2])d1mart[n2]=factor;
	  }
	  printf("joint MR, %ld %ld, %le\n",ielem,nelem,d1jpmr[ielem]);

	}
	  }
}

static void YMassRatioTET4el(  /* small strain elastic 4-noded tetrahedra  */
                        INT    nelem,  INT    iprop,
                        DBL *d1emct, DBL  *d1nmct, INT *i1elpr,
			INT **i2elto, DBL *d1mart
                    )
{ 
    INT ielem,i,*i1elto,itmp;
    DBL tm;

        for(ielem=0;ielem<nelem;ielem++)
        {
            if(i1elpr[ielem]==iprop)
            {
                i1elto = i2elto[ielem];

                d1emct[ielem] = d1mart[ielem]*d1emct[ielem]; //Z property ro - density

                for(i=0;i<4;i++)  
                {

                    tm=d1emct[ielem]/4.0;
                        d1nmct[i1elto[i]]=d1nmct[i1elto[i]]+tm;
                }
	     printf("mass %le ielem %ld	\n",d1emct[ielem],ielem);
            }
        }
    }


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief constants database initialization  
 *  \par Details:
 *  Ykd() initializes all global constant elements of the YKD structure...
 *
 */
void Yinit(YDK ydk)
{
  ydk->d2csmm = TalDBL2(NNODE, NNODE);
  ydk->d1lpmm = TalDBL1(NNODE);
  ydk->d3dsh  = TalDBL3(NGRSH, NDIME, NNODE);

  CalConMassM(ydk->d2csmm);

  CalLumMassM(ydk->d2csmm, ydk->d1lpmm);

  SHAPEGRAD(ydk->d3dsh);

  ydk->d2ka = DBL2NULL;

  /* meh */
  initcode();

  /* if ever this is reported returning false, TzXXXn() will need revised */
  if(!testIEEE())
  {
    CHRw(stderr,"Yinit: unexpected IEEE floating point implementation...\n");
    CHRw(stderr,"       contact authors for a software patch");
    CHRwcr(stderr); exit(1);
  }
}

void Ymass3TET(YDE yde, YDN ydn, YDP ydp, YDK ydk, YDJ ydj, YDJP ydjp, int nthreads, DBL dcstec)
{
    INT iprop,i,j;
    

    /* zero nodal forces and masses */
    
    for(i=0; i<yde->nelemi; i++)
    {
        for(j=0; j<(yde->nelno-1); j++)
        {
            ydn->i1noel[yde->i2elto[i][j]]=yde->i2elto[i][yde->nelno-1];
        }
    }
    
    TzDBL1(ydn->d1nmct, ydn->mnopo);
    
    for(iprop=0; iprop<ydp->nprop; iprop++)
    {
        if(ydp->i1ptyp[iprop] == YTE3TET10ELS) /* small strain elastic tetrahetra  */
        {
            Ymass3TET10ELS(
                           yde->nelemi,  iprop,
                           ydp->d1pero[iprop],
                           ydn->d2ncc,  ydn->d2nci,
                           yde->i1elpr,  yde->i2elto,
                           yde->d1emct, ydn->d1nmct,
                           ydk->d1lpmm,ydk->d3dsh
                           );
        }
        else if((ydp->i1ptyp[iprop])==(YTE3TET4ELS))
        {
            Ymass3TET4el(  /* small strain elastic 4-noded tetrahedra  */
                                   yde->nelemi,  iprop,
                                   ydp->d1pero[iprop],ydn->d2nci[0],ydn->d2nci[1],
                                   ydn->d2nci[2],yde->d1emct,ydn->d1nmct,yde->i1elpr,  
				   yde->i2elto,yde->d1mart,dcstec,ydp->d1pela[iprop],
				   ydp->d1pemu[iprop],ydp->d1pepe[iprop],ydp->d1pemr[iprop]
                         );

        }
    }

}




/*! \brief constants database configuration for current run
 *  \par Details:
 *  Ykd() initializes all per-run constant elements of the YKD structure...
 *
 */
void Yconfig(YDK ydk, YDC ydc, YDB ydb)   
{
  INT i;

  if(ydk->d2ka) FREE(ydk->d2ka);

  ydk->d2ka = TalDBL2(NDIME, ydb->mbcon);

  for(i=0; i<ydb->nbcon; i++)
  {
    ydk->d2ka[0][i] = ydb->d1bnax[i] + ydc->dcgrax;
    ydk->d2ka[1][i] = ydb->d1bnay[i] + ydc->dcgray;
    ydk->d2ka[2][i] = ydb->d1bnaz[i] + ydc->dcgraz;
  }
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/



