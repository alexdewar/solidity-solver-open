/*! \file Y3Dcd.c
 *  \brief Y contact detection
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
//#include "Ycd3TETintersection.h"
#include <omp.h>


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


static INLINE int Ycd3TETintersectionbuf( double *d1nccx, double *d1nccy, double *d1nccz,
                                         INT iele, INT jele, INT **i2elto,double buf)
{
    DBL xorig,yorig,zorig,dct[4][4];
    DBL xc[3],yc[3],zc[3],xnc,ync,znc,xtmp,ytmp,ztmp,xct[4],yct[4],zct[4];
    int ipc[3],i,j,icons,k,jn,inext,jnext,index;
    INT ii,iinext,in_index[4];
    DBL dct_tmp,ratio,xline,yline,zline,xv,yv,zv,xint,yint,zint;
    int p1[4] = { 1, 3, 3, 1 };
    int p2[4] = { 2, 2, 0, 0 };
    int p3[3]={3,2,1};

    
    for(i=0;i<4;i++)
    {ii=i2elto[jele][i];
        xct[i]=d1nccx[ii];
        yct[i]=d1nccy[ii];
        zct[i]=d1nccz[ii];
        in_index[i]=0;
    }
    
    for(icons=0;icons<4;icons++)
    {
        ipc[0]=icons;
        ipc[1]=p1[icons];
        ipc[2]=p2[icons];
        
        for(i=0;i<3;i++)
        {
            xc[i]=xct[ipc[i]];
            yc[i]=yct[ipc[i]];
            zc[i]=zct[ipc[i]];
        }
        
        xorig=xc[0]; yorig=yc[0]; zorig=zc[0];
        for(i=0;i<3;i++)
        { xc[i]=xc[i]-xorig; yc[i]=yc[i]-yorig; zc[i]=zc[i]-zorig;
        }
        /* contactor normal, e-base and target points in e-base */
        V3DCro(xnc,ync,znc,xc[1],yc[1],zc[1],xc[2],yc[2],zc[2]);
        V3DNor(xtmp,xnc,ync,znc);
        for(j=0;j<4;j++)
        {ii=i2elto[iele][j];
            xtmp=d1nccx[ii]-xorig; ytmp=d1nccy[ii]-yorig; ztmp=d1nccz[ii]-zorig;
            V3DDot(dct[icons][j],xnc,ync,znc,xtmp,ytmp,ztmp);
        }
        
        if((dct[icons][0]<-buf)&&(dct[icons][1]<-buf)&&
           (dct[icons][2]<-buf)&&(dct[icons][3]<-buf))  return 0;  //modified by JXiang
        for(i=0;i<4;i++)
        {
            if(dct[icons][i]>=-buf)in_index[i]++;
               if((icons==3)&&(in_index[i]==4))return 1;
        }
        
        
        for(i=0;i<3;i++)
        {
            jn=p3[i];
            
            for(j=0;j<jn;j++)
            { inext=i+j+1;
                		  if(((dct[icons][i]>=-buf)&&(dct[icons][inext]<=-buf))||
                             ((dct[icons][i]<=-buf)&&(dct[icons][inext]>=-buf)))

                    ///Modified by JXiang
                { ratio=DABS(dct[icons][i]-dct[icons][inext]);
                    if(ratio>EPSILON)
                    {
                        ratio=DABS((dct[icons][i]+buf)/ratio);
                        iinext=i2elto[iele][inext];
                        ii=i2elto[iele][i];
                        xint=ratio*d1nccx[iinext]+(R1-ratio)*d1nccx[ii];
                        yint=ratio*d1nccy[iinext]+(R1-ratio)*d1nccy[ii];
                        zint=ratio*d1nccz[iinext]+(R1-ratio)*d1nccz[ii];
                        index=0;
                        for(k=0;k<3;k++)
                        {
                            jnext=k+1;if(k==2)jnext=0;
                            
                            xtmp=xint-xct[ipc[k]];
                            ytmp=yint-yct[ipc[k]];
                            ztmp=zint-zct[ipc[k]];
                            
                            xline=xc[jnext]-xc[k];
                            yline=yc[jnext]-yc[k];
                            zline=zc[jnext]-zc[k];
                            V3DCro(xv,yv,zv,xnc,ync,znc,xline,yline,zline);
                            V3DNor(dct_tmp,xv,yv,zv);
                            V3DDot(dct_tmp,xv,yv,zv,xtmp,ytmp,ztmp);
                            if(dct_tmp<-buf){index=1;break;}
                            
                        }
                        if(index==0)return 1;
                        
                    } }  }
                        
        }
        
        
    } 

    return 2;
}




static INT Ycddis3TET(
                      INT nnopo,
                      DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
                      DBL **d2nck,
                      INT *i1nobf,
                      DBL diezon,INT *i1nobfse, INT nnobf,
                      int nthreads, INT *i1noremv
                      )
{
    INT ip,i,*index;
    int tid;
    DBL vx, vy, vz, v2,vx0,vy0,vz0,sdiezon;
        sdiezon=diezon*0.5;
    index=TalINT1(nthreads+1);


#pragma omp parallel default(none) shared(d1nccx, d1nccy, d1nccz,d2nck,sdiezon,index,i1nobfse, nnobf,i1noremv) \
private (i,ip,vx, vy, vz,v2,vx0,vy0,vz0,tid)
    {
        tid = omp_get_thread_num();
        index[tid]=0;
#pragma omp for schedule(guided) nowait
    for(i=0; i<nnobf; i++)
    {

            ip=i1nobfse[i];
            if(i1noremv[ip]==1)continue;
            vx = d1nccx[ip];
            vy = d1nccy[ip];
            vz = d1nccz[ip];
            vx0 = d2nck[0][ip];
            vy0 = d2nck[1][ip];
            vz0 = d2nck[2][ip];
            
            v2 = SQR(vx-vx0) + SQR(vy-vy0) + SQR(vz-vz0);
			
        if(SQRT(v2)>=sdiezon) index[tid]=1;
 
    }
    }

    for(i=0;i<nthreads;i++)
    {
    if(index[i])
    {
        FREE(index);
        return 1;
    }
    }
    FREE(index);
    return 0;
}

/* find maximum squared velocity (3D)
*/
static DBL Ycdvel3TET(
           INT nelem,
           INT iprop,
           DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
           INT *i1elpr,
           INT **i2elto
           ) 
{ 
  INT ielem, i, ip, *i1elto;
  DBL vx, vy, vz, v2;
  DBL mxv2 = R0;

  for(ielem=0; ielem<nelem; ielem++)
  { 
    if(i1elpr[ielem] != iprop) continue;

    i1elto = i2elto[ielem];

    for(i=0; i<NNODEX; i++)
    { 
      ip = i1elto[i]; 

      vx = d1nvcx[ip];
      vy = d1nvcy[ip];
      vz = d1nvcz[ip];

      v2 = SQR(vx) + SQR(vy) + SQR(vz);
      if(v2 > mxv2) mxv2 = v2;
    }
  }

  return mxv2;
}


/**********************************************************************/
/**********************************************************************/


/* find element coordinates (3D)
*/
static void Ycdcor3TET(  
            INT nelem,
            INT iprop,
            DBL *d1erad,
            DPT *dp1ecc,
            DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
            INT *i1elbe,
            INT *i1elpr,
            INT **i2elto,
            INT iwfast
            ) 
{ 
  INT ielem, *i1elto;
  INT ip0, ip1, ip2, ip3;
  DBL x0, x1, x2, x3;
  DBL y0, y1, y2, y3;
  DBL z0, z1, z2, z3;
  DBL dx0, dx1, dx2, dx3;
  DBL dy0, dy1, dy2, dy3;
  DBL dz0, dz1, dz2, dz3;
  DBL xc, yc, zc;
  DBL d2, mxd2;

  for(ielem=0; ielem<nelem; ielem++)
  { 
    if((i1elpr[ielem] != iprop) || ((i1elbe[ielem] <= 0) && (iwfast != 0))) continue;

    i1elto = i2elto[ielem];

    ip0 = i1elto[0]; 
    ip1 = i1elto[1]; 
    ip2 = i1elto[2]; 
    ip3 = i1elto[3]; 

    x0 = d1nccx[ip0]; x1 = d1nccx[ip1]; x2 = d1nccx[ip2]; x3 = d1nccx[ip3];
    y0 = d1nccy[ip0]; y1 = d1nccy[ip1]; y2 = d1nccy[ip2]; y3 = d1nccy[ip3];
    z0 = d1nccz[ip0]; z1 = d1nccz[ip1]; z2 = d1nccz[ip2]; z3 = d1nccz[ip3];

    xc = (x0 + x1 + x2 + x3)/R4;
    yc = (y0 + y1 + y2 + y3)/R4;
    zc = (z0 + z1 + z2 + z3)/R4;

    dp1ecc[ielem][0] = xc;
    dp1ecc[ielem][1] = yc;
    dp1ecc[ielem][2] = zc;

    /* use temp variables to avoid double-calculation within MAXIM() */
    dx0 = xc - x0; dx1 = xc - x1; dx2 = xc - x2; dx3 = xc - x3;
    dy0 = yc - y0; dy1 = yc - y1; dy2 = yc - y2; dy3 = yc - y3;
    dz0 = zc - z0; dz1 = zc - z1; dz2 = zc - z2; dz3 = zc - z3;

    /* faster than using MAXIM() */
    mxd2 = SQR(dx0) + SQR(dy0) + SQR(dz0);

    d2 = SQR(dx1) + SQR(dy1) + SQR(dz1);
    if(d2 > mxd2) mxd2 = d2;

    d2 = SQR(dx2) + SQR(dy2) + SQR(dz2);
    if(d2 > mxd2) mxd2 = d2;

    d2 = SQR(dx3) + SQR(dy3) + SQR(dz3);
    if(d2 > mxd2) mxd2 = d2;

    d1erad[ielem] = SQRT(mxd2);
  } 
}

/* find element coordinates (3D)
 */
static void YcdminmaxTET(
                       INT nelem,
                       DPT *dp1eccmin,DPT *dp1eccmax,DBL dizone,
                       DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
                       INT *i1elbe,INT **i2elto,INT iwfast, INT *i1elremv
                       )
{
    INT ielem;
    INT ip0, ip1, ip2, ip3;
    DBL x0, x1, x2, x3;
    DBL y0, y1, y2, y3;
    DBL z0, z1, z2, z3;
    DBL minx,miny,minz,maxx,maxy,maxz,dbuff;
    DBL minxx,minyy,minzz;

    minxx=BEPSILON; minyy=BEPSILON; minzz=BEPSILON;
    dbuff=dizone*0.5+EPSILON;

/*    #pragma omp parallel shared(nelem,iprop,dp1eccmin,dp1eccmax,dbuff,\
    d1nccx, d1nccy, d1nccz,i1elbe,i1elpr,i2elto,iwfast) private(ielem,ip0, ip1, ip2, \
ip3,x0, x1, x2, x3,y0, y1, y2, y3,z0, z1, z2, z3,minx,miny,minz,maxx,maxy,maxz)
 */
    {
		
//#pragma omp for schedule(static,1000)
	
    for(ielem=0; ielem<nelem; ielem++)
    {
        if(((i1elbe[ielem] <= 0) && (iwfast != 0))||(i1elremv[ielem]==1))
	{
		dp1eccmin[ielem][0]=R0;
		dp1eccmin[ielem][1]=R0;
		dp1eccmin[ielem][2]=R0;

		dp1eccmax[ielem][0]=R0;
		dp1eccmax[ielem][1]=R0;
		dp1eccmax[ielem][2]=R0;
	}
	else
	{
        ip0 = i2elto[ielem][0];
        ip1 = i2elto[ielem][1];
        ip2 = i2elto[ielem][2];
        ip3 = i2elto[ielem][3];
        
        x0 = d1nccx[ip0]; x1 = d1nccx[ip1]; x2 = d1nccx[ip2]; x3 = d1nccx[ip3];
        y0 = d1nccy[ip0]; y1 = d1nccy[ip1]; y2 = d1nccy[ip2]; y3 = d1nccy[ip3];
        z0 = d1nccz[ip0]; z1 = d1nccz[ip1]; z2 = d1nccz[ip2]; z3 = d1nccz[ip3];
        
        minx=x0;
        if (x1<minx)minx=x1;
        if (x2<minx)minx=x2;
        if (x3<minx)minx=x3;
        maxx=x0;
        if(x1>maxx)maxx=x1;
        if(x2>maxx)maxx=x2;
        if(x3>maxx)maxx=x3;
        
        miny=y0;
        if (y1<miny)miny=y1;
        if (y2<miny)miny=y2;
        if (y3<miny)miny=y3;
        maxy=y0;
        if(y1>maxy)maxy=y1;
        if(y2>maxy)maxy=y2;
        if(y3>maxy)maxy=y3;
        
        minz=z0;
        if (z1<minz)minz=z1;
        if (z2<minz)minz=z2;
        if (z3<minz)minz=z3;
        maxz=z0;
        if(z1>maxz)maxz=z1;
        if(z2>maxz)maxz=z2;
        if(z3>maxz)maxz=z3;
        
        dp1eccmin[ielem][0]=minx-dbuff;
        dp1eccmin[ielem][1]=miny-dbuff;
        dp1eccmin[ielem][2]=minz-dbuff;
        
        dp1eccmax[ielem][0]=maxx+dbuff;
        dp1eccmax[ielem][1]=maxy+dbuff;
        dp1eccmax[ielem][2]=maxz+dbuff;

	if(dp1eccmin[ielem][0]<minxx)minxx=dp1eccmin[ielem][0];
	if(dp1eccmin[ielem][1]<minyy)minyy=dp1eccmin[ielem][1];
	if(dp1eccmin[ielem][2]<minzz)minzz=dp1eccmin[ielem][2];
	}
    }

    for(ielem=0; ielem<nelem; ielem++)
        {
        if(((i1elbe[ielem] <= 0) && (iwfast != 0))||(i1elremv[ielem]==1))
        {
                dp1eccmin[ielem][0]=R0;
                dp1eccmin[ielem][1]=R0;
                dp1eccmin[ielem][2]=R0;

                dp1eccmax[ielem][0]=R0;
                dp1eccmax[ielem][1]=R0;
	       	dp1eccmax[ielem][2]=R0;
        }
       else
	 {

		 dp1eccmin[ielem][0]=dp1eccmin[ielem][0]-minxx+EPSILON;
		 dp1eccmin[ielem][1]=dp1eccmin[ielem][1]-minyy+EPSILON;
		 dp1eccmin[ielem][2]=dp1eccmin[ielem][2]-minzz+EPSILON;

		 dp1eccmax[ielem][0]=dp1eccmax[ielem][0]-minxx+EPSILON;
		 dp1eccmax[ielem][1]=dp1eccmax[ielem][1]-minyy+EPSILON;
		 dp1eccmax[ielem][2]=dp1eccmax[ielem][2]-minzz+EPSILON;
	}
		}
    }
}

static void Cell_sorting(INT nelemv,
            INT *i1elbe, INT *i1orig,INT *ielemnext, INT *ielemhea,
                         INT *i1nghost,INT nelem,DBL diam,INT *i1elemnbox,
            DPT *dp1eccmin,DPT *dp1eccmax, INT *i1elremv)
{
    INT ncelx,ncely,ncelz;
    INT *i1cnz;      /* contactor next  z                      */
    INT *i1cfz;      /* contactor first z                      */
    INT iz;             /* x, y, z cell                  */
    INT nelemd,ielem,nelemvv;             /* twice total number of elements  */
    INT *i1eccxmin,*i1eccxmax,*i1eccx;     /* element coordinate current ix          */
    INT *i1eccymin,*i1eccymax,*i1eccy;     /* element coordinate current iy          */
    INT *i1ecczmin,*i1ecczmax,*i1eccz;     /* element coordinate current iz          */
    INT iminx,iminy,iminz;    /* space boundaries              */
    INT imaxx,imaxy,imaxz;    /* space boundaries              */
    INT ni,nj,nk,i,j,k;
    
    iminx = INT_MAX; iminy = INT_MAX; iminz = INT_MAX;
    imaxx = INT_MIN; imaxy = INT_MIN; imaxz = INT_MIN;
    
    i1eccxmin = TalINT1(nelem);     /* element coordinate current ix  */
    i1eccymin = TalINT1(nelem);     /* element coordinate current iy  */
    i1ecczmin = TalINT1(nelem);     /* element coordinate current iz  */
    
    i1eccxmax = TalINT1(nelem);     /* element coordinate current ix  */
    i1eccymax = TalINT1(nelem);     /* element coordinate current iy  */
    i1ecczmax = TalINT1(nelem);     /* element coordinate current iz  */
    
    i1eccx = TalINT1(nelemv);     /* element coordinate current ix  */
    i1eccy = TalINT1(nelemv);     /* element coordinate current iy  */
    i1eccz = TalINT1(nelemv);     /* element coordinate current iz  */

    for(ielem=0; ielem<nelem; ielem++)
    {
        if((i1elbe[ielem] > 0)&&(i1elremv[ielem]!=1))
        {
        i1eccxmin[ielem] = (INT)(dp1eccmin[ielem][0]/diam);
        i1eccymin[ielem] = (INT)(dp1eccmin[ielem][1]/diam);
        i1ecczmin[ielem] = (INT)(dp1eccmin[ielem][2]/diam);
        
        i1eccxmax[ielem] = (INT)(dp1eccmax[ielem][0]/diam);
        i1eccymax[ielem] = (INT)(dp1eccmax[ielem][1]/diam);
        i1ecczmax[ielem] = (INT)(dp1eccmax[ielem][2]/diam);
        
        
        if(i1eccxmin[ielem] < iminx) iminx = i1eccxmin[ielem];
        if(i1eccxmax[ielem] > imaxx) imaxx = i1eccxmax[ielem];
        
        if(i1eccymin[ielem] < iminy) iminy = i1eccymin[ielem];
        if(i1eccymax[ielem] > imaxy) imaxy = i1eccymax[ielem];
        
        if(i1ecczmin[ielem] < iminz) iminz = i1ecczmin[ielem];
        if(i1ecczmax[ielem] > imaxz) imaxz = i1ecczmax[ielem];
        }
	else
	{
		i1eccxmin[ielem] = 0;
		i1eccymin[ielem] = 0;
		i1ecczmin[ielem] = 0;

		i1eccxmax[ielem] = 0;
		i1eccymax[ielem] = 0;
		i1ecczmax[ielem] = 0;
	}
    }
    
    iminx -= 1; iminy -= 1; iminz -= 1;
    imaxx += 2; imaxy += 2; imaxz += 2;
#pragma omp parallel for default(none) shared(nelem,i1elbe,i1nghost,i1orig,i1eccxmin,iminx,\
i1eccymin,iminy,i1ecczmin,iminz,i1eccxmax,i1eccymax,i1ecczmax,i1eccx,i1eccy,i1eccz,\
i1elemnbox,i1elremv) private(ielem,ni,nj,nk,i,j,k,nelemvv) schedule(dynamic,64)
    for(ielem=0; ielem<nelem; ielem++)
    {
        i1orig[ielem]=ielem;
        if((i1elbe[ielem] > 0)&&(i1elremv[ielem]!=1))
        {
            nelemvv=i1nghost[ielem];
            i1eccxmin[ielem] -= iminx;
            i1eccymin[ielem] -= iminy;
            i1ecczmin[ielem] -= iminz;
            i1eccxmax[ielem] -= iminx;
            i1eccymax[ielem] -= iminy;
            i1ecczmax[ielem] -= iminz;
            
            ni=i1eccxmax[ielem]-i1eccxmin[ielem];
            nj=i1eccymax[ielem]-i1eccymin[ielem];
            nk=i1ecczmax[ielem]-i1ecczmin[ielem];
            
            if((ni+nj+nk)==0)
            {
                i1eccx[ielem]=i1eccxmin[ielem];
                i1eccy[ielem]=i1eccymin[ielem];
                i1eccz[ielem]=i1ecczmin[ielem];
                i1elemnbox[ielem]=-1;
                
            }
            else
            {
                
                for(i=0;i<(ni+1);i++)
                {
                    for(j=0;j<(nj+1);j++)
                    {
                        for(k=0;k<(nk+1);k++)
                        {
                            if((i+j+k)==0)
                            {
                                i1eccx[ielem]=i1eccxmin[ielem];
                                i1eccy[ielem]=i1eccymin[ielem];
                                i1eccz[ielem]=i1ecczmin[ielem];
                                i1elemnbox[ielem]=nelemvv;
                            }
                            else
                            {
                                i1eccx[nelemvv]=i1eccxmin[ielem]+i;
                                i1eccy[nelemvv]=i1eccymin[ielem]+j;
                                i1eccz[nelemvv]=i1ecczmin[ielem]+k;
                                i1orig[nelemvv]=ielem;
                                i1elemnbox[nelemvv]=nelemvv+1;
                                nelemvv++;
                            }
                        }
                    }
                }
                i1elemnbox[nelemvv-1]=-1;
            }
        }
    }
    
    /* calculate cells, allocate memory  */
    ncelx = imaxx - iminx;
    ncely = imaxy - iminy;
    ncelz = imaxz - iminz;


    i1cnz    = TalINT1(nelemv); /* contactor next z                    */
    i1cfz    = TalINT1(ncelz); /* contactor first z                   */

    Tn1INT1(i1cfz,ncelz);
for(ielem=0; ielem<nelemv; ielem++)
{
    if((i1elbe[i1orig[ielem]] > 0)&&(i1elremv[i1orig[ielem]]!=1))
    {
        
        i1cnz[ielem] = i1cfz[i1eccz[ielem]];
        i1cfz[i1eccz[ielem]] = ielem;
        
    }
}

nelemd = nelemv*2;
/* scan all loaded z cells */ /* to make parallel, i1heaz[] etc. must be local to loop */
    
#pragma omp parallel default(none) shared(ncelz,i1cfz,i1cnz,ielemhea,ielemnext,\
	nelemd,nelemv,i1eccy,i1eccx,ncely,ncelx)
{
	INT *i1cnx;      /* contactor next  x                      */
    INT *i1cny;      /* contactor next  y                      */

    INT  i1heax;  /* heads of 5 connected lists for x cells */
    INT  i1heay;  /* heads of 5 connected lists for y cells */
    INT  i1heaz;  /* heads of 2 connected lists for z cells */
    INT *i1cfx;   /* contactor first x                      */
    INT *i1cfy;   /* contactor first y                      */
    INT ix,iy,iz;             /* x, y, z cell                  */
    INT ielemx,ielemy,ielemz,jelem;
    
/* assume no contactors at any cell */

    i1cnx    = TalINT1(nelemv); /* contactor next x                    */
    i1cny    = TalINT1(nelemv); /* contactor next y                    */
    i1cfy = TalINT1(ncely); /* contactor first y  (iz-1)           */
    i1cfx = TalINT1(ncelx); /* contactor first x  (iz-1,iy-1)      */

    Tn1INT1(i1cfx,ncelx);
    Tn1INT1(i1cfy,ncely);
    
	#pragma omp for schedule(dynamic,64) //nowait
    for(iz=0; iz<ncelz; iz++)
    {
        jelem=i1cfz[iz];

        if((jelem>=0))
        {
            i1heaz = i1cfz[iz];
            /* load elements from cells iz & iz-1 onto y cells */
            ielemz = i1heaz;
            while(ielemz >= 0)
            {
                i1cny[ielemz] = i1cfy[i1eccy[ielemz]];
                i1cfy[i1eccy[ielemz]] = ielemz;
                ielemz = i1cnz[ielemz];
            }
            /* scan all loaded y cells */
            ielemz = i1heaz;
            while(ielemz >= 0)
            {
                iy = i1eccy[ielemz];
                if(i1cfy[iy] < nelemv)
                {
                    i1heay = i1cfy[iy];
                    i1cfy[iy] += nelemd;
                    
                    /* load elements from y cells onto x cells */
                    ielemy = i1heay;
                    while(ielemy >= 0)
                    {
                        i1cnx[ielemy] = i1cfx[i1eccx[ielemy]];
                        i1cfx[i1eccx[ielemy]] = ielemy;
                        ielemy = i1cny[ielemy];
                    }
                    /* scan all loaded x cells */
                    ielemy = i1heay;
                    while(ielemy >= 0)
                    {
                        ix = i1eccx[ielemy];
                        if(i1cfx[ix] < nelemv)
                        {
                            i1heax = i1cfx[ix];
                            
                            i1cfx[ix] += nelemd;
                            
                            /* detect contacts for cell (ix,iy,iz) */
                            ielemx = i1heax;
                            while(ielemx >= 0)
                            {
                                /* form the neighbour list */
                                ielemnext[ielemx] = i1cnx[ielemx];
                                ielemhea[ielemx] = i1heax;
                                ielemx = i1cnx[ielemx];
                            } 
                        }
                        ielemy = i1cny[ielemy];
                    }
                    
                    /* unload elements from x cells */
                    Tn1INT1(i1cfx,ncelx);
                }
                ielemz = i1cnz[ielemz];
            }
            
            /* unload elements from y cells */
            Tn1INT1(i1cfy,ncely);
        }
	}
	FREE(i1cfx);
    FREE(i1cfy);
	FREE(i1cny);
    FREE(i1cnx);
    }
    FREE(i1cfz);
    FREE(i1cnz);
    FREE(i1eccz);
    FREE(i1eccy);
    FREE(i1eccx);
    FREE(i1ecczmax);
    FREE(i1eccymax);
    FREE(i1eccxmax);
    FREE(i1ecczmin);
    FREE(i1eccymin);
    FREE(i1eccxmin);
}

static void Search_contact(INT nelem,INT *i1elbe,INT *ielemhea,INT *i1orig,INT iwfast,INT *i1elcf,INT *i1iect,
                           DPT *dp1eccmin,DPT *dp1eccmax,INT *i1stat,INT *i1iecn,INT *i0iecff,
                           DBL *d1nccx,DBL *d1nccy,DBL *d1nccz,INT **i2elto,DBL diezon,
                           INT *i1elemnbox,INT *ielemnext,INT **i2eljp,INT *i1elpr,
                           INT *i1elremv, INT *i1elst,INT *i1fcstep,INT *i1fcjp,DBL *d1pepe,
                           INT ncstep,INT *i1eljt, INT *i1elncstep)
{
//#if defined(_OPENMP)
#pragma omp parallel default(none) shared(nelem,i1elbe,ielemhea,i1orig,iwfast,stderr,i1elcf,i1iect,\
dp1eccmin,dp1eccmax,i1stat,i1iecn,i0iecff,d1nccx,d1nccy,d1nccz,i2elto,diezon,\
i1elemnbox,ielemnext,i2eljp,i1elpr,i1elremv,i1elst,i1fcstep,i1fcjp,d1pepe,ncstep,i1eljt,i1elncstep)
{
    INT i1iect_t[500],nlist,ielem,icoup,ielemx,ielemfbox,i,jcoup;
    INT icouco;             /* couple's contactor              */
    INT icouta;             /* couple's target                 */
    int index,flag,isurf,jsurf, ijoint,ibreak;
#pragma omp for schedule(dynamic,64)
    //#endif
    
    for(ielem=1; ielem<nelem; ielem++)
    {
        if((i1elbe[ielem] > 0)&&(i1elremv[ielem]!=1))
        {
            nlist=0;
            ielemfbox=ielem;  /* ielem's all sub-element */
            while(ielemfbox>=0)
            {
                
                ielemx=ielemhea[ielemfbox];
                /*ielemhea the first element of bounding box in which ielemfbox is */
                while(ielemx >= 0)
                {
                    icouco = ielem;
                    icouta = i1orig[ielemx];
                    //                indexp=1;

                           flag=0;
                    ijoint=-1;
			if(icouco>icouta)
			{

			if((i1elbe[icouco]!=i1elbe[icouta])|| (iwfast==0))flag=1;
	                     else
        	                {
					       if((i1elbe[icouco]>=YIPROPMAX)||(i1elbe[icouta]>=YIPROPMAX)) /* has one problem, should be ||  */
						 {
						   flag=1;

						   for(isurf=0;isurf<4;isurf++)
						     {
							     ibreak=0;
						       if(i2eljp[icouco][isurf]>=0)
							 {
							   for(jsurf=0;jsurf<4;jsurf++)
							     {
							       if((i2eljp[icouco][isurf]==i2eljp[icouta][jsurf]))
								 {
								   flag=0;
								   ibreak=1;
								   break;	    
								 }
							     }

							 }
                                 else if(i2eljp[icouco][isurf]<=(-YIPROPMAX))
                                 {
			        
                                     if(i1elst[-(i2eljp[icouco][isurf]+YIPROPMAX)]==1)
                                     {
                                         ijoint=-1;
                                     for(jsurf=0;jsurf<4;jsurf++)
                                       {
                                         if((i2eljp[icouco][isurf]==i2eljp[icouta][jsurf]))
                                       {
                                        ijoint=-(i2eljp[icouco][isurf]+YIPROPMAX);
                                           i1elst[ijoint]=2;
					   flag=1;
					   ibreak=1;
                                         break;
                                       }
                                       }
                                     }

                                 }
				 if(ibreak==1)break;
					 
						     }
                                                   }
					}
				}

                    if(flag==1)
                    {
                        if(ijoint<0&&(i1eljt[icouta]==1||i1eljt[icouco]==1))
			{
                if((i1elbe[icouco]!=i1elbe[icouta])&&
		((ABS(i1elbe[icouco]-i1elbe[icouta])%YIPROPMAX)==0))
                   {
			ijoint=-2;
                }
                   else if((i1elbe[icouco]==i1elbe[icouta])&&(i1elbe[icouco]>=YIPROPMAX))
                   {
            ijoint=-2;
                }
			}
                        /* exclude old couple */
                        icoup = i1elcf[icouco];
                        jcoup = icoup;
                        while((icoup>=0) && (icouta>=0))
                        {
                            if(i1iect[icoup] == (-1 - icouta))
                            {
                                icouta = -1;
                            }
                            else if((i1iect[icoup] == icouta))
                            {

                                if((dp1eccmin[icouta][0]>dp1eccmax[icouco][0])||
                                   (dp1eccmax[icouta][0]<dp1eccmin[icouco][0])||
                                   (dp1eccmin[icouta][1]>dp1eccmax[icouco][1])||
                                   (dp1eccmax[icouta][1]<dp1eccmin[icouco][1])||
                                   (dp1eccmin[icouta][2]>dp1eccmax[icouco][2])||
                                   (dp1eccmax[icouta][2]<dp1eccmin[icouco][2]))
                                {
                                    
                                    i1stat[icoup]=5;
                                    i1iect[icoup] = (-1 - icouta);
                                    
                                    
                                    
                                }
                                else
                                {
                                    /*    if(i1nelsuf[icouta]>0)
                                     {  */
                                    index=Ycd3TETintersectionbuf(d1nccx,d1nccy,d1nccz,icouco,icouta,i2elto,diezon);
                                    if(index==2)index=Ycd3TETintersectionbuf(d1nccx,d1nccy,d1nccz,icouta,icouco,i2elto,diezon);
                                    
                                    if(index!=1)
                                    {
                                        i1stat[icoup]=5;
                                        i1iect[icoup] = (-1 - icouta);
                                        //                   i1elstat[icouco]=i1elstat[icouco]+1;
                                        
                                    }
                                    else
                                    {
                                        i1stat[icoup]=2;
                                        i1iect[icoup] = (-1 - icouta);
                                    }
                                }
                                icouta = -1;
                            }
                            jcoup = icoup;
                            icoup = i1iecn[icoup];
                        }
                        
                        if(icouta >= 0)
                        {
                            for(i=0;i<nlist;i++)
                            {
                                if(i1iect_t[i]==icouta)
                                {
                                    icouta=-1;
                                    break;
                                }
                            }
                            
                        }
                        
                        /* add new couple if close */
                        if(icouta >= 0)
                        {
                            if(!((dp1eccmin[icouta][0]>dp1eccmax[icouco][0])||
                                 (dp1eccmax[icouta][0]<dp1eccmin[icouco][0])||
                                 (dp1eccmin[icouta][1]>dp1eccmax[icouco][1])||
                                 (dp1eccmax[icouta][1]<dp1eccmin[icouco][1])||
                                 (dp1eccmin[icouta][2]>dp1eccmax[icouco][2])||
                                 (dp1eccmax[icouta][2]<dp1eccmin[icouco][2])))
                            {
                                
                                index=Ycd3TETintersectionbuf(d1nccx,d1nccy,d1nccz,icouco,icouta,i2elto,diezon);
                                if(index==2)index=Ycd3TETintersectionbuf(d1nccx,d1nccy,d1nccz,icouta,icouco,i2elto,diezon);
                                
                                
                                if(index==1)
                                {
                                 #pragma omp critical
                                    {
                                        icoup = *i0iecff;
                                        *i0iecff = i1iecn[icoup];
                                        i1iecn[icoup] = i1elcf[icouco];
                                        
                                    }
                                    if(icoup < 0)
                                    {
                                        CHRw(stderr,"Ycd: to small meicoc");
                                        CHRwcr(stderr); exit(1);
                                    }
                                    i1elcf[icouco] = icoup;
                                    i1iect[icoup] = -1 - icouta;/* mark new */
                                    i1stat[icoup]=1;
                                    i1fcjp[icoup]=ijoint;
                                    i1fcstep[icoup]=ncstep;
                                    d1pepe[icoup]=1.0;
                                    
                                }
                                else
                                {
                                    i1iect_t[nlist]=icouta;
                                    nlist++;
                                    if(nlist>=1500)
                                    {
                                    CHRw(stderr,"nlist: too small meicoc");
                                    CHRwcr(stderr); 
                                    exit(1);
                                    }
                                    
                                }
                                
                            }
                        }
                        
                    }
                    /*end check */
                    
                    ielemx = ielemnext[ielemx];
                }
                
                ielemfbox=i1elemnbox[ielemfbox];
                
            }
            
        }
        
    }
}
}


/**********************************************************************/
/**********************************************************************/


/* Y elements contact detection (Z - this can be split up using INLINE)
*/
static void Ycdprocess(
            INT micoup,  INT nelem,
            DBL diezon,  INT *i0iecff,
            INT *nicoup, DBL *d0iedi,
            DBL *d1nccx, DBL *d1nccy, DBL *d1nccz, 
            INT *i1elcf, INT *i1iecn,  INT *i1iect,
            INT *i1elbe, INT **i2elto, INT iwfast,DBL diam,
            INT *i1stat,DBL **d2nck, INT nnopo, INT *i1nobf,
            int *icheck, INT **i2eljp,INT *i1elpr,INT *i1nobfse,
            INT nnobf,int nthreads, INT *i1elremv, INT *i1noremv,
            INT *i1elst,INT *i1fcstep,INT *i1fcjp,DBL *d1pepe, 
	    INT ncstep,INT *i1eljt,INT *i1elncstep
            )
{ 
    DPT *dp1eccmin,*dp1eccmax;     /* element coordinate current xyz         */
    INT icoup,jcoup;        /* couple                          */
    INT ielem,i,ni,nj,nk;                /* element                       */
    INT index,nelemv;
    INT *i1orig,ihyis,*i1nghost,nelemvv;
    INT *ielemhea, *ielemnext, *i1elemnbox,ielemfbox;

    ihyis=0;
  if(*i0iecff == -2)
  {
	  
//  #if defined(_OPENMP)
//#pragma omp parallel shared(melem,i1elcf,micoup,i1iecn,i1stat,d2nck,d1nccx,d1nccy,d1nccz) private(ielem,icoup,i)
//#endif
//{
      
//#if defined(_OPENMP)
//#pragma omp for schedule(static,100) nowait
//#endif

    for(ielem=0; ielem<nelem; ielem++)
    { 
      i1elcf[ielem] = -1;
    }

//#if defined(_OPENMP)
//#pragma omp for schedule(static,100) nowait
//#endif

    for(icoup=0; icoup<micoup; icoup++)
    { 
      i1iecn[icoup] = icoup + 1;
        i1stat[icoup]=-1;
    }

//    #if defined(_OPENMP)
//#pragma omp for schedule(static,100) nowait
//#endif

      for(i=0; i<nnopo; i++)
      {
          d2nck[0][i]=d1nccx[i];
          d2nck[1][i]=d1nccy[i];
          d2nck[2][i]=d1nccz[i];
      }
 //     }
    i1iecn[micoup - 1] = -1;
    *i0iecff = 0;
      ihyis=1;
  }

  /* find maximum  velocity */
        
        index=Ycddis3TET(
                         nnopo,
                         d1nccx, d1nccy, d1nccz,
                         d2nck,
                         i1nobf,
                         diezon,i1nobfse, nnobf,nthreads, i1noremv
                         );

  /* if not time for contact detection, return */
    if(*d0iedi > 200.0)
    {
     ihyis=1;
    *d0iedi=0.0;
    }
    if((index==0)&&(ihyis==0))return;
//#pragma omp parallel shared(d2nck,d1nccx,d1nccy,d1nccz,i1iect,i1stat) private(icoup,i)

//    {
		
//#pragma omp  for schedule(static,1000) nowait
 
	
    for(i=0; i<nnopo; i++)
    {
        d2nck[0][i]=d1nccx[i];
        d2nck[1][i]=d1nccy[i];
        d2nck[2][i]=d1nccz[i];
    }
  /* unmark new couples */
	

//#pragma omp  for schedule(static,32) nowait


  for(icoup=0; icoup<micoup; icoup++)
  {
      if(i1stat[icoup]>0)
      {
          if(i1stat[icoup]==10)
         {
            i1iect[icoup] = -1 - i1iect[icoup];
             i1stat[icoup]=2;
         }
          else if (i1stat[icoup]>10)
          {
            i1iect[icoup] = -1 - i1iect[icoup];  
          }
      else
      {
       i1stat[icoup]=3;
      }
  }
    }
  /* find element coordinates */
  dp1eccmin = TalzDPT1(nelem);    /* element coordinate current min xyz */
  dp1eccmax = TalzDPT1(nelem);    /* element coordinate current max xyz */


        YcdminmaxTET(
                     nelem ,
                     dp1eccmin,dp1eccmax,diezon,
                     d1nccx, d1nccy, d1nccz,
                     i1elbe,i2elto,iwfast,i1elremv
                     );
    
    nelemv=0;
    i1nghost=TalINT1(nelem+1);
    nelemvv=nelem;
    for(ielem=0; ielem<nelem; ielem++)
    {
 //       if((i1elbe[ielem] > 0)&&(i1elremv[ielem]!=1))
        {
        ni=(INT)(dp1eccmax[ielem][0]/diam)-(INT)(dp1eccmin[ielem][0]/diam)+1;
        nj=(INT)(dp1eccmax[ielem][1]/diam)-(INT)(dp1eccmin[ielem][1]/diam)+1;
        nk=(INT)(dp1eccmax[ielem][2]/diam)-(INT)(dp1eccmin[ielem][2]/diam)+1;
        if((ni+nj+nk)==3)
          {
            i1nghost[ielem]=ielem;
            nelemv=nelemv+ni*nj*nk;
           }
          else
           {
            i1nghost[ielem]=nelemv+nelemvv-ielem;
            nelemv=nelemv+ni*nj*nk;

            }
        }
	}
        i1orig = TalINT1(nelemv);
        i1elemnbox=TalINT1(nelemv);
        ielemnext=TalINT1(nelemv);
        ielemhea=TalINT1(nelemv);
    
    Cell_sorting(nelemv,i1elbe,i1orig,ielemnext,ielemhea,i1nghost,
                nelem,diam,i1elemnbox,dp1eccmin,dp1eccmax,i1elremv);
    Search_contact(nelem,i1elbe,ielemhea,i1orig,iwfast,i1elcf,i1iect,
                    dp1eccmin,dp1eccmax,i1stat,i1iecn,i0iecff,
                   d1nccx,d1nccy,d1nccz,i2elto,diezon,i1elemnbox,ielemnext,
                   i2eljp,i1elpr,i1elremv,i1elst,i1fcstep,i1fcjp,d1pepe,ncstep,i1eljt,i1elncstep);
    
    /* free memory (note: in reverse order of allocation) */
    FREE(ielemhea);
    FREE(ielemnext);
    FREE(i1elemnbox);
    FREE(i1orig);
    FREE(i1nghost);

    FREE(dp1eccmax);
    FREE(dp1eccmin);


 if((*i0iecff)>(*nicoup))(*nicoup)=(*i0iecff);
    *icheck=1;
	
//#if defined(_OPENMP)
//#pragma omp for schedule(dynamic)
//#endif

    for(ielem=1; ielem<nelem; ielem++)
    {
 //      if(i1elstat[ielem]>-1)
       if(i1elbe[ielem] > 0)
        {
        icoup=i1elcf[ielem];
            jcoup=-1;
               while(icoup>=0)
            {
                           if((i1stat[icoup]==5)||(i1stat[icoup]==3))
                {
                    i1stat[icoup]=-1;

                    if(jcoup<0)
                {

                    {
                        i1elcf[ielem]=i1iecn[icoup];
                        i1iecn[icoup]=*i0iecff;
                        *i0iecff=icoup;
                        icoup=i1elcf[ielem];
                    }
                }
                else
                {
                    {
                        i1iecn[jcoup]=i1iecn[icoup];
                        i1iecn[icoup]=*i0iecff;
                        *i0iecff=icoup;
                        
                        icoup=i1iecn[jcoup];
                    }
                }
                }
                else
                { jcoup=icoup;
                    icoup=i1iecn[icoup];
                }
            }
        }
    }

}


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief contact detection
 *  \param[in]     ydc control database
 *  \param[in]     yde element database
 *  \param[in,out] ydi interaction database
 *  \param[in]     ydn nodal database
 *  \param[in]     ydp property database
 *  \par Details:
 *  Ycd() ... chapter 3 (Munjiza NBS)
 *
    \verbatim
    initialise data if needed
    find maximum  velocity
    if not time for contact detection, return
    find element coordinates - Ycdcor3TET()
    unmark new couples
    find maximum diameter
    intigerise element coordinates
    find space boundaries
    normalise coordinates
    assume no contactors at any cell
    assign all contactors to z-cells
    scan all loaded z cells (for each element)
      load elements from cells iz & iz-1  onto y cells
      scan all loaded y cells
        load elements from y cells onto x cells
        scan all loaded x cells
          detect contacts for cell (ix,iy,iz)
            for(ihx=0;ihx<5;ihx++)
              exclude old couple
              add new couple if close
        unload elements from x cells
      unload elements from y cells
 \endverbatim
 *
 */
void Ycd(YDE yde, YDI ydi, YDN ydn, INT iwfast,int nthreads,INT ncstep, INT *i1elst)
{ 
  /* Process contact detection if specified by input data */
  if(ydi->micoup > 0)
  { 
    Ycdprocess( /* Y elements contact detection */
      ydi->micoup,   yde->nelemi,
      ydi->diezon,    &ydi->iiecff,
      &ydi->nicoup, &(ydi->diedi),
      ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
      yde->i1elcf,   ydi->i1iecn,    ydi->i1iect,
      yde->i1elbe,   yde->i2elto,   iwfast, ydi->drect,
      ydi->i1stat,ydn->d2nck, ydn->nnopo,ydn->i1nobf,&(ydi->icheck),
      yde->i2eljp,yde->i1elpr,ydn->i1nobfse,ydn->nnobf,nthreads,
      yde->i1elchprm,ydn->i1nochprm,i1elst,ydi->i1fcstep,ydi->i1fcjp,
               ydi->d1pepe,ncstep,yde->i1eljt,yde->i1elncstep
      );
  }
}
/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/
