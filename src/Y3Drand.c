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
/* DEBUG                                                              */
/**********************************************************************/



static void gidmsh(INT **i2elto,DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,INT iele,INT jele)

{
    FILE *fout=FILENULL;
    INT icount,jcount;
    fout=fopen("gidmsh.msh","w");
    
    CHRw(fout,"MESH    dimension 3 ElemType Tetrahedra  Nnode ");
    INTw(fout,4, 4);
    CHRwsp(fout);
    CHRwcr(fout);
    CHRw(fout,"Coordinates ");
    CHRwcr(fout);
    
    icount=0;
        for(jcount=0;jcount<4;jcount++)
        {
            icount++;
        INTw(fout,icount,5);
        CHRwsp(fout);
        
        DBLw(fout,d1nccx[i2elto[iele][jcount]],19);
        CHRwsp(fout);
        DBLw(fout,d1nccy[i2elto[iele][jcount]],19);
        CHRwsp(fout);
        DBLw(fout,d1nccz[i2elto[iele][jcount]],19);
        CHRwsp(fout);
        CHRwcr(fout);
            
        }
 
    for(jcount=0;jcount<4;jcount++)
    {
        icount++;
        INTw(fout,icount,5);
        CHRwsp(fout);
        
        DBLw(fout,d1nccx[i2elto[jele][jcount]],19);
        CHRwsp(fout);
        DBLw(fout,d1nccy[i2elto[jele][jcount]],19);
        CHRwsp(fout);
        DBLw(fout,d1nccz[i2elto[jele][jcount]],19);
        CHRwsp(fout);
        CHRwcr(fout);
    }

    CHRw(fout,"end coordinates ");
    CHRwcr(fout);
    CHRw(fout,"Elements ");
    CHRwcr(fout);
    
    for(icount=0;icount<2;icount++)
    {
        INTw(fout,icount+1,5);
        CHRwsp(fout);
        for(jcount=0;jcount<4;jcount++)
        {
            INTw(fout,icount*4+jcount+1,10);
            CHRwsp(fout);
        }
        CHRwcr(fout);
    }
    CHRw(fout,"end elements ");
    CHRwcr(fout);
    fclose(fout);
    
}

 
/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/

static INLINE int Ycd3TETintersection( DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
                                      INT iele, INT jele, INT **i2elto)
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
        //xtmp=xt-xorig; ytmp=yt-yorig; ztmp=zt-zorig;
        for(i=0;i<3;i++)
        { xc[i]=xc[i]-xorig; yc[i]=yc[i]-yorig; zc[i]=zc[i]-zorig;
        }
        /* contactor normal, e-base and target points in e-base */
        V3DCro(xnc,ync,znc,xc[1],yc[1],zc[1],xc[2],yc[2],zc[2]);
        //        V3DNor(xtmp,xnc,ync,znc);
        for(j=0;j<4;j++)
        {ii=i2elto[iele][j];
            xtmp=d1nccx[ii]-xorig; ytmp=d1nccy[ii]-yorig; ztmp=d1nccz[ii]-zorig;
            V3DDot(dct[icons][j],xnc,ync,znc,xtmp,ytmp,ztmp);
        }
        
        if((dct[icons][0]<R0)&&(dct[icons][1]<R0)&&
           (dct[icons][2]<R0)&&(dct[icons][3]<R0))  return 0;  //modified by JXiangJXiang
        
        for(i=0;i<4;i++)
        {
            if(dct[icons][i]>EPSILON)in_index[i]++;
               if((icons==3)&&(in_index[i]==4))return 1;
        }
        
        
        for(i=0;i<3;i++)
        {
            jn=p3[i];
            
            for(j=0;j<jn;j++)
            { inext=i+j+1;
                if((dct[icons][i]*dct[icons][inext])<R0)
                    
                    ///Modified by JXiang
                { ratio=dct[icons][i]-dct[icons][inext];
                    
                    if(ratio>EPSILON||ratio<NEPSILON)
                    { ratio=dct[icons][i]/ratio;
                    }
                    else
                    {
                        ratio=RP5;
                    }
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
                            V3DDot(dct_tmp,xv,yv,zv,xtmp,ytmp,ztmp);
                            if(dct_tmp<=R0){index=1;break;}
                            
                        }
                        if(index==0)return 1;
                        
                    } }
            
        }
        
        
    } 
    
    return 2;
}




/**********************************************************************/
/**********************************************************************/


/* find element coordinates (3D)
 */
static void YcdminmaxTET(
                       INT nelem_start,INT nelem,
                       DPT *dp1eccmin,DPT *dp1eccmax,
                       DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
                       INT *i1elbe,INT **i2elto
                       )
{
    INT ielem;
    INT ip0, ip1, ip2, ip3;
    DBL x0, x1, x2, x3;
    DBL y0, y1, y2, y3;
    DBL z0, z1, z2, z3;
    DBL minx,miny,minz,maxx,maxy,maxz,dbuff;

    dbuff=R0;

/*    #pragma omp parallel shared(nelem,iprop,dp1eccmin,dp1eccmax,dbuff,\
    d1nccx, d1nccy, d1nccz,i1elbe,i1elpr,i2elto,iwfast) private(ielem,ip0, ip1, ip2, \
ip3,x0, x1, x2, x3,y0, y1, y2, y3,z0, z1, z2, z3,minx,miny,minz,maxx,maxy,maxz)
 */
    {
		
//#pragma omp for schedule(static,1000)
	
    for(ielem=0; ielem<nelem; ielem++)
    {
        if((i1elbe[ielem] <= 0)) continue;
        
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
        
    }
    }
}

static void Cell_sorting(INT nelemv,
            INT *i1elbe, INT *i1orig,INT *ielemnext, INT *ielemhea,
                         INT *i1nghost,INT nelem_start,INT nelem,DBL diam,INT *i1elemnbox,
            DPT *dp1eccmin,DPT *dp1eccmax)
{
    INT ncelx,ncely,ncelz;
    INT *i1cnz;      /* contactor next  z                      */
    INT *i1cfz;      /* contactor first z                      */
//    INT iz;             /* x, y, z cell                  */
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
    
    iminx -= 1; iminy -= 1; iminz -= 1;
    imaxx += 2; imaxy += 2; imaxz += 1;
#pragma omp parallel for default(none) shared(nelem,nelem_start,i1elbe,i1nghost,i1orig,i1eccxmin,iminx,\
i1eccymin,iminy,i1ecczmin,iminz,i1eccxmax,i1eccymax,i1ecczmax,i1eccx,i1eccy,i1eccz,\
i1elemnbox) private(ielem,ni,nj,nk,i,j,k,nelemvv) schedule(guided)
    for(ielem=0; ielem<nelem; ielem++)
    {
        i1orig[ielem]=ielem;
        if(i1elbe[ielem] > 0)
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
    if(i1elbe[i1orig[ielem]] > 0)
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
    
	#pragma omp for schedule(guided) //nowait
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

int Search_contact(INT nelem_start,INT nelem,INT *i1elbe,INT *ielemhea,INT *i1orig,
                           DPT *dp1eccmin,DPT *dp1eccmax,
                           DBL *d1nccx,DBL *d1nccy,DBL *d1nccz,INT **i2elto,
                           INT *i1elemnbox,INT *ielemnext,int nthreads,INT inew)
{
//#if defined(_OPENMP)
	    int *index_s,i;
    index_s=TalINT1(nthreads+1);
    TzINT1(index_s,nthreads+1);
#pragma omp parallel default(none) shared(nelem_start,nelem,i1elbe,ielemhea,i1orig,\
dp1eccmin,dp1eccmax,d1nccx,d1nccy,d1nccz,i2elto,\
i1elemnbox,ielemnext,inew,index_s)
{
    INT i1iect_t[2000],nlist,ielem,icoup,ielemx,ielemfbox,i,jcoup;
    INT icouco;             /* couple's contactor              */
    INT icouta;             /* couple's target                 */
    int index,tid;

	tid = omp_get_thread_num();
#pragma omp for schedule(dynamic,64)
    for(ielem=nelem_start; ielem<nelem; ielem++)
    {
        if(i1elbe[ielem] > 0)
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
                    if((i1elbe[icouco]>0)&&(i1elbe[icouta]>0))
                    {
                    if((icouco>icouta)&&((i1elbe[icouco]!=i1elbe[icouta]) &&(i1elbe[icouco]==inew))) //new to do
                    {
                                if(!((dp1eccmin[icouta][0]>dp1eccmax[icouco][0])||
                                 (dp1eccmax[icouta][0]<dp1eccmin[icouco][0])||
                                 (dp1eccmin[icouta][1]>dp1eccmax[icouco][1])||
                                 (dp1eccmax[icouta][1]<dp1eccmin[icouco][1])||
                                 (dp1eccmin[icouta][2]>dp1eccmax[icouco][2])||
                                 (dp1eccmax[icouta][2]<dp1eccmin[icouco][2]))) 
                            {
                                index=Ycd3TETintersection(d1nccx,d1nccy,d1nccz,icouco,icouta,i2elto);
                                if(index==2)index=Ycd3TETintersection(d1nccx,d1nccy,d1nccz,icouta,icouco,i2elto);
                                if(index==1||index==2)
                                 {

index_s[tid]=1;

}
                                }
                                
                            }                    
                     }
                    ielemx = ielemnext[ielemx];
                }
                
                ielemfbox=i1elemnbox[ielemfbox];
                
            }
            
        }
        
    }
}
    for(i=0;i<nthreads;i++)
    {
    if(index_s[i]==1)
    {
        FREE(index_s);
        return 1;
    }
    }
    FREE(index_s);
    return 0;




}


/**********************************************************************/
/**********************************************************************/


/* Y elements contact detection (Z - this can be split up using INLINE)
*/
int Ycdran(
            INT nelem_start,INT nelem,
            
            DBL *d1nccx, DBL *d1nccy, DBL *d1nccz, 
            INT *i1elbe, INT **i2elto, DBL diam,
            INT nnopo,
            int nthreads,INT inew
            )
{ 
    DPT *dp1eccmin,*dp1eccmax;     /* element coordinate current xyz         */
    INT icoup,jcoup;        /* couple                          */
    INT ielem,i,ni,nj,nk;                /* element                       */
    INT nelemv;
    INT *i1orig,ihyis,*i1nghost,nelemvv;
    INT *ielemhea, *ielemnext, *i1elemnbox;
    int index;
    index=0;
      /* find element coordinates */
  dp1eccmin = TalzDPT1(nelem);    /* element coordinate current min xyz */
  dp1eccmax = TalzDPT1(nelem);    /* element coordinate current max xyz */

        YcdminmaxTET(
                     nelem_start,nelem ,
                     dp1eccmin,dp1eccmax,
                     d1nccx, d1nccy, d1nccz,
                     i1elbe,i2elto
                     );
    
    nelemv=0;
    i1nghost=TalINT1(nelem+1);
    nelemvv=nelem;
    for(ielem=0; ielem<nelem; ielem++)
    {
        
        ni=(INT)(dp1eccmax[ielem][0]/diam)-(INT)(dp1eccmin[ielem][0]/diam)+1;
        nj=(INT)(dp1eccmax[ielem][1]/diam)-(INT)(dp1eccmin[ielem][1]/diam)+1;
        nk=(INT)(dp1eccmax[ielem][2]/diam)-(INT)(dp1eccmin[ielem][2]/diam)+1;
        i1nghost[ielem]=nelemv+nelemvv-ielem;
            nelemv=nelemv+ni*nj*nk;
	}
        i1orig = TalINT1(nelemv);
        i1elemnbox=TalINT1(nelemv);
        ielemnext=TalINT1(nelemv);
        ielemhea=TalINT1(nelemv);
    
    Cell_sorting(nelemv,i1elbe,i1orig,ielemnext,ielemhea,i1nghost,
                nelem_start,nelem,diam,i1elemnbox,dp1eccmin,dp1eccmax);
    index=Search_contact(nelem_start,nelem,i1elbe,ielemhea,i1orig,
                    dp1eccmin,dp1eccmax,
                   d1nccx,d1nccy,d1nccz,i2elto,i1elemnbox,ielemnext,nthreads,inew);
    
    /* free memory (note: in reverse order of allocation) */
    FREE(ielemhea);
    FREE(ielemnext);
    FREE(i1elemnbox);
    FREE(i1orig);
    FREE(i1nghost);

    FREE(dp1eccmax);
    FREE(dp1eccmin);
return index;
//	}

 
	

}
