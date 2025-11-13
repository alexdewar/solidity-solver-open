/**********************************************************************/
/** Copyright (C) 2008,                                              **/
/** Queen Mary University of London (QMUL) & Imperial College        **/
/** of Science, Technology and Medicine (ICSTM). All rights reserved.**/
/** Implemented for you by Prof Antonio Munjiza & Dr Jiansheng Xiang **
 
 
 
* This code is part of the Virtual Geoscience Workbench (VGW) developed
* jointly by ICSTM and QMUL through two related parallel projects at 
* ICSTM and QMUL respectively funded by EPSRC. 
*
* This code is provided by copyright holders under the GNU Lesser 

* General Public License (LGPL). It is open source code; you can 
* redistribute it and/or modify it under the terms of the GNU Lesser 
* General Public License version 3.  
*  
* This code is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty 
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
* the GNU Lesser General Public License for more details,
* http://www.gnu.org/licenses/lgpl.txt. 
*  
* You should have received a copy of the GNU Lesser General Public 

* License along with this code; if not, write to 

* Dr Jiansheng Xiang Prof Antonio Munjiza or Dr John-Paul Latham 

* j.xiang@imperial.ac.uk a.munjiza@qmul.ac.uk 

* or j.p.latham@imperial.ac.uk 

* ******************************************************************* */ 

/*1 File   Y3Did.c */

#include "Yproto.h"
static void gidmsh_interaction(INT **i2elto,DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,INT iele,INT jele)

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
        
        DBLw(fout,d1nccx[i2elto[iele][jcount]],15);
        CHRwsp(fout);
        DBLw(fout,d1nccy[i2elto[iele][jcount]],15);
        CHRwsp(fout);
        DBLw(fout,d1nccz[i2elto[iele][jcount]],15);
        CHRwsp(fout);
        CHRwcr(fout);
        
    }
    
    for(jcount=0;jcount<4;jcount++)
    {
        icount++;
        INTw(fout,icount,5);
        CHRwsp(fout);
        
        DBLw(fout,d1nccx[i2elto[jele][jcount]],15);
        CHRwsp(fout);
        DBLw(fout,d1nccy[i2elto[jele][jcount]],15);
        CHRwsp(fout);
        DBLw(fout,d1nccz[i2elto[jele][jcount]],15);
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

static INLINE int Ycd3TETintersection( double *d1nccx, double *d1nccy, double *d1nccz,
                                      INT iele, INT jele, INT **i2elto)
{
    DBL xorig,yorig,zorig,dct[4];
    DBL xc[4],yc[4],zc[4],xnc,ync,znc,xtmp,ytmp,ztmp,xct[4],yct[4],zct[4];
    int ipc[4],i,j,icons,k,jn,inext,jnext,index;
    INT ii,iinext;
    DBL dct_tmp,ratio,xline,yline,zline,xv,yv,zv,xint,yint,zint;
    int p1[4] = { 1, 3, 3, 1 };
    int p2[4] = { 2, 2, 0, 0 };
    int p3[3]={3,2,1};
    int p[4] = { 3, 0, 1, 2 };
    DBL dnor; 
    for(i=0;i<4;i++)
    {ii=i2elto[jele][i];
        xct[i]=d1nccx[ii];
        yct[i]=d1nccy[ii];
        zct[i]=d1nccz[ii];
    }
    for(icons=0;icons<4;icons++)
    {
        ipc[0]=icons;
        ipc[1]=p1[icons];
        ipc[2]=p2[icons];
        ipc[3]=p[icons];
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
		V3DDot(dnor,xnc,ync,znc,xc[3],yc[3],zc[3]);
		dnor=ABS(dnor);
        for(j=0;j<4;j++)
        {ii=i2elto[iele][j];
            xtmp=d1nccx[ii]-xorig; ytmp=d1nccy[ii]-yorig; ztmp=d1nccz[ii]-zorig;
            V3DDot(dct[j],xnc,ync,znc,xtmp,ytmp,ztmp);
        }
        
        if((dct[0]/dnor<1.0e-8)&&(dct[1]/dnor<1.0e-8)&&
           (dct[2]/dnor<1.0e-8)&&(dct[3]/dnor<1.0e-8))  return 0;  //modified by JXiang        
        //  if(buf<EPSILON)
        //  {
        for(i=0;i<3;i++)
        {
            jn=p3[i];
            
            for(j=0;j<jn;j++)
            { inext=i+j+1;
                if(((dct[i]>EPSILON1)&&(dct[inext]<NEPSILON1))||
                ((dct[i]<NEPSILON1)&&(dct[inext]>EPSILON1)))
                    ///Modified by JXiang
                { ratio=DABS(dct[i]-dct[inext]);
                    if(ratio/dnor>1.0e-08)
                    { ratio=DABS(dct[i]/ratio);
                if ((ratio > 1.0e-08)&&(ratio < (R1 - 1.0e-08)))
		{
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
                            if(dct_tmp/dnor<-1.0e-08){index=1;break;}
                            
                        }
                        if(index==0)return 1;
                        
                    } }  }
            			  }
            
            
        }
        
        
    } 
  
    return 2;
}



void V3DRot(afx,afy,afz,bfx,bfy,bfz,u,v,w,cosin,sine)

     DBL *afx; DBL *afy; DBL *afz; DBL bfx; DBL bfy; 

     DBL bfz; DBL u; DBL v; DBL w; DBL cosin; DBL sine;

{

  DBL a0,a1,a2,a3,a4,a5,a6;

  a1=u*bfx;

  a2=v*bfy;

  a3=w*bfz;

  a4=u*u;

  a5=v*v;

  a6=w*w;

  a0=a1+a2+a3;

  (*afx)=u*a0+(bfx*(a5+a6)-u*(a2+a3))*cosin+(-w*bfy+v*bfz)*sine;

  (*afy)=v*a0+(bfy*(a4+a6)-v*(a1+a3))*cosin+(w*bfx-u*bfz)*sine;

  (*afz)=w*a0+(bfz*(a5+a4)-w*(a1+a2))*cosin+(-v*bfx+u*bfy)*sine;



}

static void Yid2TRITRI(              /* Triangle to Triangle */
             nelem,
            dcstec,diezon, iprop, jprop,
            d1iesl,d1nccx,d1nccy,d1nfcx,d1nfcy,
            d1nvcx,d1nvcy,d1pepe,i1elcf,i1elpr,
            i1icff,i1iecn,i1iect,i2elto 
            )
  INT    nelem;
  DBL   dcstec; DBL   diezon; INT    iprop; INT    jprop;
  DBL  *d1iesl; DBL  *d1nccx; DBL  *d1nccy; DBL  *d1nfcx; DBL  *d1nfcy;
  DBL  *d1nvcx; DBL  *d1nvcy; DBL  *d1pepe; INT  *i1elcf; INT  *i1elpr;
  INT  *i1icff; INT  *i1iecn; INT  *i1iect; INT **i2elto;
{ INT kprop,icontact,ielem,jelem,icoup,jcoup,it,jt,in,jn,ie,je,ip,jp,np;
  DBL a0,a1,a2,b0,b1,b2,c0,c1,c2,n0,n1,n2,fn,fna,fnb;
  DBL pen,tmp,dmin2,smin,smax;
  DBL small=EPSILON;
  DBL nsmall=-EPSILON;
  DBL big=BEPSILON;
  DBL zone2;
  DBL p[10];
  DBL s[10];
  DBL fx[3];
  DBL fy[3];
  DBL vol[2];
  DBL rx[2][3];
  DBL ry[2][3];
  DBL nx[2][3];
  DBL ny[2][3];
  DBL d[2][3][3];
  INT i2to[2][3];

  zone2=(R4*diezon*diezon);
  pen=MINIM(d1pepe[iprop],d1pepe[jprop]);
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { kprop=jprop;
    }
    else if(i1elpr[ielem]==jprop)
    { kprop=iprop;
    }
    else
    { kprop=-1;
    }
    if(kprop>=0)
    { icoup=i1elcf[ielem];
      jcoup=-1;
      while(icoup>=0)
      { if(i1iect[icoup]<0)i1iect[icoup]=-1-i1iect[icoup];
        jelem=i1iect[icoup];
        icontact=-2;
        if(i1elpr[jelem]==kprop)
        { icontact=-1;
          jt=jelem;
          for(it=0;it<2;it++)
          { for(in=0;in<3;in++)
            { i2to[it][in]=i2elto[in][jt];
            }
            jt=ielem; 
          }
          for(it=0;it<2;it++)
          { for(in=0;in<3;in++)
            { rx[it][in]=d1nccx[i2to[it][in]];
              ry[it][in]=d1nccy[i2to[it][in]];
          } }
          for(it=0;it<2;it++)
          { vol[it]=(rx[it][1]-rx[it][0])*(ry[it][2]-ry[it][0])-
                    (ry[it][1]-ry[it][0])*(rx[it][2]-rx[it][0]);
            for(ie=0;ie<3;ie++)
            { je=ie+1; if(je>2)je=0;
              nx[it][ie]=ry[it][je]-ry[it][ie];
              ny[it][ie]=rx[it][ie]-rx[it][je];
          } }
          for(it=0;it<2;it++)
          { jt=it+1; if(jt>1)jt=0;
            for(in=0;in<3;in++)
            { for(ie=0;ie<3;ie++)
              { d[it][in][ie]=((rx[jt][ie]-rx[it][in])*nx[jt][ie]+
                (ry[jt][ie]-ry[it][in])*ny[jt][ie])/vol[jt];
          } } }
          dmin2=big;
          for(it=0;it<2;it++)
          { jt=it+1; if(jt>1)jt=0;
            for(in=0;in<3;in++)
            { fx[in]=R0; fy[in]=R0;
            }
            n0=(nx[jt][0]*nx[jt][0]+ny[jt][0]*ny[jt][0])/
               (vol[jt]*vol[jt]);
            n1=(nx[jt][1]*nx[jt][1]+ny[jt][1]*ny[jt][1])/
               (vol[jt]*vol[jt]);
            n2=(nx[jt][2]*nx[jt][2]+ny[jt][2]*ny[jt][2])/
               (vol[jt]*vol[jt]);
            for(in=0;in<3;in++)
            { jn=in+1; if(jn>2)jn=0;
              a0=d[it][in][0];
              a1=d[it][in][1];
              a2=d[it][in][2];
              b0=d[it][jn][0];
              b1=d[it][jn][1];
              b2=d[it][jn][2];
              c0=d[jt][0][in];
              c1=d[jt][1][in];
              c2=d[jt][2][in];
              /* check if contact */
              if((((c0>nsmall)&&(c1>nsmall)&&(c2>nsmall))||
                  ((c0<small)&&(c1<small)&&(c2<small)))||
                 (((a0<small)&&(b0<small))||((a1<small)&&(b1<small))||
                  ((a2<small)&&(b2<small))))
              { if((a0<=a1)&&(a0<=a2))
                { dmin2=MINIM(dmin2,(a0*a0/n0));
                }
                else if((a1<=a0)&&(a1<=a2))
                { dmin2=MINIM(dmin2,(a1*a1/n1));
                }
                else
                { dmin2=MINIM(dmin2,(a2*a2/n2));
                }
              }
              else
              { icontact=it;
                /* domain of contact */
                smin=R0; smax=R1;
                if((a0<R0)&&(b0>small))smin=MAXIM(smin,(a0/(a0-b0)));
                if((a1<R0)&&(b1>small))smin=MAXIM(smin,(a1/(a1-b1)));
                if((a2<R0)&&(b2>small))smin=MAXIM(smin,(a2/(a2-b2)));
                if((a0>small)&&(b0<R0))smax=MINIM(smax,(a0/(a0-b0)));
                if((a1>small)&&(b1<R0))smax=MINIM(smax,(a1/(a1-b1)));
                if((a2>small)&&(b2<R0))smax=MINIM(smax,(a2/(a2-b2)));
                if(smax>smin)
                { s[0]=smin;
                  p[0]=MINIM((a0+smin*(b0-a0)),(a1+smin*(b1-a1)));
                  p[0]=MINIM(p[0],(a2+smin*(b2-a2)));
                  np=1;
                  /* intermediate points */
                  tmp=b0-a0+a1-b1;
                  if((DABS(tmp))>small)
                  { tmp=(a1-a0)/tmp;
                    if((tmp>smin)&&(tmp<smax)&&
                       ((a0+tmp*(b0-a0))<(a2+tmp*(b2-a2))))
                    { s[np]=tmp;
                      p[np]=a0+tmp*(b0-a0);
                      np=np+1;
                  } }
                  tmp=b0-a0+a2-b2;
                  if((DABS(tmp))>small)
                  { tmp=(a2-a0)/tmp;
                    if((tmp>smin)&&(tmp<smax)&&
                       ((a0+tmp*(b0-a0))<(a1+tmp*(b1-a1))))
                    { s[np]=tmp;
                      p[np]=a0+tmp*(b0-a0);
                      np=np+1;
                  } }
                  tmp=b1-a1+a2-b2;
                  if((DABS(tmp))>small)
                  { tmp=(a2-a1)/tmp;
                    if((tmp>smin)&&(tmp<smax)&&
                       ((a1+tmp*(b1-a1))<(a0+tmp*(b0-a0))))
                    { s[np]=tmp;
                      p[np]=a1+tmp*(b1-a1);
                      np=np+1;
                  } }
                  s[np]=smax;
                  p[np]=MINIM((a0+smax*(b0-a0)),(a1+smax*(b1-a1)));
                  p[np]=MINIM(p[np],(a2+smax*(b2-a2)));
                  np=np+1;
                  /* order intermediate points */
                  for(ip=0;ip<(np-1);ip++)
                  { for(jp=(ip+1);jp<np;jp++)
                    { if(s[ip]>s[jp])
                      { tmp=s[jp]; s[jp]=s[ip]; s[ip]=tmp;
                        tmp=p[jp]; p[jp]=p[ip]; p[ip]=tmp;
                  } } }
                  /* integrate normal force */
                  fn=p[0]*(s[1]-s[0])+p[np-1]*(s[np-1]-s[np-2]);
                  fnb=p[0]*(s[1]-s[0])*(s[1]+R2*s[0])+
                      p[np-1]*(s[np-1]-s[np-2])*(s[np-2]+R2*s[np-1]);
                  for(ip=1;ip<(np-1);ip++)
                  { fn=fn+p[ip]*(s[ip+1]-s[ip-1]);
                    fnb=fnb+p[ip]*(
                    (s[ip]-s[ip-1])*(s[ip-1]+R2*s[ip])+
                    (s[ip+1]-s[ip])*(s[ip+1]+R2*s[ip]));
                  }
                  fnb=fnb*pen*RP5;
                  fn=fn*pen*RP15;
                  fna=fn-fnb;
                  /* update total force */
                  fx[in]=fx[in]-fna*nx[it][in];
                  fy[in]=fy[in]-fna*ny[it][in];
                  fx[jn]=fx[jn]-fnb*nx[it][in];
                  fy[jn]=fy[jn]-fnb*ny[it][in];               
            } } }
            if(icontact==it) /* update nodal forces  */
            { for(in=0;in<3;in++)
              { d1nfcx[i2to[it][in]]=d1nfcx[i2to[it][in]]+fx[in];
                d1nfcy[i2to[it][in]]=d1nfcy[i2to[it][in]]+fy[in];
                ie=in+1; if(ie>2)ie=0;
                for(jn=0;jn<3;jn++)
                { d1nfcx[i2to[jt][in]]=d1nfcx[i2to[jt][in]]-
                  fx[jn]*d[it][jn][ie];
                  d1nfcy[i2to[jt][in]]=d1nfcy[i2to[jt][in]]-
                  fy[jn]*d[it][jn][ie];
          } } } }
        }
        /* remove the couple if too far from each other  */
        if((icontact==(-1))&&(dmin2>zone2))
        { if(jcoup<0)
          { i1elcf[ielem]=i1iecn[icoup];
            i1iecn[icoup]=*i1icff;
            *i1icff=icoup;
            icoup=i1elcf[ielem];
          }
          else
          { i1iecn[jcoup]=i1iecn[icoup];
            i1iecn[icoup]=*i1icff;
            *i1icff=icoup;           
            icoup=i1iecn[jcoup];       
        } } 
        else 
        { jcoup=icoup;
          icoup=i1iecn[icoup];
  } } } }
}


static int S_C_points(DBL dct[4],DBL us[6],DBL vs[6], int inners[6], DBL ut[4], DBL vt[4])
{
    int nspoin,inext,itag,i,j,k;
    DBL factor,tmp;
    /* u,v coordinates of S-points and C-points   */
    nspoin=0;
    itag = 0;
    for(i=0;i<3;i++)
    { for(j=0;j<2;j++)
    { inext=i+1; if(inext>2)inext=0; if(j==0)inext=3;
        				  if(((dct[i]>EPSILON1)&&(dct[inext]<NEPSILON1))||
         ((dct[i]<NEPSILON1)&&(dct[inext]>EPSILON1)))
                           ///Modified by JXiang
              { factor=ABS(dct[i]-dct[inext]);
                if(factor>EPSILON1)
                { factor=ABS(dct[i]/factor);
        if ((factor > EPSILON1) && (factor < (R1 - EPSILON1)))
        {
                  us[nspoin]=factor*ut[inext]+(R1-factor)*ut[i];
                  vs[nspoin]=factor*vt[inext]+(R1-factor)*vt[i];
                  inners[nspoin]=0;
                  nspoin=nspoin+1;
                    if(nspoin==2) itag=i;
                    }
          } } } }
    if((nspoin<3)||(nspoin>4)) return nspoin;
    /* check odering of S-points */
    if(nspoin==4 && itag == 2)    /*Zhou LEI (alex) 20/05/2010 21:50:07*/
    { i=2; j=3;
        k=inners[i]; inners[i]=inners[j]; inners[j]=k;
        tmp=us[i];  us[i]=us[j];  us[j]=tmp;
        tmp=vs[i];  vs[i]=vs[j];  vs[j]=tmp;
        printf("%s", "actively.\n");
    }
    
    if(((us[1]-us[0])*(vs[2]-vs[0])-
        (vs[1]-vs[0])*(us[2]-us[0]))<R0)
    { i=0; j=nspoin-1;
        while(i<j)
        { k=inners[i]; inners[i]=inners[j]; inners[j]=k;
            tmp=us[i];   us[i]=us[j];         us[j]=tmp;
            tmp=vs[i];   vs[i]=vs[j];         vs[j]=tmp;
            i++; j--;
        } }
    
    return nspoin;
}

static int CtoSdistance(DBL dcs[3][6],DBL dsc[6][3],DBL uc[4],DBL vc[4],
                        DBL us[6],DBL vs[6],DBL ub[10],DBL vb[10],int nspoin,int inners[6])
{
    int ninerc,niners,nbpoin,innerc[3],i,j,inext,jnext;
    DBL anb[10],tmp,factor;
/* distances of C-points from S edges */

niners=0;
ninerc=0;
innerc[0]=0;
innerc[1]=0;
innerc[2]=0;
for(i=0;i<nspoin;i++)
{ inext=i+1;
    if(inext>=nspoin)inext=0;
        for(j=0;j<3;j++)
        { jnext=j+1;
            if(jnext>2)jnext=0;
            dcs[j][i]=(uc[jnext]-uc[j])*(vs[i]-vc[j])-
            (vc[jnext]-vc[j])*(us[i]-uc[j]);
            dsc[i][j]=(us[inext]-us[i])*(vc[j]-vs[i])-
            (vs[inext]-vs[i])*(uc[j]-us[i]);
            if(dsc[i][j]>=EPSILON1)
            { innerc[j]=innerc[j]+1;
                if(innerc[j]==nspoin)ninerc=ninerc+1;
            }
            if(dcs[j][i]>=EPSILON1)
            { inners[i]=inners[i]+1;
                if(inners[i]==3)niners=niners+1;
            } } }
/* B-points */
if(ninerc==3)           /* triangle inside poligon      */
{ nbpoin=3;
    for(i=0;i<nbpoin;i++)
    { ub[i]=uc[i]; vb[i]=vc[i];
    } }
else if(niners==nspoin) /* poligon inside triangle      */
{ nbpoin=nspoin;
    for(i=0;i<nbpoin;i++)
    { ub[i]=us[i]; vb[i]=vs[i];
    } }
else            /* intersection points poligon triangle */
{ nbpoin=0;
    for(i=0;i<nspoin;i++)
    { if(inners[i]==3)
    { ub[nbpoin]=us[i]; vb[nbpoin]=vs[i]; nbpoin++;
    } }
    for(i=0;i<3;i++)  /* grab inner C-points */
    { if(innerc[i]==nspoin)
    { ub[nbpoin]=uc[i]; vb[nbpoin]=vc[i]; nbpoin++;
    } }
    for(i=0;i<nspoin;i++)        /* intersection points   */
    { inext=i+1; if(inext>=nspoin)inext=0;
        for(j=0;j<3;j++)
        {jnext=j+1; if(jnext>2)jnext=0;
            if((((dsc[i][j]>EPSILON1)&&(dsc[i][jnext]<NEPSILON1))||
             ((dsc[i][j]<NEPSILON1)&&(dsc[i][jnext]>EPSILON1)))&&
             (((dcs[j][i]>EPSILON1)&&(dcs[j][inext]<NEPSILON1))||
             ((dcs[j][i]<NEPSILON1)&&(dcs[j][inext]>EPSILON1))))
             { factor=ABS(dsc[i][j]-dsc[i][jnext]);
             if(factor<EPSILON1){ factor=RP5;                  }
             else              { factor=ABS(dsc[i][j]/factor);}
             ub[nbpoin]=(R1-factor)*uc[j]+factor*uc[jnext];
             vb[nbpoin]=(R1-factor)*vc[j]+factor*vc[jnext];
             nbpoin++;
             } } }  

    for(i=1;i<nbpoin;i++)
    { if(vb[i]<vb[0])
    { tmp=vb[i]; vb[i]=vb[0]; vb[0]=tmp;
        tmp=ub[i]; ub[i]=ub[0]; ub[0]=tmp;
    } }
    for(i=1;i<nbpoin;i++)
    { tmp=ub[i]-ub[0];
        if((tmp<R0)&&(tmp>(-EPSILON1)))
        { tmp=tmp-EPSILON1;
        }
        else if((tmp>=R0)&&(tmp<EPSILON1))
        { tmp=tmp+EPSILON1;
        }
        anb[i]=(vb[i]-vb[0]+EPSILON1)/tmp;
    }
    for(i=1;i<nbpoin;i++)  /* sort B-points */
    { for(j=i+1;j<nbpoin;j++)
    { if(((anb[i]>=R0)&&(anb[j]>=R0)&&(anb[j]<anb[i]))||
         ((anb[i]<R0)&&((anb[j]>=R0)||(anb[j]<anb[i]))))
    { tmp=vb[i];  vb[i]=vb[j];   vb[j]=tmp;
        tmp=ub[i];  ub[i]=ub[j];   ub[j]=tmp;
        tmp=anb[i]; anb[i]=anb[j]; anb[j]=tmp;
    } } } }

return nbpoin;

}

static DBL PenforceN(DBL xc[4],DBL yc[4],DBL zc[4],DBL xt[4],DBL yt[4],DBL zt[4],DBL xe[4],DBL ye[4],DBL ze[4],DBL fforco[3],int nbpoin,DBL ub[10],DBL vb[10],DBL penalty)
{
    DBL theigh,penetr,peneto,penetu,penetv,xnt,ynt,znt;
    DBL fact0,facti,fact1,force,penetb[10];
    int i;
    /* Target-plain normal and penetration at B-points */
    V3DCro(xnt,ynt,znt,xt[1]-xt[0],yt[1]-yt[0],zt[1]-zt[0],
           xt[2]-xt[0],yt[2]-yt[0],zt[2]-zt[0]);
    V3DDot(theigh,xt[3]-xt[0],
           yt[3]-yt[0],zt[3]-zt[0],xnt,ynt,znt);
    /* penetration at origin of the e-base and dp/du dp/dv; */
    V3DDot(peneto,xc[0]-xt[0],yc[0]-yt[0],
           zc[0]-zt[0],xnt,ynt,znt)
    V3DDot(penetu,xe[0],ye[0],ze[0],xnt,ynt,znt);
    V3DDot(penetv,xe[1],ye[1],ze[1],xnt,ynt,znt);
    peneto=peneto/theigh;
    penetu=penetu/theigh;
    penetv=penetv/theigh;
    for(i=0;i<nbpoin;i++)
    { penetb[i]=peneto+ub[i]*penetu+vb[i]*penetv;
    }
    /* force and center of force */
    fforco[0]=R0; fforco[1]=R0; fforco[2]=R0;
    for(i=1;i<(nbpoin-1);i++)
    { penetr=penetb[0]+penetb[i]+penetb[i+1];
        if(penetr>EPSILON1){
            force=((ub[i]-ub[0])*(vb[i+1]-vb[0])-
                   (vb[i]-vb[0])*(ub[i+1]-ub[0]))*penetr*penalty;
            fact0=(RP5*penetb[0]+
                   RP25*(penetb[i]+penetb[i+1]))/penetr;
            facti=(RP5*penetb[i]+RP25*(penetb[0]+
                                       penetb[i+1]))/penetr;
            fact1=R1-fact0-facti;
            if(ABS(force+fforco[0])>EPSILON1)
            { fforco[1]=(fforco[0]*fforco[1]+force*(fact0*ub[0]+
                                        facti*ub[i]+fact1*ub[i+1]))/(fforco[0]+force);
                fforco[2]=(fforco[0]*fforco[2]+force*(fact0*vb[0]+
                                          facti*vb[i]+fact1*vb[i+1]))/(fforco[0]+force);
                fforco[0]=fforco[0]+force;
            } } }
    return theigh;
}
static void ForceN_cal(DBL ut[4],DBL vt[4],DBL uc[4],DBL vc[4],DBL fforco[3],
                       DBL fcn[4][3],DBL ftn[4][3],
                       DBL xnc,DBL ync,DBL znc,DBL factor1,DBL factor2,DBL theigh,
                       int ipcn1[4],int iptn1[4],DBL zt[4])
{
    int i,j,k,inext,jnext;
    DBL tmp,ft[4],fc[4];
    /*             resultant at C-points */
    for(i=0;i<4;i++)
    { fc[i]=R0; ft[i]=R0;
    }
    tmp=((uc[1]-uc[0])*(vc[2]-vc[0])-
         (vc[1]-vc[0])*(uc[2]-uc[0]));
    for(i=0;i<3;i++)
    { j=i+1; if(j>2)j=0; k=j+1; if(k>2)k=0;
        fc[k]=fforco[0]*(((uc[j]-uc[i])*(fforco[2]-vc[i])-
                      (vc[j]-vc[i])*(fforco[1]-uc[i]))/tmp);
    }
    /*             resultant at T-points  */
    tmp=((ut[1]-ut[0])*(vt[2]-vt[0])-
         (vt[1]-vt[0])*(ut[2]-ut[0]));
    inext=-1;
    if(ABS(tmp)<RP1*theigh)
    { inext=0; tmp=ABS(ut[1]-ut[0])+ABS(vt[1]-vt[0]);
        for(i=0;i<3;i++)
        { j=i+1; if(j>2)j=0;
            if(tmp>(ABS(ut[j]-ut[i])+ABS(vt[j]-vt[i])))
            { tmp=ABS(ut[j]-ut[i])+ABS(vt[j]-vt[i]);  inext=i;
            } }
        j=inext+1; if(j>2)j=0;
        if(ABS(zt[j])>ABS(zt[inext]))inext=j;
        j=inext+1; if(j>2)j=0;  k=j+1;  if(k>2)k=0;
        tmp=(ut[k]-ut[j])*(vt[3]-vt[j])-
        (vt[k]-vt[j])*(ut[3]-ut[j]);
    }
    for(jnext=0;jnext<3;jnext++)
    { i=jnext; j=i+1; if(j>2)j=0; k=j+1; if(k>2)k=0;
        if(i==inext)i=3; if(j==inext)j=3; if(k==inext)k=3;
        ft[k]=fforco[0]*(((ut[j]-ut[i])*(fforco[2]-vt[i])-
                      (vt[j]-vt[i])*(fforco[1]-ut[i]))/tmp);
    }
    ft[3]=RP25*ft[3];
    for(i=0;i<3;i++)
    { ft[i]=ft[i]+ft[3];
    }
    /* add forces into global vector */
    
    for(i=0;i<4;i++)
    {   fcn[ipcn1[i]][0]+=fc[i]*xnc*factor1;
        fcn[ipcn1[i]][1]+=fc[i]*ync*factor1;
        fcn[ipcn1[i]][2]+=fc[i]*znc*factor1;
        ftn[iptn1[i]][0]-=ft[i]*xnc*factor2;
        ftn[iptn1[i]][1]-=ft[i]*ync*factor2;
        ftn[iptn1[i]][2]-=ft[i]*znc*factor2;
        
    }
}


static int force_normal(DBL fcn[4][3],DBL ftn[4][3],DBL *d1nccx,DBL  *d1nccy, DBL  *d1nccz,
            INT *i1celto, INT *i1telto,DBL factor1, DBL factor2, DBL penalty)
{
    DBL xt[4],yt[4],zt[4],ut[4],vt[4],xcent,ycent,zcent;
    DBL xc[4],yc[4],zc[4],uc[4],vc[4],xcenc,ycenc,zcenc,xnc,ync,znc;
    DBL xorig,yorig,zorig,xe[2],ye[2],ze[2],dct[4];
    INT ipt[4],ipc[4];
    int itars,icons,iptn1[4],ipcn1[4],i,inners[6];
    int p1[4] = { 1, 3, 3, 1 };
    int p2[4] = { 2, 2, 0, 0 };
    int p3[4] = { 3, 0, 1, 2 };
    DBL us[6],vs[6];
    DBL dcs[3][6],dsc[6][3];
    DBL ub[10],vb[10],fforco[3];
    
    DBL theigh;
    int nspoin,nbpoin,fnonzero;

    
    /*set centres of contactor and target object */
    xcent=R0; ycent=R0; zcent=R0; xcenc=R0; ycenc=R0; zcenc=R0;
    for(i=0;i<4;i++)
    { xcenc=xcenc+RP25*d1nccx[i1celto[i]];
        ycenc=ycenc+RP25*d1nccy[i1celto[i]];
        zcenc=zcenc+RP25*d1nccz[i1celto[i]];
        xcent=xcent+RP25*d1nccx[i1telto[i]];
        ycent=ycent+RP25*d1nccy[i1telto[i]];
        zcent=zcent+RP25*d1nccz[i1telto[i]];
    }
    fnonzero=0;
    
    /*********************************************************/
    /*                loop over target surfaces              */
    /*********************************************************/
    for(itars=0;itars<4;itars++)
    {
        iptn1[0] = itars;
        ipt[0] = i1telto[iptn1[0]];
        
        iptn1[1] = p1[itars];
        ipt[1] = i1telto[iptn1[1]];
        
        iptn1[2] = p2[itars];
        ipt[2] = i1telto[iptn1[2]];
        
        iptn1[3] = p3[itars];
        ipt[3] = i1telto[iptn1[3]];
        
		/*           loop over contactor surfaces              */
        for(icons=0;icons<4;icons++)
        {
            
            ipcn1[0] = icons;
            ipc[0] = i1celto[ipcn1[0]];
            
            ipcn1[1] = p1[icons];
            ipc[1] = i1celto[ipcn1[1]];
            
            ipcn1[2] = p2[icons];
            ipc[2] = i1celto[ipcn1[2]];
            
            ipcn1[3] = p3[icons];
            ipc[3] = i1celto[ipcn1[3]];
            
            /* set nodal coordinates */
            for(i=0;i<3;i++)
            { xt[i]=d1nccx[ipt[i]];
				yt[i]=d1nccy[ipt[i]];
				zt[i]=d1nccz[ipt[i]];
				xc[i]=d1nccx[ipc[i]];
				yc[i]=d1nccy[ipc[i]];
				zc[i]=d1nccz[ipc[i]];
            }
            xt[3]=xcent; yt[3]=ycent; zt[3]=zcent;
            xc[3]=xcenc; yc[3]=ycenc; zc[3]=zcenc;
            xorig=xc[0]; yorig=yc[0]; zorig=zc[0];
            for(i=0;i<4;i++)
            { xt[i]=xt[i]-xorig; yt[i]=yt[i]-yorig; zt[i]=zt[i]-zorig;
				xc[i]=xc[i]-xorig; yc[i]=yc[i]-yorig; zc[i]=zc[i]-zorig;
            }
            /* contactor normal, e-base and target points in e-base */
            V3DCro(xnc,ync,znc,xc[1],yc[1],zc[1],xc[2],yc[2],zc[2]);
            V3DNor(xe[0],xnc,ync,znc);
            xe[0]=xc[1]; ye[0]=yc[1]; ze[0]=zc[1];
            V3DNor(xe[1],xe[0],ye[0],ze[0]);
            V3DCro(xe[1],ye[1],ze[1],xnc,ync,znc,xe[0],ye[0],ze[0]);
            for(i=0;i<4;i++)
            { V3DDot(dct[i],xnc,ync,znc,xt[i],yt[i],zt[i]);
				V3DDot(ut[i],xt[i],yt[i],zt[i],xe[0],ye[0],ze[0]);
				V3DDot(vt[i],xt[i],yt[i],zt[i],xe[1],ye[1],ze[1]);
            }
            			  if((dct[0]<=EPSILON1)&&(dct[1]<=EPSILON1)&&
            				  (dct[2]<=EPSILON1))continue;  //modified by JXiang
            
            /* u,v coordinates of S-points and C-points   */
            nspoin=S_C_points(dct,us,vs,inners,ut,vt);
 //           if(index==0)continue;
            if((nspoin<3)||(nspoin>4)) continue;
            for(i=0;i<3;i++)
            { V3DDot(uc[i],xc[i],yc[i],zc[i],xe[0],ye[0],ze[0]);
                V3DDot(vc[i],xc[i],yc[i],zc[i],xe[1],ye[1],ze[1]);
            }
            
            nbpoin=CtoSdistance(dcs,dsc,uc,vc,us,vs,ub,vb,nspoin,inners);
            if(nbpoin<3)continue;
            
            theigh=PenforceN(xc,yc,zc,xt,yt,zt,xe,ye,ze,fforco,nbpoin,ub,vb,penalty);
            ForceN_cal(ut,vt,uc,vc,fforco,fcn,ftn,xnc,ync,znc,factor1,factor2,
                       theigh,ipcn1,iptn1,zt);
            
            fnonzero=1;
            
        }
    }
    return fnonzero;
}


static void orientation(DBL  *d1deltat1, DBL *d1deltat2,DBL *d1deltan,
                        DBL  *d1nvx, DBL  *d1nvy, DBL  *d1nvz,
                        DBL  *d1t1vx, DBL  *d1t1vy, DBL  *d1t1vz,
                        DBL  *d1t2vx, DBL  *d1t2vy, DBL  *d1t2vz,
                        INT *i1stat, DBL nx,DBL ny,DBL nz)
{
    DBL pnx,pny,pnz,cosin,t1x,t1y,t1z,t2x,t2y,t2z,sine,tmp;
    if((*i1stat)==12)
    {
        if(DABS(((*d1nvx)*nx+(*d1nvy)*ny+(*d1nvz)*nz)-R1)>1.e-03)
        {
            V3DCro(pnx,pny,pnz,*d1nvx,*d1nvy,*d1nvz,
                   nx,ny,nz);
            V3DNor(tmp,pnx,pny,pnz);
            V3DDot(cosin,nx,ny,nz,*d1nvx,*d1nvy,*d1nvz);
            tmp=R1-cosin*cosin;
            if(tmp<R0)
                sine=R0;
            else
                sine=SQRT(R1-cosin*cosin);
            V3DRot(&t1x,&t1y,&t1z,*d1t1vx,*d1t1vy,*d1t1vz,
                   pnx,pny,pnz,cosin,sine);
            V3DRot(&t2x,&t2y,&t2z,*d1t2vx,*d1t2vy,*d1t2vz,
                   pnx,pny,pnz,cosin,sine);
            *d1nvx=nx;
            *d1nvy=ny;
            *d1nvz=nz;
            *d1t1vx=t1x;
            *d1t1vy=t1y;
            *d1t1vz=t1z;
            *d1t2vx=t2x;
            *d1t2vy=t2y;
            *d1t2vz=t2z;
            
        }
    }
    else
    {
        *i1stat=12;
        *d1deltat1=R0;
        *d1deltat2=R0;
        *d1deltan=R0;
        
        if( (ABS(nx)>=ABS(ny)) && (ABS(nx)>=ABS(nz)) )
        { t1y = t1z = 1;
            t1x = - (ny+nz)/nx; }
        else if( (ABS(ny)>=ABS(nx)) && (ABS(ny)>=ABS(nz)) )
        { t1x = t1z =1;
            
            t1y = - (nx+nz) / ny;
        }
        else
            
        { t1x = t1y =1; t1z = - (nx+ny) / nz;
        }
        V3DNor(tmp,t1x,t1y,t1z);
        V3DCro(t2x,t2y,t2z,nx,ny,nz,t1x,t1y,t1z);
        *d1nvx=nx;
        *d1nvy=ny;
        *d1nvz=nz;
        *d1t1vx=t1x;
        *d1t1vy=t1y;
        *d1t1vz=t1z;
        *d1t2vx=t2x;
        *d1t2vy=t2y;
        *d1t2vz=t2z;
    }
    
}

static void force_friction_cal(DBL mud,DBL mus,DBL muss,DBL s_vel,
                DBL penalty,DBL *deltat1,DBL *deltat2,DBL *deltan,
                DBL fcn[4][3],DBL ftn[4][3],DBL fn_normal,
                DBL nx, DBL ny, DBL nz, DBL t1x, DBL t1y,DBL t1z,
                DBL t2x, DBL t2y,DBL t2z, DBL ru, DBL rv,DBL rw,
                DBL factor1, DBL factor2,DBL dcstec,DBL fnc[4],
                DBL fnt[4],DBL fnctot,DBL fnttot,
                INT ncstep, INT ifcstep, INT ifcjp,INT icoup, DBL *d1pepe,
                DBL **d2elfc,int *index,DBL d_fact)
{
    DBL mudd,h_cr,ktss,vs1,vs2,vn,sqvel;
    DBL fs,fs1,fs2,mu,ftt[3],sfactor,fn_jp,factor_jp;
    int i,j;
    FILE *fout=FILENULL;


    if((ifcstep==ncstep)&&(ifcjp>=0))
    {
        V3DDot(fn_jp,nx,ny,nz,d2elfc[0][ifcjp],d2elfc[1][ifcjp],d2elfc[2][ifcjp]);
        V3DDot(fs1,t1x,t1y,t1z,d2elfc[0][ifcjp],d2elfc[1][ifcjp],d2elfc[2][ifcjp]);
        V3DDot(fs2,t2x,t2y,t2z,d2elfc[0][ifcjp],d2elfc[1][ifcjp],d2elfc[2][ifcjp]);

        if(fn_jp/fn_normal<R0)
        {
 
 if(fn_normal>2000.0||-fn_jp>2000)
 {
	     fout=fopen("gidmsh1.txt","w");
 fprintf(fout,"fn_jp is %f, %f, %f, %f \n", fn_jp, d2elfc[0][ifcjp],d2elfc[1][ifcjp],d2elfc[2][ifcjp]);
fprintf(fout," fn is %f, %f, %f, %f \n", fn_normal, nx, ny,nz);
fprintf(fout," joint no is %ld \n",ifcjp);
fclose(fout);

*index=1;
 }
            fn_jp=-fn_jp;
            factor_jp=fn_jp/fn_normal;
            if(factor_jp>3.0)factor_jp=3.0;
             else if(factor_jp<0.01)factor_jp=0.01;
            d1pepe[icoup]=penalty*factor_jp;
              
            ktss=2.0/(7.0)*d1pepe[icoup];
	    
            (*deltat1)=-fs1/ktss;
            (*deltat2)=-fs2/ktss;
	    
	    *deltat1=R0;
	    *deltat2=R0;
	    
	}
        else
        {
		if(fn_normal>2000.0||fn_jp>2000)
 
		{
			
             fout=fopen("gidmsh2.txt","w");
 fprintf(fout,"fn_jp is %f, %f, %f, %f \n", fn_jp, d2elfc[0][ifcjp],d2elfc[1][ifcjp],d2elfc[2][ifcjp]);
fprintf(fout," fn is %f, %f, %f, %f \n", fn_normal, nx, ny,nz);
fclose(fout);

 }
            factor_jp=fn_jp/fn_normal;
            if(factor_jp>3.0)factor_jp=3.0;
             else if(factor_jp<0.01)factor_jp=0.01;
            d1pepe[icoup]=penalty*factor_jp;
            ktss=2.0/(7.0)*d1pepe[icoup];

	    (*deltat1)=fs1/ktss;
            (*deltat2)=fs2/ktss;
	    
	    *deltat1=R0;
	    *deltat2=R0;

        }
        sqvel=SQRT(ru*ru+rv*rv+rw*rw);
        vs1=-rw*t1z-rv*t1y-ru*t1x;
        vs2=-rw*t2z-rv*t2y-ru*t2x;
        
        vn=-rw*nz-rv*ny-ru*nx;
        (*deltat1)+=vs1*dcstec;
        (*deltat2)+=vs2*dcstec;
        (*deltan)+=vn*dcstec;
        fs1=ktss*(*deltat1);
        fs2=ktss*(*deltat2);
        fs=fs1*fs1+fs2*fs2;
        fs=SQRT(fs);
	/*
        mudd=mud;
        if((DABS(*deltan)>h_cr)&&(DABS(muss-mud)>EPSILON))
            mudd=mud+(muss-mud)*pow(DABS(*deltan)/h_cr,-2);
        if(DABS(mus-R0)>EPSILON)
            mu=mudd+(mus-mudd)*exp(-s_vel*sqvel);
        else
            mu=mudd;
	    */
	mu=mud*d_fact;
        if((fs>(mu*fn_jp))&&(fs>EPSILON1))
        {sfactor=mu*fn_jp/(fs+EPSILON1);
            fs1=fs1*sfactor;
            fs2=fs2*sfactor;
            *deltat1=fs1/ktss;
            *deltat2=fs2/ktss;
        }
        
    
        for(i=0;i<4;i++)
        {
            for(j=0;j<3;j++)
            {
            fcn[i][j]=R0;
            ftn[i][j]=R0;
            }
        }
	
    }
    else
    {
        ktss=2.0/(7.0)*penalty;
        h_cr=mus/penalty;
        fn_jp=fn_normal;
    
    sqvel=SQRT(ru*ru+rv*rv+rw*rw);
    vs1=-rw*t1z-rv*t1y-ru*t1x;
    vs2=-rw*t2z-rv*t2y-ru*t2x;
    
    vn=-rw*nz-rv*ny-ru*nx;
    (*deltat1)+=vs1*dcstec;
    (*deltat2)+=vs2*dcstec;
    (*deltan)+=vn*dcstec;
//    h_cr=mus/penalty;
//    ktss=2.0/(7.0)*penalty;
    fs1=ktss*(*deltat1);
    fs2=ktss*(*deltat2);
    fs=fs1*fs1+fs2*fs2;
    fs=SQRT(fs);
    /*
    mudd=mud;
    if((DABS(*deltan)>h_cr)&&(DABS(muss-mud)>EPSILON))
        mudd=mud+(muss-mud)*pow(DABS(*deltan)/h_cr,-2);
    if(DABS(mus-R0)>EPSILON)
        mu=mudd+(mus-mudd)*exp(-s_vel*sqvel);
    else
        mu=mudd;
	*/
    mu=mud*d_fact;
    if((fs>(mu*fn_normal))&&(fs>EPSILON1))
    {sfactor=mu*fn_normal/(fs+EPSILON1);
        fs1=fs1*sfactor;
        fs2=fs2*sfactor;
 
 	*deltat1=fs1/ktss;
        *deltat2=fs2/ktss;
	
//	*deltat1=R0;
//	*deltat2=R0;
    }
    
//    d1deltat1[icoup]=deltat1;
//    d1deltat2[icoup]=deltat2;
    
    ftt[0]=fs1*t1x+fs2*t2x;
    ftt[1]=fs1*t1y+fs2*t2y;
    ftt[2]=fs1*t1z+fs2*t2z;
    
    for(i=0;i<4;i++)
    {
        for(j=0;j<3;j++)
        {
        fcn[i][j]=fcn[i][j]-ftt[j]*fnc[i]/fnctot*factor1;
        ftn[i][j]=ftn[i][j]+ftt[j]*fnt[i]/fnttot*factor2;
        }
    }
}
}



static int force_friction(DBL fcn[4][3],DBL ftn[4][3],int iquad,INT *i1celto, INT *i1telto,
                           DBL factor1, DBL factor2, DBL penalty,DBL  *d1deltat1, DBL *d1deltat2,
                           DBL *d1deltan,DBL *d1nvx,DBL *d1nvy,DBL *d1nvz,DBL  *d1t1vx, DBL  *d1t1vy, DBL  *d1t1vz,
                           DBL  *d1t2vx, DBL  *d1t2vy, DBL  *d1t2vz, INT *i1stat,
                           DBL dpefr,DBL  dpesf, DBL  dpevf, DBL  dpepf, DBL  dpepsf,DBL dcstec,
                           DBL *d1nvcx,DBL *d1nvcy,DBL *d1nvcz,INT ncstep, INT ifcstep, INT ifcjp,INT icoup,
                          DBL *d1pepe,DBL **d2elfc,DBL d_fact)

{
    DBL fnc[4],fnt[4],fn[3],fnctot,fnttot,fn_normal;
    DBL ru,rv,rw,vxc,vyc,vzc,vxt,vyt,vzt;
    int i,j,index;
    index=0;
    
    fn[0]=R0;fn[1]=R0;fn[2]=R0;fnctot=R0;fnttot=R0; fn_normal=R0;
    
    for(i=0;i<4;i++)
    {
        fnc[i]=R0;
        fnt[i]=R0;
        for(j=0;j<3;j++)
        {
        fnc[i]+=fcn[i][j]*fcn[i][j];
        fnt[i]+=ftn[i][j]*ftn[i][j];
        fn[j]+=fcn[i][j];
        }
        fnc[i]=SQRT(fnc[i]);
        fnt[i]=SQRT(fnt[i]);
        fnctot=fnctot+fnc[i];
        fnttot=fnttot+fnt[i];
    }
    if(fnctot>EPSILON1)
    {
        vxc=R0;vyc=R0;vzc=R0;vxt=R0;vyt=R0;vzt=R0;
        V3DNor(fn_normal,fn[0],fn[1],fn[2]);
        if(iquad)fn_normal=fn_normal*2.5;
        
        for(i=0;i<4;i++)
        {
            vxc=vxc+fnc[i]/fnctot*d1nvcx[i1celto[i]];
            vyc=vyc+fnc[i]/fnctot*d1nvcy[i1celto[i]];
            vzc=vzc+fnc[i]/fnctot*d1nvcz[i1celto[i]];
            vxt=vxt+fnt[i]/fnttot*d1nvcx[i1telto[i]];
            vyt=vyt+fnt[i]/fnttot*d1nvcy[i1telto[i]];
            vzt=vzt+fnt[i]/fnttot*d1nvcz[i1telto[i]];
        }
        ru=vxt-vxc;
        rv=vyt-vyc;
        rw=vzt-vzc;
        
        orientation(d1deltat1, d1deltat2,d1deltan,
                    d1nvx,d1nvy,d1nvz,
                    d1t1vx,d1t1vy,d1t1vz,
                    d1t2vx, d1t2vy, d1t2vz,i1stat,fn[0],fn[1],fn[2]);
        
        force_friction_cal(dpefr,dpesf,dpepf,dpevf,
                           penalty,d1deltat1,d1deltat2,d1deltan,fcn,ftn,fn_normal,
                           *d1nvx,*d1nvy,*d1nvz,*d1t1vx,*d1t1vy,*d1t1vz,
                           *d1t2vx, *d1t2vy, *d1t2vz,ru,rv,rw,
                           factor1,factor2,dcstec,fnc,fnt,fnctot,fnttot,
                           ncstep,ifcstep,ifcjp,icoup,d1pepe,d2elfc,&index,d_fact);
        /* end of calculating sliding frictional force  */
    }
    
    if(fn_normal>2.0e+3)return 1;
    else return 0;
    
}


static void Yid3TET2TET(/* tetrahedra to tetrahedra */
						nelem,i1elem,
            dcstec,
            d1nccx,d1nccy,d1nccz,
			d2nfc,
            d1nvcx,d1nvcy,d1nvcz,
			d1nvx,d1nvy,d1nvz,
			d1deltat1,d1deltat2,d1deltan,
			d1t1vx,d1t1vy,d1t1vz,
			d1t2vx,d1t2vy,d1t2vz,
			d1pepe,i1elcf,i1elpr,d1pefr,
			d1pesf,d1pevf,d1pepf,d1pepsf,
            i1iecn,i1iect,i2elto,i1stat,iquad,d1mart,d1emct,
            ncstep,i1fcstep,i1fcjp,d1jppepe,d2elfc,i1elst,d1df
//            nthreads,nblock,nnopo,nmode
            )
INT    nelem; INT *i1elem;
  DBL   dcstec;
  DBL  *d1nccx; DBL  *d1nccy; DBL  *d1nccz;
DBL  **d2nfc;
  DBL  *d1nvcx; DBL  *d1nvcy; DBL  *d1nvcz;
  DBL   *d1nvx; DBL   *d1nvy; DBL   *d1nvz;
  DBL  *d1deltat1;  DBL   *d1deltat2; DBL   *d1deltan; 
  DBL  *d1t1vx; DBL  *d1t1vy; DBL  *d1t1vz; 
  DBL  *d1t2vx; DBL  *d1t2vy; DBL  *d1t2vz; DBL *d1pefr; 
  DBL  *d1pesf; DBL  *d1pevf; DBL  *d1pepf; DBL  *d1pepsf; 
  DBL  *d1pepe; INT  *i1elcf; INT  *i1elpr;
  INT  *i1iecn; INT  *i1iect; INT **i2elto;
//INT **i2elsur; INT *i1nelsuf;
INT *i1stat; int iquad; DBL *d1mart; DBL *d1emct;
INT ncstep; INT *i1fcstep; INT *i1fcjp; DBL *d1jppepe;
DBL **d2elfc; INT *i1elst; DBL *d1df;
//int nthreads; INT nblock; INT nnopo; INT nmode;
  {

      int ij[10] = { -1, -1, -1, -1, 0, 1, 2, 0, 1, 2 };
      int ik[10] = { -1, -1, -1, -1, 1, 2, 0, 3, 3, 3 };

#pragma omp parallel num_threads(1) default(none) shared(nelem,dcstec,d1nccx,d1nccy,d1nccz,d2nfc,d1df,\
d1nvcx,d1nvcy,d1nvcz,d1nvx,d1nvy,d1nvz,d1deltat1,d1deltat2,d1deltan,d1t1vx,d1t1vy,d1t1vz,\
d1t2vx,d1t2vy,d1t2vz,d1pepe,i1elcf,i1elpr,d1pefr,d1pesf,d1pevf,d1pepf,d1pepsf,i1iecn,i1iect,\
i2elto,iquad,i1stat,ij,ik,i1elem,d1mart,d1emct,ncstep,i1fcstep,i1fcjp,d1jppepe,d2elfc,i1elst,stderr)
//private(i1celto,i1telto,fcn,ftn,i1celto,i1telto,penalty,icoup,tid,fnonzero,i,ielem,index)
//#endif
      {
          INT *i1celto, *i1telto;
          DBL fcn[4][3],ftn[4][3];
          DBL factor1,factor2, penalty,dcstec_el,penalty_el;
          INT icoup,ielem,iprop,jprop,jelem;
          int tid,fnonzero,i,index,j,k,jj,*i1nfs,ibig;
          INT ii;
	  DBL d_fact;

          
//#if defined(_OPENMP)
#pragma omp for schedule(dynamic,64) nowait
//#pragma omp for schedule(guided) nowait
//#endif
  for(ii=0;ii<nelem;ii++)
  {  ielem=i1elem[ii];
    icoup=i1elcf[ielem];
//      jcoup=-1;
      while(icoup>=0)
      { //ihycon=1;
		if(i1iect[icoup]<0)
	    {
	     i1iect[icoup]=-1-i1iect[icoup];
            
            if(i1stat[icoup]==1)
            {
//            ihycon=0;
            d1deltat1[icoup]=R0;
            d1deltat2[icoup]=R0;
            d1deltan[icoup]=R0;
            i1stat[icoup]=2;
            }
	    }
	    
        jelem=i1iect[icoup];
          index=0;
        if(i1stat[icoup]>=11)
        {
            index=1;
        }
        else
        {
            index=Ycd3TETintersection(d1nccx,d1nccy,d1nccz,ielem,jelem,i2elto);
            if(index==2)index=Ycd3TETintersection(d1nccx,d1nccy,d1nccz,jelem,ielem,i2elto);
        }

		fnonzero=0;

        if(index==1)
        {
            iprop=i1elpr[ielem];
            jprop=i1elpr[jelem];

          penalty=MINIM(d1pepe[iprop],d1pepe[jprop])/200.0;
	         penalty_el=0.25*d1emct[ielem]*d1emct[jelem]/
          ((d1emct[ielem]+d1emct[jelem])*dcstec*dcstec);
          penalty=MINIM(penalty,penalty_el);
            if((i1fcjp[icoup]>=0)&&(i1fcstep[icoup]!=ncstep))
	    {
		    if(ncstep<=(i1fcstep[icoup]+1000))
	    penalty_el=(ncstep-i1fcstep[icoup])*(penalty-d1jppepe[icoup])/1000.0+d1jppepe[icoup];
		              penalty=MINIM(penalty,penalty_el);

	    }
	    else if(i1fcjp[icoup]==-2)
	    {
             if(ncstep<=(i1fcstep[icoup]+100))
            penalty_el=(ncstep-i1fcstep[icoup])*(penalty)/100.0;
                              penalty=MINIM(penalty,penalty_el);

            }

	    d_fact=MINIM((R1-d1df[ielem]),(R1-d1df[jelem]));
	    if(d_fact<0.041)d_fact=d_fact/1000.0;
            penalty=penalty*d_fact;

          i1telto = i2elto[jelem];
          i1celto = i2elto[ielem];
            for(i=0;i<4;i++)
            {
                for(j=0;j<3;j++)
                {
                fcn[i][j]=R0;
                ftn[i][j]=R0;
            }
            }
            if(iquad)
            {
                factor2=R2/R5;
                factor1=R2/R5;

            }
            else
            {
                factor2=R1;
                factor1=R1;
            }
        fnonzero=force_normal(fcn,ftn,d1nccx,d1nccy, d1nccz,
                     i1celto,i1telto,factor1,factor2,penalty);
/* Calculating sliding frictional force  */
	 if(fnonzero==1)
	 {   if(d1pefr[iprop]>d1pefr[jprop])iprop=jprop;
         ibig=0;
         ibig=force_friction(fcn,ftn,iquad,i1celto,i1telto,factor1, factor2, penalty,
          &d1deltat1[icoup],&d1deltat2[icoup],&d1deltan[icoup],
                        &d1nvx[icoup], &d1nvy[icoup], &d1nvz[icoup],&d1t1vx[icoup],
          &d1t1vy[icoup], &d1t1vz[icoup],&d1t2vx[icoup], &d1t2vy[icoup], &d1t2vz[icoup],&i1stat[icoup],
        d1pefr[iprop],d1pesf[iprop],d1pevf[iprop], d1pepf[iprop], d1pepsf[iprop],dcstec,
            d1nvcx,d1nvcy,d1nvcz,ncstep,i1fcstep[icoup],i1fcjp[icoup],icoup,d1jppepe,d2elfc,d_fact);
         if(ibig)
         {
             gidmsh_interaction(i2elto,d1nccx, d1nccy, d1nccz,ielem,jelem);
//             CHRw(stderr,"Contact interaction");
//             exit(1);
         }
         
         for(i=0;i<4;i++)
         {
#pragma omp atomic
                 d2nfc[0][i1celto[i]]+=fcn[i][0];
#pragma omp atomic
                 d2nfc[1][i1celto[i]]+=fcn[i][1];
#pragma omp atomic                
                 d2nfc[2][i1celto[i]]+=fcn[i][2];
#pragma omp atomic
                 d2nfc[0][i1telto[i]]+=ftn[i][0];
#pragma omp atomic
                 d2nfc[1][i1telto[i]]+=ftn[i][1];
#pragma omp atomic
                 d2nfc[2][i1telto[i]]+=ftn[i][2];

         }
         if(iquad)
         {
             for(i=4;i<10;i++)
             {
                 j = ij[i]; 
                 k = ik[i];
#pragma omp atomic
                         d2nfc[0][i1celto[i]]+=(fcn[j][0]+fcn[k][0])/R2;
#pragma omp atomic
                         d2nfc[1][i1celto[i]]+=(fcn[j][1]+fcn[k][1])/R2;
#pragma omp atomic
                         d2nfc[2][i1celto[i]]+=(fcn[j][2]+fcn[k][2])/R2;
#pragma omp atomic
                         d2nfc[0][i1telto[i]]+=(ftn[j][0]+ftn[k][0])/R2;
#pragma omp atomic
                         d2nfc[1][i1telto[i]]+=(ftn[j][1]+ftn[k][1])/R2;
#pragma omp atomic
                         d2nfc[2][i1telto[i]]+=(ftn[j][2]+ftn[k][2])/R2;

                 }
                 
             
         }

     }
        }
   
          if(index==1)
          {
        if(fnonzero==0)

        {
            i1stat[icoup]=10;
            i1fcjp[icoup]=-1;
        }
        else
        {
            if(i1stat[icoup]<11)i1stat[icoup]=11;
        }
          }
          else
          {
              i1fcjp[icoup]=-1;
          }
//              jcoup=icoup;
              icoup=i1iecn[icoup];
          
	  }
    }
          
  }
      
  }


/*********************PUBLIC***********************/
void Yid(YDE  yde, YDI ydi, YDN ydn, YDP ydp, YDJ ydj,int nthreads,DBL dcstec, INT ncstep     /***  nodal forces  ***/
        )
//  YDC ydc; YDE yde; YDI ydi; YDN ydn; YDP ydp;
{
    int iquad;
    INT i;
    INT nblock,nmode;
if(ydi->nicoup>0)
{
    nblock=(INT)(ydn->nnopo/nthreads);
    nmode=ydn->nnopo-nblock*nthreads;
    if(ydi->icheck==1)
    {
        yde->nnelm=0;
        for(i=0;i<yde->nelemi;i++)
        {
            if(yde->i1elcf[i]>=0)
            {yde->i1elem[yde->nnelm]=i;
                yde->nnelm++;
            }
        }
        ydi->icheck=0;
    }

    if(yde->nnelm>0)
    {
        
      if(ydp->i1ptyp[0]==(YTE3TET10ELS))iquad=1;
      else iquad=0;
      

Yid3TET2TET(         /* tetrahedra to tetrahedra */
            yde->nnelm,yde->i1elem,
            dcstec,
            ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
			ydn->d2nfc,
            ydn->d2nvc[0],ydn->d2nvc[1],ydn->d2nvc[2],
			ydi->d2nv[0],ydi->d2nv[1],ydi->d2nv[2],
            ydi->d1deltat1,ydi->d1deltat2,ydi->d1deltan,
			ydi->d2t1v[0],ydi->d2t1v[1],ydi->d2t1v[2],
			ydi->d2t2v[0],ydi->d2t2v[1],ydi->d2t2v[2],
			ydp->d1pepe,yde->i1elcf,yde->i1elpr,ydp->d1pefr,
			ydp->d1pesf,ydp->d1pevf,ydp->d1pepf,ydp->d1pepsf,
            ydi->i1iecn,ydi->i1iect,yde->i2elto,ydi->i1stat,iquad,
	    yde->d1mart,yde->d1emct,ncstep,ydi->i1fcstep,ydi->i1fcjp,ydi->d1pepe,
            ydj->d2ejfc,ydj->i1elst,yde->d1eldam
            );

  }
}

}




