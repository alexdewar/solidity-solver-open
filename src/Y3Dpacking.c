/*! \file Y3Drd.c
 *  \brief read model parameters
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

static  CHR *cdig="0123456789";


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/

/* random float between 0 and 1 */
double mrnds(r)   
  double *r;   
//  int n;   
  { int m;   
    double s,u,v,p;   
    s=65536.0; u=2053.0; v=13849.0;   
       *r=u*(*r)+v; m=(int)(*r/s);   
        *r=*r-m*s; p=*r/s;   
         
    return p;   
  }   

/* random integer between a and b  */
int mrabs(a,b,r)
  int a,b,*r;
  { int k,l,m,i;
    int p;
    k=b-a+1; l=2;
    while (l<k) l=l+l;
    m=4*l; k=*r; i=1;
	while(i<=1)
    { k=k+k+k+k+k;
	k=k%m; l=k/4+a;
        if (l<=b) { p=l;i=i+1;}
      }
    *r=k;
    return p;
  }



void Ypacking(YD yd,int nthreads)
{ 
  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
//  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YDB ydb=&(yd->ydb);
  YPAR ypar=&(yd->ypar);
  YSP ysp=&(yd->ysp);
    YDX ydx=&(yd->ydx);
  YDK ydk=&(yd->ydk);
  YDJ ydj=&(yd->ydj);
  YDJP ydjp=&(yd->ydjp);
  
  DBL osp[3],*max_h,max,disp[3];
  DBL d1u,d1v,d1w,d1cosine,d1sine;
  DBL zelta,d1ang[3],d1rot[3],cosin,sine,tmp,tmp1;//,xmax,ymax,zmax;
  INT ipos,nnet,i,j,index,rindex,k,m,n,ind;
  const double mult[10]={1./6.,1./24.,1./24.,1./24.,1./60.,1./60.,1./60.,1./120.,1./120.,1./120.};
  DBL **inertia,cm[3],**points,pnx,pny,pnz,xe;
  DBL *intg;
    DBL *v1,*v2,*v3,ro,dx,dy,dz;
  INT tmp_node, tmp_nelem_t,nelem_start;
  DBL rcm[3],cm2ori;
  INT icount,jcount;
  INT r,tid;

    FILE *fp=FILENULL,*fout=FILENULL;

//srand(time(NULL));
  if(ypar->nelem>0)nnet=(INT)((yde->melem-yde->nelem-1)/ypar->nelem);
  else
   nnet=0;
  if(nnet>=1)
  {

	  if(nnet>ypar->nlpar)nnet=ypar->nlpar;



  tmp=-BEPSILON;
/* maximum projection of existing packs    */

    max_h=TalDBL1(nthreads+1);
#pragma omp parallel default(none) shared(ydn,yde,ypar,max_h) \
private (i,j,tmp,tid)
    {
        tid = omp_get_thread_num();
        max_h[tid]=-BEPSILON;
#pragma omp for schedule(guided)
for(i=0;i<yde->nelem;i++)
{
	  if(yde->i1elbe[i]>1)
	  { 
	  for(j=0;j<4;j++)
	  {   tmp=ydn->d2ncc[0][yde->i2elto[i][j]]*ypar->d1lv[0]+
	          ydn->d2ncc[1][yde->i2elto[i][j]]*ypar->d1lv[1]+
			  ydn->d2ncc[2][yde->i2elto[i][j]]*ypar->d1lv[2];
		  if(tmp>max_h[tid])max_h[tid]=tmp;
		  
	    }
    }
}
}
max=-BEPSILON;
    for(i=0;i<nthreads;i++)
    {
    if(max_h[i]>max)max=max_h[i];
    }
    FREE(max_h);

   //temporarily switch off
/* change coordinates according to position of exisiting packs  */
	  disp[0]=R0;
	  disp[1]=R0;
	  disp[2]=R0;
	  index=0;
      tmp1=ypar->radius;
      if((ypar->radius<ypar->depoff)&&(ypar->ifeed==1))
      {ypar->radius=ypar->depoff;
//      ypar->max=ypar->depheg;
      } 
if(max>(ypar->depheg-ypar->radius))index=1;

      ypar->radius=tmp1;
/* place particles in random grid and rotate it in random angle */
if(index==0)
{
      fp=fopen("orientation.txt","a");


        v1=TalDBL1(3);
        v2=TalDBL1(3);
        v3=TalDBL1(3);
        inertia=TalDBL2(3,3);
        points=TalDBL2(4,3);
        intg=TalDBL1(10);
	
	        nelem_start=yde->nelem;
	ydi->diedi=R2*ydi->diezon;
    //    ydi->iiecff=-2;
	for(i=0;i<(ypar->ngrid);i++)
	{
	ipos=mrabs(1,ypar->ngrid,&ypar->ir);

        if(i>=nnet)break;
//	if(ipos<=nnet)
	{
	osp[0]=(ypar->d2gxyz[0][ipos-1]+disp[0]);
	osp[1]=(ypar->d2gxyz[1][ipos-1]+disp[1]);
        osp[2]=(ypar->depheg+disp[2]);

	/*node    */

	 for(k=0;k<ypar->nnode;k++)
	 {

if(ypar->ifeed)
{

    
    ydn->d2ncc[0][ydn->nnopo]=ypar->d2xyz[0][k]-ypar->d1cg[0]+osp[0];
    ydn->d2ncc[1][ydn->nnopo]=ypar->d2xyz[1][k]-ypar->d1cg[1]+osp[1];
    ydn->d2ncc[2][ydn->nnopo]=ypar->d2xyz[2][k]-ypar->d1cg[2]+osp[2];

}
else
{
		 
    d1u=ypar->d2gxyz[3][i];
    d1v=ypar->d2gxyz[4][i];
    d1w=ypar->d2gxyz[5][i];

    
	V3DRotAxes(&ydn->d2ncc[0][ydn->nnopo],&ydn->d2ncc[1][ydn->nnopo],&ydn->d2ncc[2][ydn->nnopo],
    ypar->d2xyz[0][k]-ypar->d1cg[0],ypar->d2xyz[1][k]-ypar->d1cg[1],ypar->d2xyz[2][k]-ypar->d1cg[2],d1u,d1v,d1w);

        ydn->d2ncc[0][ydn->nnopo]=ydn->d2ncc[0][ydn->nnopo]+ypar->d2gxyz[0][i];
        ydn->d2ncc[1][ydn->nnopo]=ydn->d2ncc[1][ydn->nnopo]+ypar->d2gxyz[1][i];
        ydn->d2ncc[2][ydn->nnopo]=ydn->d2ncc[2][ydn->nnopo]+ypar->d2gxyz[2][i];
}
    
	ydn->d2nci[0][ydn->nnopo]=ydn->d2ncc[0][ydn->nnopo];
	ydn->d2nci[1][ydn->nnopo]=ydn->d2ncc[1][ydn->nnopo];
	ydn->d2nci[2][ydn->nnopo]=ydn->d2ncc[2][ydn->nnopo];
	
	ydn->nnopo=ydn->nnopo+1;

	 }

	 /*nelem  */
	 ypar->nindex=ypar->nindex+1; //particle index
	 ysp->nspd=ysp->nspd+1;

 	 ydn->nnopo=ydn->nnopo-ypar->nnode;
        tmp_node=ydn->nnopo;
        tmp_nelem_t=ydx->nelem_t;


	 for(k=0;k<ypar->nelem;k++)
	 {
		 for(m=0;m<(yde->nelno-1);m++)
		 {
			 yde->i2elto[yde->nelem][m]=ypar->i2elto[k][m]+ydn->nnopo;
			 ydn->i1nind[yde->i2elto[yde->nelem][m]]=ypar->nindex;
		 }
		 yde->i2elto[yde->nelem][yde->nelno-1]=ypar->nindex;
                 yde->i1elpr[yde->nelem]=ydp->nprop-1;
                 yde->i1elbe[yde->nelem]=0;
		 if(ypar->i1elbe[k]!=0)yde->i1elbe[yde->nelem]=ypar->nindex;
		 yde->nelem=yde->nelem+1;
	 }

        for(k=0;k<ypar->nelemb;k++)
        {
            for(m=0;m<ydx->nelno;m++)
            {
                ydx->i2elto_t[ydx->nelem_t][m]=ypar->i2eltob[k][m]+ydn->nnopo;
            }
            ydx->nelem_t=ydx->nelem_t+1;
        }

        

	 ydn->nnopo=ydn->nnopo+ypar->nnode;
        
        
        for(k=tmp_nelem_t; k<ydx->nelem_t; k++)
        {
            for(m=0; m<(ydx->nelno); m++)
            {
                ydn->i1nobf[ydx->i2elto_t[k][m]]=1;
                if(m<3)ydn->i1nobf[ydx->i2elto_t[k][m]]=2;
            }
        }
  
	 for(j=0;j<10;j++)intg[j]=R0;
	 
	 for(k=0;k<ypar->nelem;k++)
	 {	 
		 for(m=0;m<4;m++)
		 {
	     points[m][0]=ydn->d2ncc[0][yde->i2elto[yde->nelem-ypar->nelem+k][m]];
		 points[m][1]=ydn->d2ncc[1][yde->i2elto[yde->nelem-ypar->nelem+k][m]];
		 points[m][2]=ydn->d2ncc[2][yde->i2elto[yde->nelem-ypar->nelem+k][m]];
		 }
		 ro=ydp->d1pero[yde->i1elpr[yde->nelem-ypar->nelem+k]];
		 mirtichRoutine(points,intg,ro);
	 }

     for(n=0;n<10;n++)
      intg[n]=intg[n]*mult[n];

     //this is an error in the original publication!!!!!!!!!!!!!!!!!
     //they dont have the three!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 ind=ysp->nspd-1;
    ysp->d1mass[ind] = intg[0];  //density = 1.0 assumed 

     //center of mass 
    ysp->d2pp[0][ind]=intg[1]/ysp->d1mass[ind];
    ysp->d2pp[1][ind]=intg[2]/ysp->d1mass[ind];
    ysp->d2pp[2][ind]=intg[3]/ysp->d1mass[ind];
	cm[0]=ysp->d2pp[0][ind];
	cm[1]=ysp->d2pp[1][ind];
	cm[2]=ysp->d2pp[2][ind];

    
    inertia[0][0] = intg[5] + intg[6] -ysp->d1mass[ind]*(cm[1]*cm[1]+cm[2]*cm[2]);
    inertia[1][1] = intg[4] + intg[6] -ysp->d1mass[ind]*(cm[0]*cm[0]+cm[2]*cm[2]);
    inertia[2][2] = intg[4] + intg[5] -ysp->d1mass[ind]*(cm[1]*cm[1]+cm[0]*cm[0]);

    inertia[0][1] = -(intg[7] -ysp->d1mass[ind]*cm[1]*cm[0]);
    inertia[1][0]=inertia[0][1];

    inertia[1][2] = -(intg[8] -ysp->d1mass[ind]*cm[1]*cm[2]);
    inertia[2][1]=inertia[1][2];

    inertia[0][2] = -(intg[9] -ysp->d1mass[ind]*cm[2]*cm[0]);
    inertia[2][0]=inertia[0][2];

    //this is because the vertices could have been DoubleLinkedListed in the wrong direction
    if((inertia[0][0]<0.0)||(inertia[1][1]<0.0)||(inertia[2][2]<0.0))
      {
      //inertia*=(-1.0);
	  for(m=0;m<3;m++)
	  {
		  for(n=0;n<3;n++)
	  {
       inertia[m][n]=-inertia[m][n];
	  }
	  }
       ysp->d1mass[ind]*=(-1.0);
      }


	solveSymetricEigenProblem(inertia, v1,v2,v3,
						&ysp->d2prine[0][ind],&ysp->d2prine[1][ind],&ysp->d2prine[2][ind]);
		for(m=0;m<3;m++)
	  {
		ysp->d2prinn1[m][ind]=v1[m];
		ysp->d2prinn2[m][ind]=v2[m];
		ysp->d2prinn3[m][ind]=v3[m];

		ysp->d2pvt[m][ind]=ypar->d1vc[m];
		}

        
        if(ypar->ifeed)
        {
            
            d1rot[0]=R0; d1rot[1]=R0;d1rot[2]=R0;
            d1rot[0]=-R1+R2*mrnds(&ypar->r);
            d1rot[1]=-R1+R2*mrnds(&ypar->r);
            tmp=SQRT(mrnds(&ypar->r));
            V2DNor(xe,d1rot[0],d1rot[1]);
                    
            d1rot[0]=d1rot[0]*tmp;
            d1rot[1]=d1rot[1]*tmp;
            d1rot[2]=SQRT(R1-d1rot[0]*d1rot[0]-d1rot[1]*d1rot[1]);

                    V3DCro(pnx,pny,pnz,v1[0],v1[1],v1[2],
                           d1rot[0],d1rot[1],d1rot[2]);
                    V3DNor(tmp,pnx,pny,pnz);
                    V3DDot(cosin,d1rot[0],d1rot[1],d1rot[2],v1[0],v1[1],v1[2]);
                    tmp=R1-cosin*cosin;
                    if(tmp<R0)
                        sine=R0;
                    else
                        sine=SQRT(R1-cosin*cosin);
            
            V3DRot(&ysp->d2prinn1[0][ind],&ysp->d2prinn1[1][ind],&ysp->d2prinn1[2][ind],
                   ysp->d2prinn1[0][ind],ysp->d2prinn1[1][ind],ysp->d2prinn1[2][ind],
                   pnx,pny,pnz,cosin,sine);
            V3DRot(&ysp->d2prinn2[0][ind],&ysp->d2prinn2[1][ind],&ysp->d2prinn2[2][ind],
                   ysp->d2prinn2[0][ind],ysp->d2prinn2[1][ind],ysp->d2prinn2[2][ind],
                   pnx,pny,pnz,cosin,sine);
            V3DRot(&ysp->d2prinn3[0][ind],&ysp->d2prinn3[1][ind],&ysp->d2prinn3[2][ind],
                   ysp->d2prinn3[0][ind],ysp->d2prinn3[1][ind],ysp->d2prinn3[2][ind],
                   pnx,pny,pnz,cosin,sine);
            
            V3DNor(tmp,ysp->d2prinn1[0][ind],ysp->d2prinn1[1][ind],ysp->d2prinn1[2][ind]);
            V3DNor(tmp,ysp->d2prinn2[0][ind],ysp->d2prinn2[1][ind],ysp->d2prinn2[2][ind]);
            V3DNor(tmp,ysp->d2prinn3[0][ind],ysp->d2prinn3[1][ind],ysp->d2prinn3[2][ind]);
            
            for(k=tmp_node; k<ydn->nnopo; k++)
            {
                V3DRot(&ydn->d2ncc[0][k],&ydn->d2ncc[1][k],&ydn->d2ncc[2][k],
                       ydn->d2ncc[0][k]-cm[0],ydn->d2ncc[1][k]-cm[1],ydn->d2ncc[2][k]-cm[2],
                       pnx,pny,pnz,cosin,sine);
                ydn->d2ncc[0][k]=ydn->d2ncc[0][k]+cm[0];
                ydn->d2ncc[1][k]=ydn->d2ncc[1][k]+cm[1];
                ydn->d2ncc[2][k]=ydn->d2ncc[2][k]+cm[2];
                ydn->d2nci[0][k]=ydn->d2ncc[0][k];
                ydn->d2nci[1][k]=ydn->d2ncc[1][k];
                ydn->d2nci[2][k]=ydn->d2ncc[2][k];
                            }


            index=1;
			rcm[0]=cm[0];
			rcm[1]=cm[1];
			rcm[2]=cm[2];
			cm2ori=cm[2];
            while(index==1)
            {rindex=1;
				       while(rindex==1)
                       {
                       rcm[0]=ypar->depora*(-R1+R2*mrnds(&ypar->r));
                       rcm[1]=ypar->depora*(-R1+R2*mrnds(&ypar->r));
                       if(SQRT(rcm[0]*rcm[0]+rcm[1]*rcm[1])>ypar->depora)rindex=1;
					   else rindex=0;
					   rcm[2]=ypar->depoff*(mrnds(&ypar->r)-0.5)*2.0+cm2ori;
                       }
            for(k=tmp_node; k<ydn->nnopo; k++)
            {
                ydn->d2ncc[0][k]=ydn->d2ncc[0][k]+rcm[0]-cm[0];
                ydn->d2ncc[1][k]=ydn->d2ncc[1][k]+rcm[1]-cm[1];
                ydn->d2ncc[2][k]=ydn->d2ncc[2][k]+rcm[2]-cm[2];
                            }
			cm[0]=rcm[0];
				cm[1]=rcm[1];
				cm[2]=rcm[2];
				
			ysp->d2pp[0][ind]=rcm[0];
			ysp->d2pp[1][ind]=rcm[1];
			ysp->d2pp[2][ind]=rcm[2];
     index=Ycdran( /* Y elements contact detection */
      nelem_start,yde->nelem,      
      ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
      yde->i1elbe,   yde->i2elto, ydi->drect,ydn->nnopo,nthreads,ypar->nindex);
        } 

            for(k=tmp_node; k<ydn->nnopo; k++)
            {
                ydn->d2nci[0][k]=ydn->d2ncc[0][k];
                ydn->d2nci[1][k]=ydn->d2ncc[1][k];
                ydn->d2nci[2][k]=ydn->d2ncc[2][k];
 
             }
   
            dx=ypar->d1rv_max[0]-ypar->d1rv_min[0];
            dy=ypar->d1rv_max[1]-ypar->d1rv_min[1];
            dz=ypar->d1rv_max[2]-ypar->d1rv_min[2];
            if(dx>EPSILON)ysp->d2pvt[0][ind]=dx*mrnds(&ypar->r)+ypar->d1rv_min[0];
            if(dy>EPSILON)ysp->d2pvt[1][ind]=dy*mrnds(&ypar->r)+ypar->d1rv_min[1];
            if(dz>EPSILON)ysp->d2pvt[2][ind]=dz*mrnds(&ypar->r)+ypar->d1rv_min[2];
            
            for(k=0;k<3;k++)
                
            {
                DBLw(fp,ysp->d2prinn1[k][ind],18);
                CHRwsp(fp);
                
            }
            CHRwcr(fp);

            
        }
        

	}
	}
              Ymass3TET(yde, ydn, ydp, ydk,ydj,ydjp,nthreads,ydc->dcstec);
        FREE(intg);
        FREE(points);
        FREE(inertia);
        FREE(v3);
        FREE(v2);
        FREE(v1);
      fclose(fp);
}

}

	return;
}



/* run packing test   */


