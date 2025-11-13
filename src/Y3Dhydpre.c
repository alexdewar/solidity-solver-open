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
* Prof A.Munjiza or Dr J.P.Latham 
* a.munjiza@qmul.ac.uk or j.p.latham@imperial.ac.uk 
* ******************************************************************* */ 

/*1 File   Yhpd.c */

#include "Yproto.h" 
/****************** PUBLIC ************************************/
/* random float between 0 and 1 */
double mrndsf(r)   
  double *r;   
//  int n;   
  { int m;   
    double s,u,v,p;   
    s=65536.0; u=2053.0; v=13849.0;   
       *r=u*(*r)+v; m=(int)(*r/s);   
        *r=*r-m*s; p=*r/s;   
         
    return p;   
  }   

void Yhpd(ydx, ydn,ydc)
/*** integrate hydrostatic pressure over surface mesh  ***/
YDX ydx; YDN ydn; YDC ydc;
{ DBL A[3],x[3],y[3],z[3],area_total;
	DBL xnc,ync,znc,area,xc,yc,zc,nx,ny,nz;
	DBL xe,ye,ze,fx,fy,fz,x0,y0,z0,distance;
	INT i,j,k,n0,n1,n2,l,m,kc;
	DBL dredf; /* ramping reduction factor */
	DBL dredf_jet; /* ramping reduction factor for jet */
        DBL dpressure,jet_pre;
	DBL random_number;
  /* staged loading for PETRONAS project only 21/05/2021 */
  if(ydc->dctime<ydc->dcgrst)
  { dredf=R0;
  }
  else if(ydc->dctime<(ydc->dcgrst+ydc->dcrmpt))
  { dredf=(ydc->dctime-ydc->dcgrst)/ydc->dcrmpt;
  }
  else if(ydc->dctime<(ydc->dcgrst+ydc->dcrmpt+ydc->dcrmct))
  { dredf=R1;
  }
    else if(ydc->dctime<(ydc->dcgrst+ydc->dcrmpt+ydc->dcrmct+ydc->dcrmut))
  { dredf=R1-(ydc->dctime-ydc->dcgrst-ydc->dcrmpt-ydc->dcrmct)/ydc->dcrmut;
  }
  else
  { dredf=R0;
  }


    if(ydc->dctime<(ydc->dcgrst+ydc->dcrmpt))
  { dredf_jet=R0;
  }
  else if(ydc->dctime<(ydc->dcgrst+ydc->dcrmpt+ydc->dcrmjt))
  { dredf_jet=(ydc->dctime-ydc->dcrmpt-ydc->dcgrst)/ydc->dcrmjt;
  }
  else
  { dredf_jet=R1;
  }
//random_number=rand()/RAND_MAX;
//srand(ydc->ncstep);
random_number=(mrndsf(&ydc->r)+0.5);
//random_number=rand()/RAND_MAX+0.5;

area_total=R0; fx=R0; fy=R0; fz=R0;
	for (i=0; i<ydn->nnopo; i++)
	{
		for(j=0;j<3;j++)
		{
		ydn->d2nfp[j][i]=R0;
		}
        ydn->d1npre[i]=R0;
	}
	for (i=0; i<ydx->nelem_t; i++)
	{
        
		n0=ydx->i2elto_t[i][0];
		n1=ydx->i2elto_t[i][1];
		n2=ydx->i2elto_t[i][2];
        
        x0=(ydn->d2ncc[0][n0]+ydn->d2ncc[0][n1]+ydn->d2ncc[0][n2])/3.0;
        y0=(ydn->d2ncc[1][n0]+ydn->d2ncc[1][n1]+ydn->d2ncc[1][n2])/3.0;
        z0=(ydn->d2ncc[2][n0]+ydn->d2ncc[2][n1]+ydn->d2ncc[2][n2])/3.0;
        yc=R0; 
        xc=R0;
	dpressure=dredf*(ydx->d1nap[i]+ydx->d1nbp[i]);
        /*jet pressure boundary z=-141.86x^2+0.3669x+97.653*/
         distance=SQRT(((x0-xc)*(x0-xc))+(y0-yc)*(y0-yc));
        if((ydx->i1jtp[i]>0)&&(distance<=(0.85*1.0e-03)))
	{
		distance=distance*1000.0;
		jet_pre=(-141.86*distance*distance+0.3669*distance+97.653)*2.09e+06;
		/*random flutuation */
		if((ydc->dcstdn-z0)>0.006)jet_pre=jet_pre*(-6.725*(ydc->dcstdn-z0)*1000.0+137.25)/96.9;		
		if(jet_pre<R0)jet_pre=R0;
		dpressure=dpressure+jet_pre*dredf_jet;

	}
       ydn->d1npre[n0]=ydn->d1npre[n0]+dpressure;
        ydn->d1npre[n1]=ydn->d1npre[n0]+dpressure;
        ydn->d1npre[n2]=ydn->d1npre[n0]+dpressure;

		x[0]=ydn->d2ncc[0][n1]-ydn->d2ncc[0][n0];
		y[0]=ydn->d2ncc[1][n1]-ydn->d2ncc[1][n0];
		z[0]=ydn->d2ncc[2][n1]-ydn->d2ncc[2][n0];
		x[2]=ydn->d2ncc[0][n2]-ydn->d2ncc[0][n0];
		y[2]=ydn->d2ncc[1][n2]-ydn->d2ncc[1][n0];
		z[2]=ydn->d2ncc[2][n2]-ydn->d2ncc[2][n0];
		x[1]=ydn->d2ncc[0][n2]-ydn->d2ncc[0][n1];
		y[1]=ydn->d2ncc[1][n2]-ydn->d2ncc[1][n1];
		z[1]=ydn->d2ncc[2][n2]-ydn->d2ncc[2][n1];
		V3DCro(nx,ny,nz,x[0],y[0],z[0],x[2],y[2],z[2]);
		V3DNor(area,nx,ny,nz);
		area=area/R2;
        area_total+=area;

			for (j=0; j<ydx->nelno; j++)
			{  
				k=ydx->i2elto_t[i][j];

				ydn->d2nfp[0][k]=ydn->d2nfp[0][k]+dpressure*nx*area/ydx->nelno;
				ydn->d2nfp[1][k]=ydn->d2nfp[1][k]+dpressure*ny*area/ydx->nelno;
				ydn->d2nfp[2][k]=ydn->d2nfp[2][k]+dpressure*nz*area/ydx->nelno;
				ydn->d2nft[0][k]+=dpressure*nx*area/ydx->nelno;
                                ydn->d2nft[1][k]+=dpressure*ny*area/ydx->nelno;
                                ydn->d2nft[2][k]+=dpressure*nz*area/ydx->nelno;
                                fx+=dpressure*nx*area/ydx->nelno;
				fy+=dpressure*ny*area/ydx->nelno;
				fz+=dpressure*nz*area/ydx->nelno;

	
		}
		
	}

  }

