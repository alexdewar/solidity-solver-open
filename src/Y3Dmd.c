/*! \file Y3Dmd.c by LG
 *  \mesh elements
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
/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/

/* random float between 0 and 1 */

double mrnd(r)
  double *r;
//  int n;
  { int m;
    double s,u,v,p;
    s=65536.0; u=2053.0; v=13849.0;
       *r=u*(*r)+v; m=(int)(*r/s);
        *r=*r-m*s; p=*r/s;
         
    return p;
  }

#include "Yproto.h"
static void Yjd3ELEMENT4(
                         INT ielem, INT i,
                         INT **i2eljp, INT **i2ejjp)
{
// CHRw(stdout,"inside Yjd3ELEMENT4"); 
if(i2ejjp[ielem][0]>=0)
{
    i2ejjp[ielem][1]=i;
}
else
{
    i2ejjp[ielem][0]=i;
}

if(i2eljp[i][0]==ielem)
{
    i2eljp[i][0]=-1;
    printf("i2eljp[%d][0]=-1\n",i);
}
else if(i2eljp[i][1]==ielem)
{
    i2eljp[i][1]=-1;
    printf("i2eljp[%d][1]=-1\n",i);
}
else if(i2eljp[i][2]==ielem)
{
    i2eljp[i][2]=-1;
    printf("i2eljp[%d][2]=-1\n",i);
}
else if(i2eljp[i][3]==ielem)
{
    i2eljp[i][3]=-1;
    printf("i2eljp[%d][3]=-1\n",i);
}
else if(i2eljp[i][0]<0)
{
    i2eljp[i][0]=ielem;
}
else if(i2eljp[i][1]<0)
{
    i2eljp[i][1]=ielem;
}
else if(i2eljp[i][2]<0)
{
    i2eljp[i][2]=ielem;
}
else if(i2eljp[i][3]<0)
{
    i2eljp[i][3]=ielem;
}
}



/*JOINT ELEMENTS*/
static void Yjd3TET4(		/* mesh 4-node tetrahedra */
     INT nelest, INT mnopo, INT nnopst,
     INT *i1pejp,
     INT melem,
     INT *n0elem, INT *n0nopo,
     DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
     DBL *d1ncix, DBL *d1nciy, DBL *d1nciz,
     DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
     INT *i1elpr, INT *i1ejpr, INT **i2ejto,
     INT *i1jnef, INT *i1jnen, INT *i1nobf, INT *i1nopr,
  INT **i2elto, INT **i2eljp, INT *i1elty,INT *i1elbe, 
  INT *i1nold, INT **i2ejjp, INT *i1ptyp)
{
     INT nelem,nnopo;
     INT i,j,k,in,jn,kn,ln,ijnew,ielem;
     INT *i1elto,*i1eltonew,ijoint;
     INT iprop, ipropj;

     nelem=(*n0elem);
     nnopo=(*n0nopo);
     for(ielem=0;ielem<nelest;ielem++)
     {

   	     if(i1ptyp[i1elpr[ielem]]==3)
	  {
		  iprop=i1elpr[ielem];
		  ipropj=i1pejp[iprop];
	       i1elto=i2elto[ielem];
	       
          ijoint=-1;
	       for(i=0;i<4;i++)
	       {
		    j=i+1;
		    if(j>3)
		    {
			 j=0;
		    }		    
		    k=j+1;
		    if(k>3)
		    {
			 k=0;
		    }
		    in=i1elto[i];
		    jn=i1elto[j];
		    kn=i1elto[k];
		    ln=MAXIM(in,jn);
                    ln=MAXIM(ln,kn);
		    if(nnopo>=mnopo)
		    {
			 CHRw(stderr,"Yjd: MNOPO too small");
			 CHRwcr(stderr);
			 exit(1);
		    }

		    /*create new node*/
		    d1nccx[nnopo]=d1nccx[in];
		    d1nccy[nnopo]=d1nccy[in];
		    d1nccz[nnopo]=d1nccz[in];
		    d1ncix[nnopo]=d1ncix[in];
		    d1nciy[nnopo]=d1nciy[in];
		    d1nciz[nnopo]=d1nciz[in];
		    d1nvcx[nnopo]=d1nvcx[in];
		    d1nvcy[nnopo]=d1nvcy[in];
		    d1nvcz[nnopo]=d1nvcz[in];
		    i1nopr[nnopo]=i1nopr[in];
		    i1nobf[nnopo]=i1nobf[in];
                    i1nold[nnopo]=in;

		    nnopo=nnopo+1;

		    /*check if joint element already existent*/
		    ijnew=i1jnef[ln];

		    while(ijnew>=0)
		    {
                         i1eltonew=i2ejto[ijnew];
		      
			 if((i1eltonew[3]==in)&&(i1eltonew[4]==jn)&&(i1eltonew[5]==kn))
			 {
			      i1ejpr[ijnew]=MAXIM(i1ejpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo;
				   i1eltonew[5]=nnopo+1;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo;
				   i1eltonew[5]=nnopo-3;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo-4;
				   i1eltonew[5]=nnopo-3;
			      }
                 ijoint=ijnew;
			      ijnew=-100;

			 }
			 else if((i1eltonew[3]==jn)&&(i1eltonew[4]==kn)&&(i1eltonew[5]==in))
			 {
			      i1ejpr[ijnew]=MAXIM(i1ejpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo;
				   i1eltonew[4]=nnopo+1;
				   i1eltonew[5]=nnopo-1;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo;
				   i1eltonew[4]=nnopo-3;
				   i1eltonew[5]=nnopo-1;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-4;
				   i1eltonew[4]=nnopo-3;
				   i1eltonew[5]=nnopo-1;
			      }
                 ijoint=ijnew;
			      ijnew=-100;
			 }
			 else if((i1eltonew[3]==kn)&&(i1eltonew[4]==in)&&(i1eltonew[5]==jn))
			 {
			      i1ejpr[ijnew]=MAXIM(i1ejpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo+1;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo-3;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-3;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo-4;
			      }
                 ijoint=ijnew;
			      ijnew=-100;

			 }
			 else if((i1eltonew[3]==kn)&&(i1eltonew[4]==jn)&&(i1eltonew[5]==in))
			 {
			      i1ejpr[ijnew]=MAXIM(i1ejpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo+1;
				   i1eltonew[4]=nnopo;
				   i1eltonew[5]=nnopo-1;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo-3;
				   i1eltonew[4]=nnopo;
				   i1eltonew[5]=nnopo-1;				   
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-3;
				   i1eltonew[4]=nnopo-4;
				   i1eltonew[5]=nnopo-1;
			      }
                 ijoint=ijnew;
			      ijnew=-100;

			 }
			 else if((i1eltonew[3]==jn)&&(i1eltonew[4]==in)&&(i1eltonew[5]==kn))
			 {
			      i1ejpr[ijnew]=MAXIM(i1ejpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo+1;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo-3;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-4;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo-3;
			      }
                 ijoint=ijnew;
			      ijnew=-100;

			 }
			 else if((i1eltonew[3]==in)&&(i1eltonew[4]==kn)&&(i1eltonew[5]==jn))
			 {
			      i1ejpr[ijnew]=MAXIM(i1ejpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo+1;
				   i1eltonew[5]=nnopo;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo-3;
				   i1eltonew[5]=nnopo;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo-3;
				   i1eltonew[5]=nnopo-4;
			      }
                 ijoint=ijnew;
			      ijnew=-100;

			 }
			 else
			 {
			      ijnew=i1jnen[ijnew];

			 }
		    }


		    /*create new mid-edge joint node*/
		    if(ijnew>(-10))
		    {
			 if(nelem>=melem)
			 {
			      CHRw(stderr,"Yjd: MELEM too small");
			      CHRwcr(stderr);
			      exit(1);
			 }

			 i1eltonew=i2ejto[nelem];
			 
			 if(i<2)
			 {
			      i1eltonew[2]=nnopo+1;
			      i1eltonew[1]=nnopo;
			      i1eltonew[0]=nnopo-1;
			 }
			 else if(i==2)
			 {
			      i1eltonew[2]=nnopo-3;
			      i1eltonew[1]=nnopo;
			      i1eltonew[0]=nnopo-1;
			 }
			 else
			 {
			      i1eltonew[2]=nnopo-3;
			      i1eltonew[1]=nnopo-4;
			      i1eltonew[0]=nnopo-1;
			 }

			 i1eltonew[3]=kn;
			 i1eltonew[4]=jn;
			 i1eltonew[5]=in;
//			 i2eljp[nelem][0]=ielem;

				 
			 i1ejpr[nelem]=ipropj;
			 i1jnen[nelem]=i1jnef[ln];

			 i1jnef[ln]=nelem;
			 i1elty[nelem]=0;  //-ao asiri added 021017 (all joints are i1elty =0)
             ijoint=nelem;
			 nelem=nelem+1;
		    }
          if(ijoint>=0)Yjd3ELEMENT4(ijoint, ielem,i2eljp,i2ejjp);
	       }

	       /*detach element*/
	       for(i=0;i<4;i++)
	       {
		    i1elto[i]=nnopo-4+i;
	       }
	       
	  }
     }

     for(ielem=0;ielem<nelem;ielem++)
     {
               i1elto=i2ejto[ielem];
	       
	       if(i1elto[3]<nnopst)
	       {
		    i1elto[5]=i1elto[0];
		    i1elto[4]=i1elto[1];
		    i1elto[3]=i1elto[2];
               i1ejpr[ielem]=i1ejpr[ielem]-YIPROPMAX;

		    if (i1elty[ielem]<2)
		    {
			i1elty[ielem]=1; //-ao asiri added (boundary joints are i1elty = 1)
		    }

	       }
     }
    
     (*n0nopo)=nnopo;
     (*n0elem)=nelem;


}

/*DISTINGUISH JOINT ELEMENTS*/
static void Yjd3JOINTPR(	/* distinguish two kinds of joint elements */
			INT nelem,
			INT **i2ejto, INT *i1eljo, INT *i1elty,INT **i2eljp,INT **i2elto)
{
  INT ielem,i,j,inopo,jnopo,knopo;
  INT joflag;

  for(ielem=0;ielem<nelem;ielem++) /* joint element loop */
    {

      joflag=(-1);

      for(j=0;j<2;j++)	/* finite element loop */
	{
        i=i2eljp[ielem][j];
        if(i>0)
        {
	  for(inopo=0;inopo<4;inopo++)
	    {
	      jnopo=inopo+1;
	      if(jnopo>3)
		{
		  jnopo=0;
		}
	      knopo=jnopo+1;
	      if(knopo>3)
		{
		  knopo=0;
		}
	      
	      if((i2ejto[ielem][0]==i2elto[i][inopo])&&(i2ejto[ielem][1]==i2elto[i][jnopo])&&(i2ejto[ielem][2]==i2elto[i][knopo]))
		{
		  joflag=inopo;
		  break;
		}
	    }
	  if(joflag!=(-1))
            {
              break;
            }
	}
    }
      i1eljo[ielem]=joflag;
    }
}





/* LINK NODE TO ITS NEIGHBOURING NODES*/
static void Yjd3NODE2NODE(INT nelest, INT nnopo, INT **i2elto, INT **i2eljp,
                          INT *i1elbe,
                          DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
                          INT *n0neigh, INT **i2nnei, INT *i1nei
                          )
{
  INT ielem,jelem,i,j,k,m,n;
  INT *i1eltoi,*i1eltoj;
  INT nmax,jcon;

  for(ielem=0;ielem<nelest;ielem++)
    {
      i1eltoi=i2elto[ielem];
      for(i=0;i<4;i++)
        {
          i2nnei[i1eltoi[i]][0]=ielem; /* i2nnei[node No.][0] - store elem No. where this node is on */
	  nmax=1;

          for(jelem=0;jelem<nelest;jelem++)
            {
              if(ielem!=jelem)
                {
                  i1eltoj=i2elto[jelem];
                  for(j=0;j<4;j++)
                    {
                      if((i1eltoi[i]!=i1eltoj[j])&&
                         ((DABS(d1nccx[i1eltoi[i]]-d1nccx[i1eltoj[j]])<EPSILON)
                         &&(DABS(d1nccy[i1eltoi[i]]-d1nccy[i1eltoj[j]])<EPSILON)
                          &&(DABS(d1nccz[i1eltoi[i]]-d1nccz[i1eltoj[j]])<EPSILON)))
                        {
                          nmax+=1;
			  for(n=1;n<100;n++)
			    {
			      if(i2nnei[i1eltoi[i]][n]<0)
				{
				  i2nnei[i1eltoi[i]][n]=jelem;
				  break;
				}
			    }
                        }
                    }
                }
            }

	  i1nei[i1eltoi[i]]=nmax; /* number of node neighbours, including itself */

          if(nmax>(*n0neigh))
            {
              (*n0neigh)=nmax;
            }
        }
    }
}


/********************* FIND PREXISTING FRACTURE ********** added by Asiri 20/10/2017 *************/
static void Ypref (INT nelem, DBL *d1ncix, DBL *d1nciy,
       DBL *d1nciz, INT nnopi, INT *i1elty, INT **i2elto)
{
  INT i, j, k;
  INT ielem, jelem;
  INT in0, in1, in2, jn0, jn1, jn2;
  DBL SEP;

  SEP=1.0e-7;


  for (ielem = 0; ielem < nelem; ielem++) //nelemi
    {

      if (i1elty[ielem] > 0)
	{
	  // checking joint whether overlaps with another one
	  in0 = i2elto[ielem][0];
	  in1 = i2elto[ielem][1];
	  in2 = i2elto[ielem][2];

	  for (jelem = ielem + 1; jelem < nelem; jelem++) //nenlmi
	    {
	      //  if (i1elty[jelem]>0)
	      //{

	        for (i = 3; i < 6; i++)
	       {
	       j = i + 1;

	       if (j > 5)
	       j = 3;

	       k = j + 1;

	       if (k > 5)
	       k = 3;


		  jn0 = i2elto[jelem][i];
		  jn1 = i2elto[jelem][j];
		  jn2 = i2elto[jelem][k];

		  /* this does not work properly, need to identify better*/
		  if (DABS(d1ncix[in0]-d1ncix[jn2]) < SEP
		      && DABS(d1ncix[in1]-d1ncix[jn1]) < SEP
		      && DABS(d1ncix[in2]-d1ncix[jn0]) < SEP

		      && DABS(d1nciy[in0]-d1nciy[jn2]) < SEP
		      && DABS(d1nciy[in1]-d1nciy[jn1]) < SEP
		      && DABS(d1nciy[in2]-d1nciy[jn0]) < SEP

		      && DABS(d1nciz[in0]-d1nciz[jn2]) < SEP
		      && DABS(d1nciz[in1]-d1nciz[jn1]) < SEP
		      && DABS(d1nciz[in2]-d1nciz[jn0]) < SEP)
		    {
		      i1elty[ielem] = 2;
		      i1elty[jelem] = 2;
		      break;
		    }
		}
	    }
	}
    }

}


static void Yjd3boundary(    /* distinguish two kinds of joint elements */
            INT nelem, INT nelem_t,
            INT **i2elto, INT **i2elto_t, INT *i1nold, INT *i1elty)
{
  INT ielem,jelem;
  INT in,jn,kn,ln,mn,nn,index;

    for(jelem=0;jelem<nelem_t;jelem++) /* joint element loop */
      {
          in=i2elto_t[jelem][0];
          jn=i2elto_t[jelem][1];
          kn=i2elto_t[jelem][2];
          index=0;
          for(ielem=0;ielem<nelem;ielem++) /* joint element loop */
          {
              if(i1elty[ielem]==1)
              {
                  ln=i1nold[i2elto[ielem][0]];
                  mn=i1nold[i2elto[ielem][1]];
                  nn=i1nold[i2elto[ielem][2]];

                  if((in==ln||in==mn||in==nn)&&
                     (jn==ln||jn==mn||jn==nn)&&
                     (kn==ln||kn==mn||kn==nn))
                  {
                      if(in==ln)
                      {
                          i2elto_t[jelem][0]=i2elto[ielem][0];
                      }
                      else if(in==mn)
                      {
                          i2elto_t[jelem][0]=i2elto[ielem][1];
                      }
                      else
                      {
                          i2elto_t[jelem][0]=i2elto[ielem][2];
                
                      }
                      if(jn==ln)
                      {
                          i2elto_t[jelem][1]=i2elto[ielem][0];
                      }
                      else if(jn==mn)
                      {
                          i2elto_t[jelem][1]=i2elto[ielem][1];
                      }
                      else
                      {
                          i2elto_t[jelem][1]=i2elto[ielem][2];
                
                      }
                      if(kn==ln)
                      {
                          i2elto_t[jelem][2]=i2elto[ielem][0];
                      }
                      else if(kn==mn)
                      {
                          i2elto_t[jelem][2]=i2elto[ielem][1];
                      }
                      else
                      {
                          i2elto_t[jelem][2]=i2elto[ielem][2];
                
                      }
                      index=1;
                      break;
                  }
        }
	  }
              if(index==0)
	      {
//		      CHRw(stdout,"couln't find matched joint");
//		      CHRwsp(stdout); INTw(stdout, jelem,10);
	      }
    }
      }


static void YjdmaterialID( /*assign material property index to joint*/
                          INT nelno,INT nelem,INT **i2eljp, INT nprop, INT **i2elto, INT *i1elpr, INT *i1mat_id,
			  INT *i1index, INT *i1jpr, INT *i1defect,DBL dcdefect, DBL **d2ncc,INT **i2elto_t,DBL *dran)
{
    INT ielem,iprop,i,j,k,index,istat,n0,id_i,id_j;
    DBL fact;
 //   DBL dran;
 //   dran=R1;
    DBL x0,y0,z0,xc,yc,zc,distance;
    for(ielem=0;ielem<nelem;ielem++) /* joint element loop */
    {
        i=i2eljp[ielem][0];
        j=i2eljp[ielem][1];
	 id_i=i1mat_id[i1elpr[i]];
	 id_j=i1mat_id[i1elpr[j]];

        if((i!=-1)&&(j!=-1)&&(i!=j))
        {
            if(id_i<=id_j)
	    {
		    index=id_i*10+id_j;
	    }
            else 
	    {
		    index=id_j*10+id_i;
	    }
	    istat=0;
            for(iprop=0;iprop<nprop;iprop++) /* property loop */
            {
                if(i1index[iprop]==index)
                {
                    i1jpr[ielem]=iprop;
		    istat=1;
                    break;
                }
            }
	    if(istat==0)
	    {
		    CHRw(stdout,"couln't find matched material property");
                      CHRwsp(stdout); INTw(stdout, ielem,10);CHRwsp(stdout); INTw(stdout, index,10);
		      CHRwcr(stdout);
	    }
	    fact=mrnd(dran);
            if((i1jpr[ielem]<6)&&(fact<=dcdefect)&&(dcdefect>EPSILON))i1defect[ielem]=1;
        }

    }
    
    
}


/**************************************************************************/
/* PUBLIC                                                                 */
/**************************************************************************/
	     
void Ymd(YDE yde, YDI ydi, YDN ydn, YDP ydp,YDX ydx, YDJ ydj, YDJP ydjp) /* mesh elements */
{
     INT nelest,nnopst;
     INT iprop;
     INT inopo;
     INT imesh,imestyp;
     INT *i1jnef;   /*joint node element first for each old node*/
     INT *i1jnen;   /*joint node element next for each new node*/
     DBL *d1sdel;
     INT ielem,integ;

       imesh=0; iprop=0;
	  nelest=yde->nelem;
	  nnopst=ydn->nnopo;
	  i1jnef=INT1NULL;
	  i1jnen=INT1NULL;

			 if(i1jnef==INT1NULL)
			 {
		              CHRw(stdout,"allocate memory   ");
			      ydi->diedi=ydi->diezon+ydi->diezon;
			      i1jnef=TalINT1(nnopst);

			      i1jnen=TalINT1(4*nelest);

			      for(inopo=0;inopo<nnopst;inopo++)
			      {
				   i1jnef[inopo]=-1;
			      }
                          }
			      Yjd3TET4(   /*create joints*/
				   yde->nelem, ydn->mnopo, nnopst,
				   //iprop, 
				   ydp->i1pejp,
				   ydj->melem,
				   &(ydj->nelem), &(ydn->nnopo),
				   ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
				   ydn->d2nci[0], ydn->d2nci[1], ydn->d2nci[2],
				   ydn->d2nvc[0], ydn->d2nvc[1], ydn->d2nvc[2],
				   yde->i1elpr,ydj->i1elpr,ydj->i2elto,
				   i1jnef, i1jnen, ydn->i1nobf, ydn->i1nopr,
				   yde->i2elto,yde->i2eljp,ydj->i1elty,yde->i1elbe,
				   ydn->i1nold,ydj->i2ejjp,ydp->i1pemn);


printf("propoerty no %ld %ld %ld %ld\n",imesh,iprop,ydj->nelem, ydn->nnopo);

			      Yjd3JOINTPR( /* distinguish two kinds of joint elements */
					  ydj->nelem,ydj->i2elto,ydj->i1eljo, ydj->i1elty,ydj->i2ejjp,yde->i2elto);
          /* free memory */
	if(i1jnef!=INT1NULL)
	{
	  CHRw(stdout,"before free memory    ");
          FREE(i1jnef);
          FREE(i1jnen);
	}

    ydn->i1nobfse=INT1NULL;
    ydn->i1nobfse=TalINT1(ydn->nnopo);
    if(yde->ifrac)
    {
        ydn->nnobf=0;
        
        for(inopo=0;inopo<ydn->nnopo;inopo++)
        {
        if(ydn->i1nobf[inopo]>0)
        {
            ydn->i1nobfse[ydn->nnobf]=inopo;
            ydn->nnobf++;
        }
        }
    yde->ifrac=0;
        
    }
    Yjd3boundary(    /* distinguish two kinds of joint elements */
                ydj->nelem, ydx->nelem_t,
                 ydj->i2elto, ydx->i2elto_t, ydn->i1nold, ydj->i1elty);
    YjdmaterialID(yde->nelno,ydj->nelem,ydj->i2ejjp,ydjp->nprop,yde->i2elto,yde->i1elpr,ydp->i1index,
                  ydjp->i1ngrp,ydj->i1elpr,ydj->i1defect,ydj->defect,ydn->d2ncc,ydj->i2elto,&(ydj->dran));
    
    
}

/***********************************************************************************/
/* EOF                                                                             */
/***********************************************************************************/			      
			 
		    
		    
      
