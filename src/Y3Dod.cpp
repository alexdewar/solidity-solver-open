/*! \file Y3Dod.c
 *  \brief output database for Y3D model
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

/*
extern "C" {
    #include "Yproto.h"
}*/
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkZLibDataCompressor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>

//#include <vtkSystemIncludes.h>

extern "C" {
#include "Yproto.h"
#include "frame.h"
static INT  i1num[100];    /* numbers for space saving format     */
static DBL  d1num[100];    /* numbers for space saving format     */
static CHR c1code[500];    /* coded i1para in space saving format */

}
/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


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
                                   DBL **d2tcs,
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


    for(j=0;j<3;j++)
    {
        d2tcs[ielem][j]+=T[j][j]*WEIGHT;
    }
    d2tcs[ielem][3]+=T[0][1]*WEIGHT;
    d2tcs[ielem][4]+=T[0][2]*WEIGHT;
    d2tcs[ielem][5]+=T[1][2]*WEIGHT;

}



/**********************************************************************/
/**********************************************************************/


/* small strain elastic tetrahetra
 */
static void Yod3TET10(
                         INT nelem, INT nnopo,CHR *fout,
                         DBL *dpeks,    DBL *dpela,    DBL *dpemu,
                         DBL **d2ncc,  DBL **d2nci,  DBL **d2nfc,    DBL **d2nvc,
                         INT *i1elpr,  INT **i2elto,
                         DBL ***d3dsh,INT nspd, INT *i1nind
                      )
{
    DBL **dm1tn;
    DBL **d2tcs, dtnn;
    INT *i1elto, *i1tnn,i,ielem,ip,j;
    DBL *v1,*v2,*v3,a1,a2,a3,tmp1,min,max,med;
    DBL **T,*sigma1_max;

    d2tcs=TalDBL2(nelem+1,7);
    TzDBL2(d2tcs,nelem+1,7);
#pragma omp parallel default(none) shared(nelem,dpeks,dpela,dpemu,d2ncc,d2nci,\
d2nvc,i1elpr,i2elto,d3dsh,nnopo,d2tcs) private(i1elto,ielem)
    {
        DBL d2nd[NDIME][NNODE], d2finv[NDIME][NDIME], d2ftens[NDIME][NDIME];
        DBL d2xj[NDIME][NDIME], d2xjci[NDIME][NDIME], d2xjinv[NDIME][NDIME];
        DBL L[NDIME][NDIME], T[NDIME][NDIME];
        DBL d3st[NGRSH][NDIME][NNODE];
        DBL voli, detci, detj, detf, detf0;
        INT ig;
        INT iprop;


#pragma omp for schedule(guided) nowait
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
            YMATINV3(d2xj, d2xjinv, detj); //DBG_NEG(detj, ielem);
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
                YMATINV3(d2xj,d2xjinv,detj); //DBG_NEG(detj,ielem)
                /* calculate Cartesian gradients nd, finv, and velocity gradient L */
                grad(d2nd, d2xjinv, d3dsh[ig]);
                ctvgrad(d2nd, d2finv, d2nci, L, d2nvc, i1elto);
                /* calculate inverse and determinant of finv */
                YMATINV3(d2finv, d2ftens, detf);

                /* calculate left Cauchy-Green strain tensor B, deformation rate D, Cauchy stress T, net stress tcs */
                cauchy(T, L, d2ftens, d2tcs, dpemu[iprop], dpela[iprop], detf, detf0, dpeks[iprop], ielem);

            }

        }

    }

sigma1_max=TalDBL1(nspd+1);
T=TalDBL2(3,3);

  v1=TalDBL1(3);
  v2=TalDBL1(3);
  v3=TalDBL1(3);

  a1=R0;
  a2=R0;
  a3=R0;

  for(i=0;i<3;i++)
    {
      v1[i]=R0;
      v2[i]=R0;
      v3[i]=R0;
    }
  for(i=0;i<nspd;i++)
    {
     sigma1_max[i]=R0;
    }
  tmp1=R0;
  min=R0;
  max=R0;
  med=R0;


    dm1tn=TalDBL2(nnopo+1,7);
    TzDBL2(dm1tn,nnopo+1,7);
    i1tnn = TalzINT1(nnopo+1);
    TzINT1(i1tnn,nnopo+1);
//    CHRw(stdout,"fff");
    for(ielem=0; ielem<nelem; ielem++)
    {
        i1elto = i2elto[ielem];

        for(i=0; i<NNODE; i++)
        {
            ip = i1elto[i];

            for(j=0; j<6; j++)
            {
                dm1tn[ip][j] += d2tcs[ielem][j];
            }

            i1tnn[ip]++;
        }
    }

    vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
    vtkPoints *pts = vtkPoints::New();
    vtkCellArray *polys = vtkCellArray::New();

    vtkDoubleArray *vects = vtkDoubleArray::New();
    vects->SetName("VelocityVectors");
    vects->SetNumberOfComponents(3);

    vtkDoubleArray *fcontact = vtkDoubleArray::New();
    fcontact->SetName("ContactForces");
    fcontact->SetNumberOfComponents(3);

    vtkDoubleArray *tensor = vtkDoubleArray::New();
    tensor->SetName("STRESS");
    tensor->SetNumberOfComponents(9);

    vtkDoubleArray *scalar1 = vtkDoubleArray::New();
    scalar1->SetName("Sigma1");
    scalar1->SetNumberOfComponents(1);

    vtkDoubleArray *scalar2 = vtkDoubleArray::New();
    scalar2->SetName("Sigma2");
    scalar2->SetNumberOfComponents(1);

    vtkDoubleArray *scalar3 = vtkDoubleArray::New();
    scalar3->SetName("Sigma3");
    scalar3->SetNumberOfComponents(1);

    vtkDoubleArray *scalar4 = vtkDoubleArray::New();
    scalar4->SetName("DiffStress");
    scalar4->SetNumberOfComponents(1);

    vtkDoubleArray *scalar5 = vtkDoubleArray::New();
    scalar5->SetName("Sigma1_max");
    scalar5->SetNumberOfComponents(1);

    vtkIntArray *index = vtkIntArray::New();
    index->SetName("Index");
    index->SetNumberOfComponents(1);

    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();

    vtkZLibDataCompressor* myZlibCompressor  = vtkZLibDataCompressor::New();
    myZlibCompressor->SetCompressionLevel(1);


    dataSet->Allocate();

    for(i=0; i<nnopo; i++)
    {
        dtnn = (DBL)i1tnn[i];

        for(j=0; j<6; j++)
        {

            dm1tn[i][j] /= dtnn;

        }
       T[0][0]=dm1tn[i][0];
       T[0][1]=dm1tn[i][3];
       T[0][2]=dm1tn[i][4];
       T[1][0]=dm1tn[i][3];
       T[1][1]=dm1tn[i][1];
       T[1][2]=dm1tn[i][5];
       T[2][0]=dm1tn[i][4];
       T[2][1]=dm1tn[i][5];
       T[2][2]=dm1tn[i][2];

        solveSymetricEigenProblem(T,v1,v2,v3,&a1,&a2,&a3);

        max=MAXIM(a1,a2);
        max=MAXIM(max,a3);
        min=MINIM(a1,a2);
        min=MINIM(min,a3);

        if((DABS(a1-max)>EPSILON)&&(DABS(a1-min)>EPSILON))
          {
            med=a1;
          }
        if((DABS(a2-max)>EPSILON)&&(DABS(a2-min)>EPSILON))
          {
            med=a2;
          }
        if((DABS(a3-max)>EPSILON)&&(DABS(a3-min)>EPSILON))
          {
            med=a3;
          }

        tmp1=max-min;
        sigma1_max[i1nind[i]-1]=MAXIM(sigma1_max[i1nind[i]-1],max);

    }

    for(i=0; i<nnopo; i++)
    {
     
       T[0][0]=dm1tn[i][0];
       T[0][1]=dm1tn[i][3];
       T[0][2]=dm1tn[i][4];
       T[1][0]=dm1tn[i][3];
       T[1][1]=dm1tn[i][1];
       T[1][2]=dm1tn[i][5];
       T[2][0]=dm1tn[i][4];
       T[2][1]=dm1tn[i][5];
       T[2][2]=dm1tn[i][2];

        /* prepare output */
        pts->InsertNextPoint(d2ncc[0][i],d2ncc[1][i],d2ncc[2][i]);
        vects->InsertNextTuple3(d2nvc[0][i],d2nvc[1][i],d2nvc[2][i]);
        fcontact->InsertNextTuple3(d2nfc[0][i],d2nfc[1][i],d2nfc[2][i]);
        index->InsertNextTuple1(i1nind[i]);

        tensor->InsertNextTuple9(dm1tn[i][0],dm1tn[i][3],dm1tn[i][4],
                                 dm1tn[i][3],dm1tn[i][1],dm1tn[i][5],
                                 dm1tn[i][4],dm1tn[i][5],dm1tn[i][2]);
        solveSymetricEigenProblem(T,v1,v2,v3,&a1,&a2,&a3);

        max=MAXIM(a1,a2);
        max=MAXIM(max,a3);
        min=MINIM(a1,a2);
        min=MINIM(min,a3);

        if((DABS(a1-max)>EPSILON)&&(DABS(a1-min)>EPSILON))
          {
            med=a1;
          }
        if((DABS(a2-max)>EPSILON)&&(DABS(a2-min)>EPSILON))
          {
            med=a2;
          }
        if((DABS(a3-max)>EPSILON)&&(DABS(a3-min)>EPSILON))
          {
            med=a3;
          }

        tmp1=max-min;
        scalar1->InsertNextTuple1(max);
        scalar2->InsertNextTuple1(med);
        scalar3->InsertNextTuple1(min);
        scalar4->InsertNextTuple1(tmp1);
        scalar5->InsertNextTuple1(sigma1_max[i1nind[i]-1]);

    }

    for(ielem=0; ielem<nelem; ielem++)
    {
        i1elto = i2elto[ielem];
        polys->InsertNextCell(10);
        for(i=0;i<10;i++)
        {
            polys->InsertCellPoint(i1elto[i]);
        }
    }
    dataSet->SetPoints(pts);
    pts->Delete();

    dataSet->GetPointData()->AddArray(vects);
    dataSet->GetPointData()->SetActiveAttribute("VELOCITY",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(vects);
    vects->Delete();

    dataSet->GetPointData()->AddArray(fcontact);
    dataSet->GetPointData()->SetActiveAttribute("CONTACT",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(fcontact);
    fcontact->Delete();

    dataSet->GetPointData()->AddArray(tensor);
    dataSet->GetPointData()->SetActiveAttribute("STRESS",vtkDataSetAttributes::TENSORS);
    dataSet->GetPointData()->SetTensors(tensor);
    tensor->Delete();

    dataSet->GetPointData()->AddArray(scalar1);
    dataSet->GetPointData()->SetActiveAttribute("Sigma1",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar1);
    scalar1->Delete();

    dataSet->GetPointData()->AddArray(scalar2);
    dataSet->GetPointData()->SetActiveAttribute("Sigma2",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar2);
    scalar2->Delete();

    dataSet->GetPointData()->AddArray(scalar3);
    dataSet->GetPointData()->SetActiveAttribute("Sigma3",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar3);
    scalar3->Delete();

    dataSet->GetPointData()->AddArray(scalar4);
    dataSet->GetPointData()->SetActiveAttribute("DiffStress",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar4);
    scalar4->Delete();

    dataSet->GetPointData()->AddArray(scalar5);
    dataSet->GetPointData()->SetActiveAttribute("Sigma1_max",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar5);
    scalar5->Delete();
    dataSet->GetPointData()->SetActiveAttribute("Index",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(index);
    index->Delete();

    dataSet->SetCells(VTK_QUADRATIC_TETRA,polys);

    polys->Delete();
    writer->SetInputData(dataSet);
    writer->SetFileName(fout);
    writer->SetCompressor(myZlibCompressor);
    writer->Write();
    dataSet->Delete();
    myZlibCompressor->Delete();
    writer->Delete();

    FREE(i1tnn);
    FREE(dm1tn);
    FREE(v3);
    FREE(v2);
    FREE(v1);
    FREE(T);
    FREE(sigma1_max);
    FREE(d2tcs);

}

static INLINE DBL FdftensTet4(DBL F[3][3],DBL L[3][3],DBL *d1nccx,  DBL *d1nccy,  DBL *d1nccz,
                              DBL *d1ncix,  DBL  *d1nciy,
                              DBL *d1nciz,DBL *d1nvcx,  DBL *d1nvcy,  DBL *d1nvcz, INT *i1elto)
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
    detf=volc/voli;

    //    __builtin_prefetch(&L[0][0],1,1);
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

static INLINE void CauchyTet4_pore_pressure(DBL F[3][3],DBL L[3][3],DBL dpemu,DBL dpeks,DBL dpela,DBL detf,DBL **T,DBL *T_pore_pressure)
{
    DBL  B[3][3]; /* left Cauchy-Green strain tensor */
    DBL  D[3][3]; /* rate of deformation (stretching) tensor */
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
      Tm=-(detf-1.0)*2.0*dpemu*SK*(1.0 + dups_u)/(3.0*(1.0-2.0*dups_u));
      Tm=R0;
for(i=0;i<3;i++)
  {
        T[i][i]=T[i][i]-ALPHA/SK*Tm;
  }
 //   *T_pore_pressure=-ALPHA/SK*Tm;
    *T_pore_pressure=Tm;
}


static INLINE void CauchyTet4(DBL F[3][3],DBL L[3][3],DBL dpemu,DBL dpeks,DBL dpela,DBL detf,DBL **T)
{
    DBL  B[3][3]; /* left Cauchy-Green strain tensor */
    DBL  D[3][3]; /* rate of deformation (stretching) tensor */
    int i,j,k;
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
}


/* small strain softening quadratic tetrahedron output
*/
void Yod3TET4(
            INT nelem,   INT nnopi,INT nnopo,
            CHR *fout,    DBL *dpeks,    DBL *dpela,    DBL *dpemu,
            DBL *d1nccx,  DBL *d1nccy,  DBL *d1nccz,  DBL *d1ncix,  DBL  *d1nciy,
            DBL *d1nciz,  DBL **d2nfc,DBL *d1nvcx,  DBL *d1nvcy,  DBL *d1nvcz,
            INT *i1elpr,  INT **i2elto,INT *i1ptyp, INT *i1elfrg, DBL *d1mart,
            DBL *d1eldam, DBL *d1npre,DBL ***d3tst, INT *i1elfrm, DBL *d1ndisp//,INT nspd, INT *i1nind
            )
{
    DBL **d2tcs, *sigma1,*diff,*sigma3,*T_pore_pressure;
  INT *i1elto, i,ielem,ip,j;

    d2tcs=TalDBL2(nelem+1,7);
 //   TzDBL2(d2tcs,nelem+1,7);
    sigma1=TalDBL1(nelem+1);
    sigma3=TalDBL1(nelem+1);
    diff=TalDBL1(nelem+1);
    T_pore_pressure=TalDBL1(nelem+1);
#pragma omp parallel default(none) shared(nelem,i1elpr,dpeks,dpela,dpemu,d1npre,\
d1nccx,d1nccy,d1nccz,d1ncix,d1nciy,d1nciz,d1nvcx,d1nvcy,d1nvcz,i2elto,d2tcs,sigma1,sigma3,diff,T_pore_pressure,stdout)\
private(i1elto,i,j)
    {
        DBL  F[3][3]; /* deformation gradient in global base delta ux/delta x */
        DBL  L[3][3]; /* velocity gradient in global base  delta vx/delta x    */
        DBL **T;
        DBL  detf;
        INT iprop;
        DBL *v1,*v2,*v3,a1,a2,a3,tmp1,min,max,med;
T=TalDBL2(3,3);
  v1=TalDBL1(3);
  v2=TalDBL1(3);
  v3=TalDBL1(3);

#pragma omp for schedule(guided) nowait
        for(i=0;i<nelem;i++)
        {
            iprop=i1elpr[i];
            i1elto=i2elto[i];

            detf=FdftensTet4(F,L,d1nccx,d1nccy,d1nccz,d1ncix,d1nciy,d1nciz,d1nvcx,d1nvcy,d1nvcz,i1elto);
	    //            CauchyTet4(F,L,dpemu[iprop],dpeks[iprop],dpela[iprop],detf,T);
            CauchyTet4_pore_pressure(F,L,dpemu[iprop],dpeks[iprop],dpela[iprop],detf,T,&T_pore_pressure[i]);

            for(j=0;j<3;j++)
            {
                d2tcs[i][j]=T[j][j];
            }
            d2tcs[i][3]=T[0][1];
            d2tcs[i][4]=T[0][2];
            d2tcs[i][5]=T[1][2];

  a1=R0;
  a2=R0;
  a3=R0;

  for(j=0;j<3;j++)
    {
      v1[j]=R0;
      v2[j]=R0;
      v3[j]=R0;
    }
  tmp1=R0;
  min=R0;
  max=R0;
  med=R0;

            solveSymetricEigenProblem(T,v1,v2,v3,&a1,&a2,&a3);
            max=MAXIM(a1,a2);
            max=MAXIM(max,a3);
            min=MINIM(a1,a2);
            min=MINIM(min,a3);

            if((DABS(a1-max)>EPSILON)&&(DABS(a1-min)>EPSILON))
              {
                med=a1;
              }
            if((DABS(a2-max)>EPSILON)&&(DABS(a2-min)>EPSILON))
              {
                med=a2;
              }
            if((DABS(a3-max)>EPSILON)&&(DABS(a3-min)>EPSILON))
              {
                med=a3;
              }

            tmp1=max-min;
            sigma1[i]=max;
            sigma3[i]=min;
            diff[i]=max-min;
    //        sigma1_max[i1nind[i]-1]=MAXIM(sigma1_max[i1nind[i]-1],max);
    }
    FREE(v3);
    FREE(v2);
    FREE(v1);
    FREE(T);
}
    vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
    vtkPoints *pts = vtkPoints::New();
    vtkCellArray *polys = vtkCellArray::New();

    vtkDoubleArray *vects = vtkDoubleArray::New();
    vects->SetName("VelocityVectors");
    vects->SetNumberOfComponents(3);

    vtkDoubleArray *fcontact = vtkDoubleArray::New();
    fcontact->SetName("ContactForces");
    fcontact->SetNumberOfComponents(3);

    vtkDoubleArray *ploading = vtkDoubleArray::New();
    ploading->SetName("Pressure_loading");
    ploading->SetNumberOfComponents(1);

    vtkDoubleArray *pdisp = vtkDoubleArray::New();
    pdisp->SetName("Nodal_displacement");
    pdisp->SetNumberOfComponents(1);

    vtkDoubleArray *tensor = vtkDoubleArray::New();
    tensor->SetName("STRESS");
    tensor->SetNumberOfComponents(9);

    vtkDoubleArray *tensor2 = vtkDoubleArray::New();
    tensor2->SetName("StrainRate");
    tensor2->SetNumberOfComponents(9);

    vtkDoubleArray *scalar1 = vtkDoubleArray::New();
    scalar1->SetName("Sigma1");
    scalar1->SetNumberOfComponents(1);

    vtkDoubleArray *scalar3 = vtkDoubleArray::New();
    scalar3->SetName("Sigma3");
    scalar3->SetNumberOfComponents(1);


    vtkDoubleArray *scalar4 = vtkDoubleArray::New();
    scalar4->SetName("DiffStress");
    scalar4->SetNumberOfComponents(1);

    vtkDoubleArray *scalar6 = vtkDoubleArray::New();
    scalar6->SetName("MassRatio");
    scalar6->SetNumberOfComponents(1);

    vtkDoubleArray *scalar5 = vtkDoubleArray::New();
    scalar5->SetName("Pore_pressure");
    scalar5->SetNumberOfComponents(1);

    vtkIntArray *index = vtkIntArray::New();
    index->SetName("Index");
    index->SetNumberOfComponents(1);

    vtkIntArray *findex = vtkIntArray::New();
    findex->SetName("FragNumber");
    findex->SetNumberOfComponents(1);

    vtkDoubleArray *scalar7 = vtkDoubleArray::New();
    scalar7->SetName("Damage");
    scalar7->SetNumberOfComponents(1);

    vtkDoubleArray *scalar8 = vtkDoubleArray::New();
    scalar8->SetName("Chipping_removal");
    scalar8->SetNumberOfComponents(1);

    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    
    vtkZLibDataCompressor* myZlibCompressor  = vtkZLibDataCompressor::New();
    myZlibCompressor->SetCompressionLevel(1);


    dataSet->Allocate();


    for(i=0; i<nnopo; i++)
    {

    /* prepare output */
	pts->InsertNextPoint(d1nccx[i],d1nccy[i],d1nccz[i]);
      vects->InsertNextTuple3(d1nvcx[i],d1nvcy[i],d1nvcz[i]);
	  fcontact->InsertNextTuple3(d2nfc[0][i],d2nfc[1][i],d2nfc[2][i]);
        ploading->InsertNextTuple1(d1npre[i]);
        pdisp->InsertNextTuple1(d1ndisp[i]);
  }

	for(ielem=0; ielem<nelem; ielem++)
	{
    i1elto = i2elto[ielem];
    polys->InsertNextCell(4);
		for(i=0;i<4;i++)
        {
            polys->InsertCellPoint(i1elto[i]);
        }
       index->InsertNextTuple1(i1elto[4]);
        findex->InsertNextTuple1(i1elfrg[ielem]);

              tensor->InsertNextTuple9(d2tcs[ielem][0],d2tcs[ielem][3],d2tcs[ielem][4],
            d2tcs[ielem][3],d2tcs[ielem][1],d2tcs[ielem][5],
            d2tcs[ielem][4],d2tcs[ielem][5],d2tcs[ielem][2]);

            tensor2->InsertNextTuple9(d3tst[ielem][0][0],d3tst[ielem][0][1],d3tst[ielem][0][2],
            d3tst[ielem][1][0],d3tst[ielem][1][1],d3tst[ielem][1][2],
            d3tst[ielem][2][0],d3tst[ielem][2][1],d3tst[ielem][2][2]);

                scalar1->InsertNextTuple1(sigma1[ielem]);
//                scalar2->InsertNextTuple1(med);
                scalar3->InsertNextTuple1(sigma3[ielem]);
                scalar4->InsertNextTuple1(diff[ielem]);
                scalar5->InsertNextTuple1(T_pore_pressure[ielem]);
		scalar6->InsertNextTuple1(d1mart[ielem]);
                scalar7->InsertNextTuple1(d1eldam[ielem]);
                scalar8->InsertNextTuple1(i1elfrm[ielem]);

       }
        dataSet->SetPoints(pts);
        pts->Delete();

        dataSet->GetPointData()->AddArray(vects);
        dataSet->GetPointData()->SetActiveAttribute("VELOCITY",vtkDataSetAttributes::VECTORS);
        dataSet->GetPointData()->SetVectors(vects);
        vects->Delete();

        dataSet->GetPointData()->AddArray(fcontact);
        dataSet->GetPointData()->SetActiveAttribute("CONTACT",vtkDataSetAttributes::VECTORS);
        dataSet->GetPointData()->SetVectors(fcontact);
        fcontact->Delete();

    dataSet->GetPointData()->AddArray(ploading);
    dataSet->GetPointData()->SetActiveAttribute("Pressure_Lading",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(ploading);
    ploading->Delete();

    dataSet->GetPointData()->AddArray(pdisp);
    dataSet->GetPointData()->SetActiveAttribute("Nodal_displament",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(pdisp);
    pdisp->Delete();

        dataSet->GetCellData()->AddArray(tensor);
        dataSet->GetCellData()->SetActiveAttribute("STRESS",vtkDataSetAttributes::TENSORS);
        dataSet->GetCellData()->SetTensors(tensor);
        tensor->Delete();

        dataSet->GetCellData()->AddArray(tensor2);
        dataSet->GetCellData()->SetActiveAttribute("STRAINRATE",vtkDataSetAttributes::TENSORS);
        dataSet->GetCellData()->SetTensors(tensor2);
        tensor2->Delete();

    dataSet->GetCellData()->AddArray(scalar1);
    dataSet->GetCellData()->SetActiveAttribute("Sigma1",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar1);
    scalar1->Delete();

    dataSet->GetCellData()->AddArray(scalar3);
    dataSet->GetCellData()->SetActiveAttribute("Sigma3",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar3);
    scalar3->Delete();

    dataSet->GetCellData()->AddArray(scalar4);
    dataSet->GetCellData()->SetActiveAttribute("DiffStress",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar4);
    scalar4->Delete();

    dataSet->GetCellData()->AddArray(scalar5);
    dataSet->GetCellData()->SetActiveAttribute("Pore_pressure",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar5);
    scalar5->Delete();

    dataSet->GetCellData()->AddArray(scalar6);
    dataSet->GetCellData()->SetActiveAttribute("MassRatio",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar6);
    scalar6->Delete();

    dataSet->GetCellData()->AddArray(scalar7);
    dataSet->GetCellData()->SetActiveAttribute("Damage",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar7);
    scalar7->Delete();

    dataSet->GetCellData()->AddArray(scalar8);
    dataSet->GetCellData()->SetActiveAttribute("Chipping_removal",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar8);
    scalar8->Delete();

        dataSet->SetCells(VTK_TETRA,polys);
    dataSet->GetCellData()->AddArray(index);
    dataSet->GetCellData()->SetActiveAttribute("Index",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(index);
    index->Delete();

    dataSet->SetCells(VTK_TETRA,polys);
dataSet->GetCellData()->AddArray(findex);
dataSet->GetCellData()->SetActiveAttribute("FragNumber",vtkDataSetAttributes::SCALARS);
dataSet->GetCellData()->SetScalars(findex);
findex->Delete();

        polys->Delete();
        writer->SetInputData(dataSet);
        writer->SetFileName(fout);

        writer->SetCompressor(myZlibCompressor);
        writer->Write();
        dataSet->Delete();
        myZlibCompressor->Delete();
        writer->Delete();

    FREE(T_pore_pressure);
        FREE(diff);
        FREE(sigma3);
        FREE(sigma1);
  FREE(d2tcs);

}

static void Yod3TET4JOINTSINTACT(  /* 3D joint output for 4-node tetrahedra */
     INT nelem,
     FILE *fout,
     DBL dcsizc,DBL dcsizv,
     DBL dpefs, DBL dpeft, DBL dpegf, DBL dpeks, DBL dpepe,
     INT icoutp,INT iprop ,
     DBL *d1nccx,DBL *d1nccy,DBL *d1nccz,
     DBL *d1nvcx,DBL *d1nvcy,DBL *d1nvcz,
     DBL *d1sdel,
     INT *i1elpr,INT **i2elto
     )
{
}

static void Yod3TET4JOINTSBROKEN(  /* 3D joint output for 4-node tetrahedra */
     INT nelem, INT nnopo,
     CHR *fout,
     DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
     DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
//     DBL *d1sdel,
     INT *i1elpr, INT **i2elto, INT *i1elty, INT *i1defect,
     DBL **d2jst
     )
{    INT ielem;
     INT i,ibound;
     INT ipropc;
     INT *i1elto;

    vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
    vtkPoints *pts = vtkPoints::New();
    vtkCellArray *polys = vtkCellArray::New();


    vtkDoubleArray *scalar1 = vtkDoubleArray::New();
    scalar1->SetName("Fracture");
    scalar1->SetNumberOfComponents(1);

    vtkDoubleArray *scalar2 = vtkDoubleArray::New();
    scalar2->SetName("Microcrack");
    scalar2->SetNumberOfComponents(1);

    vtkDoubleArray *scalar3 = vtkDoubleArray::New();
    scalar3->SetName("StrainRateN");
    scalar3->SetNumberOfComponents(1);

    vtkDoubleArray *scalar4 = vtkDoubleArray::New();
    scalar4->SetName("StrainRateS");
    scalar4->SetNumberOfComponents(1);

    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    //    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();

    vtkZLibDataCompressor* myZlibCompressor  = vtkZLibDataCompressor::New();
    myZlibCompressor->SetCompressionLevel(1);


    dataSet->Allocate();

    for(i=0; i<nnopo; i++)
    {

	pts->InsertNextPoint(d1nccx[i],d1nccy[i],d1nccz[i]);
    }
     for(ielem=0;ielem<nelem;ielem++)
     {
          ibound=0;
          i1elto=i2elto[ielem];
	  ipropc=i1elpr[ielem];

	       if((i1elto[0]==i1elto[5])&&(i1elto[1]==i1elto[4])&&(i1elto[2]==i1elto[3]))
	       {
		    ipropc=ipropc-YIPROPMAX*2;
                    ibound=1;
	       }

            if (i1elty[ielem]>2)
        {
          d2jst[ielem][0]=0;
          d2jst[ielem][1]=0;
        }
               /* prepare output */


               if(ibound)
		 {
          polys->InsertNextCell(3);
		for(i=0;i<3;i++)
        {
            polys->InsertCellPoint(i1elto[i]);
        }
        scalar1->InsertNextTuple1(i1elty[ielem]);
             scalar2->InsertNextTuple1(i1defect[ielem]);
        scalar3->InsertNextTuple1(d2jst[ielem][0]);
        scalar4->InsertNextTuple1(d2jst[ielem][1]);

		 }
	       else		/* boundary surface */
		 {
          polys->InsertNextCell(3);
		for(i=0;i<3;i++)
        {
            polys->InsertCellPoint(i1elto[i]);
        }
        scalar1->InsertNextTuple1(i1elty[ielem]);
             scalar2->InsertNextTuple1(i1defect[ielem]);
        scalar3->InsertNextTuple1(d2jst[ielem][0]);
        scalar4->InsertNextTuple1(d2jst[ielem][1]);

        polys->InsertNextCell(3);
		for(i=0;i<3;i++)
        {
            polys->InsertCellPoint(i1elto[i+3]);
        }
        scalar1->InsertNextTuple1(i1elty[ielem]);
        scalar2->InsertNextTuple1(i1defect[ielem]);
        scalar3->InsertNextTuple1(d2jst[ielem][0]);
        scalar4->InsertNextTuple1(d2jst[ielem][1]);
		 }
    }

        dataSet->SetPoints(pts);
        pts->Delete();

        dataSet->SetCells(VTK_TRIANGLE,polys);

    dataSet->GetCellData()->AddArray(scalar1);
    dataSet->GetCellData()->SetActiveAttribute("fract",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar1);
    scalar1->Delete();

    dataSet->GetCellData()->AddArray(scalar2);
    dataSet->GetCellData()->SetActiveAttribute("microcrack",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar2);
    scalar2->Delete();

    dataSet->GetCellData()->AddArray(scalar3);
    dataSet->GetCellData()->SetActiveAttribute("StrainRateN",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar3);
    scalar3->Delete();

    dataSet->GetCellData()->AddArray(scalar4);
    dataSet->GetCellData()->SetActiveAttribute("StrainRateS",vtkDataSetAttributes::SCALARS);
    dataSet->GetCellData()->SetScalars(scalar4);
    scalar4->Delete();

        polys->Delete();
        writer->SetInputData(dataSet);
        writer->SetFileName(fout);
        writer->SetCompressor(myZlibCompressor);
        writer->Write();
        dataSet->Delete();
        myZlibCompressor->Delete();
        writer->Delete();


}

static void Yod3BROKEN_AREA( /* output broken joint area */

          FILE *fract,INT nelem,DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
          INT **i2elto,INT *i1elty,DBL dctime

          )
{

    INT ielem;
     INT i,ibound;
     INT ipropc;
     INT *i1elto;
     DBL nx,ny,nz,nx1,ny1,nz1,nx2,ny2,nz2,sum,area;//,sum1;
     sum=R0;
//     sum1=R0;
     for(ielem=0;ielem<nelem;ielem++)
     {
          i1elto=i2elto[ielem];
          if(i1elty[ielem]>2)
          {
nx1=d1nccx[i1elto[1]]-d1nccx[i1elto[0]];
ny1=d1nccy[i1elto[1]]-d1nccy[i1elto[0]];
nz1=d1nccz[i1elto[1]]-d1nccz[i1elto[0]];
nx2=d1nccx[i1elto[2]]-d1nccx[i1elto[0]];
ny2=d1nccy[i1elto[2]]-d1nccy[i1elto[0]];
nz2=d1nccz[i1elto[2]]-d1nccz[i1elto[0]];
V3DCro(nx,ny,nz,nx1,ny1,nz1,nx2,ny2,nz2);
V3DNor(area,nx,ny,nz);
sum=sum+ABS(area)/R2;


          }

     }
  DBLw(fract,dctime,10);
  CHRwsp(fract);
  DBLw(fract,sum,10);


  CHRwcr(fract);
}

static void Yod3CONTACTFORCE( /* output contact forces */

		FILE *fcnf,INT  nelem,INT nelno, DBL **d2nfcon,DBL dctime,
		INT *i1elbe, INT **i2elto, INT nnopo, INT layer_max,INT icouti

		)


{ INT ielem,i,n,inopo,ngroup;
  INT *i1nrep;
  DBL *dsumfx,*dsumfy,*dsumfz;

  ngroup=layer_max+2;
  i1nrep=TalINT1(nnopo+1);
  dsumfx=TalDBL1(ngroup);
  dsumfy=TalDBL1(ngroup);
  dsumfz=TalDBL1(ngroup);

  for(i=0;i<nnopo;i++)
  {
   i1nrep[i]=0;
  }

  for(i=0;i<ngroup;i++)
  {
   dsumfx[i]=0;
   dsumfy[i]=0;
   dsumfz[i]=0;

  }


  for(ielem=0;ielem<nelem;ielem++)

  { n=i1elbe[ielem];
	  if((n>0)&&(n<=layer_max))

        {
		for(i=0;i<nelno;i++)
		{
			inopo=i2elto[ielem][i];
			if(i1nrep[inopo]==0)
                        {
		dsumfx[n]=dsumfx[n]+DABS(d2nfcon[0][inopo]);
		dsumfy[n]=dsumfy[n]+DABS(d2nfcon[1][inopo]);
		dsumfz[n]=dsumfz[n]+DABS(d2nfcon[2][inopo]);
		i1nrep[inopo]=1;

			}
		}
  }
	  else if(n>layer_max)

        {
                for(i=0;i<nelno;i++)
                {
                        inopo=i2elto[ielem][i];
                        if(i1nrep[inopo]==0)
                        {
                dsumfx[layer_max+1]=dsumfx[layer_max+1]+DABS(d2nfcon[0][inopo]);
                dsumfy[layer_max+1]=dsumfy[layer_max+1]+DABS(d2nfcon[1][inopo]);
                dsumfz[layer_max+1]=dsumfz[layer_max+1]+DABS(d2nfcon[2][inopo]);
                i1nrep[inopo]=1;

                        }
                }
  }

  }
  if(icouti==0)
{
fprintf(fcnf,"time");
  for(i=1;i<ngroup;i++)
  {
     fprintf(fcnf,"     group_%ld_fx    group_%ld_fy    group_%ld_fz",i, i ,i);
  }

  fprintf(fcnf,"\n");
}
  fprintf(fcnf,"   %f  ", dctime);

  for(i=1;i<ngroup;i++)
  {
     fprintf(fcnf,"   %f    %f   %f",dsumfx[i], dsumfy[i] ,dsumfz[i]);
  }
  fprintf(fcnf,"\n");


FREE(i1nrep);
FREE(dsumfx);
FREE(dsumfy);
FREE(dsumfz);

}



/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief output the Y database
 *  \param[in] namep base file name
 *  \param[in] yd pointer to Y database structure
 *  \par Details:
 *  Yod() writes to the history file every 100 time steps, and outputs frames
 *  according to the model output frequency.
 *
 *  Z - needs polished
 *
 *  \see #YDC_struct.icoutf
 */
void Yod(CHR *namep, YD yd)
{
  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDK ydk=&(yd->ydk);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YSP ysp=&(yd->ysp);
  YDJ ydj=&(yd->ydj);

  INT ihys;
  CHR namef[300];
  CHR cindex[50];
  static INT ncall=0;
  static INT ndefinitions=0;
  INT ifract,i;
  FILE *fcnf=FILENULL;
  FILE *fract=FILENULL;
  FILE *ffrag=FILENULL;

  DBL tmp;

  ncall=ncall+1;
  ifract=0;
  /* output history */
  if((ydc->ncstep%100)==0)
  {
    for(ihys=0;ihys<(ydo->nohys);ihys++)
    {
      ydo->d1ohyt[ihys]=ydc->dctime;
	    if((ydo->d1ohyt[ihys])==(ydc->dctime))
      {
        if((ydo->f2ohyf[ihys])==FILENULL)
        {
          CHRcpy(namef,namep);
          SINTw(cindex,ihys,0);
          CHRcat(namef,"h");
          CHRcat(namef,cindex);
          ydo->f2ohyf[ihys]=fopen(namef,"a");
        }
        if((ydo->f2ohyf[ihys])!=FILENULL)
        {
          tmp=(ydo->d1ohyt[ihys])*(ydo->d1ohyc[ihys]);
          DBLw((ydo->f2ohyf[ihys]),tmp,17);
          CHRwsp(ydo->f2ohyf[ihys]);
          tmp=(ydo->d1ohys[ihys])*(ydo->d1ohyf[ihys]);
          DBLw((ydo->f2ohyf[ihys]),tmp,17);
          CHRwcr(ydo->f2ohyf[ihys]);
        }
      }
    }
  }

  /* note - odd to hard-code this */
  if((ncall>100)||((ydc->ncstep)>=(ydc->mcstep-2)))
  {
    ncall=0;
    for(ihys=0;ihys<(ydo->nohys);ihys++)
    {
      if((ydo->f2ohyf[ihys])!=FILENULL) fclose(ydo->f2ohyf[ihys]);
      ydo->f2ohyf[ihys]=FILENULL;
    }
  }

  /* output animation */
  if((ydc->ncstep%ydc->icoutf)==0)
  {
    CHRcpynoext(namef,namep);
    SINTw(cindex,ydc->icouti,0);
    CHRcat(namef,cindex);
    CHRcat(namef,".vtu");

        if(ydp->i1ptyp[0] == YTE3TET10ELS) /* small strain elastic tetrahetra  */
        {

            Yod3TET10(
                      yde->nelem,ydn->nnopo, namef,
                      ydp->d1peks, ydp->d1pela, ydp->d1pemu,
                      ydn->d2ncc, ydn->d2nci, ydn->d2nfc, ydn->d2nvc,
                      yde->i1elpr, yde->i2elto,
                      ydk->d3dsh,ysp->nspd,ydn->i1nind
                     );

        }
        else if((ydp->i1ptyp[0])==(YTE3TET4ELS))
        {
            Yod3TET4(
                     yde->nelemi, ydn->nnopi,ydn->nnopo,
                     namef,ydp->d1peks,ydp->d1pela,ydp->d1pemu,
                     ydn->d2ncc[0],ydn->d2ncc[1],
                     ydn->d2ncc[2],ydn->d2nci[0],ydn->d2nci[1],
                     ydn->d2nci[2],ydn->d2nfc,ydn->d2nvc[0],
                     ydn->d2nvc[1],ydn->d2nvc[2],yde->i1elpr,
		     yde->i2elto,ydp->i1ptyp,yde->i1elfrg,yde->d1mart,
                     yde->d1eldam, ydn->d1npre,yde->d3tst,yde->i1elfrm,ydn->d1ndisp
                     );

        }

 for(i=0;i<ydp->nprop;i++)
{
        if((ydp->i1ptyp[i])==(YTE3TET4JOINT))ifract=1;
}
		    if(ifract) {
   CHRcpynoext(namef,namep);
    SINTw(cindex,ydc->icouti,0);
    CHRcat(namef,"_crack");
    CHRcat(namef,cindex);
    CHRcat(namef,".vtu");


			  Yod3TET4JOINTSBROKEN(
			       ydj->nelem, ydn->nnopo,
			       namef,
			       ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
			       ydn->d2nvc[0],ydn->d2nvc[1],ydn->d2nvc[2],
			       ydj->i1elpr,ydj->i2elto,ydj->i1elty,ydj->i1defect,ydj->d2jst
			       );

}
  fcnf=FILENULL;
  fract=FILENULL;
	  if(ydc->icouti==0)

	  { fcnf=fopen("contactforce.txt","w");
            fract=fopen("fract_area.txt","w");
            CHRw(fract,"Time        Total crack area");
            CHRwcr(fract);
          ffrag=fopen("fragment.txt","w");

	  }

	  else

	  { fcnf=fopen("contactforce.txt","a");
            fract=fopen("fract_area.txt","a");
          ffrag=fopen("fragment.txt","a");


	  }


      Yod3CONTACTFORCE(

	  fcnf,

          yde->nelemi,yde->nelno,ydn->d2nfc,ydc->dctime,yde->i1elbe,yde->i2elto,ydn->nnopo,yde->layer_max,ydc->icouti

	  );
          Yod3BROKEN_AREA(

          fract,ydj->nelem,ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
          ydj->i2elto,ydj->i1elty,ydc->dctime

          );
      CHRw(ffrag,"time   "); DBLw(ffrag,ydc->dctime,18);
      CHRwcr(ffrag);
      CHRw(ffrag,"fragment no       mass"); CHRwcr(ffrag);
      for(i=0;i<yde->nfrag;i++)
      {
          INTw(ffrag,i,10); CHRwsp(ffrag);
          DBLw(ffrag,yde->d1volfrg[i],18);
	  INTw(ffrag,yde->i1frgrm[i],5);
	  CHRwsp(ffrag);
	  CHRwcr(ffrag);
      }



	  fclose(fcnf);
          fclose(fract);
      fclose(ffrag);

      ydc->icouti=ydc->icouti+1;

  }
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/
