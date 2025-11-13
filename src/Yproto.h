/*! \file Yproto.h
 *  \brief Y main prototypes
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


#ifndef YPROTO_H
#define YPROTO_H


#include "Y3Dd.h"


/**********************************************************************/
/**********************************************************************/

void Ymd(        /* mesh elements, added by LG */
 YDE yde,        /* element database */
 YDI ydi,        /* interaction database */
 YDN ydn,        /* nodal database */
 YDP ydp,         /* property database */
 YDX ydx,          /* surface element database */
 YDJ ydj,          /* joint element database*/
 YDJP ydjp
); 

void Ycd(        /*   contact detection                               */
 YDE yde,         /* = element database                               */
 YDI ydi,         /* = interaction database                           */
 YDN ydn,         /* = nodal database                                 */
 INT iwgast,
 int nthreads,
INT ncstep,
INT *i1elst
// YDP ydp         /* - property database                               */
);

void Yfd(        /*   nodal forces                                    */
 YDE yde,         /* = element database                               */
 YDN ydn,         /* = nodal database                                 */
 YDP ydp,         /* - property database                              */
 YDK ydk,          /* = constants database                             */
 YDI ydi,         /* = interaction database                           */
 YDC ydc,         /* - control database                               */
  YDX ydx,         /* - surface database                               */
  YDJ ydj,
  YDJP ydjp,
   YDB ydb,         /* - boundary condition database                    */
 int nthreads
);

void Yid(        /*   procreate                                       */
 YDE yde,         /* = element database                               */
 YDI ydi,         /* = interaction database                           */
 YDN ydn,         /* = nodal database                                 */
 YDP ydp,          /* - property database                              */
 YDJ ydj,
 int nthreads,
 DBL dcstec,
INT ncstep

);

void Ysd(        /*  solve equations                                  */
 YDC ydc,         /* - control database                               */
 YDE yde,         /* = element database                               */ 
 YDN ydn,         /* = nodal database                                 */
 YDO ydo,         /* = output database                                */
 YDB ydb,         /* - boundary condition database                    */
 YDK ydk          /* = constants database                             */
);

int Ycdran(INT nelem_start,
            INT nelem,
            
            DBL *d1nccx, DBL *d1nccy, DBL *d1nccz, 
            INT *i1elbe, INT **i2elto, DBL diam,
            INT nnopo,
            int nthreads,
			INT inew
            );

INT Yrd(         /*  read input                                       */
 CHR *namep,      /* - name of the problem i.e. input file            */
 CHR *name_grid,
 CHR *name_mesh,
 CHR *name_mat,
 YD  yd           /* = database                                       */
);
 
void Yod(         /* output results in space-save format              */
 CHR *namep,      /* - name of the problem i.e. input file            */
 YD  yd           /* = database                                       */
);

void Yinit(       /* constants database initialization                */
 YDK ydk          /* = database                                       */
);
void Ymass3TET(
                YDE yde,
                YDN ydn,
                YDP ydp,
                YDK ydk,
		YDJ ydj,
		YDJP ydjp,
                int nthreads,
		DBL dcstec
                );
void Yconfig(     /* constants database configuration                 */
 YDK ydk,         /* = database                                       */
 YDC ydc,         /* = control database                               */
 YDB ydb          /* = boundary condition database                    */
);

void Ywrs(         /*  read input                                       */
 CHR *namep,      /* - name of the problem i.e. input file            */
 YD  yd           /* = database                                       */
);

void Ywrs(         /*  read input                                       */
          CHR *namep,      /* - name of the problem i.e. input file            */
          YD  yd           /* = database                                       */
);


void  mirtichRoutine(DBL **points,
                      DBL *intg, 
                      DBL dro);
/*
int Ycd3TETintersection( double *d1nccx, double *d1nccy, double *d1nccz, 
						INT iele, INT jele, INT **i2elto,double buf);
*/
void Ypacking(         /* packing              */
  YD  yd, int nthreads           /* = database                                       */
);

void V3DRot(
     DBL *afx, DBL *afy, DBL *afz, DBL bfx, DBL bfy,
     DBL bfz, DBL u, DBL v, DBL w, DBL cosin, DBL sine
);

void solveSymetricEigenProblem(DBL **m, DBL  *v1, DBL *v2,
							   DBL *v3, DBL *a1,DBL *a2,DBL *a3);

void V3DRotAxes(DBL *afx,DBL *afy,DBL *afz,DBL bfx,DBL bfy,DBL bfz,DBL u,DBL v,DBL w);
void Yhpd(YDX ydx, YDN ydn,YDC ydc);




/**********************************************************************/
/**********************************************************************/


#endif /* YPROTO_H */


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/
