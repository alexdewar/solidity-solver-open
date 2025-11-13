/*! \file Y3Dd.h
 *  \brief Y data base description
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


#ifndef Y3DD_H
#define Y3DD_H


#include "Ytypes.h"
#include "frame.h"
#include <omp.h>

/**********************************************************************/
/**********************************************************************/


/*! constants structure pointer
*/
typedef struct YDK_struct *YDK; 

/*! constants structure
*/
struct YDK_struct    
{ 
  DBL **d2csmm;  /*!< [NNODE][NNODE] consistent mass matrix - Ysd()  */
  DBL *d1lpmm;   /*!< [NNODE] lumped mass matrix - Yfd()             */
  DBL ***d3dsh;  /*!< [NGRSH][NDIME][NNODE] shape function - Yfd()   */
  DBL **d2ka;    /*!< [mbcon] = ydb->d1bnaXYZ[ibcon] + ydc->dcgraXYZ */
};


/**********************************************************************/
/**********************************************************************/


/*! control structure pointer
*/
typedef struct YDC_struct *YDC; 

/*! control structure
*/
struct YDC_struct    
{ 
  INT  mcstep;  /*!< maximum number of time steps */
  INT  ncstep;  /*!< current number of time steps */
  FILE *finp;
  FILE *fcheck;
  
  DBL  dcgrax;  /*!< gravity x */
  DBL  dcgray;  /*!< gravity y */
  DBL  dcgraz;  /*!< gravity z */

  DBL  dcgrst;				/* gravity settling time				  */
  DBL  dcrmpt;				/* ramping time for loading application   */
  DBL  dcrmjt;                          /* ramping time for loading for jet drilling   */
  DBL  dcstdn;  /*jet stand-off-distance */
  DBL  dcrmct; /* constant time for loading application   */
  DBL  dcrmut; /* unloading ramping time for loading application   */
  DBL  dcstec;  /*!< current time step size */
  DBL  dctime;  /*!< current time */
  DBL  dcurelx; /*!< under relaxation for mass matrix */
  DBL  dcfact1;  /*!< intergrain factor */
  DBL  dcfact2;  /*!< intergrain factor */
  DBL  dcfact3;  /*!< intergrain factor */
  DBL  dcfact4;  /*!< intergrain factor */
  DBL  dcfact5;  /*!< intergrain factor */
  DBL  dcimvelx; /*!< piston impact velocity */
  DBL  dcimvely; /*!< piston impact velocity */
  DBL  dcimvelz; /*!< piston impact velocity */ 


  INT  initer;  /*!< number of iterations for multi-pass algorithm */
  INT  icoutf;  /*!< write output frequency */
  INT  icoutf2; /*!< write output frequency */
  INT  icouti;  /*!< current write output No */
  INT  iwfast;  /*!< fast mode */
  INT  isave;   /*!< restart number of loop */
  INT  ifrag;  /* 1, activate fragment group as new crack found */
  INT  ichprm;  /* 1, activate chipping removal */


//  DBL dcmart; /* mass ratio coefficient  */////////
  DBL r;

  INT icfrm ; /* activate fragment removal */
  INT ncstepfrm; /*timestep fragment removal */
  DBL dcrmgrax;
  DBL dcrmgray;
  DBL dcrmgraz;
  DBL dcfrmca;
  DBL dcrmvx;
  DBL dcrmvy;
  DBL dcrmvz;
};


/**********************************************************************/
/**********************************************************************/


/*! element description structure pointer
*/
typedef struct YDE_struct *YDE;  

/*! element description structure
*/
struct YDE_struct
{ 
  INT melem; /*!< maximum number of elements                     */
  INT nelem; /*!< actual number of elements                      */
  INT melst; /*!< maximum number of elemen. states var.          */
  INT nelst; /*!< actual number of elemen. states var.           */
  INT melno; /*!< maximum number of elemen. nodes                */
  INT nelno; /*!< actual number of elemen. nodes                 */
  INT *i1elcf;   /*!< [melem] contacting couple first            */ 
  INT *i1elpr;   /*!< [melem] element property                   */
  INT layer_max;  /*    maximum laye number                     */
//  INT *i1elty;   /*[melem]	  element type						 */
				 /*3-5:normal, shear, shear rotation; 
					2:pre-existing; 1:boundary;		 			*/
				 /*0:unbroken; -1:non-joint elemnt				 */
  DBL **d2elst;  /*!< [melst][melem] element state               */
  INT **i2elto;  /*!< [melem][melno] element topology            */
    DBL **d2elfc; /*elemental force */
   DBL ***d3tst;       /*!<[melem][ndime][ndime]  element strain rate            */
  DBL *d1emct;   /*!< [melem] total elemental mass               */
  INT *i1elbe;   /*!< [melem] boundary element                   */

  INT **i2elfr;   /*!< [melem][4] element fracture - contact label                   */

//  INT *i1eljo;	 /*[melem] joint element property                */
  INT **i2eljp;  /*!< [melem][melno] element topology            */

    INT nelemi;
    INT nnelm;
    INT *i1elem;
 INT *i1elcft;			/* [melem] fem elem failure timestep due to compression */
INT njoint_intact;
INT *i1joint;
INT ifrac;
    
    INT *i1elfrg;
    INT nfrag;
    DBL *d1volfrg;

    INT *i1frgrm;
    DBL *d1elvol;
    DBL *d1mart;
    INT *i1elchprm;  /* 1, activate element chipping removal */

    INT *i1eljt;     /* -1 no new fracture joint broken, 1 new fracture joint brokn  */
    INT *i1elncstep;  /* timestep joint broken*/

    DBL *d1pem0n; /*damage model m0n for each element */
    DBL *d1eldam; /*damage model damage value */

    INT *i1elfrm; /*identify fragment removal, 0 remains 1 remove */
};


/**********************************************************************/
/**********************************************************************/


typedef struct YDJ_struct *YDJ;

/*! interaction structure
*/

struct YDJ_struct
{
    INT melem; /*!< maximum number of elements                     */
    INT nelem; /*!< actual number of elements                      */

    INT melno; /*!< maximum number of elemen. nodes                */
    INT nelno; /*!< actual number of elemen. nodes                 */
    INT *i1elpr;   /*!< [melem] element property                   */
    INT *i1elty;   /*[melem]      element type                         */
                   /*3-5:normal, shear, shear rotation;
                      2:pre-existing; 1:boundary;                     */
                   /*0:unbroken; -1:non-joint elemnt                 */
    INT *i1elst;  /* state of joint, -1 intact, 0 just broken, 1 broken*/
    INT **i2elto;  /*!< [melem][melno] element topology            */
  INT **i2ejjp;  /*!< [melem][melno] element topology            */
  INT *i1eljo;   /*[melem] joint element property                */
    
    DBL **d2ejfc; /* Joint element force  */

    DBL *d1ejmr;  /*mass ratio to adjust spring stffness for joint */  
    
    
  INT   njoint;		  /* number of joint elements				*/

  INT  *i1jtid;		  /* joint element index in element topology list 	*/
  DBL  *d1jkni;		  /* [mjelem] joint initial normal stiffness	  	*/
  DBL  *d1jknc;		  /* [mjelem] joint current normal stiffness	  	*/

  DBL  *d1jnst;		  /* [mjelem] joint normal stress			*/
  DBL  *d1jsst;		  /* [mjelem] joint shear stress			*/

  DBL  *d1japi;		  /* [mjelem] joint initial normal aperture		*/
  DBL  *d1japc;		  /* [mjelem] joint current normal aperture		*/
  DBL  *d1japh;		  /* [mjelem] joint hydraulic aperture			*/

  DBL  *d1jefl;		  /* [mjelem] joint effective fracture length	  	*/

  DBL  *d1jfmd;		  /* [mjelem] joint (if new crack) failure mode   	*/
  DBL **d2jnrm;		  /* [mjelem][3] joint element unit normal		*/

  DBL  *d1jtar;		  /* [mjelem] joint element area			*/
    
    INT *i1defect;   /* indicate if the joint is microcrack, 0 intact, 1 microcrack*/
    DBL defect;     /* percentage of defects*/
    DBL dran;      /*random number seeds*/

 DBL  **d2jst;     /*!<[mjelem][2]  joint strain rate          */
}; 

/**********************************************************************/
/**********************************************************************/


/*! interaction structure pointer
*/
typedef struct YDI_struct *YDI; 

/*! interaction structure
*/
struct YDI_struct
{ 
  INT micoup;  /*!< maximum possible number of contacting couples  */
  INT nicoup;  /*!< actual number of contacting couples            */
  INT iiecff;  /*!< interaction element contact. couple free first */  
  DBL  diedi;  /*!< travel since last detection                    */
  DBL  diezon; /*!< buffer zone size                               */
  DBL *d1iesl; /*!< [mcoup] contact sliding                        */
  INT *i1iecn; /*!< [mcoup] couple next                            */
  INT *i1iect; /*!< [mcoup] couple target                          */
  DBL   *d1deltat1;   /*!<[mcontn] tangential overlap between two elements */
  DBL   *d1deltat2;   /*!<[mcontn] tangential overlap between two elements */
  DBL   *d1deltan;   /*!<[mcontn] normal overlap between two elements */
  DBL   **d2nv;
  DBL   **d2t1v;
  DBL   **d2t2v;
    DBL drect;
    INT *i1stat;
    int icheck;
    /* -1 no contact pair  1 new contact pair 2 old contact pair
     3 old contact pair may not in contact 5 need remove later 
     10 in contact once, but now separate, still in buff zone
     11 in conact now
     */

  INT *i1fcstep; 		/* [mcoup] first timestep couple contact due to fracture*/
    INT *i1fcjp;        /*  * [mcoup] joint no link      */
    DBL *d1pepe;        /* [mcoup] penalty number for contact from broken joint*/
}; 


/**********************************************************************/
/**********************************************************************/


/*! node description structure pointer
*/
typedef struct YDN_struct *YDN;   

/*! node description structure
*/
struct YDN_struct
{ 
  INT mnodim;  /*!< maximum nodal dimensions number                    */
  INT nnodim;  /*!< actual nodal dimensions number                     */
  INT mnopo;   /*!< maximum number of nodal points                     */
  INT nnopo;   /*!< actual number of nodal points                      */
  INT nnopi;   /*!< actual number of nodal points                      */
  DBL *d1nmct; /*!< [mnopo] nodal mass current translation             */
  DBL **d2ncc; /*!< [mnodim][mnopo] nodal coordinate current           */
  DBL **d2nci; /*!< [mnodim][mnopo] nodal coordinate initial           */
  DBL **d2nfc; /*!< [mnodim][mnopo] nodal force current due to contact */ 
  DBL **d2nft; /*!< [mnodim][mnopo] nodal total current force          */
  DBL **d2nvc; /*!< [mnodim][mnopo] nodal velocity current             */
  INT *i1nobf; /*!< [mnopo] nodal boundary >0 is boundary              */
  INT *i1nopr; /*!< [mnopo] nodal property                             */
  INT *i1nind;
    
  INT *i1noel; /*!< [mnopo] nodal property  to elem                         */
    
    DBL **d2nck; /*!< [mnodim][mnopo] nodal coordinate last check contact   */
    DBL **d2ndf;

  DBL *d1nvct; /*2 d1nvct [mnopo] nodal volume                     */


  INT *i1ntoC2D; /* CON - DISC topologu */
    DBL *d1nvol;
    INT nnobf;
    INT *i1nobfse;

    INT *i1nold; /* link discountious node to old index*/
    DBL **d2nfp; /* < [mnodim][mnopo] nodal force current due to pressure condition*/
    INT *i1nochprm;  /* 1, activate noal chipping removal */
    DBL *d1npre;   /*pressure loading*/

    DBL **d2nfrm; /*!< [mnodim][mnopo] nodal coordinate before activate fragement removal   */
    DBL *d1ndisp;
    INT *i1rmgrav;
};


/**********************************************************************/
/**********************************************************************/


/*! output description structure pointer
*/
typedef struct YDO_struct *YDO;   

/*! output description structure
*/
struct YDO_struct
{ 
  INT mohys;     /*!<  maximum number of hystory variables             */
  INT nohys;     /*!< actual number of hystory variables               */
  DBL dohyp;     /*!< output hystory accuracy                          */ 
  DBL *d1ohyf;   /*!< [mohys] output hystory factor to scale state     */
  DBL *d1ohyc;   /*!< [mohys] output hystory factor to scale time      */
  DBL *d1ohys;   /*!< [mohys] output hystory state                     */
  DBL *d1ohyt;   /*!< [mohys] output hystory time                      */
  DBL *d1ohyx;   /*!< [mohys] output history x coordinate of the point */
  DBL *d1ohyy;   /*!< [mohys] output history y coordinate of the point */
  DBL *d1ohyz;   /*!< [mohys] output history z coordinate of the point */
  FILE **f2ohyf; /*!< [mohys] output history files                     */
  INT *i1ohyt;   /*!< [mohys] output hystory type, i.e. which variable */
};




/**********************************************************************/
/**********************************************************************/


/*! Super particle description structure pointer
*/

typedef struct YSP_struct *YSP;
struct YSP_struct        /*Y_Super-Particle_Dbs_Structure                     */
{ INT   mspd,nspd; /*maximum (actual) number of particle*/


  DBL   *d1mass;
  DBL   **d2pp;    /*[nnodim][msd] sphere_position_x                     */
  DBL   **d2poldp;    /*[nnodim][msd] sphere_position_x                     */
  DBL   **d2pvt;     /*[nnodim][msd] sphere_velocity_x              */


  DBL   **d2ptfc;    /*[nnodim][msd] sphere_total_contact force_x */
  DBL   **d2ptf;    /*[nnodim][msd] sphere_total_force_x */
  DBL   **d2ptfd;
  DBL   **d2pa;     /*[nnodim][msd] sphere_rotation angel_x      */
  DBL   **d2paold;     /*[nnodim][msd] sphere_rotation angel_x      */
  DBL   **d2pvr;     /*[nnodim][msd] sphere_rotation velocity_x      */
  DBL   **d2pm;     /*[nnodim][msd] sphere_moment_x */
  DBL   **d2pmc;     /*[nnodim][msd] sphere_moment_x due to contact */
  DBL   **d2prine;  /*[nnodim] pricinple inertia                         */
  DBL   **d2prinn1;  /*[nnodim] pricinple axe                         */
  DBL   **d2prinn2;  /*[nnodim] pricinple axe                         */
  DBL   **d2prinn3;  /*[nnodim] pricinple axe                         */
  DBL   **d2poldvr;     /*[nnodim][msd] sphere_rotation velocity_x      */
  INT   *i1sppr;    /* condition number */
  INT   *i1con;
  DBL *d1impv;
};



/**********************************************************************/
/**********************************************************************/


/*! property description structure pointer
*/
typedef struct YDP_struct *YDP;  

/*! property description structure
*/
struct YDP_struct
{ 
  INT mprop;   /*!< maximum number of properties                         */
  INT nprop;   /*!< actual number of properties                          */
  DBL *d1pefr; /*!< [mprop] property friction coefficient                */
  DBL   *d1pesf; /*!< d1pesf [mprop] property static friction coefficient            */
  DBL   *d1pepsf; /*!< d1pepsf [mprop] property pressure friction coefficient         */
  DBL   *d1pevf;  /*!< d1pefr [mprop] property velcocity dependent friction sfactor   */
  DBL   *d1pepf;  /*!< d1pefr [mprop] property pressure dependent friction sfactor   */

  DBL   *d1pefs;  /*[mprop] ultimate shear strength at joint               */
  DBL   *d1peft;  /*[mprop] ultimate tensile strength at joint             */
  DBL   *d1pegfn;  /*[mprop] ultimate fracture energy mode I at joint            */
  DBL   *d1pegfs;  /*[mprop] ultimate fracture energy mode II at joint            */

  DBL *d1peks; /*!< [mprop] dpeks=2hbeta*sqrt(E*ro) in 2D or 3D,0<beta<1 */
  DBL *d1pela; /*!< [mprop] property lamda - Lame elastic constant       */
  DBL *d1pemu; /*!< [mprop] property mu    - Lame elastic constant       */
  DBL *d1pepe; /*!< [mprop] property penalty parameter                   */
  DBL *d1pero; /*!< [mprop] property ro    - density                     */
  INT *i1ptyp; /*!< [mprop] property type                                */

  INT   *i1pejp;  /*[mprop] joint property; if<0, no joints              */
  INT   *i1pemn;  /*[mprop] number of mesh refinements                   */
  INT   *i1psde;  /*[mprop] state damage elastic id  */



//-PY changed it for 3D_fracture_coupling_with_multiphase

  DBL   *d1picf;  /*[mprop] internal friction angle                      */
  DBL   *d1pcoh;  /*[mprop] cohesion                        */
  DBL   *d1pemr;  /* mass ratio coefficient       */
  INT *i1index; /*[mprop] material number index                    */

//  DBL *d1pem0n;  /*[mprop] linear damage new delm0 */
  DBL *d1pem0;   /*[mprop] linear damage delm0 */
  DBL *d1pemf;   /*[mprop] linear damage delmf */
  DBL *d1pedm;   /*[mprop] linear damage threshold */

};

/**********************************************************************/
/**********************************************************************/


/*! property description structure pointer
*/
typedef struct YDJP_struct *YDJP;

/*! property description structure
*/
struct YDJP_struct
{
  INT mprop;   /*!< maximum number of properties                         */
  INT nprop;   /*!< actual number of properties                          */
    INT *i1ngrp;
  

 // DBL   *d1pefs;  /*[mprop] ultimate shear strength at joint               */
  DBL   *d1peft;  /*[mprop] ultimate tensile strength at joint             */
  DBL   *d1pegfn;  /*[mprop] ultimate fracture energy mode I at joint            */
  DBL   *d1pegfs;  /*[mprop] ultimate fracture energy mode II at joint            */

  
  DBL *d1pepe; /*!< [mprop] property penalty parameter                   */


  INT   *i1pejp;  /*[mprop] joint property; if<0, no joints              */


  DBL   *d1picf;  /*[mprop] internal friction angle                      */
  DBL   *d1pcoh;  /*[mprop] cohesion                        */

};

/**********************************************************************/
/**********************************************************************/


/*! boundary description structure pointer
*/
typedef struct YDB_struct *YDB;

/*! boundary description structure
*/
struct YDB_struct
{ 
  INT mbcon;   /*!< maximum number of properties            */
  INT nbcon;   /*!< actual number of properties             */
  DBL *d1bnax; /*!< [mprop] amplitude of acceleration x     */
  DBL *d1bnay; /*!< [mprop] amplitude of acceleration y     */
  DBL *d1bnaz; /*!< [mprop] amplitude of acceleration z     */
  DBL *d1bnfx; /*!< [mprop] amplitude of force x            */
  DBL *d1bnfy; /*!< [mprop] amplitude of force y            */
  DBL *d1bnfz; /*!< [mprop] amplitude of force z            */
  DBL *d1bnvx; /*!< [mprop] amplitude of velocity x         */
  DBL *d1bnvy; /*!< [mprop] amplitude of velocity y         */
  DBL *d1bnvz; /*!< [mprop] amplitude of velocity z         */
  INT *i1bnvx; /*!< [mprop] fixity x direction >0 vel.      */ //Z are these hack-ish?
  INT *i1bnvy; /*!< [mprop] fixity y direction >0 vel.      */ //Z are these hack-ish?
  INT *i1bnvz; /*!< [mprop] fixity z direction >0 vel.      */ //Z are these hack-ish? 
};

typedef struct YDX_struct *YDX;

struct YDX_struct
{  INT nelem_t;
   INT melem_t;
   INT nelno;
   INT **i2elto_t;
    DBL *d1nap;  /*pressure loading */
    DBL *d1nbp;  /* back pressure */
    DBL dbkp;
    INT *i1jtp; /*HPWJ pressure boundary*/
};

typedef struct YPAR_struct *YPAR;

/*! Particle description structure
*/
struct YPAR_struct
{ 
  INT ntpar;
  INT ngrid;
  INT nlpar;
  INT nindex;
  INT *i1elbe;
    INT ifeed;
    
    DBL depora;  /*deposition_radius */
    DBL depheg;  /*deposition_height  */
    DBL depoff;  /*deposition_offset  */


  INT nelem;   /*!< maximum number of properties            */
    INT nelemb;   /*!< maximum number of properties            */
  INT nnode;   /*!< actual number of properties             */
  INT **i2elto;
    INT **i2eltob;
  DBL **d2xyz; /*!< [mprop] amplitude of acceleration x     */
  DBL **d2gxyz; /*!< [mprop] amplitude of acceleration x     */
  DBL *d1cg;
  DBL *d1rv_max;
  DBL *d1rv_min;
  DBL *d1dst;
  DBL radius;
  DBL max;
  DBL *d1lv;
  INT ir;
  DBL r;
  DBL *d1vc; /*!< [mprop] current velocity     */
  
};


/**********************************************************************/
/**********************************************************************/


/*! Y3D database structure pointer
*/
typedef struct YD_struct *YD;

/*! Y3D database structure
*/
struct YD_struct
{ 
  struct YDK_struct ydk; /*!< constants structure */
  struct YDC_struct ydc; /*!< control structure */
  struct YDE_struct yde; /*!< element structure */
  struct YDI_struct ydi; /*!< interaction structure */
  struct YDN_struct ydn; /*!< node structure */
  struct YDO_struct ydo; /*!< output structure */
  struct YDP_struct ydp; /*!< properties structure */
  struct YDB_struct ydb; /*!< boundary structure */
  struct YPAR_struct ypar; /*!< boundary structure */
  struct YSP_struct ysp; /*!< boundary structure */
  struct YDX_struct ydx; /*!< surface structure */
  struct YDJ_struct ydj; /*!< joint element structure */
  struct YDJP_struct ydjp; /*!< Joint properties structure */

};


/**********************************************************************/
/**********************************************************************/


#endif /* Y3DD_H */


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/



