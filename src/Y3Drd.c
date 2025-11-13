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


/* read control data    
*/
static void Yrdc(YDC ydc, FILE *finp, CHR *name)
{ 
  CHR *namep;
  
  namep=name+8;
  if(CHRcmp(namep,"MCSTEP",6)==0) 
  { 
    INTr(finp,&(ydc->mcstep));
  }
  else if(CHRcmp(namep,"NCSTEP",6)==0) 
  { 
    INTr(finp,&(ydc->ncstep));
  }
  else if(CHRcmp(namep,"DCGRAX",6)==0) 
  { 
    DBLr(finp,&(ydc->dcgrax));
  }
  else if(CHRcmp(namep,"DCGRAY",6)==0) 
  { 
    DBLr(finp,&(ydc->dcgray));
  }
  else if(CHRcmp(namep,"DCGRAZ",6)==0) 
  { 
    DBLr(finp,&(ydc->dcgraz));
  }
  else if(CHRcmp(namep,"DCSTEC",6)==0) 
  { 
    DBLr(finp,&(ydc->dcstec));
  }
  else if(CHRcmp(namep,"DCTIME",6)==0) 
  { 
    DBLr(finp,&(ydc->dctime));
  }
  else if(CHRcmp(namep,"DCGRST",6)==0)
  { DBLr(finp,&(ydc->dcgrst));
  }
  else if(CHRcmp(namep,"DCRMPT",6)==0)
  { DBLr(finp,&(ydc->dcrmpt));
  }
  else if(CHRcmp(namep,"DCRMCT",6)==0)
  { DBLr(finp,&(ydc->dcrmct));
  }
  else if(CHRcmp(namep,"DCRMUT",6)==0)
  { DBLr(finp,&(ydc->dcrmut));
  }
  else if(CHRcmp(namep,"DCRMJT",6)==0)
  { DBLr(finp,&(ydc->dcrmjt));
  }
  else if(CHRcmp(namep,"DCSTDN",6)==0)
  { DBLr(finp,&(ydc->dcstdn));
  }
  else if(CHRcmp(namep,"DCURELX",7)==0) 
  { 
    DBLr(finp,&(ydc->dcurelx));
  }
  else if(CHRcmp(namep,"DCFACT1",7)==0)
  {
    DBLr(finp,&(ydc->dcfact1));
  }
  else if(CHRcmp(namep,"DCFACT2",7)==0)
  {
    DBLr(finp,&(ydc->dcfact2));
  }
  else if(CHRcmp(namep,"DCFACT3",7)==0)
  {
    DBLr(finp,&(ydc->dcfact3));
  }
  else if(CHRcmp(namep,"DCFACT4",7)==0)
  {
    DBLr(finp,&(ydc->dcfact4));
  }
  else if(CHRcmp(namep,"DCFACT5",7)==0)
  {
    DBLr(finp,&(ydc->dcfact5));
  }

  else if(CHRcmp(namep,"INITER",6)==0) 
  { 
    INTr(finp,&(ydc->initer));
  }
  else if(CHRcmp(namep,"ICOUTF",6)==0) 
  { 
    INTr(finp,&(ydc->icoutf));
    ydc->icoutf2 = ydc->icoutf;
  }
  else if(CHRcmp(namep,"ICOUTI",6)==0) 
  { 
    INTr(finp,&(ydc->icouti));
  }
  else if(CHRcmp(namep,"IWFAST",6)==0)
  { 
    INTr(finp,&(ydc->iwfast));
  }
  else if(CHRcmp(namep,"ISAVE",5)==0) 
  { 
   INTr(finp,&(ydc->isave));
  }
  else if(CHRcmp(namep,"ICHPRM",6)==0)
  {
    INTr(finp,&(ydc->ichprm));
  }
  else if(CHRcmp(namep,"DCIMVELX",8)==0)
  {
    DBLr(finp,&(ydc->dcimvelx));
  }
  else if(CHRcmp(namep,"DCIMVELY",8)==0)
  {
    DBLr(finp,&(ydc->dcimvely));
  }
  else if(CHRcmp(namep,"DCIMVELZ",8)==0)
  {
    DBLr(finp,&(ydc->dcimvelz));
  }

  else
  { 
    CHRw(stderr,"Yrdc: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(103);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* default values 
*/
static void Yrdd(YD yd)
{ 
  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YDB ydb=&(yd->ydb);
  YPAR ypar=&(yd->ypar);
  YSP ysp=&(yd->ysp);
  YDX ydx=&(yd->ydx);
  YDJ ydj=&(yd->ydj);
  YDJP ydjp=&(yd->ydjp);
  /* Set Control to default   */

  ydc->mcstep=0; ydc->ncstep=0;
  ydc->finp=(FILE*)NULL;   /*Z FILENULL */
  ydc->fcheck=(FILE*)NULL;
  ydc->dcgrax=R0;
  ydc->dcgray=R0;
  ydc->dcgraz=R0;
  ydc->dcstec=R1;
  ydc->dctime=R0;
  ydc->dcurelx=R1;
  ydc->initer=2;
  ydc->icoutf=0;
  ydc->icoutf2=0;
  ydc->icouti=0;
  ydc->iwfast=1;
  ydc->isave=1000;
  ydc->dcfact1=1.0;
  ydc->dcfact2=1.0;
  ydc->dcfact3=1.0;
  ydc->dcfact4=1.0;
  ydc->dcfact5=1.0;
  ydc->dcrmpt=R0;
  ydc->dcrmjt=R0;
  ydc->dcrmct=R0;
  ydc->dcrmut=R0;
  ydc->dcgrst=R0;
  ydc->dcstdn=R0;
  ydc->ifrag=1;
    ydc->ichprm=0;
//  ydc->dcmart=R1;
  ydc->r=R1;

  ydc->icfrm=1;
  ydc->ncstepfrm=20000;
  ydc->dcrmgrax=0;
  ydc->dcrmgray=0;
  ydc->dcrmgraz=98.10*2.0;
  ydc->dcfrmca=R0;
  ydc->dcrmvx=R0;
  ydc->dcrmvy=R0;
  ydc->dcrmvz=0.01;
  ydc->dcimvelx=R0;
  ydc->dcimvely=R0;
  ydc->dcimvelz=R0;
  /* Set Elements to default    */
  yde->melem=0; yde->nelem=0;
  yde->melst=0; yde->nelst=0;
  yde->melno=0; yde->nelno=0;
  yde->i1elcf=INT1NULL; 
  yde->i1elpr=INT1NULL;
  yde->d2elst=DBL2NULL;
  yde->i2elto=INT2NULL;
//  yde->d3tcs=DBL3NULL;
  yde->d3tst=DBL3NULL;
  yde->d1emct=DBL1NULL;
  yde->i1elbe=INT1NULL;
//  yde->i2elsur=INT2NULL;
//  yde->i1nelsur=INT1NULL;
//  yde->i1stat=INT1NULL;
    yde->nnelm=0;
    yde->i1elem=INT1NULL;
    yde->i2elfr=INT2NULL;

//  yde->i1eljo=INT1NULL;
  yde->i2eljp=INT2NULL;
//  yde->i2elbnp=INT2NULL;

//  yde->i1elty=INT1NULL;
  yde->i1elcft=INT1NULL;
    yde->d2elfc=DBL2NULL;
 //   yde->d2elfsx=DBL2NULL;
 //   yde->d2elfsy=DBL2NULL;
 //   yde->d2elfsz=DBL2NULL;
    yde->i1joint=INT1NULL;
    yde->i1elfrg=INT1NULL;
    yde->nfrag=0;
    yde->d1volfrg=DBL1NULL;
    yde->d1elvol=DBL1NULL;
    yde->layer_max=1;
    yde->d1mart=DBL1NULL;
    yde->i1elchprm=INT1NULL;
    yde->i1eljt=INT1NULL;
    yde->i1elncstep=INT1NULL;
    yde->d1pem0n=DBL1NULL;
    yde->d1eldam=DBL1NULL;
    yde->i1frgrm=INT1NULL;

    yde->i1elfrm=INT1NULL;
  /* Set Interaction to default */
  ydi->micoup=0; ydi->nicoup=0;
    ydi->icheck=1;
  ydi->iiecff=-2;
  ydi->diedi=BEPSILON;
  ydi->diezon=R0;
  ydi->d1iesl=DBL1NULL;
  ydi->i1iecn=INT1NULL;
  ydi->i1stat=INT1NULL;
  ydi->i1iect=INT1NULL;
  ydi->d1deltat1=DBL1NULL;
  ydi->d1deltat2=DBL1NULL;
  ydi->d1deltan=DBL1NULL;
  ydi->d2nv=DBL2NULL;
  ydi->d2t1v=DBL2NULL;
  ydi->d2t2v=DBL2NULL;
  ydi->drect=R0;
  ydi->i1fcstep=INT1NULL;
    ydi->i1fcjp=INT1NULL;
    ydi->d1pepe=DBL1NULL;

  /* Set Nodes to default  */
  ydn->mnodim=0;  ydn->nnodim=0;
  ydn->mnopo=0;  ydn->nnopo=0;
  ydn->nnopi=0;
  ydn->d1nmct=DBL1NULL;
  //ydn->d1nmcf=DBL1NULL; //Z?
  ydn->d2ncc=DBL2NULL;
  ydn->d2nci=DBL2NULL;
  ydn->d2nck=DBL2NULL;
  ydn->d2nfc=DBL2NULL;
  ydn->d2nft=DBL2NULL;
  ydn->d2nfp=DBL2NULL;
  ydn->i1nold=INT1NULL;

  ydn->d2nvc=DBL2NULL;
  ydn->i1nobf=INT1NULL;
  ydn->i1nopr=INT1NULL;
  ydn->i1noel=INT1NULL;
  ydn->d2ndf=DBL2NULL;
  ydn->i1nind=INT1NULL;
  ydn->d1npre=DBL1NULL;

//  ydn->nneigh=0;
//  ydn->i2nnei=INT2NULL;
//  ydn->i1nei=INT1NULL;
  ydn->i1ntoC2D=INT1NULL;
    ydn->i1nochprm=INT1NULL;

ydn->d2nfrm=DBL2NULL;
ydn->d1ndisp=DBL1NULL;
ydn->i1rmgrav=INT1NULL;
  /* Set Output to default  */
  ydo->mohys=0;  ydo->nohys=0;
  ydo->dohyp=0.05;  /* 5% accuracy */
  ydo->d1ohyc=DBL1NULL;
  ydo->d1ohyf=DBL1NULL;
  ydo->d1ohys=DBL1NULL;
  ydo->d1ohyt=DBL1NULL;
  ydo->d1ohyx=DBL1NULL;
  ydo->d1ohyy=DBL1NULL;
  ydo->i1ohyt=INT1NULL;

  /* Set Properties to default  */
  ydp->mprop=0; ydp->nprop=0;
  ydp->d1peks=DBL1NULL;
  ydp->d1pela=DBL1NULL;
  ydp->d1pemu=DBL1NULL;
  ydp->d1pepe=DBL1NULL;
  ydp->d1pero=DBL1NULL;
  ydp->d1pefr=DBL1NULL;
  ydp->i1ptyp=INT1NULL;
  ydp->d1pesf=DBL1NULL;
  ydp->d1pepsf=DBL1NULL;
  ydp->d1pevf=DBL1NULL;
  ydp->d1pepf=DBL1NULL;
  ydp->d1pefs=DBL1NULL;
  ydp->d1peft=DBL1NULL;
  ydp->d1pegfn=DBL1NULL;
  ydp->d1pegfs=DBL1NULL;

  ydp->i1pejp=INT1NULL;
  ydp->i1pemn=INT1NULL;
  ydp->i1psde=INT1NULL;
 ydp->i1index=INT1NULL;
//-PY changed it for 3D_fracture_coupling_with_multiphase
  ydp->d1picf=DBL1NULL;
  ydp->d1pcoh=DBL1NULL;
  ydp->d1pemr=DBL1NULL;  

//  ydp->d1pem0n=DBL1NULL;
  ydp->d1pem0=DBL1NULL;
  ydp->d1pemf=DBL1NULL;
  ydp->d1pedm=DBL1NULL;

    /* Set Elements to default    */
    ydj->melem=0; ydj->nelem=0;
    ydj->melno=7; ydj->nelno=6;
    ydj->i1elpr=INT1NULL;
    ydj->i2elto=INT2NULL;
    ydj->i1elty=INT1NULL;
    ydj->i2ejjp=INT2NULL;
    ydj->i1eljo=INT1NULL;
    ydj->i1elst=INT1NULL;
    ydj->d2ejfc=DBL2NULL;
//    ydj->d2ejfjx=DBL2NULL;
//    ydj->d2ejfjy=DBL2NULL;
//    ydj->d2ejfjz=DBL2NULL;
    ydj->i1defect=INT1NULL;
    ydj->defect=0.0;
    ydj->dran=R1;
    ydj->d1ejmr=DBL1NULL;
    ydj->d2jst=DBL2NULL;

    
    ydjp->mprop=0; ydjp->nprop=0;
    ydjp->i1ngrp=INT1NULL;
    ydjp->d1pepe=DBL1NULL;
    ydjp->d1peft=DBL1NULL;
    ydjp->d1pegfn=DBL1NULL;
    ydjp->d1pegfs=DBL1NULL;
    ydjp->i1pejp=INT1NULL;
    ydjp->d1picf=DBL1NULL;
    ydjp->d1pcoh=DBL1NULL;


  /* Set Boundary Condition to default  */
  ydb->mbcon=0; ydb->nbcon=0;
  ydb->d1bnax=DBL1NULL;
  ydb->d1bnay=DBL1NULL;
  ydb->d1bnaz=DBL1NULL;
  ydb->d1bnfx=DBL1NULL;
  ydb->d1bnfy=DBL1NULL;
  ydb->d1bnfz=DBL1NULL;
  ydb->d1bnvx=DBL1NULL;
  ydb->d1bnvy=DBL1NULL;
  ydb->d1bnvz=DBL1NULL;
  ydb->i1bnvx=INT1NULL;
  ydb->i1bnvy=INT1NULL;
  ydb->i1bnvz=INT1NULL;

  // ydb->d1bnap=DBL1NULL;


//-PY changed it for 3D_fracture_coupling_with_multiphase
  ypar->nelem=0;   /*!< maximum number of properties            */
  ypar->nelemb=0;   /*!< maximum number of properties            */
  ypar->nnode=0;   /*!< actual number of properties             */
  ypar->i2elto=INT2NULL;
  ypar->d2xyz=DBL2NULL; /*!< [mprop] amplitude of acceleration x     */
  ypar->d2gxyz=DBL2NULL;
  ypar->d1cg=DBL1NULL;
  ypar->d1dst=DBL1NULL;
  ypar->radius=R0;
  ypar->max=R0;
  ypar->d1lv=DBL1NULL;
  ypar->nindex=0;
  ypar->i1elbe=INT1NULL;
  ypar->r=5.0;
  ypar->ir=1;
  ypar->d1vc=DBL1NULL;
  //ypar->zmax=R0;
  ypar->d1rv_min=DBL1NULL;
  ypar->d1rv_max=DBL1NULL;
    ypar->depora=R0;
    ypar->depheg=R0;
    ypar->depoff=R0;



  ysp->mspd=0;
  ysp->nspd=0; /*maximum (actual) number of particle*/

  ysp->d1mass=DBL1NULL;
  ysp->d2pp=DBL2NULL;    /*[nnodim][msd] sphere_position_x                     */
  ysp->d2poldp=DBL2NULL;    /*[nnodim][msd] sphere_position_x                     */
  ysp->d2pvt=DBL2NULL;     /*[nnodim][msd] sphere_velocity_x              */
  ysp->d2ptfd=DBL2NULL;    /*[nnodim][msd] sphere_normal force_x         */
//  ysp->d2pft1=DBL2NULL;   /*[nnodim][msd] sphere_tangential force_x1    */
//  ysp->d2pft2=DBL2NULL;   /*[nnodim][msd] sphere_tangential force_x2    */

  ysp->d2ptfc=DBL2NULL;    /*[nnodim][msd] sphere_total_contact force_x */
  ysp->d2ptf=DBL2NULL;    /*[nnodim][msd] sphere_total_force_x */
  ysp->d2pa=DBL2NULL;     /*[nnodim][msd] sphere_rotation angel_x      */
      ysp->d2paold=DBL2NULL;     /*[nnodim][msd] sphere_rotation angel_x      */
  ysp->d2pvr=DBL2NULL;     /*[nnodim][msd] sphere_rotation velocity_x      */
  ysp->d2pm=DBL2NULL;     /*[nnodim][msd] sphere_moment_x */
    ysp->d2pmc=DBL2NULL;     /*[nnodim][msd] sphere_moment_x due to contact*/

    //  DBL    dsengy;
  ysp->d2prine=DBL2NULL;  /*[nnodim] pricinple inertia                         */
  ysp->d2prinn1=DBL2NULL;  /*[nnodim] pricinple axe                         */
  ysp->d2prinn2=DBL2NULL;  /*[nnodim] pricinple axe                         */
  ysp->d2prinn3=DBL2NULL;  /*[nnodim] pricinple axe                         */
  ysp->d2poldvr=DBL2NULL;     /*[nnodim][msd] sphere_rotation velocity_x      */
  ysp->i1sppr=INT1NULL;
//    ysp->i1spmat=INT1NULL;
  ysp->i1con=INT1NULL;
    ysp->d1impv=DBL1NULL;

    ydx->nelem_t=0;
	ydx->melem_t=0;
    ydx->nelno=0;
    ydx->i2elto_t=INT2NULL;
    ydx->d1nap=DBL1NULL;
    ydx->i1jtp=INT1NULL;
    ydx->d1nbp=DBL1NULL;
    ydx->dbkp=R0;
}


/**********************************************************************/
/**********************************************************************/


/* read data elements   */
static void Yrde(YDE yde, FILE *finp, CHR *name, INT ntelem)
{ 
    CHR *namep; 
    INT i,j,m;

  namep=name+8;
  if(CHRcmp(namep,"MELEM",5)==0) 
  {  
    INTr(finp,&(yde->melem));

yde->melem=yde->melem+ntelem;
      yde->i1elem=TalINT1(yde->melem);
      
      yde->d2elfc=TalDBL2(3,yde->melem);
 
      yde->d3tst=TalDBL3(yde->melem,3,3);
      for(i=0;i<yde->melem;i++)  
        {

          for(j=0;j<3;j++)
            {
              for(m=0;m<3;m++)
              yde->d3tst[i][j][m]=R0;
            }
        }

  yde->i1elcft=TalINT1(yde->melem);
  yde->i1joint=TalINT1(yde->melem);
  yde->i1eljt=TalINT1(yde->melem);
 yde->ifrac=1;
yde->njoint_intact=0;

    for(i=0;i<yde->melem;i++)
      {
	yde->i1elcft[i]=-1;
        yde->i1eljt[i]=-1;

      }
  }
  else if(CHRcmp(namep,"NELEM",5)==0) 
  {  
    INTr(finp,&(yde->nelem));
    yde->nelemi=yde->nelem;
  }
   else if(CHRcmp(namep,"LAYMX",5)==0)
  {
    INTr(finp,&(yde->layer_max));
  }
  else if(CHRcmp(namep,"MELST",5)==0) 
  {  
    INTr(finp,&(yde->melst));
  }
  else if(CHRcmp(namep,"NELST",5)==0) 
  {  
    INTr(finp,&(yde->nelst));
  }
  else if(CHRcmp(namep,"MELNO",5)==0) 
  {  
    INTr(finp,&(yde->melno));
  }
  else if(CHRcmp(namep,"NELNO",5)==0) 
  {  
    INTr(finp,&(yde->nelno));
  }
  else if(CHRcmp(namep,"I1ELCF",6)==0) 
  { 
    TformINT1(finp,-1,yde->melem,&(yde->i1elcf));
  }
  else if(CHRcmp(namep,"I1ELBE",6)==0) 
  { 
    TformINT1(finp,-1,yde->melem,&(yde->i1elbe));
  }
  else if(CHRcmp(namep,"I1ELPR",6)==0) 
  { 
    TformINT1(finp,0,yde->melem,&(yde->i1elpr));

  }
  else if(CHRcmp(namep,"D2ELST",6)==0) 
  { 
    TformDBL2(finp,R0,yde->melst,yde->melem,&(yde->d2elst));
  }
  else if(CHRcmp(namep,"I2ELTO",6)==0) 
  {
    // transformed array order
	  TformINT2_inv(finp,-1,yde->melem,yde->melno,&(yde->i2elto));

          yde->i2eljp=TalINT2(yde->melem,yde->melno);
      yde->i1elfrg=TalINT1(yde->melem);
      yde->d1volfrg=TalDBL1(yde->melem);
      yde->d1elvol=TalDBL1(yde->melem);
      yde->d1mart=TalDBL1(yde->melem);
      yde->i1elchprm=TalINT1(yde->melem);
      yde->i1elncstep=TalINT1(yde->melem);
      yde->i1elfrm=TalINT1(yde->melem);
      yde->i1frgrm=TalINT1(yde->melem);

          for(i=0;i<yde->melem;i++)
            {
                yde->d1volfrg[i]=R0;
                yde->d1elvol[i]=R0;
                yde->i1elfrg[i]=-1;
                yde->i1elchprm[i]=-1;
		yde->d1mart[i]=R1;
                yde->i1elncstep[i]=-1;
		yde->i1elfrm[i]=0;
		yde->i1frgrm[i]=0;
              for(j=0;j<yde->melno;j++)
                {
                  yde->i2eljp[i][j]=-1;
                }
            }

  }
//-PY changed it for 3D_fracture_coupling_with_multiphase
    else if(CHRcmp(namep,"I2ELBNP",6)==0) 
  {

  }

  else if(CHRcmp(namep,"D1EMCT",6)==0) 
  { 
    TformDBL1(finp,R0,yde->melem,&(yde->d1emct));

  }
  else
  { 
    CHRw(stderr,"Yrde: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
  } 
}

/**********************************************************************/
/**********************************************************************/


/* read data joint elements
*/
static void Yrdj(YDJ ydj, FILE *finp, CHR *name)
{
    CHR *namep;
    INT i,j;

  namep=name+8;
  if(CHRcmp(namep,"MELEM",5)==0)
  {
    INTr(finp,&(ydj->melem));

      ydj->i1elty=TalINT1(ydj->melem);
      ydj->i2ejjp=TalINT2(ydj->melem,2);
      ydj->i1elpr=TalINT1(ydj->melem);
      ydj->i1eljo=TalINT1(ydj->melem);
      ydj->i1elst=TalINT1(ydj->melem);
      ydj->d2ejfc=TalDBL2(4,ydj->melem);
/*      ydj->d2ejfjx=TalDBL2(ydj->melem,7);
      ydj->d2ejfjy=TalDBL2(ydj->melem,7);
      ydj->d2ejfjz=TalDBL2(ydj->melem,7);
      */
      ydj->d1ejmr=TalDBL1(ydj->melem);
      ydj->i1defect=TalINT1(ydj->melem);
      ydj->d2jst=TalDBL2(ydj->melem,2);
      
    for(i=0;i<ydj->melem;i++)
      {
	/*      for(j=0;j<6;j++)
	      {
		      ydj->d2ejfjx[i][j]=R0;
		      ydj->d2ejfjy[i][j]=R0;
		      ydj->d2ejfjz[i][j]=R0;
	      }
	      */
          ydj->i1elty[i]=-1;
          ydj->i1elpr[i]=-1;
          ydj->i1eljo[i]=-1;
          ydj->i2ejjp[i][0]=-1;
          ydj->i2ejjp[i][1]=-1;
          ydj->i1elst[i]=-1;
          ydj->d2ejfc[0][i]=R0;
          ydj->d2ejfc[1][i]=R0;
          ydj->d2ejfc[2][i]=R0;
          ydj->d1ejmr[i]=R1;
          ydj->i1defect[i]=-1;
          ydj->d2jst[i][0]=R0;
          ydj->d2jst[i][1]=R0;

//	  INTw(stdout,ydj->i1elst[i],5);

      }
  }

  else if(CHRcmp(namep,"MELNO",5)==0)
  {
    INTr(finp,&(ydj->melno));
  }
  else if(CHRcmp(namep,"NELNO",5)==0)
  {
    INTr(finp,&(ydj->nelno));
    ydj->i2elto=TalINT2(ydj->melem,ydj->nelno+1);

  }
  else if(CHRcmp(namep,"DEFECT",6)==0)
  {
    DBLr(finp,&(ydj->defect));
  }
  else if(CHRcmp(namep,"DRAN",6)==0)
  {
    DBLr(finp,&(ydj->dran));
  }

  else
  {
    CHRw(stderr,"Yrdj: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
  }
}
/**********************************************************************/
/**********************************************************************/


/* read data interaction 
*/
static void Yrdi(YDI ydi, FILE *finp, CHR *name)
{ 
    CHR *namep; 
    INT i,j;

  namep=name+8;
  if(CHRcmp(namep,"MICOUP",6)==0) 
  {  
	  INTr(finp,&(ydi->micoup));
/*
	  ydi->i1fcstep=TalINT1(ydi->micoup);
	  for(i=0;i<ydi->micoup;i++)
	    {
	      ydi->i1fcstep[i]=-1;
	    }
	    */
  }
  else if(CHRcmp(namep,"NICOUP",6)==0) 
  {  
	  INTr(finp,&(ydi->nicoup));
  }
  else if(CHRcmp(namep,"IIECFF",6)==0) 
  {  
	  INTr(finp,&(ydi->iiecff));
  }
  else if(CHRcmp(namep,"DIEDI",5)==0) 
  {  
	  DBLr(finp,&(ydi->diedi));
  }
  else if(CHRcmp(namep,"DIEZON",6)==0) 
  {  
	  DBLr(finp,&(ydi->diezon));
  }
  else if(CHRcmp(namep,"D1IESL",6)==0) 
  { 
	  TformDBL1(finp,R0,ydi->micoup,&(ydi->d1iesl)); 
  }
  else if(CHRcmp(namep,"DRECT",6)==0)
  {
	  DBLr(finp,&(ydi->drect));
  }
  else if(CHRcmp(namep,"I1IECN",6)==0)
  { 
	  TformINT1(finp,-1,ydi->micoup,&(ydi->i1iecn)); 
  }
  else if(CHRcmp(namep,"I1STAT",6)==0)
  {
	  TformINT1(finp,-1,ydi->micoup,&(ydi->i1stat));
  }
  else if(CHRcmp(namep,"I1IECT",6)==0)
  { 
	  TformINT1(finp,-1,ydi->micoup,&(ydi->i1iect)); 
  }
  else if(CHRcmp(namep,"D1DELTAT1",9)==0) 
  { 
	  TformDBL1(finp,R0,ydi->micoup,&(ydi->d1deltat1)); 
  }
  else if(CHRcmp(namep,"D1DELTAT2",9)==0) 
  { 
	  TformDBL1(finp,R0,ydi->micoup,&(ydi->d1deltat2)); 
  }
  else if(CHRcmp(namep,"D1DELTAN",8)==0) 
  { 
	  TformDBL1(finp,R0,ydi->micoup,&(ydi->d1deltan));
      ydi->i1fcstep=TalINT1(ydi->micoup);
      ydi->i1fcjp=TalINT1(ydi->micoup);
      ydi->d1pepe=TalDBL1(ydi->micoup);
      for(i=0;i<ydi->micoup;i++)
      {
          ydi->i1fcstep[i]=-1;
          ydi->i1fcjp[i]=-1;
          ydi->d1pepe[i]=R0;
      }

/*
      ydi->d2nv=TalDBL2(4,ydi->micoup);
      ydi->d2t1v=TalDBL2(4,ydi->micoup);
      ydi->d2t2v=TalDBL2(4,ydi->micoup);
      for(i=0;i<ydi->micoup;i++)
      {
          for(j=0;j<3;j++)
          {
          ydi->d2nv[j][i]=R0;
          ydi->d2t1v[j][i]=R0;
          ydi->d2t2v[j][i]=R0;
          }
      } */
  }
  else if(CHRcmp(namep,"D2NV",4)==0) 
  { 
	  TformDBL2(finp,R0,3,ydi->micoup,&(ydi->d2nv));
  }
  else if(CHRcmp(namep,"D2T1V",5)==0) 
  { 
	  TformDBL2(finp,R0,3,ydi->micoup,&(ydi->d2t1v));
  }
  else if(CHRcmp(namep,"D2T2V",5)==0) 
  { 
	  TformDBL2(finp,R0,3,ydi->micoup,&(ydi->d2t2v));
  }

  else
  { 
    CHRw(stderr,"Yrdi: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* read data nodes       
*/
static void Yrdn(YDN ydn, FILE *finp, CHR *name, INT ntnode, DBL dcimvelx,DBL dcimvely,DBL dcimvelz)
{ 
  CHR *namep;
    INT i,j;

  namep=name+8;
  if(CHRcmp(namep,"MNODIM",5)==0) 
  {  
    INTr(finp,&(ydn->mnodim));
  }
  else if(CHRcmp(namep,"NNODIM",5)==0) 
  {  
    INTr(finp,&(ydn->nnodim));
  }
  else if(CHRcmp(namep,"MNOPO",5)==0) 
  {  
    INTr(finp,&(ydn->mnopo));
    ydn->mnopo=ydn->mnopo+ntnode;

/*    ydn->i2nnei=TalINT2(ydn->mnopo,100);
    for(i=0;i<ydn->mnopo;i++)
      {
        for(j=0;j<100;j++)
          {
            ydn->i2nnei[i][j]=-1;
          }
      }
    ydn->i1nei=TalINT1(ydn->mnopo);
    for(i=0;i<ydn->mnopo;i++)
      {
	ydn->i1nei[i]=0;
      }
*/
//-PY changed it for 3D_fracture_coupling_with_multiphase
    ydn->mnopo=ydn->mnopo+ntnode; //-ao there could be a memory issue with this 
	ydn->i1nind=TalINT1(ydn->mnopo);

	ydn->d1nvct=TalDBL1(ydn->mnopo);
  }
  else if(CHRcmp(namep,"NNOPO",5)==0) 
  {  
    INTr(finp,&(ydn->nnopo));
ydn->nnopi=ydn->nnopo;
  }
    /*
  else if(CHRcmp(namep,"D1NMCT",6)==0) 
  { 
    TformDBL1(finp,R0,ydn->mnopo,&(ydn->d1nmct));
  }*/
  else if(CHRcmp(namep,"D2NCC",5)==0) 
  { 
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2ncc)); 
  }
  else if(CHRcmp(namep,"D2NCK",5)==0)
  {
      TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nck));
  }
  else if(CHRcmp(namep,"D2NCI",5)==0)
  { 
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nci));
      ydn->d2nfrm=TalDBL2(ydn->mnodim,ydn->mnopo);
      ydn->d2nfc=TalDBL2(ydn->mnodim,ydn->mnopo);
      ydn->d2nft=TalDBL2(ydn->mnodim,ydn->mnopo);
      ydn->d2nfp=TalDBL2(ydn->mnodim,ydn->mnopo);

      ydn->d1nmct=TalDBL1(ydn->mnopo);
      ydn->i1noel=TalINT1(ydn->mnopo);
      ydn->i1nold=TalINT1(ydn->mnopo);
      ydn->i1nochprm=TalINT1(ydn->mnopo);

      ydn->d1npre=TalDBL1(ydn->mnopo);
      ydn->d1ndisp=TalDBL1(ydn->mnopo);
      ydn->i1rmgrav=TalINT1(ydn->mnopo);
      
      for(i=0;i<ydn->mnopo;i++)
      {
          ydn->d1nmct[i]=R0;
          ydn->i1noel[i]=-1;
          ydn->i1nold[i]=i;
          ydn->i1nochprm[i]=-1;
          ydn->d1npre[i]=R0;
	  ydn->d1ndisp[i]=R0;
	  ydn->i1rmgrav[i]=0;
      for(j=0;j<3;j++)
      {
          ydn->d2nfrm[j][i]=R0;
          ydn->d2nfc[j][i]=R0;
          ydn->d2nft[j][i]=R0;
          ydn->d2nfp[j][i]=R0;

          
      }
      }

  }
  else if(CHRcmp(namep,"D2NDF",5)==0)
  {
      TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2ndf));
  }

    /*
  else if(CHRcmp(namep,"D2NCK",5)==0)
  {
      TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nck));
  }
  else if(CHRcmp(namep,"D2NFC",5)==0)
  {
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nfc)); 
  }  
  else if(CHRcmp(namep,"D2NFT",5)==0) 
  { 
//    TformDBL2_inv(finp,R0,ydn->mnopo,ydn->mnodim,&(ydn->d2nft)); 
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nft)); 
  }
     */
    
  else if(CHRcmp(namep,"D2NVC",5)==0) 
  { 
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nvc)); 
/*      for(i=0;i<ydn->nnopo;i++)
      {
	      if((ydn->d2nvc[0][i]>1.0e-05)||(ydn->d2nvc[0][i]<-1.0e-05))ydn->d2nvc[0][i]=dcimvelx;
	      if((ydn->d2nvc[1][i]>1.0e-05)||(ydn->d2nvc[1][i]<-1.0e-05))ydn->d2nvc[1][i]=dcimvely;
	      if((ydn->d2nvc[2][i]>1.0e-05)||(ydn->d2nvc[2][i]<-1.0e-05))ydn->d2nvc[2][i]=dcimvelz;
      }
      */
  }
  else if(CHRcmp(namep,"I1NOBF",6)==0) 
  { 
    TformINT1(finp,0,ydn->mnopo,&(ydn->i1nobf));
  }
  else if(CHRcmp(namep,"I1NOPR",6)==0) 
  { 
    TformINT1(finp,0,ydn->mnopo,&(ydn->i1nopr));
  }
  else
  { 
    CHRw(stderr,"Yrdn: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* read output specification     
*/
static void Yrdo(YDO ydo, FILE *finp, CHR *name)
{ 
  INT i; CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MOHYS",5)==0) 
  {  
    INTr(finp,&(ydo->mohys))
    i=(ydo->mohys)*sizeof(FILE*);
    if(i>0)ydo->f2ohyf=(FILE**)MALLOC(i);
    for(i=0;i<(ydo->mohys);i++)
    { 
      ydo->f2ohyf[i]=FILENULL;
    }
  }
  else if(CHRcmp(namep,"NOHYS",5)==0) 
  {  
    
    INTr(finp,&(ydo->nohys));
  }
  else if(CHRcmp(namep,"DOHYP",5)==0) 
  { 
    DBLr(finp,&(ydo->dohyp));
  }
  else if(CHRcmp(namep,"D1OHYS",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohys));
  }
  else if(CHRcmp(namep,"D1OHYC",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyc));
  }
  else if(CHRcmp(namep,"D1OHYF",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyf));
  }
  else if(CHRcmp(namep,"D1OHYT",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyt));
  }
  else if(CHRcmp(namep,"D1OHYX",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyx));
  }
  else if(CHRcmp(namep,"D1OHYY",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyy));
  }
  else if(CHRcmp(namep,"D1OHYZ",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyz));
  }
  else if(CHRcmp(namep,"I1OHYT",6)==0) 
  { 
    TformINT1(finp,-1,ydo->mohys,&(ydo->i1ohyt));
  }
  else
  { 
    CHRw(stderr,"Yrdo: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* read data properties  
*/
static void Yrdp(YDP ydp, FILE *finp, CHR *name, YDE yde)
{ 
  CHR *namep;
  INT i;
  namep=name+8;
  if(CHRcmp(namep,"MPROP",5)==0) 
  {  
    INTr(finp,&(ydp->mprop));
      ydp->mprop++;
  }
  else if(CHRcmp(namep,"NPROP",5)==0) 
  {  
    INTr(finp,&(ydp->nprop));
  }
  else if(CHRcmp(namep,"I1INDEX",7)==0)
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1index));
  }
  else if(CHRcmp(namep,"D1PEKS",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1peks)); 
  }
  else if(CHRcmp(namep,"D1PELA",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pela)); 
  }
  else if(CHRcmp(namep,"D1PEMU",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pemu)); 
  }
  else if(CHRcmp(namep,"D1PEPE",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepe)); 
  }
  else if(CHRcmp(namep,"D1PERO",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pero)); 
  }
  else if(CHRcmp(namep,"D1PEFR",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pefr)); 
  }
  else if(CHRcmp(namep,"D1PESF",6)==0) 
  { 
	TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pesf)); 
  }
  else if(CHRcmp(namep,"D1PEPSF",7)==0) 
  { 
	TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepsf)); 
  }
  else if(CHRcmp(namep,"D1PEVF",6)==0) 
  { 
	TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pevf)); 
  }
  else if(CHRcmp(namep,"D1PEPF",6)==0) 
  { 
	TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepf)); 
  }
  else if(CHRcmp(namep,"D1PEFS",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pefs));
  }
  else if(CHRcmp(namep,"D1PEFT",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1peft));
		ydp->d1pefs=TalDBL1(ydp->mprop);
  }
  else if(CHRcmp(namep,"D1PEGFN",7)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pegfn));
  }
  else if(CHRcmp(namep,"D1PEGFS",7)==0)
  {
      TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pegfs));
  }
  else if(CHRcmp(namep,"D1PCOH",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pcoh));
  }
  else if(CHRcmp(namep,"D1PICF",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1picf));
  }
  else if(CHRcmp(namep,"D1PEMR",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pemr));
  }
  else if(CHRcmp(namep,"I1PTYP",6)==0) 
  { 
        TformINT1(finp,-1,ydp->mprop,&(ydp->i1ptyp));
  }
  else if(CHRcmp(namep,"I1PEJP",6)==0)
  {
        TformINT1(finp,0,ydp->mprop,&(ydp->i1pejp));
  }
  else if(CHRcmp(namep,"I1PEMN",6)==0)
  {
        TformINT1(finp,0,ydp->mprop,&(ydp->i1pemn));
  }
  else if(CHRcmp(namep,"I1PSDE",6)==0)
  {
        TformINT1(finp,0,ydp->mprop,&(ydp->i1psde));
  }
  else if(CHRcmp(namep,"D1PEM0",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pem0));
        yde->d1pem0n=TalDBL1(yde->melem);
	yde->d1eldam=TalDBL1(yde->melem);
        for(i=0;i<(yde->nelem);i++)
	{
		yde->d1pem0n[i]=ydp->d1pem0[yde->i1elpr[i]];
		yde->d1eldam[i]=R0;
//		printf("d1pe0n %f\n",yde->d1pem0n[i]);
	}
  }
  else if(CHRcmp(namep,"D1PEMF",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pemf));
  }
  else if(CHRcmp(namep,"D1PEDM",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pedm));
   }
  else
  { 
    CHRw(stderr,"Yrdp: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* read data properties  
*/
static void Yrdb(YDB ydb, FILE *finp, CHR *name)
{ 
  CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MBCON",5)==0) 
  {  
    INTr(finp,&(ydb->mbcon));
	ydb->mbcon+=1;  //packing only in case initial input (.Y3D) has only one condition
  }
  else if(CHRcmp(namep,"NBCON",5)==0) 
  {  
    INTr(finp,&(ydb->nbcon));
  }
  else if(CHRcmp(namep,"D1BNAX",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnax)); 
  }
  else if(CHRcmp(namep,"D1BNAY",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnay)); 
  }
  else if(CHRcmp(namep,"D1BNAZ",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnaz)); 
  }
  else if(CHRcmp(namep,"D1BNFX",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnfx));
  }
  else if(CHRcmp(namep,"D1BNFY",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnfy));
  }
  else if(CHRcmp(namep,"D1BNFZ",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnfz));
  }
  else if(CHRcmp(namep,"I1BNVX",6)==0) 
  { 
    TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnvx));
  }
  else if(CHRcmp(namep,"I1BNVY",6)==0) 
  { 
    TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnvy));
  }
  else if(CHRcmp(namep,"I1BNVZ",6)==0) 
  { 
    TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnvz));
  }
  else if(CHRcmp(namep,"D1BNVX",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnvx));
  }
  else if(CHRcmp(namep,"D1BNVY",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnvy));
  }
  else if(CHRcmp(namep,"D1BNVZ",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnvz));
  }
  else
  { 
    CHRw(stderr,"Yrdp: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
  } 
}
/* read data nodes       
 */
static void Yrsp(YSP ysp, FILE *finp, CHR *name, INT ntpar)
{ 
	CHR *namep;
    INT nspd_tmp,i;
    INT mnodim;
	mnodim=4;
	namep=name+8;
	if(CHRcmp(namep,"MSPD",4)==0) 
	{  
		INTr(finp,&(ysp->mspd));
		ysp->mspd=ysp->mspd+ntpar+1;
	}
	else if(CHRcmp(namep,"NSPD",4)==0) 
	{  
		INTr(finp,&(ysp->nspd));
	}
        else if(CHRcmp(namep,"I1CON",5)==0)
        {
			/*
                INTr(finp,&nspd_tmp);
                for(i=0;i<nspd_tmp;i++)
                {
                INTr(finp,&(ysp->i1con));
                }*/
        TformINT1(finp,-1,ysp->mspd,&(ysp->i1con));
        }
	else if(CHRcmp(namep,"D1MASS",6)==0) 
	{ 
		TformDBL1(finp,R0,ysp->mspd,&(ysp->d1mass));
	}
	else if(CHRcmp(namep,"D2PP",4)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2pp)); 
	}
	/*
	else if(CHRcmp(namep,"D2   DBLr(fgrid,&d1tmp[0]);POLDP",7)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2poldp)); 
	}*/
	else if(CHRcmp(namep,"D2PVT",5)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2pvt)); 
	}
	else if(CHRcmp(namep,"D2PA",4)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2pa));
        ysp->d2paold=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation angel_x      */
        
        for(i=0;i<ysp->mspd;i++)
        {
            
                ysp->d2paold[0][i]=ysp->d2pa[0][i];
            ysp->d2paold[1][i]=ysp->d2pa[1][i];
            ysp->d2paold[2][i]=ysp->d2pa[2][i];
        }

	}
	else if(CHRcmp(namep,"D2PVR",5)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2pvr)); 
	}

	else if(CHRcmp(namep,"D2PRINE",7)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2prine)); 
	}
	else if(CHRcmp(namep,"D2PRINN1",8)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2prinn1)); 
	}
	else if(CHRcmp(namep,"D2PRINN2",8)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2prinn2)); 
	}
	else if(CHRcmp(namep,"D2PRINN3",8)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2prinn3));
        
        ysp->d1impv=TalDBL1(3);
                ysp->d1impv[0]=R0;
                ysp->d1impv[1]=R0;
                ysp->d1impv[2]=R0;
	}
	else
	{ 
		CHRw(stderr,"Yrsp: unknown name: ");
		CHRw(stderr,name); 
		CHRwcr(stderr);
		exit(105);
	} 
}


/**********************************************************************/
/**********************************************************************/


/* read data elements
 */
static void Yrdx(YDX ydx, YDN ydn, FILE *finp, CHR *name, INT ntelem, INT nelem_tet)
{
    CHR *namep; INT i,j;
    
    namep=name+8;
	if(CHRcmp(namep,"MELEM",5)==0)
    {
        INTr(finp,&(ydx->melem_t));
        /*considering fracture model, maximum number surface mesh is $*nelem_tet*/
        ydx->melem_t=nelem_tet*4;
        ydx->melem_t=ydx->melem_t+ntelem+1;
    }
    else if(CHRcmp(namep,"NELEM",5)==0)
    {
        INTr(finp,&(ydx->nelem_t));
        //ydx->nelem_t=ydx->nelem_t+ntelem;
    }
    else if(CHRcmp(namep,"NELNO",5)==0)
    {
        INTr(finp,&(ydx->nelno));
    }
    else if(CHRcmp(namep,"I2ELTO",6)==0)
    {
        //Z transformed array order
        //Z TformINT2(finp,-1,yde->melno,yde->melem,&(yde->i2elto));
        TformINT2_inv(finp,-1,ydx->melem_t,ydx->nelno+1,&(ydx->i2elto_t));

    }
    else if(CHRcmp(namep,"D1NAP",5)==0)
    {
        TformDBL1(finp,R0,ydx->melem_t,&(ydx->d1nap));
/*        ydx->i1bd=TalINT1(ydx->melem_t);
        for(i=0;i<ydx->melem_t;i++)
        {
            ydx->i1bd[i]=-1;
            if(ydx->d1nap[i]>EPSILON)ydx->i1bd[i]=1;
        }
	*/
        
    }
    else if(CHRcmp(namep,"D1NBP",5)==0)
    {
        TformDBL1(finp,R0,ydx->melem_t,&(ydx->d1nbp));
/*	for(i=0;i<ydx->melem_t;i++)
        
	{
		ydx->d1nbp[i]=20.0e+06;  //manually set back pressure in the code
            if(ydx->d1nbp[i]>EPSILON)ydx->dbkp=ydx->d1nbp[i];
        }*/
    }
    else if(CHRcmp(namep,"I1JTP",5)==0)
    {
        TformINT1(finp,0,ydx->melem_t,&(ydx->i1jtp));
/*
	for(i=0;i<ydx->melem_t;i++)

        {
                if(ydn->d2ncc[2][ydx->i2elto_t[i][0]]>-0.003)ydx->i1jtp[i]=1;  //manually set back pressure in the code
        }
*/
    }
    else
    { 
        CHRw(stderr,"Yrdx: unknown name: ");
        CHRw(stderr,name); 
        CHRwcr(stderr);
        exit(106);
    } 
}




/* allocate memory for variable of super particle  */
  void Yspini(YSP ysp,INT ncindex,INT ntpar,INT icon)
  {INT i,j,mnodim;
  mnodim=4;
  if (ysp->mspd==0) 
  {

  ysp->mspd=ntpar+ncindex+1;
  ysp->nspd=ncindex; /*maximum (actual) number of particle*/

  ysp->d1mass=TalDBL1(ysp->mspd);
  ysp->d2pp=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_position_x                     */
  ysp->d2poldp=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_position_x                     */
  ysp->d2pvt=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_velocity_x              */
  //ysp->d2pfn=TalDBL2(3,ysp->mspd);    /*[nnodim][msd] sphere_normal force_x         */
  //ysp->d2pft1=TalDBL2(3,ysp->mspd);   /*[nnodim][msd] sphere_tangential force_x1    */
  //ysp->d2pft2=TalDBL2(3,ysp->mspd);   /*[nnodim][msd] sphere_tangential force_x2    */

  ysp->d2ptfc=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_total_contact force_x */
  ysp->d2ptf=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_total_force_x */
	    ysp->d2ptfd=TalDBL2(mnodim,ysp->mspd);  
  ysp->d2pa=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation angel_x      */
        ysp->d2paold=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation angel_x      */
  ysp->d2pvr=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation velocity_x      */
  ysp->d2pm=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_moment_x */
      ysp->d2pmc=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_moment_x due to contact */

  ysp->d2prine=TalDBL2(mnodim,ysp->mspd);  /*[nnodim] pricinple inertia                         */
  ysp->d2prinn1=TalDBL2(mnodim,ysp->mspd);  /*[nnodim] pricinple axe                         */
  ysp->d2prinn2=TalDBL2(mnodim,ysp->mspd);  /*[nnodim] pricinple axe                         */
  ysp->d2prinn3=TalDBL2(mnodim,ysp->mspd);  /*[nnodim] pricinple axe                         */
  ysp->d2poldvr=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation velocity_x    */
  ysp->i1sppr=TalINT1(ysp->mspd);
 //     ysp->i1spmat=TalINT1(ysp->mspd);
  ysp->i1con=TalINT1(ysp->mspd);

  ysp->d1impv=TalDBL1(3);
  ysp->d1impv[0]=R0;
  ysp->d1impv[1]=R0;
  ysp->d1impv[2]=R0;

  for(i=0;i<ysp->mspd;i++)
  {ysp->d1mass[i]=R0;
   ysp->i1sppr[i]=-1;
  //    ysp->i1spmat[i]=-1;
   if(icon==1){ysp->i1con[i]=0;}
   else {ysp->i1con[i]=1;}

	  for(j=0;j<mnodim-1;j++)
	  {ysp->d2pp[j][i]=R0;
	  ysp->d2poldp[j][i]=R0;
	  ysp->d2pvt[j][i]=R0;
	  //ysp->d2pfn
	  //ysp->d2pft1
      //ysp->d2pft2
	  ysp->d2ptfc[j][i]=R0;
	  ysp->d2ptf[j][i]=R0;
	  ysp->d2pa[j][i]=R0;
          	  ysp->d2paold[j][i]=R0;
	  ysp->d2pvr[j][i]=R0;
	  ysp->d2pm[j][i]=R0;
          ysp->d2pmc[j][i]=R0;
	  ysp->d2prine[j][i]=R0;
	  ysp->d2prinn1[j][i]=R0;
	  ysp->d2prinn2[j][i]=R0;
	  ysp->d2prinn3[j][i]=R0;

	  }
  }
  }
	  else
	  {
		  
		  ysp->d2ptfc=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_total_contact force_x */
		  ysp->d2ptf=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_total_force_x */
		  ysp->d2ptfd=TalDBL2(mnodim,ysp->mspd);
		  ysp->d2pm=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_moment_x */
          ysp->d2pmc=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_moment_x due to contact */
          ysp->d2poldp=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_position_x                     */
                  ysp->d1impv=TalDBL1(3);

		  ysp->i1sppr=TalINT1(ysp->mspd);
  //        ysp->i1spmat=TalINT1(ysp->mspd);

          //ysp->i1con=TalINT1(ysp->mspd);
                  ysp->d2poldvr=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation velocity_x    */		  

		  for(i=ysp->nspd;i<ysp->mspd;i++)
		  {
	
			  if(icon==1){ysp->i1con[i]=0;}
			  else {ysp->i1con[i]=1;}
		  }

		  for(i=0;i<ysp->mspd;i++)
		  {
			  ysp->i1sppr[i]=-1;
 //             ysp->i1spmat[i]=-1;

			  if(icon==1){ysp->i1con[i]=0;}
			  else {ysp->i1con[i]=1;}
			  
			  for(j=0;j<mnodim-1;j++)
			  {
				  ysp->d2ptfc[j][i]=R0;
				  ysp->d2ptf[j][i]=R0;
				  ysp->d2pm[j][i]=R0;
                  ysp->d2pmc[j][i]=R0;
                                  ysp->d1impv[0]=R0;
                                  ysp->d1impv[1]=R0;
                                  ysp->d1impv[2]=R0;



				  
			  }
		  }  
	  }



}

void static read_mat_file(YDJP ydjp,CHR *name_mat)
{
    FILE *fmat;
    CHR name[300];
    
    INT ngroup;
    ngroup=0;
    fmat=fopen(name_mat,"r");
    CHRr(fmat,name); 
  //  CHRw(stdout,name);
    while(FILEND(fmat)==0)
    {
      if(CHRcmp(name,"Material_number",15)==0)
      {
          INTr(fmat,&ydjp->nprop);
          ydjp->i1ngrp=TalINT1(ydjp->nprop);
	  ydjp->d1peft=TalDBL1(ydjp->nprop);
          ydjp->d1pepe=TalDBL1(ydjp->nprop);
          ydjp->d1pegfn=TalDBL1(ydjp->nprop);
          ydjp->d1pegfs=TalDBL1(ydjp->nprop);
          ydjp->d1picf=TalDBL1(ydjp->nprop);
          ydjp->d1pcoh=TalDBL1(ydjp->nprop);
      }
        else if(CHRcmp(name,"Material_Group_index",20)==0)
        {
            INTr(fmat,&ydjp->i1ngrp[ngroup]);
            ngroup++;
        }
        else if(CHRcmp(name,"Tensile_Strength",16)==0)
        {
            DBLr(fmat,&ydjp->d1peft[ngroup-1]);
//	    fprintf(stdout,"ngroup is %ld, d1peft is %f",ngroup, ydjp->d1peft[ngroup-1]);
        }
        else if(CHRcmp(name,"Internal_Friction_Angle",23)==0)
        {
            DBLr(fmat,&ydjp->d1picf[ngroup-1]);
        }
        else if(CHRcmp(name,"Cohesion",8)==0)
        {
            DBLr(fmat,&ydjp->d1pcoh[ngroup-1]);
        }
        else if(CHRcmp(name,"GI_Energy_Release_Rate",22)==0)
        {
            DBLr(fmat,&ydjp->d1pegfn[ngroup-1]);
        }
        else if(CHRcmp(name,"GII_Energy_Release_Rate",23)==0)
        {
            DBLr(fmat,&ydjp->d1pegfs[ngroup-1]);
        }
        else if(CHRcmp(name,"Spring_Stiffness",16)==0)
        {
            DBLr(fmat,&ydjp->d1pepe[ngroup-1]);
        }
	else
       {
        CHRw(stderr,"Yrdjp: unknown name: ");
        CHRw(stderr,name);
        CHRwcr(stderr);
        exit(107);
        }

        CHRr(fmat,name);

    }
    
    fclose(fmat);
    
}


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief read input
 *  \param[in]     namep file name of input
 *  \param[in,out] yd Y database
 *  \par Details:
 *  Yrd() ...
 *
 */
INT Yrd( CHR *namep, CHR *name_grid,CHR *name_mesh, CHR *name_mat, YD yd)
{ 
  INT icount;
  CHR name[300];
  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YDB ydb=&(yd->ydb);
  YPAR ypar=&(yd->ypar);
  YSP ysp=&(yd->ysp);
  YDJ ydj=&(yd->ydj);
  YDJP ydjp=&(yd->ydjp);

  YDX ydx=&(yd->ydx);
    
  FILE *fgrid, *fgid;
    INT i,j,index;//nelemb,**i2eltob,
    INT ntyp,icon;
  DBL dis,tmp,**points;
  INT ntelem,ntelemb,ntnode,*i1nobf,k,m,n,index_sp;
  const double mult[10]={1./6.,1./24.,1./24.,1./24.,1./60.,1./60.,1./60.,1./120.,1./120.,1./120.};
  DBL *v1,*v2,*v3;
  DBL **inertia,ro,cm[3],*intg,d1impv[3];
  DBL d1tmp[11];
  INT inewmat;
INT l,kk,ll,jj;
//  nelemb=0;
  if((ydc->finp)==(FILENULL))
  { 
    Yrdd(yd);
    ydc->finp=fopen(namep,"r");
    ydc->fcheck=fopen("Ytmp","w");
  }
  if(((ydc->finp)==(FILENULL))||((ydc->fcheck)==(FILENULL)))
  { 
    CHRw(stderr,"Yrd: Could not open input file - usage -i inputfile"); 
    CHRwcr(stderr);   
    return 0;
  }


/* read grid file to define particle position also include total particle number,
 initial particle number, and minimum height distance */
ypar->d1dst=TalDBL1(4);
ypar->d1lv=TalDBL1(4);
ypar->d1vc=TalDBL1(4);
    ypar->d1rv_max=TalDBL1(4);
    ypar->d1rv_min=TalDBL1(4);

v1=TalDBL1(4);
v2=TalDBL1(4);
v3=TalDBL1(4);
inertia=TalDBL2(4,4);
points=TalDBL2(5,4);
intg=TalDBL1(11);

    read_mat_file(ydjp,name_mat);
    
inewmat=0;
  fgrid=fopen(name_grid,"r");
  fgid=fopen(name_mesh,"r");
  ypar->ifeed=0;
  while(FILEND(fgrid)==0) 
  { 
    if(CHRcmp(name,"total_number",12)==0)    
    {    
		INTr(fgrid,&ypar->ntpar); 
	}
	else if(CHRcmp(name,"particle_number_per_layer",15)==0)    
    {    
		INTr(fgrid,&ypar->nlpar); 
	}
	else if(CHRcmp(name,"grid_number",11)==0)    
    {    
		INTr(fgrid,&ypar->ngrid); 
    }
		else if(CHRcmp(name,"maximum_dimension",11)==0)    
    {    
		DBLr(fgrid,&ypar->d1dst[0]); 
		DBLr(fgrid,&ypar->d1dst[1]); 
		DBLr(fgrid,&ypar->d1dst[2]); 
    }
	else if(CHRcmp(name,"normal_direction",11)==0)    
    {    
		DBLr(fgrid,&ypar->d1lv[0]); 
		DBLr(fgrid,&ypar->d1lv[1]); 
		DBLr(fgrid,&ypar->d1lv[2]); 
    }
   else if(CHRcmp(name,"particle_velocity",11)==0)    
    {    
		DBLr(fgrid,&ypar->d1vc[0]); 
		DBLr(fgrid,&ypar->d1vc[1]); 
		DBLr(fgrid,&ypar->d1vc[2]); 
    }
	   else if(CHRcmp(name,"control_velocity",11)==0)    
    {   icon=0; 
		CHRr(fgrid,name);
		if(CHRcmp(name,"yes",3)==0)
		{
			icon=1;
		}
    }

       else if(CHRcmp(name,"particle_feedin",15)==0)
       {
           CHRr(fgrid,name);
           if(CHRcmp(name,"yes",3)==0)
           {
               ypar->ifeed=1;
           }
       }
       else if(CHRcmp(name,"new_material",12)==0)
       {
           CHRr(fgrid,name);
           if(CHRcmp(name,"yes",3)==0)
           {
               inewmat=1;
           }
       }

  else if(CHRcmp(name,"D1PEKS",6)==0)
  {
   DBLr(fgrid,&d1tmp[0]);
  }
  else if(CHRcmp(name,"D1PELA",6)==0)
  {
   DBLr(fgrid,&d1tmp[1]);
  }
  else if(CHRcmp(name,"D1PEMU",6)==0)
  {
   DBLr(fgrid,&d1tmp[2]);
  }
  else if(CHRcmp(name,"D1PEPE",6)==0)
  {
   DBLr(fgrid,&d1tmp[3]);
  }
  else if(CHRcmp(name,"D1PERO",6)==0)
  {
   DBLr(fgrid,&d1tmp[4]);
  }
  else if(CHRcmp(name,"D1PEFR",6)==0)
  {
   DBLr(fgrid,&d1tmp[5]);
  }
  else if(CHRcmp(name,"D1PESF",6)==0)
  {
   DBLr(fgrid,&d1tmp[6]);
  }
  else if(CHRcmp(name,"D1PEPSF",7)==0)
  {
   DBLr(fgrid,&d1tmp[7]);
  }
  else if(CHRcmp(name,"D1PEVF",6)==0)
  {
   DBLr(fgrid,&d1tmp[8]);
  }
  else if(CHRcmp(name,"D1PEPF",6)==0)
  {
   DBLr(fgrid,&d1tmp[9]);
  }


       else if(CHRcmp(name,"deposition_radius",17)==0)
       {
           DBLr(fgrid,&ypar->depora);
       }
       else if(CHRcmp(name,"deposition_height",17)==0)
       {
           DBLr(fgrid,&ypar->depheg);
       }
       else if(CHRcmp(name,"deposition_offset",17)==0)
       {
           DBLr(fgrid,&ypar->depoff);
       }
   else if(CHRcmp(name,"impact_velocity",11)==0)
    {
                DBLr(fgrid,&d1impv[0]);
                DBLr(fgrid,&d1impv[1]);
                DBLr(fgrid,&d1impv[2]);
    }
   else if(CHRcmp(name,"random_max_velocity",19)==0)
    {
                DBLr(fgrid,&ypar->d1rv_max[0]);
                DBLr(fgrid,&ypar->d1rv_max[1]);
                DBLr(fgrid,&ypar->d1rv_max[2]);
    }
   else if(CHRcmp(name,"random_min_velocity",19)==0)
   {
       DBLr(fgrid,&ypar->d1rv_min[0]);
       DBLr(fgrid,&ypar->d1rv_min[1]);
       DBLr(fgrid,&ypar->d1rv_min[2]);
   }

		CHRr(fgrid,name);
  }

/*
	else if(CHRcmp(name,"distance_to_bottom",11)==0)
    {    
 		DBLr(fgrid,&ypar->d1dxyz[0]);
		DBLr(fgrid,&ypar->d1dxyz[1]);
		DBLr(fgrid,&ypar->d1dxyz[2]);
  }*/
  fclose(fgrid);

  V3DNor(tmp,ypar->d1lv[0],ypar->d1lv[1],ypar->d1lv[2]);

  ypar->d2gxyz=TalDBL2(7,ypar->ngrid+1);

  fgrid=fopen(name_grid,"r");

  while(FILEND(fgrid)==0) 
  { 
	if(CHRcmp(name,"grid_number",11)==0)    
    {    
		CHRr(fgrid,name); 

		for(i=0;i<ypar->ngrid;i++)
		{
			for(j=0;j<6;j++)
			{
			DBLr(fgrid,&ypar->d2gxyz[j][i]);
		}
                  //   ypar->d2gxyz[1][i]=ypar->d2gxyz[1][i]+0.2;
//ypar->d2gxyz[2][i]=ypar->depheg;

		}
	}
	CHRr(fgrid,name);
  }


/* read mesh of particle  */

CHRr(fgid,name);
index=0;

while(FILEND(fgid)==0) 
{ if(CHRcmp(name, "Tetrahedra",10)==0)
  {
   	  CHRr(fgid,name); 
	  INTr(fgid,&ntyp);
  }
else if(CHRcmp(name, "Coordinates",11)==0)
{    
	CHRr(fgid,name); 
while((CHRcmp(name, "end",3)!=0))//&&(CHRcmp(name, "End",3)!=0))
{
 ypar->nnode++;
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
}

}
else if(CHRcmp(name, "Elements",8)==0)
{    if(index==0)
{index=index+1;
	CHRr(fgid,name); 
while((CHRcmp(name, "end",3)!=0))//&&(CHRcmp(name, "End",3)!=0))
{
 ypar->nelem++;
 if(ntyp==4){
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 }
 else if(ntyp==10)
 { CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 

 }
}
}
else 
{
	CHRr(fgid,name); 
while((CHRcmp(name, "end",3)!=0))//&&(CHRcmp(name, "End",3)!=0))
{
 ypar->nelemb++;
 if(ntyp==4){
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
  CHRr(fgid,name); 

  }
 else if(ntyp==10)
 { CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
  CHRr(fgid,name); 

 }
}
}
}
    CHRr(fgid,name); 

}
fclose(fgid);


  fgid=fopen(name_mesh,"r");

ypar->i2elto=TalINT2(ypar->nelem+1,ntyp+1);
ypar->i2eltob=TalINT2(ypar->nelemb+1,ntyp/2+1+1);
i1nobf=TalINT1(ypar->nnode+1);
ypar->d2xyz=TalDBL2(4,ypar->nnode+1);
ypar->i1elbe=TalINT1(ypar->nelem+1);


index=0;
while(FILEND(fgid)==0) 
{ 
	
if(CHRcmp(name, "Elements",8)==0)
{
	if(index==0)
   {
	index=index+1;
	for(i=0;i<ypar->nelem;i++)
	{    
 CHRr(fgid,name); 
 for(j=0;j<ntyp;j++)
 {
 INTr(fgid,&ypar->i2elto[i][j]);
 ypar->i2elto[i][j]=ypar->i2elto[i][j]-1;
 }
  CHRr(fgid,name); 

 }
}
else
{	for(i=0;i<ypar->nelemb;i++)
	{    
 CHRr(fgid,name); 
 for(j=0;j<(ntyp/2+1);j++)
 {
 INTr(fgid,&ypar->i2eltob[i][j]);
 ypar->i2eltob[i][j]=ypar->i2eltob[i][j]-1;
 }
  CHRr(fgid,name); 

 }
}
}
else if(CHRcmp(name, "Coordinates",11)==0)
{    	if(index==0)
   {

		for(i=0;i<ypar->nnode;i++)
		{
			  CHRr(fgid,name); 

		for(j=0;j<3;j++)
		{
			DBLr(fgid,&ypar->d2xyz[j][i]);
		}
		}
}
}

    CHRr(fgid,name); 
}

if(ntyp==10)
{

for(i=0;i<ypar->nelem;i++)
        {
     for(j=4;j<10;j++)
         {
         k=j-4;
         l=k+1;
         if(k==2)l=0;
         if(k>2)
         {
          k=j-7;
          l=3;
         }
         jj=ypar->i2elto[i][j];
         kk=ypar->i2elto[i][k];
         ll=ypar->i2elto[i][l];


                          for(m=0;m<3;m++)
                {
                        ypar->d2xyz[m][jj]=0.5*(ypar->d2xyz[m][ll]+ypar->d2xyz[m][kk]);

                }


          }

}
}
/*assign layer number as facenode  */
//if(index2==0)
//{

for(i=0;i<(ypar->nnode);i++)
{ 
	i1nobf[i]=0;   //mark 1 on boundary node 

}

for(i=0;i<(ypar->nelemb);i++)
{ for(j=0;j<(ntyp/2+1);j++)
{ 
	i1nobf[ypar->i2eltob[i][j]]=1;   //mark 1 on boundary node 
} 
}

for(i=0;i<(ypar->nelem);i++)
{
ypar->i1elbe[i]=0;
  for(j=0;j<ntyp;j++)
  {
	  if(i1nobf[ypar->i2elto[i][j]]>0)
		 ypar->i1elbe[i]=1;
  }

}

ypar->d1cg=TalDBL1(3);


ypar->d1cg[0]=R0;
ypar->d1cg[1]=R0;
ypar->d1cg[2]=R0;
for(i=0;i<(ypar->nnode);i++)
{
	 	for(j=0;j<3;j++)
		{
	     ypar->d1cg[j]=ypar->d1cg[j]+ypar->d2xyz[j][i];
		}
}

ypar->d1cg[0]=ypar->d1cg[0]/ypar->nnode;
ypar->d1cg[1]=ypar->d1cg[1]/ypar->nnode;
ypar->d1cg[2]=ypar->d1cg[2]/ypar->nnode;


for(i=0;i<(ypar->nnode);i++)
{dis=R0;
	 	for(j=0;j<3;j++)
		{
	     dis=dis+(ypar->d1cg[j]-ypar->d2xyz[j][i])*(ypar->d1cg[j]-ypar->d2xyz[j][i]);
		}
		dis=SQRT(dis);
		if(dis>ypar->radius)ypar->radius=dis;
}

/*generate random number for postion and rotate angle  */

  ypar->max=BEPSILON;
  //ypar->ymax=EPSILON;
  //ypar->zmax=EPSILON;
  
  for(i=0;i<ypar->ngrid;i++)
  {
	  tmp=ypar->d2gxyz[0][i]*ypar->d1lv[0]+ypar->d2gxyz[1][i]*ypar->d1lv[1]+ypar->d2gxyz[2][i]*ypar->d1lv[2];
	  if(tmp<ypar->max)ypar->max=tmp;
	  //if(ypar->d2xyz[1][i]<ypar->xmax)ypar->xmax=ypar->d2xyz[1][i];
	  //if(ypar->d2xyz[2][i]<ypar->xmax)ypar->xmax=ypar->d2xyz[2][i];
	  
  }

ntelem=ypar->ntpar*ypar->nelem;
ntnode=ypar->ntpar*ypar->nnode;
    ntelemb=ypar->ntpar*ypar->nelemb;

/* Read container mesh or .Y3D */

  
 
  SETLINEBUF(ydc->fcheck);  
  CHRr(ydc->finp,name);
  while(FILEND(ydc->finp)==0) 
  { 
    if(CHRcmp(name, "$YSTOP",6)==0)
    { 
      CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
      return 0;
    }
    else if(CHRcmp(name, "$YDOIT",6)==0)
    {
       //        printf("YDOIT, OLD element nodes:%ld %ld %ld %ld\n",yde->i2elto[0][0],yde->i2elto[0][1],yde->i2elto[0][2],yde->i2elto[0][3]);

        index_sp=0;
      CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
		  for(i=0;i<yde->nelemi;i++)
          {
	        if(yde->i2elto[i][yde->nelno-1]>ypar->nindex)ypar->nindex=yde->i2elto[i][yde->nelno-1];
			for(j=0;j<(yde->nelno-1);j++)
			{
				ydn->i1nind[yde->i2elto[i][j]]=yde->i2elto[i][yde->nelno-1];
			}
            }
if(index_sp==1)
{
if(ysp->mspd==0)index_sp=1;
		  Yspini(ysp,ypar->nindex,ypar->ntpar,icon);

        ysp->d1impv[0]=d1impv[0];
        ysp->d1impv[1]=d1impv[1];
        ysp->d1impv[2]=d1impv[2];
        
      for(i=0;i<ysp->nspd;i++)
	  {

		  /* to be modified properly   */
		  for(n=0;n<ydn->nnopo;n++)
    	  {   
			 // ysp->i1con[i]=1;
			  if(ydn->i1nind[n]==(i+1))
			  {
				  ysp->d2pvt[0][i]=ydn->d2nvc[0][n];
				  ysp->d2pvt[1][i]=ydn->d2nvc[1][n];
				  ysp->d2pvt[2][i]=ydn->d2nvc[2][n];
				  if(ydn->i1nopr[n]!=(-1))ysp->i1sppr[i]=ydn->i1nopr[n];}
			  
	       }
		  /* end */
		  
          if((ydb->i1bnvx[ysp->i1sppr[i]]!=0)||
			  (ydb->i1bnvy[ysp->i1sppr[i]]!=0)||
			  (ydb->i1bnvz[ysp->i1sppr[i]]!=0)) ysp->i1con[i]=1;


	if (index_sp==1)
	{
				  
		  
		  for(n=0;n<10;n++)intg[n]=R0;
	     
	  	 for(k=0;k<yde->nelemi;k++)
		 {	 if(yde->i2elto[k][yde->nelno-1]==(i+1))

		 {
		 for(m=0;m<4;m++)
		 {
	     points[m][0]=ydn->d2ncc[0][yde->i2elto[k][m]];
		 points[m][1]=ydn->d2ncc[1][yde->i2elto[k][m]];
		 points[m][2]=ydn->d2ncc[2][yde->i2elto[k][m]];
		 }
		 ro=ydp->d1pero[yde->i1elpr[k]];
		 mirtichRoutine(points,intg,ro);
	 }
		 }

     for(n=0;n<10;n++)
      intg[n]=intg[n]*mult[n];

     //this is an error in the original publication!!!!!!!!!!!!!!!!!
     //they dont have the three!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ysp->d1mass[i] = intg[0];  //density = 1.0 assumed 

     //center of mass 
    ysp->d2pp[0][i]=intg[1]/ysp->d1mass[i];
    ysp->d2pp[1][i]=intg[2]/ysp->d1mass[i];
    ysp->d2pp[2][i]=intg[3]/ysp->d1mass[i];
	cm[0]=ysp->d2pp[0][i];
	cm[1]=ysp->d2pp[1][i];
	cm[2]=ysp->d2pp[2][i];

    
    inertia[0][0] = intg[5] + intg[6] -ysp->d1mass[i]*(cm[1]*cm[1]+cm[2]*cm[2]);
    inertia[1][1] = intg[4] + intg[6] -ysp->d1mass[i]*(cm[0]*cm[0]+cm[2]*cm[2]);
    inertia[2][2] = intg[4] + intg[5] -ysp->d1mass[i]*(cm[1]*cm[1]+cm[0]*cm[0]);

    inertia[0][1] = -(intg[7] -ysp->d1mass[i]*cm[1]*cm[0]);
    inertia[1][0]=inertia[0][1];

    inertia[1][2] = -(intg[8] -ysp->d1mass[i]*cm[1]*cm[2]);
    inertia[2][1]=inertia[1][2];

    inertia[0][2] = -(intg[9] -ysp->d1mass[i]*cm[2]*cm[0]);
    inertia[2][0]=inertia[0][2];

//    inertia*=((REAL)(1.0/3.0));

	for(m=0;m<3;m++)
	  {
		v1[m]=R0;
		v2[m]=R0;
		v3[m]=R0;

		  for(n=0;n<3;n++)
	  {
//       inertia[m][n]=inertia[m][n]*(R1/R3);
          inertia[m][n]=inertia[m][n];	     
	  }
	  }

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
       ysp->d1mass[i]*=(-1.0);
      }
    
  
	solveSymetricEigenProblem(inertia, v1,v2,v3,
						&ysp->d2prine[0][i],&ysp->d2prine[1][i],&ysp->d2prine[2][i]);
		for(m=0;m<3;m++)
	  {
		ysp->d2prinn1[m][i]=v1[m];
		ysp->d2prinn2[m][i]=v2[m];
		ysp->d2prinn3[m][i]=v3[m];
		}


	  }
	}
}
	  /*check condition  */
	  index=0;
	  for(i=0;i<ydb->nbcon;i++)
	  {
		  if((ydb->i1bnvx[i]==0)&&(ydb->i1bnvy[i]==0)&&(ydb->i1bnvz[i]==0))
            {
              index=i+1;
            }

	  }
	  if(index==0)
	  {   //ydb->nbcon+=1;
  		  ydb->d1bnvx[ydb->nbcon]=ydb->d1bnvx[0];
		  ydb->d1bnvy[ydb->nbcon]=ydb->d1bnvy[0];
		  ydb->d1bnvz[ydb->nbcon]=ydb->d1bnvz[0];
		  ydb->d1bnax[ydb->nbcon]=ydb->d1bnax[0];
		  ydb->d1bnay[ydb->nbcon]=ydb->d1bnay[0];
		  ydb->d1bnaz[ydb->nbcon]=ydb->d1bnaz[0];
		  ydb->d1bnfx[ydb->nbcon]=ydb->d1bnfx[0];
		  ydb->d1bnfy[ydb->nbcon]=ydb->d1bnfy[0];
		  ydb->d1bnfz[ydb->nbcon]=ydb->d1bnfz[0];
		  ydb->i1bnvx[ydb->nbcon]=ydb->i1bnvx[0];
		  ydb->i1bnvy[ydb->nbcon]=ydb->i1bnvy[0];
		  ydb->i1bnvz[ydb->nbcon]=ydb->i1bnvz[0];
        FREE(intg);
        FREE(points);
        FREE(inertia);
        FREE(v3);
        FREE(v2);
        FREE(v1);
  		  ydb->i1bnvx[0]=0;
		  ydb->i1bnvy[0]=0;
		  ydb->i1bnvz[0]=0;

		  for(i=0;i<ydn->nnopo;i++)
		  {
			  if(ydn->i1nopr[i]==0)ydn->i1nopr[i]=ydb->nbcon;
		  }
                  for(i=0;i<ydn->mnopo;i++)
                  {
                          if(ydn->i1nopr[i]==(-1))ydn->i1nopr[i]=0;
                  }

                  for(i=0;i<ysp->mspd;i++)
                   {
                    if(ysp->i1sppr[i]==0) ysp->i1sppr[i]=ydb->nbcon;
                   }


                  for(i=0;i<ysp->mspd;i++)
                   {
                    if(ysp->i1sppr[i]==(-1)) ysp->i1sppr[i]=0;
                   }

                ydb->nbcon+=1;

	  }
          else
            {
                  for(i=0;i<ydn->mnopo;i++)
                  {
                          if(ydn->i1nopr[i]==(-1))ydn->i1nopr[i]=index-1;
                  }
                  for(i=0;i<ysp->mspd;i++)
                   {
                    if(ysp->i1sppr[i]==(-1)) ysp->i1sppr[i]=index-1;
                   }

             }
        
        /*check material property  */
/*        
        ydp->d1peks[ydp->nprop]=ydp->d1peks[0];
        ydp->d1pefr[ydp->nprop]=ydp->d1pefr[0];
        ydp->d1pela[ydp->nprop]=ydp->d1pela[0];
        ydp->d1pemu[ydp->nprop]=ydp->d1pemu[0];
        ydp->d1pepe[ydp->nprop]=ydp->d1pepe[0];
        ydp->d1pero[ydp->nprop]=ydp->d1pero[0];
        ydp->d1pesf[ydp->nprop]=ydp->d1pesf[0];
        ydp->d1pepsf[ydp->nprop]=ydp->d1pepsf[0];
        ydp->d1pevf[ydp->nprop]=ydp->d1pevf[0];
        ydp->d1pepf[ydp->nprop]=ydp->d1pepf[0];
        ydp->i1ptyp[ydp->nprop]=ydp->i1ptyp[0];
        
            for(i=0;i<yde->nelem;i++)
            {
                if(yde->i1elpr[i]==0)yde->i1elpr[i]=ydp->nprop;
            }
            for(i=0;i<yde->melem;i++)
            {
                if(yde->i1elpr[i]==(-1))yde->i1elpr[i]=0;
            }
            
        
            ydp->nprop+=1;
            
*/


if(inewmat==1)
{

        ydp->d1peks[ydp->nprop]=d1tmp[0];
        ydp->d1pefr[ydp->nprop]=d1tmp[5];
        ydp->d1pela[ydp->nprop]=d1tmp[1];
        ydp->d1pemu[ydp->nprop]=d1tmp[2];
        ydp->d1pepe[ydp->nprop]=d1tmp[3];
        ydp->d1pero[ydp->nprop]=d1tmp[4];
        ydp->d1pesf[ydp->nprop]=d1tmp[6];
        ydp->d1pepsf[ydp->nprop]=d1tmp[7];
        ydp->d1pevf[ydp->nprop]=d1tmp[8];
        ydp->d1pepf[ydp->nprop]=d1tmp[9];
        ydp->i1ptyp[ydp->nprop]=ydp->i1ptyp[0];
            ydp->nprop+=1;

}



//    FREE(i2eltob);
    FREE(i1nobf);
/*    FREE(v1);
	FREE(v2);
	FREE(v3);
    FREE(inertia);
	FREE(points);
	FREE(intg);
*/

/*
        FREE(intg);
        FREE(points);
        FREE(inertia);
        FREE(v3);
        FREE(v2);
        FREE(v1);
*/
      return 1;
    }
    else if(CHRcmp(name,"/*",2)==0)    /* read and ignore comments */
    { 
      icount=0;
      do
      { 
        CHRr(ydc->finp,name); icount++;
        if(icount>100)
        { 
          CHRw(stderr,"Yrd: too long comment near - ");
          CHRw(stderr,name);
          CHRwcr(stderr);      
          return 0;
        } 
      } while((FILEND(ydc->finp)==0)&&(CHRcmp(name,"*/",2)!=0));
    }
    else if(CHRcmp(name, "/YD/",4)==0)
    { 
      CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
      if(CHRcmp(name,"/YD/YDC/",8)==0)        /* read control data    */
      { 
        Yrdc(ydc,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDE/",8)==0)   /* read data elements   */
      { 
        Yrde(yde,ydc->finp,name,ntelem);

      }
      else if(CHRcmp(name,"/YD/YDJ/",8)==0)   /* read data elements   */
      {
        Yrdj(ydj,ydc->finp,name);

      }
      else if(CHRcmp(name,"/YD/YDI/",8)==0)  /* read data interaction */
      { 
        Yrdi(ydi,ydc->finp,name);

      }
      else if(CHRcmp(name,"/YD/YDN/",8)==0)  /* read data nodes       */
      { 
        Yrdn(ydn,ydc->finp,name,ntnode,ydc->dcimvelx,ydc->dcimvely,ydc->dcimvelz);

      }
      else if(CHRcmp(name,"/YD/YDO/",8)==0)  /* read output spec      */
      { 
        Yrdo(ydo,ydc->finp,name);

      }
      else if(CHRcmp(name,"/YD/YDP/",8)==0)  /* read data properties  */
      { 
        Yrdp(ydp,ydc->finp,name,yde);

      }
      else if(CHRcmp(name,"/YD/YDB/",8)==0)  /* read data properties  */
      { 
		    Yrdb(ydb,ydc->finp,name);

      }
      else if(CHRcmp(name,"/YD/YSP/",8)==0)  /* read data properties  */
      { 
		  Yrsp(ysp,ydc->finp,name,ypar->ntpar);
      }
      else if(CHRcmp(name,"/YD/YDX/",8)==0)  /* read data properties  */
      {
		  Yrdx(ydx,ydn,ydc->finp,name,ntelemb,yde->nelemi);

      }
      else
      { 
        CHRw(stderr,"Yrd: unknown name: ");
        CHRw(stderr,name); 
        CHRwcr(stderr);   
        return 0;
      } 
    }
    else
    { 
      CHRw(stderr,"Yrd: unknown name: ");
      CHRw(stderr,name); 
      CHRwcr(stderr);
      return 0;
    }
    CHRr(ydc->finp,name); 
  }
  fclose(ydc->finp);
  fclose(ydc->fcheck);
  fclose(fgrid);
  fclose(fgid);
//  FREE(i2eltob);
  FREE(i1nobf);

  return 0;
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/


