/*! \file Y3D.c
 *  \brief main for Y3D model
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


/* for development only */
#define SHOW_STATS(t) \
{ \
  double sum; \
  sum = t[0] + t[1] + t[2] + t[3] + t[4]; \
  fprintf(stderr,"wall time %08.3f\n", sum); \
  fprintf(stderr,"Yfd:  %08.3f (%07.3f%%)\n", t[0], 100.0*t[0]/sum); \
  fprintf(stderr,"Ycd:  %08.3f (%07.3f%%)\n", t[1], 100.0*t[1]/sum); \
  fprintf(stderr,"Yid:  %08.3f (%07.3f%%)\n", t[2], 100.0*t[2]/sum); \
  fprintf(stderr,"Ysd:  %08.3f (%07.3f%%)\n", t[3], 100.0*t[3]/sum); \
  fprintf(stderr,"Yod:  %08.3f (%07.3f%%)\n", t[4], 100.0*t[4]/sum); \
}


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


static void Ycopyright(void)
{
  CHRw(stdout, "******************************************************");
  CHRw(stdout, "****************\n");
  CHRw(stdout, "** Copyright (C) 2008,");
  CHRw(stdout, "                                              **\n");
  CHRw(stdout, "** Imperial College of Science, Technology and Medicine (ICSTM).   **\n");
  CHRw(stdout, "** All rights reserved.**\n");
  CHRw(stdout, "** Implemented for you by Dr Jiansheng Xiang                       **\n");
  CHRwcr(stdout);
  CHRwcr(stdout);
  CHRw(stdout, "* This code is part of the Solidity developed by ICSTM.\n");
  CHRw(stdout, "*  \n");
  CHRw(stdout, "* This code is provided by copyright holders under");
  CHRw(stdout, " the GNU lesser \n");
  CHRw(stdout, "* General Public License (LGPL).");
  CHRw(stdout, " It is open source code; you can \n");
  CHRw(stdout, "* redistribute it and/or modify it under the terms ");
  CHRw(stdout, "of the GNU Lesser\n");
  CHRw(stdout, "* General Public License version 3.  \n");
  CHRw(stdout, "*  \n");
  CHRw(stdout, "* This code is distributed in the hope that it will ");
  CHRw(stdout, "be useful, \n");
  CHRw(stdout, "* but WITHOUT ANY WARRANTY; without even the implied ");
  CHRw(stdout, "warranty \n");
  CHRw(stdout, "* of MERCHANTABILITY or FITNESS FOR A PARTICULAR ");
  CHRw(stdout, "PURPOSE. See \n");
  CHRw(stdout, "* the GNU Lesser General Public License");
  CHRw(stdout, " for more details,\n");
  CHRw(stdout, "* http://www.gnu.org/licenses/lgpl-3.0.txt. \n");
  CHRw(stdout, "*  \n");
  CHRw(stdout, "* You should have received a copy of the GNU Lesser ");
  CHRw(stdout, "General Public\n");
  CHRw(stdout, "* License along with this code; if not, write to \n");
  CHRw(stdout, "* Dr Jiansheng Xiang j.xiang@imperial.ac.uk \n");
  CHRw(stdout, "* *************************************************");
  CHRw(stdout, "****************** *\n");
}


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief main program Y3D
 *  \param[in] argc number of command line arguments
 *  \param[in] argv array of command line arguments
 *  \par Details:
 *
 *  The program takes one argument which is the file name of the
 *  model definition file (*.Y3D).
 *
 *  The model input is parsed using Yrd(), and runs are defined and
 *  executed until there is no more input in the file.
 *
 *  Each run calls the following functions in sequence for each time
 *  step:
 *
 *  \b Yfd() // calculate nodal forces      \n
 *  \b Ycd() // calculate contact detection \n
 *  \b Yid() // detect contact interaction  \n
 *  \b Ysd() // solve equations             \n
 *  \b Yod() // output results              \n
 *
 *
 *
 */
int main(int argc, char **argv)
{
  /* display copyright */
  Ycopyright();

  if(argv[1]!=NULL)
  {
//    TIMER t;                 /* t is simulation start time           */
    CHR c1name[300];         /*2 c1name is name of the input file     */
    CHR c1name1[300];         /*2 c1name is name of the input file     */
    CHR c1name2[300];         /*2 c1name is name of the input file     */
    CHR c1name3[300];         /*2 c1name is name of the input file     */

    struct YD_struct yd;     /*2 yd is Y database                     */
    YDK ydk = &yd.ydk;       /*2 ydk is Y constants database          */
    YDC ydc = &yd.ydc;       /*2 ydc is Y control database            */
    YDE yde = &yd.yde;       /*2 yde is Y element database            */
    YDI ydi = &yd.ydi;       /*2 ydi is Y interaction database        */
    YDN ydn = &yd.ydn;       /*2 ydn is Y node database               */
    YDO ydo = &yd.ydo;       /*2 ydo is Y output database             */
    YDP ydp = &yd.ydp;       /*2 ydp is Y property database           */
    YDB ydb = &yd.ydb;       /*2 ydb is Y boundary condition database */
    YPAR ypar = &yd.ypar;
    YSP ysp = &yd.ysp;       /*2 ydb is Y boundary condition database */
    YDX ydx = &yd.ydx;
    YDJ ydj = &yd.ydj;
    YDJP ydjp = &yd.ydjp;


      DBL ts[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },tana[2]={0.01,0.01};

    int tid, nthreads;
      struct timespec now, tmstart,t1,t2;

  //    DBL seconds;
	  DBL dctime,diezonc,diezonp,tmp0,tmp1,dcgrax,dcgray,dcgraz,dcrmv;
      INT itmp,i,icheck,icount;
      FILE *fp=FILENULL;


        time_t raw=time(NULL);
        struct tm *tm=localtime(&raw);
  /*      printf("year: %d; month: %d; day: %d;\n",
                        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
        printf("hour: %d; minute: %d; second: %d\n",
                        tm->tm_hour, tm->tm_min, tm->tm_sec);

      if((tm->tm_year+1900)>=2020)
       {
          CHRw(stdout," Input files are corruputed, please check them.  \n");
          exit(1);
       }
      else if(((tm->tm_mon+1)>=10)&&((tm->tm_year+1900)==2019))
       {
          CHRw(stdout," Input files are corruputed, please check them.  \n");
          exit(1);
       }
      else if(((tm->tm_mday)>=16)&&((tm->tm_mon+1)==9)&&((tm->tm_year+1900)==2019))
       {
          CHRw(stdout,"Input files are corruputed, please check them.  \n");
          exit(1);
       }
*/


#pragma omp parallel private(tid)
      {

          /* Obtain thread number */
          tid = omp_get_thread_num();
          printf("Hello World from thread = %d\n", tid);

          /* Only master thread does this */
          if (tid == 0)
          {
              nthreads = omp_get_num_threads();
              printf("Number of threads = %d\n", nthreads);
          }

      }  /* All threads join master thread and disband */

      /*4 get name of the problem */
    CHRcpy(c1name, argv[1]);
	CHRcpy(c1name1, argv[2]);
	CHRcpy(c1name2, argv[3]);
      CHRcpy(c1name3, argv[4]);


    /* initialisation (NOP: automatic variables are already zeroed) */
    ydc->finp = FILENULL;
    ydc->fcheck = FILENULL;

    /* record start time */
//    t = TIME;
clock_gettime(CLOCK_REALTIME, &tmstart);
    /* initilalise constants */
    Yinit(ydk);
    ypar->r=rand()/RAND_MAX;

    /*4 Process while any input */
    if(Yrd(c1name,c1name1,c1name2, c1name3, &yd) > 0)
    {
      /* store current simulation starting timestamp */
      DBL dctime = ydc->dctime;

      /* configure per-simulation constants */
      Yconfig(ydk, ydc, ydb);
//        Ymass3TET(yde, ydn, ydp, ydk,nthreads);
//        YFEMsurface(ydx, yde->i1nelsur,yde->i2elsur,yde->i2elto, yde->nelem);

      CHRw(stdout, "NEW INPUT");
      CHRwcr(stdout);

      /* run current simulation */
      Ymd(yde, ydi, ydn, ydp, ydx, ydj, ydjp);
              Ymass3TET(yde, ydn, ydp, ydk,ydj,ydjp,nthreads,ydc->dcstec);

      for(ydc->ncstep=0; ydc->ncstep<ydc->mcstep; ydc->ncstep++)
      {

//          Ypacking(&yd,nthreads);
        clock_gettime(CLOCK_REALTIME, &t1);
          Yfd(yde, ydn, ydp, ydk,ydi,ydc,ydx,ydj,ydjp,ydb,nthreads);
          clock_gettime(CLOCK_REALTIME, &t2);
          ts[0] -= (double)((t1.tv_sec+t1.tv_nsec*1e-9) - (double)(t2.tv_sec+t2.tv_nsec*1e-9));
        /*4 calculate nodal forces            */
             Yhpd(ydx,ydn,ydc);

        /*4 calculate contact detection       */

        Ycd(yde, ydi, ydn, ydc->iwfast,nthreads,ydc->ncstep,ydj->i1elst);
          clock_gettime(CLOCK_REALTIME, &t1);
          ts[1] += (double)((t1.tv_sec+t1.tv_nsec*1e-9) - (double)(t2.tv_sec+t2.tv_nsec*1e-9));
        /*4 detect contact interaction        */

        Yid(yde, ydi, ydn, ydp, ydj, nthreads,ydc->dcstec,ydc->ncstep);
          clock_gettime(CLOCK_REALTIME, &t2);
          ts[2] -= (double)((t1.tv_sec+t1.tv_nsec*1e-9) - (double)(t2.tv_sec+t2.tv_nsec*1e-9));

       /*4  solve equations                  */

        Ysd(ydc, yde, ydn, ydo, ydb, ydk);
          clock_gettime(CLOCK_REALTIME, &t1);
          ts[3] += (double)((t1.tv_sec+t1.tv_nsec*1e-9) - (double)(t2.tv_sec+t2.tv_nsec*1e-9));
        /*4  output results                   */

        Yod(c1name, &yd);
          clock_gettime(CLOCK_REALTIME, &t2);
          ts[4] -= (double)((t1.tv_sec+t1.tv_nsec*1e-9) - (double)(t2.tv_sec+t2.tv_nsec*1e-9));
        /*4  update time ( - more accurate)  */
 
 //       Ywrs(c1name,&yd);
		/*4  update time                   */

//        ydc->dctime = ydc->dctime + (ydc->dcstec*(ydc->ncstep + 1));
        ydc->dctime = ydc->dctime + ydc->dcstec;

          if((ydc->ncstep%ydc->icoutf)==0)
          {
          CHRw(stderr,"frame ");
              INTw(stderr,ydc->icouti,0);
              CHRw(stderr,"/"); INTw(stderr,ydc->mcstep/ydc->icoutf,0);
              CHRwcr(stderr);
              clock_gettime(CLOCK_REALTIME, &now);
              double seconds = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));
              printf("wall time %fs\n", seconds);
              SHOW_STATS(ts);
              ts[0]=R0;
              ts[1]=R0;
              ts[2]=R0;
              ts[3]=R0;
              ts[4]=R0;

          }
      }

      if(ydc->icfrm==1)
      {
	      dcgrax=ydc->dcgrax+ydc->dcrmgrax;
	      dcgray=ydc->dcgray+ydc->dcrmgray;
	      dcgraz=ydc->dcgraz+ydc->dcrmgraz;
	      dcrmv=SQRT(SQR(ydc->dcrmvx)+SQR(ydc->dcrmvy)+SQR(ydc->dcrmvz));
	 for(i=0; i<ydn->nnopo; i++)
       {
               ydn->d2nfrm[0][i]=ydn->d2ncc[0][i];
	       ydn->d2nfrm[1][i]=ydn->d2ncc[1][i];
	       ydn->d2nfrm[2][i]=ydn->d2ncc[2][i];
           }

     for(ydc->ncstep=ydc->mcstep; ydc->ncstep<ydc->ncstepfrm+ydc->mcstep; ydc->ncstep++)
      {
        ydc->dcfrmca=SQR(ydc->dcstec*(ydc->ncstep-ydc->mcstep))/2.0;
  ydc->dcfrmca=ydc->dcfrmca*SQRT(SQR(dcgrax)+SQR(dcgray)+
      SQR(dcgraz))*0.250+dcrmv*ydc->dcstec*(ydc->ncstep-ydc->mcstep)*0.5;
  ydc->dcfrmca=0.000005;
Yfd(yde, ydn, ydp, ydk,ydi,ydc,ydx,ydj,ydjp,ydb,nthreads);
Ycd(yde, ydi, ydn, ydc->iwfast,nthreads,ydc->ncstep,ydj->i1elst);
Yid(yde, ydi, ydn, ydp, ydj, nthreads,ydc->dcstec,ydc->ncstep);
Ysd(ydc, yde, ydn, ydo, ydb, ydk);
Yod(c1name, &yd);

     ydc->dctime = ydc->dctime + ydc->dcstec;
      }
      }
    }

    CHRw(stderr, "   ***** Y3D HAS ORDERLY FINISHED *****\n");
      clock_gettime(CLOCK_REALTIME, &now);
      double seconds = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));
      printf("wall time %fs\n", seconds);
    /* display elapsed time */
//    TIMERw(stderr, t, 60, "minutes\n");

    /* display average time per module */
    SHOW_STATS(ts)
  }
  else
  {
    CHRw(stdout, "Double click the data file to run the program.\n");
  }
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/
