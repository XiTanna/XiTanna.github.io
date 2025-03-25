#include <gromacs/fileio/confio.h>
#include <gromacs/fileio/xvgr.h>
#include <gromacs/utility/arraysize.h>
#include <gromacs/commandline/cmdlineinit.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <gromacs/utility/smalloc.h>
#include <gromacs/topology/topology.h>
#include <gromacs/topology/index.h>
#include <gromacs/utility/fatalerror.h>
#include <gromacs/mdtypes/inputrec.h>
#include <gromacs/fileio/tpxio.h>
#include <gromacs/fileio/trxio.h>
#include <gromacs/commandline/pargs.h>
#ifdef HAVE_CONFIG_H
#endif


/****************************************************************************/
/* This program calculates the local density in 3-D.                  */
/* fengjie peng july 2016                                                 */
/* xitan add the xtc not match promble                           */
/* xitan add 2 dynamic output
/* attention the differenrce of load position center and no center, 
/*it can load no continue result now use center always not uese the no center version function*/ 
/****************************************************************************/

static int num2ave=20;
static real zmin  = 0;
static real zmax  = 0;
static real rcut  = 1.0;
static rvec nslice= {1, 1, 1};
static int  rmind = 0;  /* this is for mismatch in index when using xtc file while xtc did not storage all groups*/
static rvec min   = {0.0, 0.0, 0.0};
static rvec max   = {0.0, 0.0, 0.0};
static real kB    = 1.3806505;/* Boltzmann constant without 10^-23 */
static real NA    = 6.02214129;  /* Avogadro constant without 10^23 */
static real fcomb = 138.935485; /* kJ mol^−1 nm e^−2 factor in coulomb interaction */

double * CreateVector(int rows);
double ** CreateMatrix(int rows,int cols);
double *** CreateMatrix_3d(int m,int n, int p);
double pbc(double dx,double box);
int atomindex2mol(int *molindex,t_topology *top);
int load_positionCenter(rvec *x0,int nm,int napm,int *index,double **pos,double Lbox[3],double *deloc, double totMNC);
int load_position(rvec *x0,int nm,int napm,int nStart,double ***posmolecule,double Lbox[3]);
static void corr_print(const char *fn,const char *title,const char *yaxis,double dbinx,double dbiny,double dbinz,int nx,int ny,int nz,
         double ***numden,const char **dens_opt,char *grpname[],const gmx_output_env_t* oenv);
void do_corr(const char *trx_file, const char *ndx_file, const char *msd_file, int *tot_atom,char *grpname[],gmx_bool bMOL,const char **dens_opt,
         t_topology *top,const gmx_output_env_t* oenv)
{
  // t_trxframe  fr;
  rvec        *x0;      /* coordinates without pbc */
  matrix       box;     /* box (3x3) */
  t_trxstatus *status;
  real        t,kden;
  const char  *ylabel = NULL;
  int         natoms;  /* nr. atoms in trj */
  int         i,j,k,nx,ny,nz,mex,mey,mez,step,nm,nstart,napm;
  // int         flags = TRX_READ_X;   /* read only position */
  double      Lbox[3],dt,halfLbox[3];
  int     **index;
  int         molindex[top->atoms.nr];
  double      **posCen,dbinx,dbiny,dbinz,*mass,tmass,***numden;                /* some 2d-matrix to be used */
  FILE *outputfile;
  // read_first_frame(oenv,&status,trx_file,&fr,flags);
  atomindex2mol(molindex,top);
  natoms=read_first_x(oenv, &status, trx_file, &t, &x0, box);
  if (natoms!=top->atoms.nr)//lzeng add in 31.08.2021
  {
	  printf("\n|* warning! *|\n");
	  printf("\n|* warning! *|\n");
	  printf("\n|* warning! *|\n");
	  printf("\n|* warning! *|\n");
	  printf("\n|* warning! *|\n");
	  printf("\n|* warning! *|\n");
	  printf("\n|* warning! *|\n");
	  printf("number of atoms is %d in tracjectory file, %d in topfile\n",natoms,top->atoms.nr);
	  printf("preseted number of atoms not in tracjectory file is %d\n",rmind);
	  if (rmind!=(top->atoms.nr-natoms))
	  {
		  rmind=top->atoms.nr-natoms;
	  	  printf("now the number of atoms not in tracjectory file is rested as: %d\n",rmind);
	  }
  }
  Lbox[0]   = box[XX][XX]; halfLbox[0] = 0.5*Lbox[0];
  Lbox[1]   = box[YY][YY]; halfLbox[1] = 0.5*Lbox[1];
  Lbox[2]   = box[ZZ][ZZ]; halfLbox[2] = 0.5*Lbox[2];
  
  if(zmax == 0){ zmax = Lbox[2]; }
  fprintf(stdout," \n ");
  snew(index,1);
  get_index(&top->atoms,ndx_file,1,tot_atom,index,grpname);

  if(max[0] == 0){
    max[0] = Lbox[0];
  }
  if(max[1] == 0){
    max[1] = Lbox[1];
  }
  if(max[2] == 0){
    max[2] = Lbox[2];
  }
  nx   = (int)(nslice[0]+0.5);
  ny   = (int)(nslice[1]+0.5);
  nz   = (int)(nslice[2]+0.5);  

  nstart = index[0][0];
  napm   = top->mols.index[1+molindex[nstart]] - top->mols.index[molindex[nstart]];
  if(!bMOL){napm = 1;}
  nm     = tot_atom[0]/napm;
  dbinx  = (max[0] - min[0])/nx;
  dbiny  = (max[1] - min[1])/ny;
  dbinz  = (max[2] - min[2])/nz;
  posCen = CreateMatrix(nm,3);
  mass   = CreateVector(napm);
  numden = CreateMatrix_3d(nx,ny,nz);
  tmass  = 0.0;
  for(i=0;i<napm;i++){
	  mass[i] = top->atoms.atom[nstart+i].m;
	  tmass += mass[i];
  }
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
       for(k=0;k<nz;k++){
	  numden[i][j][k] = 0.0;
      }
    }
  }
  switch (dens_opt[0][0])
    {
        case 'n': ylabel = "Number density (nm\\S-3\\N)"; break;
		case 'm': ylabel = "Density (kg m\\S-3\\N)"; break;
        case 'c': ylabel = "Charge density (e nm\\S-3\\N)"; break;
    }
  outputfile=xvgropen(msd_file,"density.",output_env_get_xvgr_tlabel(oenv),ylabel,oenv);
  fprintf(outputfile,"# density of %s in 3-D\n",grpname[0]);
  if(dens_opt[0][0] == 'n'){
		fprintf(outputfile,"# location x,y,z (nm), number density (#/nm^3)\n");
  }
  if(dens_opt[0][0] == 'm'){
		fprintf(outputfile,"# location x,y,z (nm), mass density (kg/m^3)\n");
  }
  if(dens_opt[0][0] == 'c'){
		fprintf(outputfile,"# location x,y,z (nm), charge density (e/nm^3)\n");
  }
  for(i=0; i<nx; i++)
		for(j=0; j<ny; j++)
			for(k=0; k<nz; k++)
				fprintf(outputfile," %4g %4g %4g ",min[0]+(i+0.5)*dbinx,min[1]+(j+0.5)*dbiny,min[2]+(k+0.5)*dbinz);
  fprintf(outputfile,"\n");
	
  kden = 1; // for number density ( unit #/nm^3);
  step = 0;
  do{
	  
	  load_positionCenter(x0,nm,napm,index[0],posCen,Lbox,mass,tmass);
	  
	  for(i=0;i<nm;i++){
		  if(dens_opt[0][0] == 'm'){
			  kden = 0;
			  for(k=0;k<napm;k++){
				  kden += top->atoms.atom[nstart+i*napm+k].m;
			  }
			  kden /= (NA/10); // mass density has a unit of kg/m^3
		  }
		  if(dens_opt[0][0] == 'c'){
			  kden = 0;
			  for(k=0;k<napm;k++){
				  kden += top->atoms.atom[nstart+i*napm+k].q;
			  }// charge density has a unit of e/nm^3
		  }
		  mex = (int)((posCen[i][0]-min[0])/dbinx + 1) -1;
		  mey = (int)((posCen[i][1]-min[1])/dbiny + 1) -1;
		  mez = (int)((posCen[i][2]-min[2])/dbinz + 1) -1;
                  // by +1 and -1 to avoid (-1,0) -> 0 , another way is using floor() function;
		  if(mex<nx && mex>-1){
			  if(mey<ny && mey>-1){
				  if(mez<nz && mez>-1){
					  numden[mex][mey][mez] += kden/(dbinx*dbiny*dbinz);
				  }
			  }
		  }
	  }
	  if(step%num2ave==(num2ave-1))//clean data and printf
		{
			for(i=0; i<nx; i++)
				for(j=0; j<ny; j++)
					for(k=0; k<nz; k++)
					{
							  numden[i][j][k] /= num2ave;
						      fprintf(outputfile,"%10g  ",numden[i][j][k]);
							  numden[i][j][k] =0.0;
					}
			fprintf(outputfile,"\n");
		}
	  step++;
	  // dt = fr.time - fr.tpf;
	  // }while(read_next_frame(oenv,status,&fr));
  }while(read_next_x(oenv, status, &t, x0, box));

	fclose(outputfile);



	// fprintf(stdout," \n");
  //corr_print(msd_file,"density.",ylabel,dbinx,dbiny,dbinz,nx,ny,nz,numden,dens_opt,grpname,oenv);
}

int main_func(int argc,char *argv[])
{
  const char *desc[] = {
    "This program calculates the density.",
  };
  static const char *dens_opt[] =
    { NULL, "number", "mass", "charge", NULL };
  static gmx_bool bMOL       = TRUE;
  static gmx_bool bRmCOMM    = FALSE;
  t_pargs pa[] = {
	{ "-mol",    FALSE, etBOOL, {&bMOL},
      "count by molecules or atoms, default is treated as molecules" },
    { "-dens",    FALSE, etENUM, {dens_opt},
      "Density"},
	{ "-rm",FALSE, etINT, {&rmind},
      "number of atoms didn't storage in xtc file before selected group" },
    { "-sl",FALSE, etRVEC, {nslice},
      "Divide the box in #nr slices in x,y,z direction" },
	{ "-num2ave",FALSE, etINT, {&num2ave},
      "for dynamic output, number of step to output" },
    { "-min",FALSE, etRVEC, {min},
      "min (x,y,z) to calculate (nm),set to 0 when no limit"},
    { "-max",FALSE, etRVEC, {max},
      "max (x,y,z) to calculate (nm),set to 0 when no limit"}
  };

  t_filenm fnm[] = { 
    { efTRX, NULL, "run",  ffREAD },
    { efTPS, NULL, "run",  ffREAD }, 
    { efNDX, NULL, NULL,  ffOPTRD },
    { efXVG, NULL, "density", ffWRITE },
  };
#define NFILE asize(fnm)

  t_topology  *top;
  int         ePBC;
  const char  *trx_file, *tps_file, *ndx_file, *msd_file;
  int         *gnx; /* the selected groups' sizes */
  char        **grpname;
  gmx_output_env_t* oenv;


  if(!parse_common_args(&argc,argv,
                    PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_TIME_UNIT | PCA_CAN_TIME,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv))
{
	exit(0);
}
  trx_file = ftp2fn_null(efTRX,NFILE,fnm);
  tps_file = ftp2fn_null(efTPS,NFILE,fnm);
  ndx_file = ftp2fn_null(efNDX,NFILE,fnm);
  msd_file = ftp2fn_null(efXVG,NFILE,fnm);


  snew(gnx,1); 
  snew(grpname,1);
  top = read_top(tps_file,&ePBC);     /* read topology file */ 
 // get_index(&top->atoms,ndx_file,1,gnx,index,grpname);
/* the function you are going to use */
 do_corr(trx_file,ndx_file,msd_file,gnx,grpname,bMOL,dens_opt,
	  top,oenv); 
/*  the function end  */
  
  return 0;
}

double * CreateVector(int rows)
{
   double  *m;

   m = (double *) calloc((unsigned int) rows,sizeof(double ));

   return m;
}
double ** CreateMatrix(int rows,int cols)
{
   int  i;
   double  **m;
 
   m = (double **) calloc((unsigned int) rows,sizeof(double *));
   for (i=0; i < rows; i++) {
      m[i] =(double *)  calloc((unsigned int) cols,sizeof(double ));
   }
 
   return m;
}
double *** CreateMatrix_3d(int m,int n, int p)
{
  int i;
  double ***m2;
  m2 = (double ***) calloc((unsigned int) m, sizeof(double **));
  for(i=0;i<m;i++)
    m2[i] = CreateMatrix(n,p);
  return m2;
}
double pbc(double dx,double box)
{
  while (dx > box/2.0)
    dx -=  box;
  while (dx < -box/2.0) 
    dx +=  box;

  return dx;
}
int atomindex2mol(int *molindex,t_topology *top)
{
	int i,j;
	for(i=0;i<top->mols.nr;i++){
		for(j=top->mols.index[i]; j<top->mols.index[i+1]; j++ ){
			molindex[j]=i;
		}
	}
	return 0;
}
int load_positionCenter(rvec *x0,int nm,int napm,int *oneindex,double **pos,double Lbox[3],double *deloc, double totMNC)
{
  int i,j,k,molN,m;
  double center[3], tmpPos,tmpCenter,halfLbox[3];
  for(k=0;k<3;k++){
	  halfLbox[k] = Lbox[k]/2;
  }

  for(i=0;i<nm;i++){
    molN = i*napm;
    for(k=0;k<3;k++){  
      tmpCenter = 0.0;  // initiate for each molecule
      for(j=0;j<napm;j++){
	tmpPos = x0[oneindex[molN+j]-rmind][k];// xitan add to for x0 is not all group for xtc not all group
	// first unmap the atoms of big molecule into one whole molecule (use 1st atom--x0[molN][k] as absolute location)
	if(Lbox[k]>0)
	{ // when pbc=xy Lz=0;
    // fprintf(stdout,"halfLbox %.4f \n",Lbox[k]);
	while(tmpPos-x0[oneindex[molN]-rmind][k]> halfLbox[k])
	  tmpPos -= Lbox[k];
	while(tmpPos-x0[oneindex[molN]-rmind][k]<-halfLbox[k])
	  tmpPos += Lbox[k];
	}
	tmpCenter += tmpPos*deloc[j];
      }
      pos[i][k] = tmpCenter/totMNC;
	  // move molecules into box.
	  if(Lbox[k]>0)
	{ // when pbc=xy Lz=0;
	  while(pos[i][k]>Lbox[k])
		  pos[i][k] -= Lbox[k];
	  while(pos[i][k]<0)
		  pos[i][k] += Lbox[k];
	 }
    }
  }

  return 0;
}
int load_position(rvec *x0,int nm,int napm,int nStart,double ***posmolecule,double Lbox[3])
{
  int i,j,k,molN;
  double center[3], tmpPos[3],halfLbox[3];
  for(k=0;k<3;k++)
    halfLbox[k] = Lbox[k]/2;
  
  for(i=0;i<nm;i++){
    molN = nStart+i*napm;
    for(j=0;j<napm;j++)
      for(k=0;k<3;k++){
	tmpPos[k] = x0[molN+j][k];
	// first unmap the atoms of big molecule into one whole molecule (use 1st atom--x0[molN][k] as absolute location)
	while(tmpPos[k]-x0[molN][k]> halfLbox[k])
	  tmpPos[k] -= Lbox[k];
	while(tmpPos[k]-x0[molN][k]<-halfLbox[k])
	  tmpPos[k] += Lbox[k];
	
	posmolecule[i][j][k] = tmpPos[k];
      }
  }

  return 0;
}

static void corr_print(const char *fn,const char *title,const char *yaxis,double dbinx,double dbiny,double dbinz,int nx,int ny,int nz,
         double ***numden,const char **dens_opt,char *grpname[],const gmx_output_env_t* oenv)
{
  FILE *out;
  int  i,j,k;
  
  out=xvgropen(fn,title,output_env_get_xvgr_tlabel(oenv),yaxis,oenv);
    fprintf(out,"# density of %s in 3-D\n",grpname[0]);
	if(dens_opt[0][0] == 'n'){
		fprintf(out,"# location x,y,z (nm), number density (#/nm^3)\n");
	}
	if(dens_opt[0][0] == 'm'){
		fprintf(out,"# location x,y,z (nm), mass density (kg/m^3)\n");
	}
	if(dens_opt[0][0] == 'c'){
		fprintf(out,"# location x,y,z (nm), charge density (e/nm^3)\n");
	}
    
    for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
        for(k=0;k<nz;k++){
          fprintf(out,"  %10g    ",min[0]+(i+0.5)*dbinx);
          fprintf(out,"%10g    ",min[1]+(j+0.5)*dbiny);
          fprintf(out,"%10g    ",min[2]+(k+0.5)*dbinz);
          fprintf(out,"%10g",numden[i][j][k]);
          fprintf(out,"\n");
        }
      }
    }
  fclose(out);
}

 

int main(int argc, char** argv)
{
	return gmx_run_cmain(argc, argv, &main_func);
}

