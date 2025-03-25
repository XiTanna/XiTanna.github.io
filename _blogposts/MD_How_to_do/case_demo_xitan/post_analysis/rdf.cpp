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
#define _USE_MATH_DEFINES
#include <math.h>


int load_positionCenter(t_trxframe *fr, int nm, int atoms, int *index, double **posc, double *deloc);
int load_position(t_trxframe *fr, int nm, int atoms, int *index, double ***pos);
double periodicity(double dx, double box);
double *CreateVector(int cols);
double **CreateMatrix(int rows, int cols);
double ***CreateMatrix_3d(int m, int n, int p);
void getMolInfo(t_topology *top, int *ngx, int **index, int grp, int *nm, int *atoms, double **mass, double **charge);
void getLJParameter(t_topology *top, int **index, int grp1, int grp2, int atoms1, int atoms2, double ***c6, double ***c12);
static void corr_print(const char* fn, const char* title, const char* xaxis, const char* yaxis, int nsize, double* xx, double* yy, const gmx_output_env_t* oenv);

int main_func(int argc, char *argv[])
{
    const char *desc[] = {
        "this is a small test program meant to serve as a template "};

    int nbin = 100;
	real cutoff = 1.0;
    real lowPos = 0.0;
    real upPos = 30.0;
    t_pargs pa[] = {
        {"-nbin", FALSE, etINT, {&nbin}, "number of bin sets for number"},
		{"-cutoff", FALSE, etREAL, {&cutoff}, "cutoff(nm)"},
        {"-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)"},
        {"-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)"}};

    t_filenm fnm[] = {
        {efTRX, "-f", NULL, ffREAD},
        {efTPR, NULL, NULL, ffREAD},
        {efNDX, NULL, NULL, ffOPTRD},
        {efXVG, NULL, "density", ffWRITE}};

    /* GMX-2020 support data structure */
    t_topology *top;
    int ePBC;
    gmx_output_env_t *oenv;
    t_trxframe fr;
    t_trxstatus *status;
    int flags = TRX_READ_X; /* read only position     */
    char **grpname;         /* group names            */
    int ngrps = 2;          /* nr. of group(s)        */
    int **index;            /* indices for all groups */
    int *ngx;               /* sizes of groups        */
	const char* RDF_file;
	


#define NFILE asize(fnm)
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_CAN_TIME,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        exit(0);
    }
	
	RDF_file = ftp2fn_null(efXVG, NFILE, fnm);
    top = read_top(ftp2fn(efTPR, NFILE, fnm), &ePBC); /* read topology file */
    snew(grpname, ngrps);
    snew(index, ngrps);
    snew(ngx, ngrps);
    printf("\nSpecify %d group%s to analysis:\n", ngrps, (ngrps > 1) ? "s" : "");
    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), ngrps, ngx, index, grpname);
    /* ======================= Main body of code ========================= */

    int nm1, nm2, atoms1, atoms2;
    double ***pos1, ***pos2;                   /* size: [nm, atoms, 3]     */
    double **posc1, **posc2;                   /* size: [nm, 3]            */
    double *mass1, *mass2, *charge1, *charge2; /* size: [atoms]            */
    // double **c6, **c12;                           /* LJ parameter. size: [atoms1, atoms2] */

    getMolInfo(top, ngx, index, 0, &nm1, &atoms1, &mass1, &charge1);
    getMolInfo(top, ngx, index, 1, &nm2, &atoms2, &mass2, &charge2);
    // getLJParameter(top, index, 0, 1, atoms1, atoms2, &c6, &c12);

    read_first_frame(oenv, &status, ftp2fn_null(efTRX, NFILE, fnm), &fr, flags);

    pos1 = CreateMatrix_3d(nm1, atoms1, 3);
    pos2 = CreateMatrix_3d(nm2, atoms2, 3);
    posc1 = CreateMatrix(nm1, 3);
    posc2 = CreateMatrix(nm2, 3);
    int count[nbin];
	///Lz_left  Lz_right	how to denfine
    double 	dr,r,dis,Vol,Voldr;
    double 	gr[nbin];
    double	Lbox[3];
    int	i,j,k,m,me,step;
    double dx,dy,dz;
    dr = cutoff / nbin;
    memset(gr,0,sizeof(gr));
    memset(count,0,sizeof(count));
    int cN1, cN2;
    step=0;

    do
    {

    memset(count,0,sizeof(count));
    load_position(&fr, nm1, atoms1, index[0], pos1);
    load_position(&fr, nm2, atoms2, index[1], pos2);
    load_positionCenter(&fr, nm1, atoms1, index[0], posc1, mass1);
    load_positionCenter(&fr, nm2, atoms2, index[1], posc2, mass2);
    Lbox[0]= fr.box[XX][XX];        //当前帧盒子坐标dx,dy,dz
    Lbox[1]= fr.box[YY][YY];        //当前帧盒子坐标dx,dy,dz
    Lbox[2]= fr.box[ZZ][ZZ];        //当前帧盒子坐标dx,dy,dz
    Vol=Lbox[0]*Lbox[1]*Lbox[2];

   for (i=0;i<nm1;i++){
	if (posc1[i][2]>=lowPos && posc1[i][2]<=upPos ){
		cN1++;
		for (k=0;k<nm2;k++){
			dx=posc1[i][0]-posc2[k][0];
			dy=posc1[i][1]-posc2[k][1];
			dz=posc1[i][2]-posc2[k][2];
			if(dx > Lbox[0]/2){dx -= Lbox[0];} //periodicity
			if(dx < -Lbox[0]/2){dx += Lbox[0];}
			if(dy > Lbox[1]/2){dy -= Lbox[1];}
			if(dy < -Lbox[1]/2){dy += Lbox[1];}
			if(dz > Lbox[2]/2){dz -= Lbox[2];}
			if(dz < -Lbox[2]/2){dz += Lbox[2];}
			dis= sqrt(dx*dx+dy*dy+dz*dz);
			if (dis>0.00000001 && dis<=cutoff)
			{
				me =dis/dr;
				if(me>=0 && me<nbin){
				count[me] = count[me] + 1;
				cN2++;
				}
			}
	    }
	
	}
   }

for(i=0;i<sizeof(count)/sizeof(int);i++)
{
	count[i] = count[i]/cN1;
}

  
for (i=0;i<nbin;i++){
 	r=(i+1)*dr;
	Voldr=4/3*M_PI*(pow((i+1)*dr,3)-pow(i*dr,3));
	gr[i]=count[i]/Voldr/(nm2/Vol);
//printf("count=%d	nm1=%d	Vdr=%lf	V=%lf	%lf\n",count[i], nm1, Voldr, Vol,gr[i]);
  }
    step++;
	gr[nbin] += gr[nbin];
    }while (read_next_frame(oenv, status, &fr));

	for (i=0;i<nbin;i++){
		gr[i] =gr[i]/step;
	}
	

	FILE *output;
	if((output = fopen("rdf.xvg","w"))==NULL)
	{
		printf("\n\nfailed to open file.\n");
		exit(0);
	}
	for(i=0;i<nbin;i++)
	{
		fprintf(output, "%lf	%lf\n", (i+1)*dr, gr[i]);
	}
	//corr_print(RDF_file, "RDF", "r (nm)", "g(r)", length, r, gr, oenv);

    return 0;
}

//reference function
int load_positionCenter(t_trxframe *fr, int nm, int atoms, int *index, double **pos, double *deloc)
{
    int i, j, k, molN;
    double tmpPos, tmpCenter, halfLbox[3];
    double totNCM = 0.0;
    for (i = 0; i != atoms; ++i)
    {
        totNCM += deloc[i];
    }

    for (k = 0; k < 3; k++)
        halfLbox[k] = fr->box[k][k] / 2;

    for (i = 0; i < nm; i++)
    {
        molN = index[i * atoms];

        for (k = 0; k < 3; k++)
        {
            tmpCenter = 0.0; // initiate for each molecule
            for (j = 0; j < atoms; j++)
            {
                tmpPos = fr->x[molN + j][k];
                // first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
                while (tmpPos - fr->x[molN][k] > halfLbox[k])
                    tmpPos -= fr->box[k][k];
                while (tmpPos - fr->x[molN][k] < -halfLbox[k])
                    tmpPos += fr->box[k][k];

                tmpCenter += tmpPos * deloc[j];
            }
            pos[i][k] = tmpCenter / totNCM;
        }
    }

    return 0;
}

int load_position(t_trxframe *fr, int nm, int atoms, int *index, double ***posmolecule)
{
    int i, j, k, molN;
    double tmpPos[3], halfLbox[3];
    for (k = 0; k < 3; k++)
        halfLbox[k] = fr->box[k][k] / 2;

    for (i = 0; i < nm; i++)
    {
        molN = index[i * atoms];
        for (j = 0; j < atoms; j++)
            for (k = 0; k < 3; k++)
            {
                tmpPos[k] = fr->x[molN + j][k];
                // first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
                while (tmpPos[k] - fr->x[molN][k] > halfLbox[k])
                    tmpPos[k] -= fr->box[k][k];
                while (tmpPos[k] - fr->x[molN][k] < -halfLbox[k])
                    tmpPos[k] += fr->box[k][k];

                posmolecule[i][j][k] = tmpPos[k];
            }
    }
    return 0;
}

double periodicity(double dx, double box)
{
    while (dx > box / 2.0)
        dx -= box;
    while (dx < -box / 2.0)
        dx += box;
    return dx;
}

double *CreateVector(int rows)
{
    return (double *)malloc(rows * sizeof(double));
}

double **CreateMatrix(int rows, int cols)
{
    double **m = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        m[i] = (double *)malloc(cols * sizeof(double));
    }
    return m;
}

double ***CreateMatrix_3d(int m, int n, int p)
{
    double ***m2 = (double ***)malloc(m * sizeof(double **));
    for (int i = 0; i < m; i++)
        m2[i] = CreateMatrix(n, p);
    return m2;
}

void getMolInfo(t_topology *top, int *ngx, int **index, int grp, int *nm, int *atoms, double **mass, double **charge)
{
    int i;
    *atoms = 0;
    for (;;)
    {
        if (top->atoms.atom[index[grp][*atoms]].resind == top->atoms.atom[index[grp][0]].resind)
            (*atoms)++;
        else
            break;
    }
    *nm = ngx[grp] / *atoms;

    printf("\n### Group %d: \n", grp);
    printf("# Nr. of atoms in molecule : %d\n", *atoms);
    printf("# Nr. of molecules in group: %d\n", *nm);
    printf("# nr.   name     mass   charge\n");
    printf("------------------------------\n");
    double sumMass = 0.0, sumCharge = 0.0;
    snew(*mass, *atoms);
    snew(*charge, *atoms);
    for (i = 0; i != *atoms; ++i)
    {
        (*mass)[i] = top->atoms.atom[index[grp][i]].m;
        (*charge)[i] = top->atoms.atom[index[grp][i]].q;
        sumMass += (*mass)[i];
        sumCharge += (*charge)[i];
        printf("%5d%7s%9.3f%9.3f\n", i, *top->atoms.atomname[index[grp][i]], (*mass)[i], (*charge)[i]);
    }
    printf("------------------------------\n");
    printf("%12s%9.3f%9.3f\n", "Total", sumMass, sumCharge);
}

void getLJParameter(t_topology *top, int **index, int grp1, int grp2, int atoms1, int atoms2, double ***c6, double ***c12)
{
    int n1, n2, i, j;
    int ntypes = top->atomtypes.nr;
    int type;
    *c6 = CreateMatrix(atoms1, atoms2);
    *c12 = CreateMatrix(atoms1, atoms2);
    printf("\n### LJ Parameters(Group %d and %d):\n", grp1, grp2);
    printf("[  at1  at2]: %12s %12s\n", "c6", "c12");
    printf("---------------------------------------\n");
    int n = 0;
    for (i = 0; i != atoms1; ++i)
    {
        for (j = 0; j != atoms2; ++j)
        {
            n1 = index[grp1][i];
            n2 = index[grp2][j];
            type = ntypes * (top->atoms.atom[n1].type) + (top->atoms.atom[n2].type);
            (*c6)[i][j] = top->idef.iparams[type].lj.c6;
            (*c12)[i][j] = top->idef.iparams[type].lj.c12;
            if (n++ <= 12)
            {
                printf("[%5s%5s]: %12.4e %12.4e\n", *top->atoms.atomname[n1], *top->atoms.atomname[n2], (*c6)[i][j], (*c12)[i][j]);
            }
        }
    }
    if (atoms1 * atoms2 > 12)
    {
        printf("(.........)");
    }
    printf("\n\n");
}

int main(int argc, char **argv)
{
    return gmx_run_cmain(argc, argv, &main_func);
}
