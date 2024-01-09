#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <direct.h>

#define NMAX 996 /*number of explicit levels*/

#define NMIN 2 /*1 if Case A; 2 if Case B*/
#define SIZE 200
#define ATOMZ 1
#define ATOMN 1/*number of nucleons*/
#define NC 2300 /*number of pseudo levels*/

#define NE 100 /*electron density in cm^-3*/
#define TE 300 /*electron temp*/

#define TR 10000 /*radiation temp*/
#define RAD 1 /*turn spontaneous radiative transitions on or off*/

#define CMB 0 /*turn CMBR on or off*/
#define DUST 0 /*turn radiation from dust on or off*/
#define FF 0 /*turn free-free radiation from the gas on or off*/
#define CONSTRADF 0

double hypgeo ( const int, const int, const int, const double );
void aval ( const int, const int, double [], double [], const double );
double anlmll ( const int, const int, const int, const double );
double anlmlu ( const int, const int, const int, const double );
double absorp ( const int, const int, const int, const int, const double, const double );
double avalaprx (const int, const int);
double bbradmat(const int, const int, const int, const int);
double hyprege(const int , const int , const int , const double );
double abnm ( const int, const int, const double [], const double [] );
double bfn ( const int, const double [] );
double bfn2 ( const int, const double [] );
double bbfactu ( const int, const int, const int);
double bbfactl ( const int, const int, const int);
double fact ( const int );
void bval ( const int, const int, const double [], const double [], double [], double [],
		double [], double [], const double );
double freq( const int, const int, const double );
double energynm ( const int, const int, const double);
void colratecoef( const int, const int, const double, const double, const int, double []);
double bbcolex(const int, const int, const double, const double, const int);
double bbcoldeex(const int, const int, const double, const int, const double);
double oscstrength( const int, const int, const double [], const double [], const double );
double oscstrengthnm( const int, const int, const double, const double );
void colionratecoef( const int, const double, const int, double [], const double );
void colmomcoef( const int, double [], double [], const int, const double, const double, const double [], const double, const int );
void colmomcoefbrock( const int, double [], double [], const int, const double, const double, const double [], const int, const double, const int );
void gullc( const int, double [], double [], double );
void radrecomb( const int, const double, double [], double [], double [] );
void photoion( const int, double [], double [], double [], double const, double const );
void stimrecomb( const int, const double, const int, double [], double [], double [], double const, double const );
double infsumnm(const int, const double, const double, const int, const double );
double infsumnlf(const int, const int, const double, const double, const int, const double, const int, double [] );
double infsumnle(const int, const int, const double, const double, const int, const double, const int, double [] );
double infsumnlfacc(const int, const int, const double, const double, const int, const double, const int,  double []);
double infsumnleacc(const int, const int, const double, const double, const int, const double, const int, double [] );
void rowvec(const int, const double, const double [], const double [], const double [], const int, double [], const double, const double);
void elgs ( double [][SIZE], int[], int );
void matinv( double [][SIZE], double [][SIZE], int);
void matmult(double [][SIZE], double [], double []);
void lagintermat(const int [], double [][NMAX]);
void stepgen(int []);
void stepgenold(int []);
void expvec(const double [], double [][NMAX], double []);
void condmat(  double [][NMAX+1],  double [][NMAX], double [][SIZE], int[]);
int fminn(const int, const int);
int fmaxx(const int, const int);
double jnu(const int, const int, const double, const double);
double mintens(const int, const double, const double, const double);
double uppernsumin(const int, const int, const double, const double, const int, const double, const int,  double []);
double uppernsumout(const int, const int, const double, const double, const int, const double, const int,  double []);
void iterativerefinement(double matrix[][SIZE], double [], double [], int);
double infnormvec(double [], int, int);
/*void gausselim(double [][SIZE+1], double [], double [], int );*/
void gausselim(double [][SIZE], double [], double [], int );
double matnorm( double matrix[][SIZE]);
double vecnorm (const double [], const int);
void matmultn(double [][SIZE],double [][SIZE], double [][SIZE], int);
int idx(const int, const int);
int idx2(const int, const int, const int);
double normvec1(double [], int, int);
int amidx(const int, const int, const int );
int apidx(const int, const int, const int );
double nround(double , int );
int findlten(double []);
double freefree(const int, const int, const double);
double freefreebf(const int, const double, const double);
double cmbr(const int, const int, const double);
double cmbrbf(const int, const double, const double);
double dust(const int, const int, const double);
double dustbf(const int, const double, const double);
int round_up_5(int);
char *makepath(const char[], const char[] );



int main()
{
	int i, n, m, temp, seek, maxn, j, nlimit, itconv,minmat, np,nlmax, n1=5,n2=10,n3=40,n4=150;
	double inff, dw, sum; /*dilution factor for rad field*/
	double rydberg, rydinf = 109737.31568508, elecmass, protmass, y, freq3, _jnu;
	double avalu[ NMAX ] = { 0 }, avall[ NMAX ] = { 0 }, anm, anl[ NMAX ], einsteina;
	double babsu[ NMAX ] = { 0 }, babsl[ NMAX ] = { 0 }, babsnm, bstimnm, bbexl, bbdexl;
	double bstimu[ NMAX ] = { 0 }, bstiml[ NMAX ] = { 0 };
	double cmomlp[ NMAX] = { 0 }, cmomup[ NMAX] = { 0 }, cmomle[ NMAX ] = { 0 }, cmomue[NMAX] = { 0 },cmomlh[ NMAX ] = { 0 }, cmomuh[NMAX] = { 0 };
	double bbex, bbdex, bfcolrate[ 2 ] = {0};
	double thetau[ NMAX ] = {0}, thetal[ NMAX ] = {0};
	double radrec[ NMAX ] = {0}, phion[NMAX] = {0}, stimrec[NMAX] = {0};
	double  veca[/*SIZE*/NMAX] = {0}, cmat[SIZE][SIZE]= {0};

	double (*lagmat)[NMAX] = calloc(1, sizeof(double[SIZE][NMAX]));
	double imat[SIZE][SIZE]={0}, cbn[SIZE]={0}, bn[NMAX] = {0},ident[SIZE][SIZE]={0};

	double (*mat)[NMAX+1] = calloc(1,sizeof(double[NMAX+1][NMAX+1]));

	double a2, a3, a4, a5, a6, a7, stim, b, aur, alr;
	int step[SIZE+1] = {0}, conv, ncor, nlcount, depoptest;
	double emptyrate, effrad,al, colionnl, upsumout;


	double* a1;
	double* emp1;
	double newbnl, error, lerror=1, inffin, inffout, bnladd, rowsum1, rowsum2,maxdiff,condno;
	clock_t tic, toc, ticinv, tocinv, ticit, tocit,tic1,toc1,tic2,toc2;
	int calcbn, depop;
	double inground=0, outground = 0, scale=1, phiong =0, radg=0, ccbn[SIZE]={0}, ingnl, outgnl, norma, residual,heionmass;
	double tolerance = 0.00005, invnorm, bnlnorm,nlinvnorm,lsingest,nlmatnorm2,resone,reszero,a567g,bnladdd = 0;;
	int lten;
	double** rate;
	char caseab;
	double tman, nman;
	int texp, nexp;
	char folder_name[50];

	FILE *anm500Ptr;
	FILE *bnPtr;
	FILE *bnlPtr;
	FILE *bnlSHPtr;
	FILE *bnlsumPtr;
	FILE *bnlsumSHPtr;
	FILE *matPtr;
	FILE *vecPtr;
	FILE *coefmatPtr;
	FILE *nlvecPtr;
	FILE *nlempPtr;
	FILE *bnlcleanPtr;
	FILE *it3Ptr;
	FILE *it10Ptr;
	FILE *it20Ptr;
	FILE *it30Ptr;
	FILE *it40Ptr;
	FILE *it50Ptr;
	FILE *ratePtr;
	FILE *invnormPtr;
	FILE *empPtr;
	FILE *sparsePtr;
	FILE *rangeinfoPtr;
	FILE *aminPtr;
	FILE *aplusPtr;
	FILE *colionPtr;
	FILE *colionnlPtr;
	FILE *radrecPtr;
	FILE *threecomPtr;
	FILE *threecomnlPtr;


	temp = TE;
	depoptest = 1; //1 to have same depop as norad, 0 to add extra level when rad is present


	tic = clock();
	nlmax = round_up_5(350 - 15*log(NE));
	if(nlmax<50){
		nlmax=50;
	}


	printf("nlmax = %d\n",nlmax);

	nlcount = nlmax*(nlmax + 1)/2;

	depop = 0; /*start off with NOT assuming ground state(s) are depopulated*/
	if(NMIN == 2){ /*always depopulate level 1 for Case B*/
		depop = 1;
	}

	dw =0*pow(10,0);

	//CREATE OUTPUT FOLDER NAME

	if(NMIN==1){caseab='A';}
	if(NMIN==2){caseab='B';}

	nexp = floor(log10(fabs(NE)));
	nman = NE / pow(10, nexp);

	texp = floor(log10(fabs(TE)));
	tman = TE / pow(10, texp);

	while(nman - (int) nman != 0){
		nman = nman*10;
	}

	while(tman - (int) tman != 0){
		tman = tman*10;
	}


	snprintf(folder_name, sizeof(folder_name), "%d%d_%d%d_%c_", (int) nman, nexp, (int) tman, texp, caseab);
	if(FF != 0 && CMB == 0) {
		snprintf(folder_name + strlen(folder_name), sizeof(folder_name) - strlen(folder_name), "F/");
	}
	if(CMB != 0 && FF == 0){
		snprintf(folder_name + strlen(folder_name), sizeof(folder_name) - strlen(folder_name), "C/");
	}
	if(CMB != 0 && FF != 0 && DUST == 0 && dw == 0){
		snprintf(folder_name + strlen(folder_name), sizeof(folder_name) - strlen(folder_name), "C_F/");
	}
	if(CMB != 0 && FF != 0 && DUST != 0){
		snprintf(folder_name + strlen(folder_name), sizeof(folder_name) - strlen(folder_name), "C_F_D/");
	}
	if(FF == 0 && dw == 0 && CMB == 0 && DUST == 0){
		snprintf(folder_name + strlen(folder_name), sizeof(folder_name) - strlen(folder_name), "0_nc2300/");
	}
	if(FF != 0 && dw != 0 && CMB != 0 && DUST == 0){
		snprintf(folder_name + strlen(folder_name), sizeof(folder_name) - strlen(folder_name), "C_F_54_-1/");
	}
	printf(folder_name);
	printf("\n");

	mkdir(folder_name);

	/*zero matrix mat*/
	for(n=1; n<=NMAX; n++){
		for(m=1; m<=NMAX; m++){
			mat[n-1][m-1] = 0;

		}
	}
	if(depoptest == 1){
		depop = depop + 1;
	}else{
	if(dw == 0){/*depopulated ground state(s) if no radiation field. */
		depop = depop + 1;
	}
	}
	/* depop = 2;
	minmat =3; */

	np = NE*0.909; /*proton density in cm^-3*/
	elecmass = 9.10938356*pow(10, -28);
	protmass = 1.672621898*pow(10, -24);
	heionmass = 6.6446572*pow(10, -24);
	rydberg = rydinf / (1 + (elecmass / (ATOMN*protmass)));
	y = 1.43877735*ATOMZ*ATOMZ*rydberg/(temp); /*y=Z^2*R*h*c/(kT)*/
	printf("NE = %d\nNP = %d\nTE = %d\nTR = %d\nW = %e\nCMB  %d\nFF   %d\nDUST %d\n",NE,np,temp,TR,dw,CMB,FF,DUST);

	maxn = NMAX;
	/*printf("Enter 1 to read bn's from file: ");
	scanf("%d", &calcbn );*/
	calcbn=1;

/* 	n=4;
	gullc( n, thetau, thetal, 0.0025);
	for(i=0;i<n;i++){
		printf("%d -> %d  %e\n",i,i+1,thetau[i]);
		printf("%d -> %d  %e\n",i,i-1,thetal[i]);
	}  */


	if(calcbn==0){
		bnPtr = fopen("bnval.dat", "r");
		for ( n = 1; n <= NMAX; n++){
			fscanf(bnPtr, "%le", &b);
			bn[n-1] = b;
			//bn[n-1] = 1;

		}
	}else{

	printf("\nCalculating bn's\n");
	stepgen(step);


	lagintermat(step, lagmat);

/* 	n=4;
	radrecomb( n, y, thetau, thetal, radrec);
	printf("%d  %e\n",n,bfn(n, radrec));
	printf("y  %e\n",y);
	printf("RyH  %e\n",rydberg);
	for(i=0;i<n;i++){
		printf("%d   %e\n",i,radrec[i]);

	} */
	/* n=10;
	radrecomb( n, y, thetau, thetal, radrec);
	printf("%d  %e\n",n,bfn(n, radrec));
	n=50;
	radrecomb( n, y, thetau, thetal, radrec);
	printf("%d  %e\n",n,bfn(n, radrec)); */
	
	
	colionPtr = fopen(makepath(folder_name,"colion.dat"),"w"); 
	fprintf(colionPtr,"COLLISIONAL IONIZATION RATES \nNE = %d  NP = %d  TE = %d  NMIN = %d  NLMAX = %d \n\n",NE,np,temp,NMIN,nlmax);
	for ( n = 1; n <= 250; n++){
		colionratecoef(n, rydberg, temp, bfcolrate, y);
		fprintf(colionPtr,"%d  %e \n",n, bfcolrate[0]);
	}
	fclose(colionPtr);
	
	
	threecomPtr = fopen(makepath(folder_name,"3comb.dat"),"w"); 
	fprintf(threecomPtr,"THREE BODY RECOMBINATION RATES \nNE = %d  NP = %d  TE = %d  NMIN = %d  NLMAX = %d \n\n",NE,np,temp,NMIN,nlmax);
	for ( n = 1; n <= 250; n++){
		colionratecoef(n, rydberg, temp, bfcolrate, y);
		bfcolrate[ 1 ] = n*n*exp(y/n/n)*bfcolrate[ 0 ];
		fprintf(threecomPtr,"%d  %e \n",n, bfcolrate[1]);
	}
	fclose(threecomPtr);
	
	
	
	radrecPtr = fopen(makepath(folder_name,"radrec_n.dat"),"w"); 
	fprintf(radrecPtr,"RADIATIVE RECOMBINATION RATES \nNE = %d  NP = %d  TE = %d  NMIN = %d  NLMAX = %d \n\n",NE,np,temp,NMIN,nlmax);
	for ( n = 1; n < 250; n++){
		radrecomb( n, y, thetau, thetal, radrec);
		fprintf(radrecPtr,"%d  %e \n",n, bfn(n, radrec));
	}
			
	fclose(radrecPtr);
	
	

	/*CALCULATION OF ROW VECTOR*/
	 for ( n = NMIN; n <= SIZE; n++){
		
		colionratecoef(step[n], rydberg, temp, bfcolrate, y);
		bfcolrate[ 1 ] = step[n]*step[n]*exp(y/step[n]/step[n])*bfcolrate[ 0 ];

		stimrecomb( step[n], y, temp, thetau, thetal, stimrec, rydberg, dw );
		if (n==depop+1){
			inground = bfn(step[n], radrec) + bfn(step[n], stimrec) + NE*bfcolrate[1] ;
		}
		if(n < 0){
			inff = 0;
		}else{
			inff = infsumnm(step[n], y, rydberg, temp, dw);
		}
		if (n==NMIN){
			radg = bfn(step[n], radrec);
		}
		
		veca[n-1] = 2.41468184*pow(10, 15)*pow(temp,3.0/2)*( RAD*bfn(step[n], radrec) + bfn(step[n], stimrec)  ) + NE*bfcolrate[1] /*+ inff*/;

		


	}
	
	
/* 	 vecPtr = fopen("veca.dat","w");
	for ( n = NMIN; n <= NMAX; n++){

		radrecomb( n, y, thetau, thetal, radrec);
		colionratecoef(n, rydberg, temp, bfcolrate, y);
		bfcolrate[ 1 ] = n*n*exp(y/n/n)*bfcolrate[ 0 ];

		stimrecomb( n, y, temp, thetau, thetal, stimrec, rydberg, dw );
		if (n==depop+1){
			inground = bfn(n, radrec) + bfn(n, stimrec) + NE*bfcolrate[1] ;
		}
		if(n < 0){
			inff = 0;
		}else{
			inff = infsumnm(n, y, rydberg, temp, dw);
		}
		if (n==NMIN){
			radg = bfn(n, radrec);
		}

		fprintf(vecPtr,"%e", 2.41468184*pow(10, 15)*pow(temp,3.0/2)*( RAD*bfn(n, radrec) + bfn(n, stimrec)  ) + NE*bfcolrate[1]
				+ inff);
		if(n != NMAX){fprintf(vecPtr," ");  }


	}




	fclose(vecPtr);   */


	printf("after vec\n");
	anm500Ptr = fopen("anm500.dat","rb");

	/*POPULATE MATRIX*/
	for ( n = NMIN; n <= NMAX; n ++){/*rows*/
		for ( m = NMIN; m <= NMAX; m ++){/*columns*/
			if(n == m){
				colionratecoef(n, rydberg, temp, bfcolrate,y);
				photoion( n, thetau, thetal, phion, rydberg, dw );

				mat[n-1][n-1] = NE*bfcolrate[0] + bfn2(n, phion);

				if(n==NMIN){

					phiong = bfn2(n, phion);
					printf("phio %e %e\n",phiong,NE*bfcolrate[0]);}
			}
			if(n > m){
				/*below diagonal*/
				if(n <= 500){
					seek = (n-1)*(n-2)/2 + m;
					fseek( anm500Ptr, (seek-1)*sizeof(double), SEEK_SET);
					fread(&anm, sizeof(double), 1, anm500Ptr);
				}
				else{
					anm = avalaprx (n, m);
				}

				bbex = bbcolex(n, m, anm, rydberg, temp);
				freq3 = pow(freq(n, m, rydberg), 3);
				bstimnm = 6.78196256*pow(10, 46)*anm / freq3;


				mat[n-1][m-1] = (RAD*anm + NE*bbcoldeex(n, m, rydberg, temp, bbex) + jnu(n, m, rydberg, dw)*bstimnm);
				//printf("mat[%d][%d] = %e\n",n,m,mat[n-1][m-1]);
				

			}

			if(n < m){
				/*above diagonal*/
				if(m <= 500){
						seek = (m-1)*(m-2)/2 + n; /*posistion of A_nm in anm.dat*/
						fseek( anm500Ptr, (seek-1)*sizeof(double), SEEK_SET);
						fread(&anm, sizeof(double), 1, anm500Ptr);
				}
				else{
					anm = avalaprx (m, n);
				}

				bbex = bbcolex(m, n, anm, rydberg, temp);

				freq3 = pow(freq(m, n, rydberg), 3);
				babsnm = (float)m*m/n/n*6.78196256*pow(10, 46)*anm / freq3;

				mat[n-1][m-1] = (NE*bbex + jnu(m, n, rydberg, dw)*babsnm);
				//printf("mat[%d][%d] = %e\n",n,m,mat[n-1][m-1]);
				//printf("%d -> %d   %e\n",n,m,bbex);
				

			}

		}
	}

	fclose(anm500Ptr);

	

	/*DIAGONAL*/ /*total depop rate for each level*/
	for(n=NMIN; n<=NMAX; n++){
		for(m=NMIN; m<=NMAX; m++){
			if(n != m){
				mat[n-1][n-1] = mat[n-1][n-1] + mat[n-1][m-1];
			}
		}
	}


	outground =   mat[depop][depop];


	n=depop+1;
	for(m=n+1; m<=NMAX; m++){
		inground = inground + mat[m-1][n-1];
	}
	printf("ground = %d\n",depop+1);
	printf("+++++inground = %e\n",inground);
	printf("-----outground = %e\n\n",outground);
	printf("radg = %e\n",radg);
	 if(dw!=0){
		/* if(phiong==0 || (floor(log10(phiong)) < floor(log10(outground)) && floor(log10(phiong)) < floor(log10(NE*np*radg)))){ */
		if(log10(inground/outground) > 10){
			printf("login-10 %e\n",log10(inground)-10);
			mat[depop][depop]=pow(10,log10(inground)-10);
			printf("mat[depop][depop] %e\n",mat[depop][depop]);
		}

	}



	/*multiply with exponential factors*/
	for(n=1; n<=NMAX; n++){
		for(m=1; m<=NMAX; m++){
			if(m == n){
				mat[n-1][m-1] = m*m*exp(y/(m*m))*mat[n-1][m-1];
			}
			if(m < n){
				mat[n-1][m-1] = -n*n*exp(y/(n*n))*mat[n-1][m-1];
			}
			if(m > n){
				mat[n-1][m-1] = -n*n*exp(y/(n*n))*mat[n-1][m-1];
			}
		}
	}


	condmat( mat, lagmat, cmat, step);
	printf("after condmat\n");

/* 	 matinv(cmat, imat, depop);

	printf("after matinv\n");

	condno = matnorm(imat)*matnorm(cmat);
	printf("CONDITION NUMBER: %e\n",condno);

	/*scale departure coefficient of ground state if ground state is populated to reduce condition number*/
	/* if(log10(condno) > 10){
		scale = (veca[depop] - cmat[depop +1][depop ])/cmat[depop][depop];
		printf("scale  %e\n",scale);
		printf("scale = (%e - %e)/%e\n",veca[depop],cmat[depop+1][depop ],cmat[depop][depop]);
		for(n=depop+1; n<=SIZE; n++){
			cmat[depop][n-1] = scale*cmat[depop][n-1];
		}
		matinv(cmat, imat, depop);

	condno = matnorm(imat)*matnorm(cmat);
	printf("NEW CONDITION NUMBER: %e\n",condno);
	}
	printf("scale = %e\n",scale);  */

	//matmult(imat, veca, ccbn);
	//bn3Ptr = fopen("invbn.dat","w");
	//for ( n = depop+1; n <= SIZE; n++){
		/*printf("%e\n", bn[n-1]);  */
	//	fprintf(bn3Ptr,"%d %.10e\n", n,ccbn[n-1]);
	//}
	//fclose(bn3Ptr);


	printf("depop = %d\n",depop);

	gausselim(cmat,veca,cbn, depop);
	printf("after gausselim\n");
//	cbn[depop] =cbn[depop]*scale;


	/*EXPAND BN*/
	expvec(cbn, lagmat, bn);
	free(lagmat);




	printf("BN GROUND  %e\n",bn[NMIN-1+depop]);
	bnPtr = fopen(makepath(folder_name,"bnval.dat"),"w"); /* open for writing */
	printf("NE = %d  NP = %d  TE = %d  TR = %d  W = %e   NMIN = %d  TOL = %e\n\n",NE,np,temp,TR,dw,NMIN,tolerance);
	/*fprintf(bnPtr,"NE = %d   TE = %d  TR = %d  W = %e   NMIN = %d\n\n",NE,temp,TR,dw,NMIN);*/
	/*fprintf(bnPtr,"%d\n", 0);*/
	if(NMIN == 2){
			fprintf(bnPtr,"%d\n", 0);
	}
	for ( n = NMIN; n <= NMAX; n++){
		/*printf("%e\n", bn[n-1]);  */
		/* fprintf(bnPtr,"%d %.10e\n", n,bn[n-1]);   */
		fprintf(bnPtr,"%.8e\n", bn[n-1]);
	}
	fclose(bnPtr);

	/*iterativerefinement(cmat, veca, cbn, depop);*/

	printf("\a");
	printf("bn's calculated\n");

	}


	lten = findlten(bn);


/*COMMENT OUT FROM HERE TO GET JUST N-MODEL*/


	for ( n = 0; n < NMAX; n++){
		for ( m = 0; m < NMAX; m++){ /*l values*/
			mat[n][m]= 0;
		}
	}

	for ( n = NMIN; n <= NMAX; n++){
		for ( m = 0; m < n; m++){ /*l values*/
			mat[n][m]= bn[n-1] ;
			//mat[n][m]= 1 ;
			if(bn[n-1] <0){mat[n][m]=0;}
			/*mat[n][m]= 1 ;
			bn[n-1] = 1;*/
		}
	}


	rate = calloc(nlmax*(nlmax + 1)+1, sizeof(double *));
	for(n = 0; n < nlmax*(nlmax + 1)+1; n++) {
		rate[n] = calloc(nlmax, sizeof(double));
	}



	if(rate == NULL/*  || nlmat == NULL  */){
		printf("OUT OF MEMORY!\n");
		exit(1);
	}


	printf("rate allocated\n");

	a1 = calloc(nlmax*(nlmax + 1)/2, sizeof(double));
	emp1 = calloc(nlmax*(nlmax + 1)/2, sizeof(double));

	//if(dw==0&&NE==0){depop = depop+1;}
	anm500Ptr = fopen("anm500.dat","rb");
	aminPtr = fopen("avalmin300.dat","rb");
	aplusPtr = fopen("avalplus300.dat","rb");

	for ( n = NMIN; n <= nlmax; n++){

		for(i = 0; i < n; i++){
			anl[ i ] = 0;
		}

		for(m = NMIN; m < n; m++){
			for(i = 0; i <= m; i++){
				if(i > 0 && i <= m){
					if(n<=300){
						fseek(aminPtr, (amidx(n,m,i))*sizeof(double), SEEK_SET);
						fread(&alr, sizeof(double), 1, aminPtr);
					}else{
						alr = anlmll(n, m, i, rydberg);
					}
					anl[i] = anl[i] + alr;
				}
				if(i < m-1){
					if(n<=300){
						fseek(aplusPtr, (apidx(n,m,i))*sizeof(double), SEEK_SET);
						fread(&aur, sizeof(double), 1, aplusPtr);
					}else{
						aur = anlmlu(n, m, i, rydberg);
					}
					anl[i] = anl[i] + aur;
				}
			}
		}
		//colmomcoefbrock( n, cmomlp, cmomup, temp, elecmass, protmass, anl, nlmax, protmass, 1);
	//	colmomcoefbrock( n, cmomle, cmomue, temp, elecmass, protmass, anl, nlmax, elecmass, -1);
	//	colmomcoefbrock( n, cmomlh, cmomuh, temp, elecmass, protmass, anl, nlmax, 4*protmass+elecmass, 1);

		colmomcoef( n, cmomlp, cmomup, temp, elecmass, protmass, anl, protmass, 1);
		colmomcoef( n, cmomle, cmomue, temp, elecmass, protmass, anl, elecmass, -1);
		colmomcoef( n, cmomlh, cmomuh, temp, elecmass, protmass, anl, heionmass, 1);


		for(i = 0; i < n; i++){

			for(m = n + 1; m <= nlmax; m++){


				if(m != n){
					seek = (m-1)*(m-2)/2 + n; /*posistion of A_mn in anm.dat*/
					fseek( anm500Ptr, (seek-1)*sizeof(double), SEEK_SET);
					fread(&anm, sizeof(double), 1, anm500Ptr);
					_jnu = jnu(m, n, rydberg, dw);
					bbex = bbcolex(m, n, anm, rydberg, temp);
					bbdex = bbcoldeex(m, n, rydberg, temp, bbex);
					freq3 = pow(freq(m, n, rydberg), 3);

					if(i != 0 && i < n){ /* l' = l-1 */
						if(m<=300){
							fseek(aplusPtr, (apidx(m,n,i-1))*sizeof(double), SEEK_SET);
							fread(&einsteina, sizeof(double), 1, aplusPtr);
						}else{
							einsteina = anlmlu(m, n, i-1, rydberg);
						}
						stim = 6.78196256*pow(10, 46)*einsteina / freq3;
						//bbexl = einsteina/anm*bbex;
						bbexl = einsteina/anm*(float)(2*(i-1) + 1)/(2*i + 1)*(float)n*n/m/m*bbex;
						
						//bbdexl = (float)(2*i + 1)/(2*(i - 1) + 1)*exp(y/n/n-y/m/m)*bbexl;
						
						bbdexl = einsteina/anm*(float)n*n/m/m*exp(y/n/n-y/m/m)*bbex;
						/*m i-1 -> n i*/
						rate[n*(n-1) + 2*i][m-1] = RAD*einsteina + NE*bbdexl + stim*_jnu;/*a2*/

						/*n i-> m i-1*/
						rate[m*(m-1) +2*i -1][n-1] = NE*bbexl + absorp(n, i, m, i-1, rydberg, stim )*_jnu;/*a6*/

					}

					if(i <= n-1){/* l' = l+1 */
						if(m<=300){
							fseek(aminPtr, (amidx(m,n,i+1))*sizeof(double), SEEK_SET);
							fread(&einsteina, sizeof(double), 1, aminPtr);
						}else{
							einsteina = anlmll(m, n, i+1, rydberg);
						}
						stim = 6.78196256*pow(10, 46)*einsteina / freq3;
						//bbexl = einsteina/anm*bbex;
						bbexl = einsteina/anm*(float)(2*(i+1) + 1)/(2*i + 1)*(float)n*n/m/m*bbex;
						
					
						//bbdexl = (float)(2*i + 1)/(2*(i + 1) + 1)*exp(y/n/n-y/m/m)*bbexl;
						bbdexl = einsteina/anm*(float)n*n/m/m*exp(y/n/n-y/m/m)*bbex;
						
						/*m i+1-> n i*/
						rate[n*(n-1) + 2*i + 1][m-1] = RAD*einsteina + NE*bbdexl + stim*_jnu;/*a2*/

						/*n i-> m i+1*/
						rate[m*(m-1) + 2*i +2][n-1] = NE*bbexl + absorp(n, i, m, i+1, rydberg, stim )*_jnu;/*a6*/

					}
				}
			}

			for(m = NMIN; m < n; m++){

				if(m != n){
					seek = (n-1)*(n-2)/2 + m;
					fseek( anm500Ptr, (seek-1)*sizeof(double), SEEK_SET);
					fread(&anm, sizeof(double), 1, anm500Ptr);
					_jnu = jnu(n, m, rydberg, dw);
					bbex = bbcolex(n, m, anm, rydberg, temp);
					bbdex = bbcoldeex(n, m, rydberg, temp, bbex);
					freq3 = pow(freq(n, m, rydberg), 3);

					if( i != 0 && i <= m ){ /* l' = l-1 */
						if(n<=300){
							fseek(aminPtr, (amidx(n,m,i))*sizeof(double), SEEK_SET);
							fread(&einsteina, sizeof(double), 1, aminPtr);
						}else{
							einsteina = anlmll(n, m, i, rydberg);
						}
						stim = 6.78196256*pow(10, 46)*einsteina / freq3;
						//bbexl = einsteina/anm*bbex;
						bbexl = einsteina/anm*(float)(2*i+1)/(2*(i-1)+1)*(float)m*m/n/n*bbex;
						
						//bbdexl = (float)(2*(i - 1) + 1)/(2*i + 1)*exp(y/m/m-y/n/n)*bbexl;
						bbdexl = einsteina/anm*(float)m*m/n/n*exp(y/m/m-y/n/n)*bbex;
						
						/*m i-1-> n i*/
						rate[n*(n-1) + 2*i][m-1] = NE*bbexl + absorp(m, i-1, n, i, rydberg, stim )*_jnu;/*a3*/

						/*n i-> m i-1*/
						rate[m*(m-1) + 2*i -1][n-1] = RAD*einsteina + NE*bbdexl + stim*_jnu;/*a5*/

					}

					if(i <= m-2){/** l' = l+1 */
						if(n<=300){
							fseek(aplusPtr, (apidx(n,m,i))*sizeof(double), SEEK_SET);
							fread(&einsteina, sizeof(double), 1, aplusPtr);
						}else{
							einsteina = anlmlu(n, m, i, rydberg);
						}
						stim = 6.78196256*pow(10, 46)*einsteina / freq3;
						//bbexl = einsteina/anm*bbex;
						bbexl = einsteina/anm*(float)(2*i+1)/(2*(i+1)+1)*(float)m*m/n/n*bbex;
						//bbdexl = (float)(2*(i + 1) + 1)/(2*i + 1)*exp(y/m/m-y/n/n)*bbexl;
						
						bbdexl = einsteina/anm*(float)m*m/n/n*exp(y/m/m-y/n/n)*bbex;
						
						/*m i+1-> n i*/
						rate[n*(n-1) + 2*i + 1][m-1] = NE*bbexl + absorp(m, i+1, n, i, rydberg, stim )*_jnu;/*a3*/

						/*n i-> m i+1*/
						rate[m*(m-1) + 2*i +2 ][n-1] = RAD*einsteina + NE*bbdexl + stim*_jnu; /*a5*/

					}
				}
			}
			if(n!=1){
				if( i != 0 ){ /* l' = l-1 */
					/*n i-1 -> n i*/
					rate[n*(n-1) + 2*i][n-1] = np*cmomup[i-1] + NE*0.09*cmomuh[i-1] + NE*cmomue[i-1]; /*a4*/

					/*n i -> n i-1*/
					rate[n*(n-1) +2*i -1][n-1] = np*cmomlp[i]  + NE*0.09*cmomlh[i] +  NE*cmomle[i]; /*a7*/

				}

				if(i < n-1){/** l' = l+1 */
					/*n i+1 -> n i*/
					rate[n*(n-1) + 2*i + 1][n-1] = (np*cmomlp[i+1]  + NE*0.09*cmomlh[i+1] + NE*cmomle[i+1]);/*a4*/

					/*n i -> n i+1*/
					rate[n*(n-1) + 2*i +2][n-1] = np*cmomup[i] + NE*0.09*cmomuh[i] + NE*cmomue[i]; /*a7*/

				}
			}
		}
	}

	ratePtr = fopen(makepath(folder_name,"rate.dat"),"w");
	fprintf(ratePtr,"NE = %d  NP = %d  TE = %d  TR = %d  W = %e  NMIN = %d  NLMAX = %d \n\n",NE,np,temp,TR,dw,NMIN,nlmax);
	for ( n = 0; n < nlmax*(nlmax + 1)+1; n++){
		for(m = 0; m < nlmax; m++){
			fprintf(ratePtr,"%.8e ",rate[n][m]);
		}
		fprintf(ratePtr,"\n");
	}
	fclose(ratePtr);

	printf("matrix filled\n");

	ingnl=0;
	outgnl=0;
	
	colionnlPtr = fopen(makepath(folder_name,"colionnl.dat"),"w");
	threecomnlPtr = fopen(makepath(folder_name,"3bodynl.dat"),"w");	
	fprintf(colionnlPtr,"COLLISIONAL IONIZATION RATES \nNE = %d  NP = %d  TE = %d  NMIN = %d  NLMAX = %d \n\n",NE,np,temp,NMIN,nlmax);
	fprintf(threecomnlPtr,"THREE BODY RECOMBINATIN RATES \nNE = %d  NP = %d  TE = %d  NMIN = %d  NLMAX = %d \n\n",NE,np,temp,NMIN,nlmax);
	
	n=10;
	colionratecoef(n, rydberg, temp, bfcolrate, y);
	for ( i = 1; i < n; i++){
		colionnl =  bfcolrate[ 0 ];
		bfcolrate[1] = 4.14133233*pow(10, -16)*pow(temp,-3.0/2)*(2*i + 1)*exp(y/n/n)*colionnl;
		fprintf(colionnlPtr,"%d %d %e \n",n, i, colionnl);
		fprintf(threecomnlPtr,"%d %d %e \n",n, i, bfcolrate[1]);
	}
	fprintf(colionnlPtr,"\n");
	fprintf(threecomnlPtr,"\n");
	
	n=40;
	colionratecoef(n, rydberg, temp, bfcolrate, y);
	for ( i = 1; i < n; i++){
		colionnl =  bfcolrate[ 0 ];
		bfcolrate[1] = 4.14133233*pow(10, -16)*pow(temp,-3.0/2)*(2*i + 1)*exp(y/n/n)*colionnl;
		fprintf(colionnlPtr,"%d %d %e \n",n, i, colionnl);
		fprintf(threecomnlPtr,"%d %d %e \n",n, i, bfcolrate[1]);
	}
	fprintf(colionnlPtr,"\n");
	fprintf(threecomnlPtr,"\n");
	
	n=100;
	colionratecoef(n, rydberg, temp, bfcolrate, y);
	for ( i = 1; i < n; i++){
		colionnl =  bfcolrate[ 0 ];
		bfcolrate[1] = 4.14133233*pow(10, -16)*pow(temp,-3.0/2)*(2*i + 1)*exp(y/n/n)*colionnl;
		fprintf(colionnlPtr,"%d %d %e \n",n, i, colionnl);
		fprintf(threecomnlPtr,"%d %d %e \n",n, i, bfcolrate[1]);
	}
	fprintf(colionnlPtr,"\n");
	fprintf(threecomnlPtr,"\n");
	
	n=200;
	colionratecoef(n, rydberg, temp, bfcolrate, y);
	for ( i = 1; i < n; i++){
		colionnl =  bfcolrate[ 0 ];
		bfcolrate[1] = 4.14133233*pow(10, -16)*pow(temp,-3.0/2)*(2*i + 1)*exp(y/n/n)*colionnl;
		fprintf(colionnlPtr,"%d %d %e \n",n, i, colionnl);
		fprintf(threecomnlPtr,"%d %d %e \n",n, i, bfcolrate[1]);
	}
	fprintf(colionnlPtr,"\n");
	fprintf(threecomnlPtr,"\n");
	
	
	fclose(colionnlPtr);
	fclose(threecomnlPtr);


	for ( n = nlmax; n >= NMIN; n--){
		tic1 = clock();
		radrecomb( n, y, thetau, thetal, radrec);
		stimrecomb( n, y, temp, thetau, thetal, stimrec, rydberg, dw );
		photoion( n, thetau, thetal, phion, rydberg, dw );
		colionratecoef(n, rydberg, temp, bfcolrate, y);
		for(i = 0; i < n; i++){
			colionnl =  bfcolrate[ 0 ];
			bfcolrate[1] = 4.14133233*pow(10, -16)*pow(temp,-3.0/2)*(2*i + 1)*exp(y/n/n)*colionnl;

			if(n < 0){
				inffin = 0;
				inffout = 0;
			}else{
				inffin = infsumnlfacc(n, i ,y, rydberg, temp, dw, nlmax, bn);
				inffout = infsumnleacc(n, i,y, rydberg, temp, dw, nlmax, bn);
				if(isinf(inffout)){
					printf("inffout inf \n n = %d \n i = %d",n,i);
					return 0;
				}
					 
			}
			
			

			a1[n*(n-1)/2 + i] = 2.41468184*pow(10, 15)*pow(temp,3.0/2)*( RAD*radrec[i] + stimrec[i] + NE*bfcolrate[1])
				+ uppernsumin(n, i ,y, rydberg, temp, dw, nlmax, bn) + inffin;
			upsumout = uppernsumout(n, i,y, rydberg, temp, dw, nlmax, bn);
			if(isinf(upsumout)){
				printf("upsumout inf \n n = %d \n i = %d",n,i);
				return 0;
			}
			if(isinf(colionnl)){
				printf("colionnl inf \n n = %d \n i = %d",n,i);
				return 0;
			}
			
			emp1[n*(n-1)/2 + i] = phion[i] + NE*colionnl + upsumout + inffout;
			if(isinf(emp1[n*(n-1)/2 + i])){
				printf("emp1 inf \n n = %d \n i = %d",n,i);
				return 0;
			}
			
			
			if(n==NMIN){
				ingnl= ingnl+RAD*radrec[i] + stimrec[i] + NE*bfcolrate[i];
				outgnl=outgnl+phion[i] + NE*bfcolrate[i];

			}

		}
		toc1 = clock();
		printf("%d   %f min\n",n,(double)(toc1 - tic1) / CLOCKS_PER_SEC /60);
	}


	minmat=depop;
	if(depoptest == 1){
		if(NMIN==2){
			minmat = 3;
		}
	}else{
		if(NMIN==2 && dw==0){
			minmat = 3;
		}
	}
	printf("minmat = %d\n",minmat);
	/* depop = 2;
	minmat =3; */

	// NL SPARSE COEFF MATRIX BEGIN

	sparsePtr = fopen(makepath(folder_name,"sparsemat.dat"),"w");
	fprintf(sparsePtr,"NE = %d  NP = %d  TE = %d  TR = %d  W = %e  NMIN = %d  NLMAX = %d \n\n",NE,np,temp,TR,dw,NMIN,nlmax);
	for ( n = nlmax; n >= NMIN; n--){
		for(i = 0; i < n; i++){

			a2 = 0;
			a3 = 0;
			a5 = 0;
			a4 = 0;
			a6 = 0;
			a7 = 0;

			for(m = n + 1; m <= nlmax; m++){
				if(m != n){

					if(i != 0 && i < n){ /* l' = l-1 */
						/*m i-1 -> n i*/
						//nlmat[idx(n,i)][idx(m,i-1)] =  -(2*i - 1)*exp(y/m/m)*(rate[n*(n-1) + 2*i][m-1]);
						if(n >= depop+1){
							fprintf(sparsePtr,"%d %d %.8e\n",idx2(n,i,minmat),idx2(m,i-1,minmat),-(2*i - 1)*exp(y/m/m)*(rate[n*(n-1) + 2*i][m-1]));
						}

						/*n i -> m i-1*/
						a6 = a6 + rate[m*(m-1) +2*i -1][n-1];

					}

					if(i <= n-1){/* l' = l+1 */
						/*m i+1 -> n i*/
						//nlmat[idx(n,i)][idx(m,i+1)] = -(2*i + 3)*exp(y/m/m)*(rate[n*(n-1) + 2*i + 1][m-1]);
						if(n >= depop+1){
							fprintf(sparsePtr,"%d %d %.8e\n",idx2(n,i,minmat),idx2(m,i+1,minmat),-(2*i + 3)*exp(y/m/m)*(rate[n*(n-1) + 2*i + 1][m-1]));
						}

						/*n i -> m i+1*/
						a6 = a6 + rate[m*(m-1) + 2*i +2][n-1];

					}
				}
			}


			for(m = NMIN; m < n; m++){
				if(m != n){

					if( i != 0 && i <= m ){ /* l' = l-1 */

						/*m i-1 -> n i*/
						//nlmat[idx(n,i)][idx(m,i-1)] = - (2*i - 1)*exp(y/m/m)*(rate[n*(n-1) + 2*i][m-1]);
						if(m>=depop+1){
							fprintf(sparsePtr,"%d %d %.8e\n",idx2(n,i,minmat),idx2(m,i-1,minmat),- (2*i - 1)*exp(y/m/m)*(rate[n*(n-1) + 2*i][m-1]));
						}

						/*n i-> m i-1*/
						a5 = a5 + rate[m*(m-1) + 2*i -1][n-1];

					}

					if(i <= m-2){/* l' = l+1 */
						if(m >= depop+1){
							/*m i+1 -> n i*/
							//nlmat[idx(n,i)][idx(m,i+1)] = -(2*i + 3)*exp(y/m/m)*(rate[n*(n-1) + 2*i + 1][m-1]);
							fprintf(sparsePtr,"%d %d %.8e\n",idx2(n,i,minmat),idx2(m,i+1,minmat),-(2*i + 3)*exp(y/m/m)*(rate[n*(n-1) + 2*i + 1][m-1]));
						}
						/*n i -> m i+1*/
						a5 = a5 + rate[m*(m-1) + 2*i +2 ][n-1];
					}
				}
			}

			if( i != 0 ){ /* l' = l-1 */
				/*n i-1-> n i*/
				//nlmat[idx(n,i)][idx(n,i-1)] = -(2*i - 1)*exp(y/n/n)*(rate[n*(n-1) + 2*i][n-1]);
				if(n>=depop+1){
					fprintf(sparsePtr,"%d %d %.8e\n",idx2(n,i,minmat),idx2(n,i-1,minmat),-(2*i - 1)*exp(y/n/n)*(rate[n*(n-1) + 2*i][n-1]));
				}
				/*n i -> n i-1*/
				a7 = a7 + rate[n*(n-1) +2*i -1][n-1];

			}

			if(i < n-1){/** l' = l+1 */
				/*n i+1 -> n i*/
				//nlmat[idx(n,i)][idx(n,i+1)] = -(2*i + 3)*exp(y/n/n)*(rate[n*(n-1) + 2*i + 1][n-1]);
				if(n>=depop+1){
					fprintf(sparsePtr,"%d %d %.8e\n",idx2(n,i,minmat),idx2(n,i+1,minmat),-(2*i + 3)*exp(y/n/n)*(rate[n*(n-1) + 2*i + 1][n-1]));
				}
				/*n i -> n i+1*/
				a7 = a7 +  rate[n*(n-1) + 2*i +2][n-1];



			}

			if(depoptest != 1){
			if(NMIN == 1 && dw !=0 && n == 1){

				//nlmat[idx(1,0)][idx(2,0)] = -RAD*exp(y/4)*8.2206;
				printf("blarg\n");
				fprintf(sparsePtr,"%d %d %.8e\n",idx2(1,0,minmat),idx2(2,0,minmat),-RAD*exp(y/4)*8.2206);
			}
			}
			if(depoptest != 1){
			if(dw !=0 && n == 2 && i == 0){
				a5 = a5 + RAD*8.2206;

			}}
			if(n>=depop+1){
				fprintf(sparsePtr,"%d %d %.8e\n",idx2(n,i,minmat),idx2(n,i,minmat),(exp(y/n/n)*(2*i + 1)*(a5 + a6 + a7 + emp1[n*(n-1)/2 + i])));
				//nlmat[idx(n,i)][idx(n,i)] = (exp(y/n/n)*(2*i + 1)*(a5 + a6 + a7 + emp1[n*(n-1)/2 + i]));
			}


		}
	}
	fclose(sparsePtr);

	// NL SPARSE COEFF MATRIX END

	for( n = 0; n < nlmax; n++){
		free(rate[n]);
	}
	free(rate);

	nlvecPtr = fopen(makepath(folder_name,"nlvec.dat"),"w"); /* open for writing */
	fprintf(nlvecPtr,"NE = %d  NP = %d  TE = %d  TR = %d  W = %e  NMIN = %d  NLMAX = %d \n\n",NE,np,temp,TR,dw,NMIN,nlmax);
	for ( n = depop+1; n <= nlmax; n++){
		for ( i = 0; i < n; i++){

			fprintf(nlvecPtr,"%.8e ", a1[n*(n-1)/2 + i]);
		}
	}

	fclose(nlvecPtr);

/*	for(j=1;j<=550;j++){
	/*printf("\nj = %d\n", j);*/
	printf("a1/emp1 = %e\n",a1[depop]/emp1[depop]);

	/*if(a1[depop]/emp1[depop] > pow(10,20)){
		depop = depop + 1;
	}*/
	/*depop = 1;*/
	printf("********depop nl = %d\n",depop);
	if(depop > 2){printf("********** DEPOP > 2\n");depop =2;}

	nlimit = depop +1;
	conv = 0;
	itconv = 0;
	j = 0;



	empPtr = fopen(makepath(folder_name,"emp.dat"),"w"); /* open for writing */
	fprintf(empPtr,"NE = %d  NP = %d  TE = %d  TR = %d  W = %e  NMIN = %d  NLMAX = %d \n\n",NE,np,temp,TR,dw,NMIN,nlmax);
	for ( n = depop+1; n <= nlmax; n++){
		for ( i = 0; i < n; i++){
			fprintf(empPtr,"%.8e ", emp1[n*(n-1)/2 + i]);
		}
	}

	fclose(empPtr);


/* 	if(NMIN==2 && dw==0){minmat = 3;}
	coefmatPtr = fopen("coefmat.dat","w");
	for ( n = minmat; n < nlmax*(nlmax + 1)/2; n++){
		for(m = minmat; m < nlmax*(nlmax + 1)/2; m++){
			fprintf(coefmatPtr,"%e   ",nlmat[n][m]);
		}
		fprintf(coefmatPtr,"\n");
	}
	fclose(coefmatPtr); */


	fclose(anm500Ptr);
	fclose(aplusPtr);
	fclose(aminPtr);


	toc = clock();
	printf("Elapsed: %f minutes\n", (double)(toc - tic) / CLOCKS_PER_SEC /60);
	printf("Elapsed: %f hours\n", (double)(toc - tic) / CLOCKS_PER_SEC /60/60);


	rangeinfoPtr = fopen(makepath(folder_name,"rangeinfo.dat"),"w"); /* open for writing */
	fprintf(rangeinfoPtr,"%d\n",minmat);
	fprintf(rangeinfoPtr,"%d\n",depop+1);
	fprintf(rangeinfoPtr,"%d\n",nlmax);
	fprintf(rangeinfoPtr,"%f\n",(double)(toc - tic) / CLOCKS_PER_SEC /60/60);

	fprintf(rangeinfoPtr,"%d\n",NE);
	fprintf(rangeinfoPtr,"%d\n",np);
	fprintf(rangeinfoPtr,"%d\n",temp);
	fprintf(rangeinfoPtr,"%d\n",TR);
	fprintf(rangeinfoPtr,"%e\n",dw);
	fprintf(rangeinfoPtr,"%d\n",NMIN);

	fclose(rangeinfoPtr);
	printf(folder_name);
	printf("\n");

	printf("\a");


	return 0;
}



//finds n after which bn=1
int findlten(double bn[])
{
	int stop, lten;

	lten=NMAX+1;
	stop=1;
	if(nround(bn[NMAX-1],4)>1.0001 || nround(bn[NMAX-1],4)<0.99999){
		lten = NMAX;
	}else{
		while(stop!=0){
			lten=lten-1;
			if(nround(bn[lten-1],4)<1){
				stop=0;
			}
		}
	}
	if(lten<=NMIN){
		stop=1;
		while(stop!=0){
			lten=lten+1;
			if(nround(bn[lten-1],4)<=1){
				stop=0;
			}
		}
	}

	return lten;
}

//rounds num to j decimal digits
double nround(double num, int j)
{
	double a,b,c;

    a = pow(10, j);
    b = num*a;
    c = round(b)/a;

    return c;
}

//gives location of A(m,i -> n,i+1) in the file avalplus300.dat generated by makealmlupper.nb
int apidx(const int m, const int n, const int l)
{
	return (m*m*m - 6*m*m + 11*m - 6)/6 + (2*l*m - l*l - 3*l)/2 + m - n - 1;
}

//gives location of A(m,i -> n,i-1) in the file avalmin300.dat generated by makealmllower.nb
int amidx(const int m, const int n, const int l)
{
	return m*(m - 1)*(m - 2)/6 + (2*l*m - l*l - l)/2 + n - m;
}

/*calculates intensity of free-free emission at transition frequency from level n to level m*/
double freefree(const int n, const int m, const double ryd)
{
	double freq3,f, gauntff, kc, tc, bte, jc, em, iff;

	f=freq(n, m, ryd);
	freq3 = pow(f, 3);

	em = (float)NE*NE;
	//em = 3.4*pow(10,58)*pow(f*pow(10,-9),2.0/3.0);
	//optical depth from 2003MNRAS.341..369D
	tc = -0.03014*pow(TE,-1.5)*pow(f/1000000000,-2)*(log(0.04955/f) + 1.5*log(TE))*em;
	bte = (1.47449944*pow(10,-47)*freq3)/(exp(4.79924466*pow(10,-11)*f/TE)-1);
	iff=(1 - exp(-tc))*bte;

	if(iff<0){iff=0;}
	return iff;
}

/*calculates intensity of cmbr at transition frequency from level n to level m*/
double cmbr(const int n, const int m, const double ryd)
{
	double f;
	f=freq(n, m, ryd);
	return 1.47449944*pow(10, -47)*f*f*f*1.0/(exp(4.79924466*pow(10, -11)*f/2.7) - 1);
}

/*calculates intensity of dust radiation at transition frequency from level n to level m*/
double dust(const int n, const int m, const double ryd)
{
	double td,beta,f0,f,t,dusti;
	//Based on Planck 2013 results.XI.All-sky model of thermal dust emission, A&A
	td=100;//20.3;
	f=freq(n, m, ryd);
	//Optical depth from Draine book p121
	t=9.23626*pow(10,-25)*pow(f,1.7);

	dusti = t*1.47449944*pow(10, -47)*f*f*f*1.0/(exp(4.79924466*pow(10, -11)*f/td) - 1);
	if(dusti<0){dusti=0;}
	return dusti;
}

/*calculates intensity of black body radiation at transition frequency from level n to level m*/
double jnu(const int n, const int m, const double ryd, const double dw)
{
	double freq3, brem, cmbri, dusti, star,consti;

	freq3 = pow(freq(n, m, ryd), 3);

	star = 0;
	if(dw!=0){
		star = dw*1.47449944*pow(10, -47)*freq3*1.0/(exp(4.79924466*pow(10, -11)*freq(n, m, ryd)/TR) - 1);
	}

	brem = 0;
	if(FF!=0){
		brem = freefree(n, m, ryd);
	}

	cmbri = 0;
	if(CMB!=0){
		cmbri = cmbr(n, m, ryd);
	}

	dusti = 0;
	if(DUST!=0){
		dusti = dust(n, m, ryd);
	}

	consti=0;
	/* if(CONSTRADF !=0){
		consti = 1.1*pow(10,-12);
	} */

	return brem + cmbri + dusti + star + consti;
}

/*calculates specific intensity of free-free emission for bound-free transitions*/
double freefreebf(const int n, const double r, const double ryd)
{
	double f, tc, bte, iff;

	f=29979245800*ryd*ATOMZ*ATOMZ*(1.0/(n*n) + r);

	//optical depth from 2003MNRAS.341..369D
	tc = -0.03014*pow(TE,-1.5)*pow(f/1000000000,-2)*(log(0.04955/f) + 1.5*log(TE))*NE*NE;
	bte = (1.47449944*pow(10,-47)*f*f*f)/(exp(4.79924466*pow(10,-11)*f/TE)-1);

	iff=(1 - exp(-tc))*bte;

	if(iff<0){iff=0;}
	return iff;
}

/*calculates (4pi*I_nu)/(h*nu) of cmbr for bound-free transitions*/
double cmbrbf(const int n, const double r, const double ryd)
{
	double f;
	f = 29979245800*ryd*ATOMZ*ATOMZ*(1.0/(n*n) + r);

	return 1.47449944*pow(10, -47)*f*f*f*1.0/(exp(4.79924466*pow(10, -11)*f/2.7) - 1);
}

/*calculates intensity of dust radiation at transition frequency from level n to level m*/
double dustbf(const int n, const double r, const double ryd)
{
	double td,beta,f0,f,t,dusti;
	//Based on Planck 2013 results.XI.All-sky model of thermal dust emission, A&A
	td=100;//20.3;
	f = 29979245800*ryd*ATOMZ*ATOMZ*(1.0/(n*n) + r);
	//Optical depth from Draine book p121
	t=9.23626*pow(10,-25)*pow(f,1.7);

	dusti = t*1.47449944*pow(10, -47)*f*f*f*1.0/(exp(4.79924466*pow(10, -11)*f/td) - 1);
	if(dusti<0){dusti=0;}
	return dusti;
}

/*calculates (4pi*I_nu)/(h*nu) for bound-free transisions*/
double mintens(const int n, const double r, const double ryd, const double dw)
{
	double freq, brem, cmbri, dusti, star, consti;

	freq = 29979245800*ryd*ATOMZ*ATOMZ*(1.0/(n*n) + r);

	star = 0;
	if(dw!=0){
		star = dw*1.47449944*pow(10, -47)*freq*freq*freq*1.0/(exp(4.79924466*pow(10, -11)*freq/TR) - 1);
	}

	brem = 0;
	if(FF!=0){
		brem = freefreebf(n, r, ryd);
	}

	cmbri = 0;
	if(CMB!=0){
		cmbri = cmbrbf(n, r, ryd);
	}

	dusti = 0;
	if(DUST!=0){
		dusti = dustbf(n, r, ryd);
	}

	consti=0;
	/* if(CONSTRADF !=0){
		consti = 1.1*pow(10,-12);
	} */

	return 1.896504344*pow(10,27)/freq*(brem + cmbri + dusti + star + consti);
}


/*matrix condensation. Multiplies full matrix mat with
lagrangian interpolation matrix*/
void condmat(  double mat[][NMAX+1],  double lmat[][NMAX], double cmat[][SIZE], int step[])
{
	int i, j, k;


	for ( i = 1; i <= SIZE; i++){
		for ( j = 1; j <= SIZE; j++){
			for ( k = 1; k <= NMAX; k++){
				cmat[i-1][j-1] = cmat[i-1][j-1] + lmat[i-1][k-1]*mat[k-1][step[j]-1];
			}
		}
	}

}



/*returns minimum of two integer values*/
int fminn(const int i, const int j)
{
	int min;

	min = i;
	if(min > j){
		min = j;
	}

	return min;
}

/*returns maximum of two integer values*/
int fmaxx(const int i, const int j)
{
	int max;

	max = i;
	if(max < j){
		max = j;
	}

	return max;
}


/*expands condensed bn vector by lagrange interpolation*/
void expvec(const double cbn[], double lmat[][NMAX], double bn[])
{
	int i, j, k;

	for(j=0; j < NMAX; j++){
		bn[j] = 0;
	}

	for(j=0; j < NMAX; j++){
		for(k=0; k < SIZE; k++){
			bn[j] = bn[j] + cbn[k]*lmat[k][j];
		}
	}
}


/*calculates pivot points for interpolation*/
void stepgen(int step[])
{
	int i = 0, m, n, l = 0;

	step[0] = 0;

	/* for(n = 1; n <= SIZE; n ++){
		step[n] = n;
	}  */

	for(n = 1; n <= 9; n += 2){
		for(m = 1; m <= 40; m ++){
			i = i+1;
			l = l + n;
			step[i] = l;
		}
		l = l-1;
	}
}





/*calculates lagrange interpolating matrix used for matrix condensation*/
void lagintermat(const int step[], double lmat[][NMAX])
{
	int t, j, l, i, t1, t2,n, m, min, max;

	for(j = 0; j < SIZE; j++){
		for(i = 0; i < SIZE; i++){
			lmat[i][j] = 0;
		}
	}

	lmat[0][0] = 1;

	for(l = 2; l <= SIZE; l++){
		j = floor((step[l] - step[l-1])/2);
		t1 = fminn(SIZE, l+j);
		t2 = fmaxx(1, l-j);
		for(n = t2; n <= t1; n++){
			min = fminn(NMAX, step[l] + j);
			max = fmaxx(1, step[l] - j);
			for(m = max; m <= min; m++){
				lmat[n-1][m-1] = 1;
				for(i = t2; i <= t1; i++){
					if(i != n){
						lmat[n-1][m-1] = lmat[n-1][m-1]*(float)(m - step[i])/(step[n] - step[i]);

					}
				}
			}
		}
	}
}


/*calculates Sum_(t=nc)^\infty t^2*exp(\Chi_t/kT)*(A_tn + NE*C_tn) by converting
sum to integral using trapizoidal rule and the evaluating integral using
20 point Gaussian integration from NMAX to NC*/
double infsumnm(const int n, const double y, const double ryd, const int te, const double dw)
{
	double sum = 0, p, q, r, a1, a2, c1, c2, f1, f2, f, colex;
	int j, t, s;
	double bs1, ba1, bs2, ba2, intens, freq3;
	double x[20] = {-0.9931285992, -0.9639719273, -0.9122344283, -0.8391169718, -0.7463319065,
					-0.6360536807, -0.5108670020, -0.3737060887, -0.2277858511, -0.0765265211,
					0.0765265211, 0.2277858511, 0.3737060887, 0.5108670020, 0.6360536807,
					0.7463319065, 0.8391169718, 0.9122344283, 0.9639719273, 0.9931285992};
	double w[20] = {0.0176140071, 0.0406014298, 0.0626720483, 0.0832767416, 0.1019301198,
					0.1181945320, 0.1316886384, 0.1420961093, 0.1491729865, 0.1527533871,
					0.1527533871, 0.1491729865, 0.1420961093, 0.1316886384, 0.1181945320,
					0.1019301198, 0.0832767416, 0.0626720483,0.0406014298, 0.0176140071};
	FILE *anm500Ptr;


	anm500Ptr = fopen("anm500.dat","rb");
	p = (NC+1 - NMAX)/2;
	q = (NC+1 + NMAX)/2;

	for(j = 0; j <= 19; j++){
		r = p*x[j] + q;
		t = floor(r);

		if(t != n){
			if(t <= 500){
				s = (t-1)*(t-2)/2 + n;
				fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
				fread(&a1, sizeof(double), 1, anm500Ptr);

			}
			else{
				a1 = avalaprx (t, n);
			}
			colex = bbcolex(t, n, a1, ryd, te);
			c1 = bbcoldeex(t, n, ryd, te, colex);

			freq3 = pow(freq(t, n, ryd), 3);
			intens = jnu(t, n, ryd, dw);
			bs1 = 6.78196256*pow(10, 46)*a1 / freq3;
			ba1 = (float)t*t/n/n*6.78196256*pow(10, 46)*a1 / freq3;


			f1 = (float)t*t*exp(y/t/t)*(RAD*a1 + NE*c1 + intens*bs1) - (float)n*n*exp(y/n/n)*(NE*colex + intens*ba1);

		}

		if(t+1 != n){
			if(t+1 <= 500){
				s = (t)*(t-1)/2 + n; /*posistion of A_(t+1)n in anm.dat*/
				fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
				fread(&a2, sizeof(double), 1, anm500Ptr);
			}
			else{
				a2 = avalaprx (t+1, n);
			}

			colex = bbcolex(t+1, n, a2, ryd, te);
			c2 = bbcoldeex(t+1, n, ryd, te, colex);

			freq3 = pow(freq(t+1, n, ryd), 3);
			intens = jnu(t+1, n, ryd, dw);
			bs2 = 6.78196256*pow(10, 46)*a2 / freq3;

			ba2 = (float)(t+1)*(t+1)/n/n*bs2;

			f2 = (float)(t + 1)*(t + 1)*exp(y/(float)(t + 1)/(float)(t + 1))*(RAD*a2 + NE*c2 + intens*bs2) - (float)n*n*exp(y/(float)n/(float)n)*(NE*colex + intens*ba2);

		}

		if(t == n){
			f = f2;
		}
		else{
			if(t+1 == n){
				f = f1;
			}
			else{
				f = (f2 - f1)*(r - t) + f1; /*linear interpolation*/
			}
		}

		sum = sum+p*w[j]*f;

	}
	fclose(anm500Ptr);

	return sum;
}


/*calculates Sum_(t=nc)^\infty t^2*exp(\Chi_t/kT)*(A_tn + NE*C_tn) by converting
sum to integral using trapizoidal rule and the evaluating integral using
20 point Gaussian integration from NMAX to NC*/
double infsumnlf(const int n, const int l, const double y, const double ryd, const int te, const double dw, const int _min,  double bval[])
{
	double sum = 0, p, q, r, a1, a2, cdex, f1, f2, f, colex;
	double al1, au1, colexl1, colexu1, colexl2, colexu2, bal1, bal2, bau1, bau2;
	double bsu1, bsu2, bsl1, bsl2, cdexl1, cdexl2, cdexu1, cdexu2, al2, au2, b;
	int j, t, s;
	double bs1, ba1, bs2, ba2, intens, freq3;
	double x[20] = {-0.9931285992, -0.9639719273, -0.9122344283, -0.8391169718, -0.7463319065,
					-0.6360536807, -0.5108670020, -0.3737060887, -0.2277858511, -0.0765265211,
					0.0765265211, 0.2277858511, 0.3737060887, 0.5108670020, 0.6360536807,
					0.7463319065, 0.8391169718, 0.9122344283, 0.9639719273, 0.9931285992};
	double w[20] = {0.0176140071, 0.0406014298, 0.0626720483, 0.0832767416, 0.1019301198,
					0.1181945320, 0.1316886384, 0.1420961093, 0.1491729865, 0.1527533871,
					0.1527533871, 0.1491729865, 0.1420961093, 0.1316886384, 0.1181945320,
					0.1019301198, 0.0832767416, 0.0626720483,0.0406014298, 0.0176140071};
	FILE *anm500Ptr;


	anm500Ptr = fopen("anm500.dat","rb");
	p = (NC - _min)/2;
	q = (NC + _min)/2;

	for(j = 0; j <= 19; j++){
		r = p*x[j] + q;
		t = floor(r);

		if(t != n){
			if(t <= 500){
				s = (t-1)*(t-2)/2 + n;
				fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
				fread(&a1, sizeof(double), 1, anm500Ptr);
			}
			else{
				a1 = avalaprx (t, n);
			}

			colex = bbcolex(t, n, a1, ryd, te);
			cdex = bbcoldeex(t, n, ryd, te, colex);

			freq3 = pow(freq(t, n, ryd), 3);
			intens = jnu(t, n, ryd, dw);

			au1 = 0;
			al1 = 0;
			bsl1 = 0;
			bsu1 = 0;
			bal1 = 0;
			bau1 = 0;
			colexl1 = 0;
			colexu1 = 0;
			cdexl1 = 0;
			cdexu1 = 0;

			if(l != 0 && l < n){ /* l' = l-1 */
				/*t l-1 <-> n l*/
				au1 = anlmlu(t, n, l-1, ryd);
				bsu1 = 6.78196256*pow(10, 46)*au1 / freq3;
				bal1 = absorp(n, l, t, l-1, ryd, bau1 );
				colexl1 = au1/a1*colex;
				cdexu1 = au1/a1*cdex;
			}

			if(l <= n-1){/* l' = l+1 */
				/*t l+1 <-> n l*/
				al1 = anlmll(t, n, l+1, ryd);
				bsl1 = 6.78196256*pow(10, 46)*al1 / freq3;
				bau1 = absorp(n, l, n, l+1, ryd, bsl1 );
				colexu1 = al1/a1*colex;
				cdexl1 = al1/a1*cdex;
			}

			if(t<NMAX){
				b=bval[t-1];}
			else{
				b=1.0;}

			f1 = exp(y/(float)t/(float)t)*(b*(2*l + 3)*(al1 + NE*cdexl1 +bsl1*intens) + (float)b*(2*l - 1)*(au1 + NE*cdexu1 + bsu1*intens));
		}


		if(t+1 != n){
			if(t+1 <= 500){
				s = (t)*(t-1)/2 + n; /*posistion of A_(t+1)n in anm.dat*/
				fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
				fread(&a2, sizeof(double), 1, anm500Ptr);
			}
			else{
				a2 = avalaprx (t+1, n);
			}

			colex = bbcolex(t+1, n, a1, ryd, te);
			cdex = bbcoldeex(t+1, n, ryd, te, colex);

			freq3 = pow(freq(t+1, n, ryd), 3);
			intens = jnu(t+1, n, ryd, dw);

			au2 = 0;
			al2 = 0;
			bsl2 = 0;
			bsu2 = 0;
			bal2 = 0;
			bau2 = 0;
			colexl2 = 0;
			colexu2 = 0;
			cdexl2 = 0;
			cdexu2 = 0;

			if(l != 0 && l < n){ /* l' = l-1 */
				au2 = anlmlu(t+1, n, l-1, ryd);
				bsu2 = 6.78196256*pow(10, 46)*au2 / freq3;
				bal2 = absorp(n, l, t+1, l-1, ryd, bau2 );
				colexl2 = au2/a2*colex;
				cdexu2 = au2/a2*cdex;
			}

			if(l <= n-1){/* l' = l+1 */
				al2 = anlmll(t+1, n, l+1, ryd);
				bsl2 = 6.78196256*pow(10, 46)*al2 / freq3;
				bau2 = absorp(n, l, t+1, l+1, ryd, bsl2 );
				colexu2 = al2/a2*colex;
				cdexl2 = al2/a2*cdex;
			}

			if(t<NMAX){
				b=bval[t];}
			else{
				b=1.0;}


			f2 = exp(y/(float)(t+1)/(float)(t+1))*(b*(2*l + 3)*(al2 + NE*cdexl2 + bsl2*intens) + (float)b*(2*l - 1)*(au2 + NE*cdexu2 + bsu2*intens));
		}


		if(t == n){
			f = f2;
		}
		else{
			if(t+1 == n){
				f = f1;
			}
			else{
				f = (f2 - f1)*(r - t) + f1; /*linear interpolation*/
			}
		}

		sum = sum+p*w[j]*f;

	}
	fclose(anm500Ptr);

	return sum;
}



/*calculates Sum_(t=nc)^\infty t^2*exp(\Chi_t/kT)*(A_tn + NE*C_tn) by converting
sum to integral using trapizoidal rule and the evaluating integral using
20 point Gaussian integration from NMAX to NC*/
double infsumnle(const int n, const int l, const double y, const double ryd, const int te, const double dw, const int _min,  double bval[])
{
	double sum = 0, p, q, r, a1, a2, cdex, f1, f2, f, colex;
	double al1, au1, colexl1, colexu1, colexl2, colexu2, bal1, bal2, bau1, bau2;
	double bsu1, bsu2, bsl1, bsl2, cdexl1, cdexl2, cdexu1, cdexu2, al2, au2, b;
	int j, t, s;
	double bs1, ba1, bs2, ba2, intens, freq3;
	double x[20] = {-0.9931285992, -0.9639719273, -0.9122344283, -0.8391169718, -0.7463319065,
					-0.6360536807, -0.5108670020, -0.3737060887, -0.2277858511, -0.0765265211,
					0.0765265211, 0.2277858511, 0.3737060887, 0.5108670020, 0.6360536807,
					0.7463319065, 0.8391169718, 0.9122344283, 0.9639719273, 0.9931285992};
	double w[20] = {0.0176140071, 0.0406014298, 0.0626720483, 0.0832767416, 0.1019301198,
					0.1181945320, 0.1316886384, 0.1420961093, 0.1491729865, 0.1527533871,
					0.1527533871, 0.1491729865, 0.1420961093, 0.1316886384, 0.1181945320,
					0.1019301198, 0.0832767416, 0.0626720483,0.0406014298, 0.0176140071};
	FILE *anm500Ptr;


	anm500Ptr = fopen("anm500.dat","rb");
	p = (NC - _min)/2;
	q = (NC + _min)/2;

	for(j = 0; j <= 19; j++){
		r = p*x[j] + q;
		t = floor(r);

		if(t != n){
			if(t <= 500){
				s = (t-1)*(t-2)/2 + n;
				fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
				fread(&a1, sizeof(double), 1, anm500Ptr);
			}
			else{
				a1 = avalaprx (t, n);
			}

			colex = bbcolex(t, n, a1, ryd, te);
			cdex = bbcoldeex(t, n, ryd, te, colex);

			freq3 = pow(freq(t, n, ryd), 3);
			intens = jnu(t, n, ryd, dw);

			au1 = 0;
			al1 = 0;
			bsl1 = 0;
			bsu1 = 0;
			bal1 = 0;
			bau1 = 0;
			colexl1 = 0;
			colexu1 = 0;
			cdexl1 = 0;
			cdexu1 = 0;

			if(l != 0 && l <= n){ /* l' = l-1 */
				au1 = anlmlu(t, n, l-1, ryd);
				bsu1 = 6.78196256*pow(10, 46)*au1 / freq3;
				bal1 = absorp(n, l, t, l-1, ryd, bau1 );
				colexl1 = au1/a1*colex;
				cdexu1 = au1/a1*cdex;
			}

			if(l <= n-1){/* l' = l+1 */
				al1 = anlmll(t, n, l+1, ryd);
				bsl1 = 6.78196256*pow(10, 46)*al1 / freq3;
				bau1 = absorp(n, l, t, l+1, ryd, bsl1 );
				colexu1 = al1/a1*colex;
				cdexl1 = al1/a1*cdex;
			}

			if(t<NMAX){
				b=bval[t-1];}
			else{
				b=1.0;}

			f1 = NE*colexu1 + bau1*intens + NE*colexl1 + bal1*intens;
		}


		if(t+1 != n){
			if(t+1 <= 500){
				s = (t)*(t-1)/2 + n; /*posistion of A_(t+1)n in anm.dat*/
				fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
				fread(&a2, sizeof(double), 1, anm500Ptr);
			}
			else{
				a2 = avalaprx (t+1, n);
			}

			colex = bbcolex(t+1, n, a1, ryd, te);
			cdex = bbcoldeex(t+1, n, ryd, te, colex);

			freq3 = pow(freq(t+1, n, ryd), 3);
			intens = jnu(t+1, n, ryd, dw);

			au2 = 0;
			al2 = 0;
			bsl2 = 0;
			bsu2 = 0;
			bal2 = 0;
			bau2 = 0;
			colexl2 = 0;
			colexu2 = 0;
			cdexl2 = 0;
			cdexu2 = 0;

			if(l != 0 && l <= n){ /* l' = l-1 */
				au2 = anlmlu(t+1, n, l-1, ryd);
				bsu2 = 6.78196256*pow(10, 46)*au2 / freq3;
				bal2 = absorp(n, l, t+1, l-1, ryd, bau2 );
				colexl2 = au2/a2*colex;
				cdexu2 = au2/a2*cdex;
			}

			if(l <= n-1){/* l' = l+1 */
				al2 = anlmll(t+1, n, l+1, ryd);
				bsl2 = 6.78196256*pow(10, 46)*al2 / freq3;
				bau2 = absorp(n, l, t+1, l+1, ryd, bsl2 );
				colexu2 = al2/a2*colex;
				cdexl2 = al2/a2*cdex;
			}

			if(t<NMAX){
				b=bval[t];}
			else{
				b=1.0;}


			f2 = NE*colexu2 + bau2*intens + NE*colexl2 + bal2*intens;
		}


		if(t == n){
			f = f2;
		}
		else{
			if(t+1 == n){
				f = f1;
			}
			else{
				f = (f2 - f1)*(r - t) + f1; /*linear interpolation*/
			}
		}

		sum = sum+p*w[j]*f;

	}
	fclose(anm500Ptr);

	return sum;
}





/*adds contribution of processes populating level nl from levels with principle quantum number NLMAX < m < NLMAX+20*/
double uppernsumin(const int n, const int l, const double y, const double ryd, const int te, const double dw, const int _min,  double bval[])
{
	double sum, p, q, r, a1, a2, cdex, f1, f2, f, colex;
	double al1, au1, colexl1, colexu1, colexl2, colexu2, bal1, bal2, bau1, bau2;
	double bsu1, bsu2, bsl1, bsl2, cdexl1, cdexl2, cdexu1, cdexu2, al2, au2, b;
	int j, t, s;
	double bs1, ba1, bs2, ba2, intens, freq3;

	FILE *anm500Ptr;
	FILE *aminPtr;
	FILE *aplusPtr;

	sum = 0;
	anm500Ptr = fopen("anm500.dat","rb");
	aminPtr = fopen("avalmin300.dat","rb");
	aplusPtr = fopen("avalplus300.dat","rb");

	for(t=_min + 1;t<=_min+20;t+=1){
		if(t <= 500){
			s = (t-1)*(t-2)/2 + n;
			fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
			fread(&a1, sizeof(double), 1, anm500Ptr);
		}
		else{
			a1 = avalaprx (t, n);
		}
		/*a1 = calcaval(t, n, ryd);*/

		colex = bbcolex(t, n, a1, ryd, te);
		cdex = bbcoldeex(t, n, ryd, te, colex);

		freq3 = pow(freq(t, n, ryd), 3);
		intens = jnu(t, n, ryd, dw);

		au1 = 0;
		al1 = 0;
		bsl1 = 0;
		bsu1 = 0;
		bal1 = 0;
		bau1 = 0;
		colexl1 = 0;
		colexu1 = 0;
		cdexl1 = 0;
		cdexu1 = 0;

		if(l != 0 && l < n){ /* l' = l-1 */
			if(t<=300){
				fseek(aplusPtr, (apidx(t,n,l-1))*sizeof(double), SEEK_SET);
				fread(&au1, sizeof(double), 1, aplusPtr);
			}else{
				au1 = anlmlu(t, n, l-1, ryd);
			}
			bsu1 = 6.78196256*pow(10, 46)*au1 / freq3;

			//colexl1 = au1/a1*colex;
			colexl1 = au1/a1*(float)(2*(l-1) + 1)/(2*l + 1)*(float)n*n/t/t*colex;
			
			//cdexu1 = (float)(2*l + 1)/(float)(2*(l - 1) + 1)*exp(y/(float)n/(float)n-y/(float)t/(float)t)*colexl1;
			cdexu1 = au1/a1*(float)n*n/t/t*exp(y/n/n-y/t/t)*colex;

		}

		if(l <= n-1){/* l' = l+1 */
			if(t<=300){
				fseek(aminPtr, (amidx(t,n,l+1))*sizeof(double), SEEK_SET);
				fread(&al1, sizeof(double), 1, aminPtr);
			}else{
				al1 = anlmll(t, n, l+1, ryd);
			}
			bsl1 = 6.78196256*pow(10, 46)*al1 / freq3;

			//colexu1 = al1/a1*colex;
			colexu1 = al1/a1*(float)(2*(l+1) + 1)/(2*l + 1)*(float)n*n/t/t*colex;
			//cdexl1 = (float)(2*l + 1)/(float)(2*(l + 1) + 1)*exp(y/(float)n/(float)n-y/(float)t/(float)t)*colexu1;
			cdexl1 = al1/a1*(float)n*n/t/t*exp(y/n/n-y/t/t)*colex;
		}

		if(t<=NMAX){
			b=bval[t-1];}
		else{
			b=1.0;
		}

		f = exp(y/(float)t/(float)t)*(b*(float)(2*l + 3)*(RAD*al1 + NE*cdexl1 +bsl1*intens) + b*(float)(2*l - 1)*(RAD*au1 + NE*cdexu1 + bsu1*intens));
		sum = sum + f;

	}
	fclose(aplusPtr);
	fclose(anm500Ptr);
	fclose(aminPtr);
	return sum;
}

/*adds contribution of processes depopulating level nl from levels with principle quantum number NLMAX < m < NLMAX+20*/
double uppernsumout(const int n, const int l, const double y, const double ryd, const int te, const double dw, const int _min,  double bval[])
{
	double sum, sump, _error, p, q, r, a1, a2, cdex, f1, f2, f, colex;
	double al1, au1, colexl1, colexu1, colexl2, colexu2, bal1, bal2, bau1, bau2;
	double bsu1, bsu2, bsl1, bsl2, cdexl1, cdexl2, cdexu1, cdexu2, al2, au2, b;
	int j, t, s;
	double bs1, ba1, bs2, ba2, intens, freq3;

	FILE *anm500Ptr;
	FILE *aminPtr;
	FILE *aplusPtr;

	anm500Ptr = fopen("anm500.dat","rb");
	aminPtr = fopen("avalmin300.dat","rb");
	aplusPtr = fopen("avalplus300.dat","rb");

	sum = 0;
	for(t=_min + 1;t<=_min+20;t+=1){
		if(t <= 500){
				s = (t-1)*(t-2)/2 + n;
				fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
				fread(&a1, sizeof(double), 1, anm500Ptr);
			}
			else{
				a1 = avalaprx (t, n);
			}

			colex = bbcolex(t, n, a1, ryd, te);
			cdex = bbcoldeex(t, n, ryd, te, colex);

			freq3 = pow(freq(t, n, ryd), 3);
			intens = jnu(t, n, ryd, dw);

			au1 = 0;
			al1 = 0;
			bsl1 = 0;
			bsu1 = 0;
			bal1 = 0;
			bau1 = 0;
			colexl1 = 0;
			colexu1 = 0;
			cdexl1 = 0;
			cdexu1 = 0;

			if(l != 0 && l <= n){ /* l' = l-1 */
				if(t<=300){
					fseek(aplusPtr, (apidx(t,n,l-1))*sizeof(double), SEEK_SET);
					fread(&au1, sizeof(double), 1, aplusPtr);
				}else{
					au1 = anlmlu(t, n, l-1, ryd);
				}
				bsu1 = 6.78196256*pow(10, 46)*au1 / freq3;
				bau1 = absorp(n, l, t, l-1, ryd, bsu1 );
				//colexl1 = au1/a1*colex;
				colexl1 = au1/a1*(float)(2*(l-1) + 1)/(2*l + 1)*(float)n*n/t/t*colex;
			}

			if(l <= n-1){/* l' = l+1 */
				if(t<=300){
					fseek(aminPtr, (amidx(t,n,l+1))*sizeof(double), SEEK_SET);
					fread(&al1, sizeof(double), 1, aminPtr);
				}else{
					al1 = anlmll(t, n, l+1, ryd);
				}
				bsl1 = 6.78196256*pow(10, 46)*al1 / freq3;
				bal1 = absorp(n, l, t, l+1, ryd, bsl1 );
				//colexu1 = al1/a1*colex;
				colexu1 = al1/a1*(float)(2*(l+1) + 1)/(2*l + 1)*(float)n*n/t/t*colex;
			}

			f = NE*colexu1 + bau1*intens + NE*colexl1 + bal1*intens;
			sum = sum + f;
	}

	fclose(aplusPtr);
	fclose(anm500Ptr);
	fclose(aminPtr);
	return sum;

}

double infsumnlfacc(const int n, const int l, const double y, const double ryd, const int te, const double dw, const int _min,  double bval[])
{
	double sum = 0, sump=100000, p, q, r, a1, a2, cdex, f1, f2, f, colex, _error;
	double al1, au1, colexl1, colexu1, colexl2, colexu2, bal1, bal2, bau1, bau2;
	double bsu1, bsu2, bsl1, bsl2, cdexl1, cdexl2, cdexu1, cdexu2, al2, au2, b;
	int j, t, s;
	double bs1, ba1, bs2, ba2, intens, freq3;
	double x[20] = {-0.9931285992, -0.9639719273, -0.9122344283, -0.8391169718, -0.7463319065,
					-0.6360536807, -0.5108670020, -0.3737060887, -0.2277858511, -0.0765265211,
					0.0765265211, 0.2277858511, 0.3737060887, 0.5108670020, 0.6360536807,
					0.7463319065, 0.8391169718, 0.9122344283, 0.9639719273, 0.9931285992};
	double w[20] = {0.0176140071, 0.0406014298, 0.0626720483, 0.0832767416, 0.1019301198,
					0.1181945320, 0.1316886384, 0.1420961093, 0.1491729865, 0.1527533871,
					0.1527533871, 0.1491729865, 0.1420961093, 0.1316886384, 0.1181945320,
					0.1019301198, 0.0832767416, 0.0626720483,0.0406014298, 0.0176140071};

	FILE *anm500Ptr;
	FILE *aminPtr;
	FILE *aplusPtr;

	anm500Ptr = fopen("anm500.dat","rb");
	aminPtr = fopen("avalmin300.dat","rb");
	aplusPtr = fopen("avalplus300.dat","rb");
	t=_min+21;
	/*integration to level n = NC*/


		p = (NC - t)/2;
		q = (NC + t)/2;

		for(j = 0; j <= 19; j++){
			r = p*x[j] + q;
			t = floor(r);

			if(t != n){
				if(t <= 500){
					s = (t-1)*(t-2)/2 + n;
					fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
					fread(&a1, sizeof(double), 1, anm500Ptr);
				}
				else{
					a1 = avalaprx (t, n);
				}
				/*a1 = calcaval(t, n, ryd);*/

				colex = bbcolex(t, n, a1, ryd, te);
				cdex = bbcoldeex(t, n, ryd, te, colex);

				freq3 = pow(freq(t, n, ryd), 3);
				intens = jnu(t, n, ryd, dw);

				au1 = 0;
				al1 = 0;
				bsl1 = 0;
				bsu1 = 0;
				bal1 = 0;
				bau1 = 0;
				colexl1 = 0;
				colexu1 = 0;
				cdexl1 = 0;
				cdexu1 = 0;

				if(l != 0 && l < n){ /* l' = l-1 */
					if(t<=300){
						fseek(aplusPtr, (apidx(t,n,l-1))*sizeof(double), SEEK_SET);
						fread(&au1, sizeof(double), 1, aplusPtr);
					}else{
						au1 = anlmlu(t, n, l-1, ryd);
					}
					bsu1 = 6.78196256*pow(10, 46)*au1 / freq3;
					bal1 = absorp(n, l, t, l-1, ryd, bau1 );
					
					//colexl1 = au1/a1*colex;
					colexl1 = au1/a1*(float)(2*(l-1) + 1)/(2*l + 1)*(float)n*n/t/t*colex;
					//cdexu1 = (float)(2*l + 1)/(2*(l - 1) + 1)*exp(y/n/n-y/t/t)*colexl1;
					cdexu1 = au1/a1*(float)n*n/t/t*exp(y/n/n-y/t/t)*colex;
				}

				if(l <= n-1){/* l' = l+1 */
					if(t<=300){
						fseek(aminPtr, (amidx(t,n,l+1))*sizeof(double), SEEK_SET);
						fread(&al1, sizeof(double), 1, aminPtr);
					}else{
						al1 = anlmll(t, n, l+1, ryd);
					}
					bsl1 = 6.78196256*pow(10, 46)*al1 / freq3;
					bau1 = absorp(n, l, t, l+1, ryd, bsl1 );
					//colexu1 = al1/a1*colex;
					colexu1 = al1/a1*(float)(2*(l+1) + 1)/(2*l + 1)*(float)n*n/t/t*colex;
					
					//cdexl1 = (float)(2*l + 1)/(float)(2*(l + 1) + 1)*exp(y/(float)n/(float)n-y/(float)t/(float)t)*colexu1;
					cdexl1 = al1/a1*(float)n*n/t/t*exp(y/n/n-y/t/t)*colex;
				}

				if(t<NMAX){
					b=bval[t-1];}
				else{
					b=1.0;}

				f1 = exp(y/(float)t/(float)t)*(b*(float)(2*l + 3)*(RAD*al1 + NE*cdexl1 +bsl1*intens) + b*(float)(2*l - 1)*(RAD*au1 + NE*cdexu1 + bsu1*intens));
			}

			if(t+1 != n){
				if(t+1 <= 500){
					s = (t)*(t-1)/2 + n; /*posistion of A_(t+1)n in anm.dat*/
					fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
					fread(&a2, sizeof(double), 1, anm500Ptr);
				}
				else{
					a2 = avalaprx (t+1, n);
				}
				/*a2 = calcaval(t+1, n, ryd);*/
				colex = bbcolex(t+1, n, a1, ryd, te);
				cdex = bbcoldeex(t+1, n, ryd, te, colex);

				freq3 = pow(freq(t+1, n, ryd), 3);
				intens = jnu(t+1, n, ryd, dw);

				au2 = 0;
				al2 = 0;
				bsl2 = 0;
				bsu2 = 0;
				bal2 = 0;
				bau2 = 0;
				colexl2 = 0;
				colexu2 = 0;
				cdexl2 = 0;
				cdexu2 = 0;

				if(l != 0 && l < n){ /* l' = l-1 */
					if(t+1<=300){
						fseek(aplusPtr, (apidx(t+1,n,l-1))*sizeof(double), SEEK_SET);
						fread(&au2, sizeof(double), 1, aplusPtr);
					}else{
						au2 = anlmlu(t+1, n, l-1, ryd);
					}
					bsu2 = 6.78196256*pow(10, 46)*au2 / freq3;
					bal2 = absorp(n, l, t+1, l-1, ryd, bau2 );
					//colexl2 = au2/a2*colex;
					colexl2 = au2/a2*(float)(2*(l-1) + 1)/(2*l + 1)*(float)n*n/(t+1)/(t+1)*colex;
					//cdexu2 = (float)(2*l + 1)/(2*(l - 1) + 1)*exp(y/n/n-y/(t+1)/(t+1))*colexl2;
					cdexu2 = au2/a2*(float)n*n/(t+1)/(t+1)*exp(y/n/n-y/(t+1)/(t+1))*colex;

				}

				if(l <= n-1){/* l' = l+1 */
					if(t+1<=300){
						fseek(aminPtr, (amidx(t+1,n,l+1))*sizeof(double), SEEK_SET);
						fread(&al2, sizeof(double), 1, aminPtr);
					}else{
						al2 = anlmll(t+1, n, l+1, ryd);
					}
					bsl2 = 6.78196256*pow(10, 46)*al2 / freq3;
					bau2 = absorp(n, l, t+1, l+1, ryd, bsl2 );
					//colexu2 = al2/a2*colex;
					colexu2 = al2/a2*(float)(2*(l+1) + 1)/(2*l + 1)*(float)n*n/(t+1)/(t+1)*colex;
					
					//cdexl2 = (float)(2*l + 1)/(2*(l + 1) + 1)*exp(y/n/n-y/(t+1)/(t+1))*colexu2;
					cdexl2 = al2/a2*(float)n*n/(t+1)/(t+1)*exp(y/n/n-y/(t+1)/(t+1))*colex;
					
				}

				if(t<NMAX){
					b=bval[t];}
				else{
					b=1.0;}


				f2 = exp(y/(float)(t+1)/(float)(t+1))*(b*(float)(2*l + 3)*(RAD*al2 + NE*cdexl2 +bsl2*intens) + b*(float)(2*l - 1)*(RAD*au2 + NE*cdexu2 + bsu2*intens));
			}

			if(t == n){
				f = f2;
			}
			else{
				if(t+1 == n){
					f = f1;
				}
				else{
					f = (f2 - f1)*(r - t) + f1; /*linear interpolation*/
				}
			}

			sum = sum+p*w[j]*f;

		}


	fclose(aplusPtr);
	fclose(anm500Ptr);
	fclose(aminPtr);
	return sum;
}

double infsumnleacc(const int n, const int l, const double y, const double ryd, const int te, const double dw, const int _min,  double bval[])
{
	double sum = 0, sump, _error, p, q, r, a1, a2, cdex, f1, f2, f, colex;
	double al1, au1, colexl1, colexu1, colexl2, colexu2, bal1, bal2, bau1, bau2;
	double bsu1, bsu2, bsl1, bsl2, cdexl1, cdexl2, cdexu1, cdexu2, al2, au2;
	int j, t, s;
	double bs1, ba1, bs2, ba2, intens, freq3;
	double x[20] = {-0.9931285992, -0.9639719273, -0.9122344283, -0.8391169718, -0.7463319065,
					-0.6360536807, -0.5108670020, -0.3737060887, -0.2277858511, -0.0765265211,
					0.0765265211, 0.2277858511, 0.3737060887, 0.5108670020, 0.6360536807,
					0.7463319065, 0.8391169718, 0.9122344283, 0.9639719273, 0.9931285992};
	double w[20] = {0.0176140071, 0.0406014298, 0.0626720483, 0.0832767416, 0.1019301198,
					0.1181945320, 0.1316886384, 0.1420961093, 0.1491729865, 0.1527533871,
					0.1527533871, 0.1491729865, 0.1420961093, 0.1316886384, 0.1181945320,
					0.1019301198, 0.0832767416, 0.0626720483,0.0406014298, 0.0176140071};
	FILE *anm500Ptr;
	FILE *aminPtr;
	FILE *aplusPtr;

	anm500Ptr = fopen("anm500.dat","rb");
	aminPtr = fopen("avalmin300.dat","rb");
	aplusPtr = fopen("avalplus300.dat","rb");

	/*integration to level n = NC*/

	t=_min+21;
		p = (NC - t)/2;
		q = (NC + t)/2;

		for(j = 0; j <= 19; j++){
			r = p*x[j] + q;
			t = floor(r);

			if(t != n){
				if(t <= 500){
					s = (t-1)*(t-2)/2 + n;
					fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
					fread(&a1, sizeof(double), 1, anm500Ptr);
				}
				else{
					a1 = avalaprx (t, n);
				}

				colex = bbcolex(t, n, a1, ryd, te);
				cdex = bbcoldeex(t, n, ryd, te, colex);

				freq3 = pow(freq(t, n, ryd), 3);
				intens = jnu(t, n, ryd, dw);

				au1 = 0;
				al1 = 0;
				bsl1 = 0;
				bsu1 = 0;
				bal1 = 0;
				bau1 = 0;
				colexl1 = 0;
				colexu1 = 0;
				cdexl1 = 0;
				cdexu1 = 0;

				if(l != 0 && l <= n){ /* l' = l-1 */
					if(t<=300){
						fseek(aplusPtr, (apidx(t,n,l-1))*sizeof(double), SEEK_SET);
						fread(&au1, sizeof(double), 1, aplusPtr);
					}else{
						au1 = anlmlu(t, n, l-1, ryd);
					}
					bsu1 = 6.78196256*pow(10, 46)*au1 / freq3;
					bal1 = absorp(n, l, t, l-1, ryd, bsu1 );
					//colexl1 = au1/a1*colex;
					colexl1 = au1/a1*(float)(2*(l-1) + 1)/(2*l + 1)*(float)n*n/t/t*colex;

				}

				if(l <= n-1){/* l' = l+1 */
					if(t<=300){
						fseek(aminPtr, (amidx(t,n,l+1))*sizeof(double), SEEK_SET);
						fread(&al1, sizeof(double), 1, aminPtr);
					}else{
						al1 = anlmll(t, n, l+1, ryd);
					}
					bsl1 = 6.78196256*pow(10, 46)*al1 / freq3;
					bau1 = absorp(n, l, t, l+1, ryd, bsl1 );
					//colexu1 = al1/a1*colex;
					colexu1 = al1/a1*(float)(2*(l+1) + 1)/(2*l + 1)*(float)n*n/t/t*colex;
				}

				f1 = NE*colexu1 + bau1*intens + NE*colexl1 + bal1*intens;
			}


			if(t+1 != n){
				if(t+1 <= 500){
					s = (t)*(t-1)/2 + n; /*posistion of A_(t+1)n in anm.dat*/
					fseek( anm500Ptr, (s-1)*sizeof(double), SEEK_SET);
					fread(&a2, sizeof(double), 1, anm500Ptr);
				}
				else{
					a2 = avalaprx (t+1, n);
				}

				colex = bbcolex(t+1, n, a1, ryd, te);
				cdex = bbcoldeex(t+1, n, ryd, te, colex);

				freq3 = pow(freq(t+1, n, ryd), 3);
				intens = jnu(t+1, n, ryd, dw);

				au2 = 0;
				al2 = 0;
				bsl2 = 0;
				bsu2 = 0;
				bal2 = 0;
				bau2 = 0;
				colexl2 = 0;
				colexu2 = 0;
				cdexl2 = 0;
				cdexu2 = 0;

				if(l != 0 && l <= n){ /* l' = l-1 */
					if(t+1<=300){
						fseek(aplusPtr, (apidx(t+1,n,l-1))*sizeof(double), SEEK_SET);
						fread(&au2, sizeof(double), 1, aplusPtr);
					}else{
						au2 = anlmlu(t+1, n, l-1, ryd);
					}
					bsu2 = 6.78196256*pow(10, 46)*au2 / freq3;
					bal2 = absorp(n, l, t+1, l-1, ryd, bsu2 );
					//colexl2 = au2/a2*colex;
					colexl2 = au2/a2*(float)(2*(l-1) + 1)/(2*l + 1)*(float)n*n/(t+1)/(t+1)*colex;

				}

				if(l <= n-1){/* l' = l+1 */
					if(t+1<=300){
						fseek(aminPtr, (amidx(t+1,n,l+1))*sizeof(double), SEEK_SET);
						fread(&al2, sizeof(double), 1, aminPtr);
					}else{
						al2 = anlmll(t+1, n, l+1, ryd);
					}
					bsl2 = 6.78196256*pow(10, 46)*al2 / freq3;
					bau2 = absorp(n, l, t+1, l+1, ryd, bsl2 );
					//colexu2 = al2/a2*colex;
					colexu2 = al2/a2*(float)(2*(l+1) + 1)/(2*l + 1)*(float)n*n/(t+1)/(t+1)*colex;

				}

				f2 =NE*colexu2 + bau2*intens + NE*colexl2 + bal2*intens;
			}

			if(t == n){
				f = f2;
			}
			else{
				if(t+1 == n){
					f = f1;
				}
				else{
					f = (f2 - f1)*(r - t) + f1; /*linear interpolation*/
				}
			}

			sum = sum+p*w[j]*f;

		}

	fclose(aplusPtr);
	fclose(anm500Ptr);
	fclose(aminPtr);

	return sum;
}

int round_up_5(int number)
{
	int remainder;

	remainder = number % 5;

	if(remainder == 0)
	{
		return number;
	}else{
		return number + 5 - remainder;
	}
}



/*calculates approx Einstein A value for transisions between high
energy levels*/
double avalaprx (const int n, const int m)
{
	double gaunt, a;

		/*nu = freq(n, m, ryd);

		/*Eq 1.42 in 1935MNRAS..96...77M*/
		/*gaunt = 1 - 0.1728*pow( nu / ryd*ATOMZ, 1.0/3 )*(2.0*ryd*ATOMZ*ATOMZ/(n*n*nu) + 1);*/
		/*Baker et al 1938*/
		gaunt = 1 - 0.1728*(1 + m*m/n/n)*pow(1 - m*m/n/n, -2.0/3)*pow(n*n, -2.0/3);

		/*Eqs 3.1 & 3.2 Brocklehurst 1970*/
		a = 15745959000*ATOMZ*ATOMZ*ATOMZ*ATOMZ*gaunt / ((float)n*n*n*m*(n*n - m*m));

	return a;
}


/*calculates hypergeometric function 2F1(a, b, g, x) using
method described by 1979JPhB...12.1053V*/
double hyprege(const int a, const int b, const int g, const double x)
{
	int j, aa;
	double hypa[NC*2];

	hypa[0] = 1;
	hypa[1] = 1 - (float)b/g*x;

	for(j = 2; j <= -a; j++){
		aa = 1 - j;
		hypa[j] = -(g - aa - b + (b - aa)*(1 - x))/(aa - g)*hypa[j-1] - aa*(1 - x)/(aa - g)*hypa[j-2];
	}
	return hypa[-a];
}

/*calculates b-b radial matrix element for transision nl -> mk (n > m)
using eq 14 of 1979JPhB...12.1053V
hyp[j] is hypergeometric function in eq 15 of 1979JPhB...12.1053V
where -j is b in 2F1(a, b, c, z)*/
double bbradmat(const int n, const int l, const int m, const int k)
{
	int t, a, b, q, j, g, aa;
	double sum = 0, it, hyp[NC*2], x, term, bt = 0;


	a = -n + l + 1;
	g = 2*l + 2;
	x = 2.0*m/(n + m);
	hyp[0] = hyprege(a, m + l + 3, g, x);  /*hyp[1] = 2F1(a, b-0, g, x)*/

	hyp[1] = hyprege(a, m + l + 2, g, x);  /*hyp[1] = 2F1(a, b-1, g, x)*/


	t=0;
	it = -log(m*n)-fact(2*l + 1) + 0.5*(fact(n + l) - fact(m + k) - fact(m - k - 1) - fact(n - l - 1))
		+ fact(m + l + 2 - t) + (m + l + 3 - t)*log((float)n*m/(n + m)) + (l + 1)*log(2.0/n) + m*log(2.0/m);

	if(hyp[0] < 0){
			term = -exp(it+log(-hyp[0]));
		}
		else{
			if(hyp[0] == 0){
				term = 0;
			}
			else{
				term = exp(it+log(hyp[0]));
			}
		}

	if((n + k + 1)%2 != 0){
		term = -term;
	}

	sum = sum + term;

	for(t = 1; t <= m - k - 1; t++){
		b = m + l + 3 - t;
		hyp[t+1] = - (2*b-g-b*x +a*x)/(g-b)*hyp[t] - b*(x-1)/(g-b)*hyp[t-1];
	}

	for(t = 1; t <= m - k - 1; t++){

		bt = log((float)m/(2*t)) + log((m - t)*(m - t + 1) - (k + 1)*k) + bt; /*bt = ln(|b_t|) from van regemorter 1979*/

		it = -log(m*n)-fact(2*l + 1) + 0.5*(fact(n + l) - fact(m + k) - fact(m - k - 1) - fact(n - l - 1))
				+ fact(m + l + 2 - t) + (m + l + 3 - t)*log((float)n*m/(n + m)) + (l + 1)*log(2.0/n) + m*log(2.0/m);

		if(hyp[t] < 0){
			term = -exp(bt+it+log(-hyp[t]));

		}
		else{
			if(hyp[t] == 0){
				term = 0;
			}
			else{
				term = exp(bt+it+log(hyp[t]));
			}
		}

		/*b_t in van regemorter negative for odd t*/
		if(t%2 != 0){
			term = -term;
		}

		if((n + k + 1)%2 != 0){
			term = -term;
		}

		sum = sum + term;

	}

	return sum;
}

/*Calculates Einstein A value for transision nl -> ml+1*/
double anlmlu ( const int nn, const int mm, const int l, const double ryd )
{
	double au, freq3;

	freq3 = pow(freq(nn, mm, ryd), 3);

	if (l >= mm-1){
		au = 0;
	}
	else{
		au = bbradmat(nn, l, mm, l+1);
		au = 7.51976932*pow(10, -38)*freq3/(ATOMZ*ATOMZ)*(float)(l + 1)/(2*l + 1)*au*au;

	}

	return au;
}


/*Calculates Einstein A value for transision nl -> nl-1*/
double anlmll ( const int nn, const int mm, const int l, const double ryd )
{
	double al, freq3;

	freq3 = pow(freq(nn, mm, ryd), 3);

	if (l == 0 || l > mm){
		al = 0;
	}
	else{
		al = bbradmat(nn, l, mm, l-1);
		al = 7.51976932*pow(10, -38)*freq3/(ATOMZ*ATOMZ)*(float)l/(2*l + 1)*al*al;
	}

	return al;
}


/*Calculates absorbtion rate coefficient for transision ml -> nk
where bs is stim em rate for same transition*/
double absorp1 ( const int nn, const int mm, const int l, const int k, const double ryd, const double bs )
{
	double ba, freq3;
	freq3 = pow(freq(nn, mm, ryd), 3);

	if( k == l-1 ){
		ba = (float)(2*l + 1)/(2*l - 1)*bs;
		printf("ba = %d/%d\n",(2*l + 1),(2*l - 1));
	}
	else if( k == l+1){
		ba = (float)(2*l + 1)/(2*l + 3)*bs;
		printf("ba = %d/%d\n",(2*l + 1),(2*l + 1));
	}
	else{
		ba = 0;
	}

	return ba;

}

/*Calculates absorbtion rate coefficient for transision ml -> nk
where bs is stim em rate for same transition*/
double absorp ( const int mm, const int l, const int nn, const int k, const double ryd, const double bs )
{
	double ba;


	ba = (float)(2*k + 1)/(2*l + 1)*bs;

	return ba;

}

/*Calculates the hydrogenic bound-bound matrix element for a transition*/
/*nl -> ml' n > m using method of Brocklehurst 1971MNRAS.153..471B*/
/*Output is Einstein A-values A(nl -> ml')*/
void aval ( const int nn, const int mm, double _au[], double _al[], const double ryd )
{
	int l;
	double freq3;

	freq3 = pow(freq(nn, mm, ryd), 3);

	for(l = 0; l < nn; l++){
		_au[ l ] = 0;
		_al[ l ] = 0;
	}

	for( l = 1; l <= mm; l++ ){
		/* from Gordon 1929AnP...394.1031G*/
		/*tested with tables of 1957ApJS....3...37G*/
		/*al[ l ] = (float)pow(-1, mm - l)/4.0*bbfactl(nn, mm, l)*(hypgeo(-nn + l + 1, -mm + l, 2*l, (float)(-4*nn*mm)/((nn - mm)*(nn - mm)))
			- pow((float)(nn - mm)/(nn +mm), 2)*hypgeo(-nn + l - 1, -mm + l, 2*l, (float)(-4*nn*mm)/((nn - mm)*(nn - mm))));*/

		_al[ l ] = bbradmat(nn, l, mm, l-1);
		_al[ l ] = 7.51976932*pow(10, -38)*freq3/(ATOMZ*ATOMZ)*(float)l/(2*l + 1)*_al[ l ]*_al[ l ];
	}

	for( l = 0; l <= mm - 2; l++ ){
		/* from Gordon 1929AnP...394.1031G with n <-> m and l -> l+1*/
		/*tested with tables of 1957ApJS....3...37G*/
		/*au[ l ] = (float)pow(-1, nn - l - 1)/4.0*bbfactu(nn, mm, l)*pow(mm - nn, nn + mm - 2*l - 4)*(hypgeo(-mm + l + 2, -nn + l + 1, 2*l + 2, (float)(-4*nn*mm)/((mm - nn)*(mm - nn)))
				- pow((float)(mm - nn)/(nn +mm), 2)*hypgeo(-mm + l, -nn + l + 1, 2*l + 2, (float)(-4*nn*mm)/((mm - nn)*(mm - nn))));*/

		_au[ l ] = bbradmat(nn, l, mm, l+1);
		_au[ l ] = 7.51976932*pow(10, -38)*freq3/(ATOMZ*ATOMZ)*(float)(l + 1)/(2*l + 1)*_au[ l ]*_au[ l ];

	}
}

/*Calculates the hypergeometric function 2F1(a,b;c;z) for the case*/
/*where either a <= 0 or b <= 0*/
double hypgeo(const int a, const int b, const int c, const double z)
{
	int j = 0, stop = 0;
	double sum = 0, term = 1;

	do {
		sum = sum + term;
		term = term*((a + j)*(b + j)*z) / ((c + j)*(j + 1));
		j++;
		if (term == 0){
			stop = 1;
		}
	} while(stop == 0);

	return sum;
}


/*calculates Einstein B coefficients from A coefficients*/
void bval ( const int nn, const int mm, const double au[], const double al[], double bau[],
			double bal[], double bsu[], double bsl[], const double ryd )
{
	int l;
	double freq3;

	for(l = 0; l < nn; l++){
		bau[ l ] = 0;
		bal[ l ] = 0;
		bsu[ l ] = 0;
		bal[ l ] = 0;
	}


	freq3 = pow(freq(nn, mm, ryd), 3);

	for( l = 1; l <= mm; l++ ){
		bsl[ l ] = 6.78196256*pow(10, 46)*al[ l ] / freq3;
		bal[ l ] = (float)(2*l + 1)/(2*l - 1)*bsl[ l ];
	}

	for( l = 0; l <= mm - 2; l++ ){
		bsu[ l ] = 6.78196256*pow(10, 46)*au[ l ] / freq3;
		bau[ l ] = (float)(2*l + 1)/(2*l + 3)*bsu[ l ];
	}
}


/*calculates Einstein coef between energy levels from ang mom coef
e.g. A(n -> m) form A(nl -> ml')*/
double abnm ( const int nn, const int mm, const double u[], const double l[] )
{
	int j;
	double abnmsum;

	abnmsum = 0;
	for( j = 0; j <= mm; j++ ){
		abnmsum = abnmsum + (float)(2*j + 1)/(nn*nn)*(l[ j ] + u[ j ]);
	}

	return abnmsum;
}

/*calculates bound-free coef from level n from ang mom coef
e.g. radrecomb(n -> E) form radrecomb(nl -> E)*/
double bfn ( const int nn, const double bf[] )
{
	int l;
	double nsum;

	nsum = 0;
	for( l = 0; l <= nn; l++ ){
		nsum = nsum + bf[ l ];
	}

	return nsum;
}

/*calculates bound-free coef from level n from ang mom coef
e.g. radrecomb(n -> E) form radrecomb(nl -> E)*/
double bfn2 ( const int nn, const double bf[] )
{
	int l;
	double nsum;

	nsum = 0;
	for( l = 0; l <= nn; l++ ){
		nsum = nsum + (float)(2*l + 1)/(nn*nn)*bf[ l ];
	}

	return nsum;
}

/*calculates 1/(2l+1)!*Sqrt[((m+l+1)!(n+l)!)/((m-l-2)!(n-l-1)!)]*/
double bbfactu ( const int fn, const int fm, const int fl)
{
	double factu;

	factu = - fact(2*fl + 1) + 0.5*(fact(fm + fl + 1) + fact(fn + fl) -
		fact(fm - fl - 2) - fact(fn - fl - 1)) + (fl + 2)*log(4*fn*fm) -
		(fm + fn)*log(fn + fm);

	return exp(factu);
}


/*calculates 1/(2l-1)!*Sqrt[((n+l)!(m+l-1)!)/((n-l-1)!(m-l)!)]*/
double bbfactl ( const int fn, const int fm, const int fl)
{
	double factl;

	factl = - fact(2*fl - 1) + 0.5*(fact(fn + fl) + fact(fm + fl - 1) -
		fact(fn - fl - 1) - fact(fm - fl)) + (fl + 1)*log(4*fn*fm) +
		(fn + fm - 2*fl - 2)*log(fn - fm) - (fn + fm)*log(fn + fm);


	return exp(factl);
}




/*calculates the natural log of a! */
double fact( const int a )
{
	int j;
	double x = 0;

	if( a != 0 && a != 1 ){
		for( j=1; j <= a; j++ ){
			x = x + log(j);
		}
	}
	return x;
}


/*Calculates energy difference E_nm in Joule between bound levels n and m of
Bohr atom in cgs*/
double energynm( const int nn, const int mm, const double ryd )
{
	return ATOMZ*ATOMZ*ryd*1.98644582*pow(10, -16)*(1.0/(mm*mm) - 1.0/(nn*nn));
}


/*Calculates frequency of photon due to transition between levels n and m*/
double freq(const int nn, const int mm, const double ryd)
{
	return ATOMZ*ATOMZ*ryd*29979245800*(1.0/(mm*mm) - 1.0/(nn*nn));
}


/*Calculates oscillator strength f_mn from Einstein A-value*/
/*F(nl -> n'l') in 1968ApJS...17..445G**/
/*Draine 2011 eq 6.20*/
/*Agrees with 1935MNRAS..96...77M up to about 10% - due to approx in his calc?*/
double oscstrength( const int nn, const int mm, const double au[], const double al[], const double ryd )
{
	double eps0, elemcharge, asum, f, osc;

	asum = abnm(nn, mm, au, al);
	f = freq(nn, mm, ryd);
	osc = -24544.3253*asum/(1.8215954*pow(10, -17)*f*f);

	if (nn > mm){
		osc = -osc*(float)(nn*nn)/(mm*mm); /*f12*/
	}

	return osc;
}

/*Calculates oscillator strength f_mn from Einstein A-value*/
/*F(nl -> n'l') in 1968ApJS...17..445G**/
/*Agrees with 1935MNRAS..96...77M up to about 10% - due to approx in his calc?*/
double oscstrengthnm( const int nn, const int mm, const double anm, const double ryd )
{
	double f, osc;

	f = freq(nn, mm, ryd);

	osc = -24544.3253*anm/(1.82159534*pow(10, -17)*f*f);

	if (nn > mm){
		osc = -osc*(float)(nn*nn)/(mm*mm); /*f12*/
	}

	return osc;
}

/*Calculates rate coefficients of b-b collisional processes*/
/*excitation rate coef: Eq 17 of 1980PhRvA..22..940V*/
/*de-excitation rate coef: detailed balace*/
void colratecoef( const int nn, const int mm, const double anm, const double ryd,
					const int te, double bbcol[] )
{
	int s;
	double apn, bpn, epi, emn, ryde, k, kte, delta, gamma, p, bp, x;
	/*bbcol = (ex rate coef, de-ex rate coef)*/

	emn = energynm(nn, mm, ryd)*6.24150764*pow(10, 11); /*energy level difference in eV*/
	ryde = 1.98644582*pow(10, -16)*ryd*6.24150764*pow(10, 11);  /*Rydberg energy in eV*/
	epi = ryde*ATOMZ*ATOMZ/(mm*mm);

	/*dimensionless coef A_pn of Eq 11 1980PhRvA..22..940V; p = m*/
	apn = oscstrengthnm(nn, mm, anm, ryd)*(2*ryde)/(emn);

	p = (float) mm;

	bp = 1.4*log(p)/p - 0.7/p - 0.51/(p*p) + 1.16/(p*p*p) - 0.55/(p*p*p*p);


	/*dimensionless coef B_pn of Eq 12 1980PhRvA..22..940V; p = m*/
	bpn = (4*ryde*ryde)/(nn*nn*nn)*(1/(emn*emn) + 4*epi/(3*emn*emn*emn) + bp*(epi*epi)/(emn*emn*emn*emn));

	s = abs(nn - mm);

	if (s > 0.6)	/*abs() only valid for integers*/
		x = s - 0.6;
	else
		x = 0.6 - s;

	k = 8.6173303*pow(10, -5); /*boltzmann const in eV/K*/
	delta = exp(-bpn/apn) + (0.06*s*s)/(nn*mm*mm);
	kte = k*te;

	gamma = ryde*log(1 + (float)(mm*mm*mm*kte)/ryde)*(3 + 11*(float)(s*s)/(mm*mm))*
			1.0/(6 + 1.6*nn*s + 0.3/(s*s) + 0.8*pow(nn, 1.5)/pow(s, 0.5)*x);

	/*excitation rate in cm^3/s*/
	bbcol[ 0 ] = 1.6*pow(10, -7)*sqrt(kte)/(kte + gamma)*exp(-emn/(kte))*(apn*log(0.3*kte/ryde + delta) + bpn);

	/*de-excitation rate in cm^3/s*/
	bbcol[ 1 ] = (bbcol[ 0 ]*mm*mm)/(nn*nn)*exp(emn/kte);

}

/*Calculates rate coefficient of b-b collisional excitation*/
/*excitation rate coef: Eq 17 of 1980PhRvA..22..940V*/
double bbcolex(const int nn, const int mm, const double anm, const double ryd, const int te)
{
	int s;
	double apn, bpn, epi, emn, ryde, k, kte, delta, gamma, p, bp, x, exrate;

	emn = energynm(nn, mm, ryd)*6.24150764*pow(10, 11); /*energy level difference in eV*/
	ryde = 1.98644582*pow(10, -16)*ryd*6.24150764*pow(10, 11);  /*Rydberg energy in eV*/
	epi = ryde*ATOMZ*ATOMZ/(mm*mm);

	/*dimensionless coef A_pn of Eq 11 1980PhRvA..22..940V; p = m*/
	apn = oscstrengthnm(nn, mm, anm, ryd)*(2*ryde)/(emn);
	
	p = (float) mm;
	bp = 1.4*log(p)/p - 0.7/p - 0.51/(p*p) + 1.16/(p*p*p) - 0.55/(p*p*p*p);

	/*dimensionless coef B_pn of Eq 12 1980PhRvA..22..940V; p = m*/
	bpn = (4*ryde*ryde)/(nn*nn*nn)*(1/(emn*emn) + 4*epi/(3*emn*emn*emn) + bp*(epi*epi)/(emn*emn*emn*emn));
	
	s = abs(nn - mm);
	
	if (s > 0.6)	/*abs() only valid for integers*/
		x = s - 0.6;
	else
		x = 0.6 - s;

	k = 8.6173303*pow(10, -5); /*boltzmann const in eV/K*/
	delta = exp(-bpn/apn) + (0.06*s*s)/(nn*mm*mm);
	kte = k*te;
	gamma = ryde*log(1 + (float)(mm*mm*mm*kte)/ryde)*(3 + 11*(float)(s*s)/(mm*mm))*
		1.0/(6 + 1.6*nn*s + 0.3/(s*s) + 0.8*pow(nn, 1.5)/pow(s, 0.5)*x);
	
	/*excitation rate in cm^3/s*/
	exrate = 1.6*pow(10, -7)*sqrt(kte)/(kte + gamma)*exp(-emn/(kte))*(apn*log(0.3*kte/ryde + delta) + bpn);
	
	if(isinf(exrate) || exrate < 0){
		exrate = 0;
	}
	
	return exrate;
}

/*calculates b-b de-excitation rate coefficient from level n to m
by using principle of detailed balance*/
double bbcoldeex(const int n, const int m, const double ryd, const int te, const double colex)
{

	/*de-excitation rate in cm^3/s*/
	return (colex*m*m)/(n*n)*exp(energynm(n, m, ryd)/(1.380658*pow(10, -16)*te));
}


/*Calculates rate coefficients of b-f collisional processes*/
/*col ionization: Eq 8 of 1980PhRvA..22..940V*/
/*3body recombination: Eq 9 of 1980PhRvA..22..940V*/
void colionratecoef( const int nn, const double ryd, const int te, double bfcol[], const double y )
{
	double ryde, eni, k, kte, x;

	ryde = 1.98644582*pow(10, -16)*ryd*6.24150764*pow(10, 11);  /*Rydberg energy in eV*/
	eni = ryde*ATOMZ*ATOMZ/(nn*nn);
	k = 8.6173303*pow(10, -5); /*boltzmann const in eV/K*/
	kte = k*te;
	eni = eni/kte;
	x = pow(eni, 2.33) + 4.38*pow(eni, 1.72) +1.32*eni;

	/*coll ionization rate in cm^3/s*/
	/*agrees with fig 1 of 1980PhRvA..22..940V*/
	bfcol[ 0 ] = (9.56*pow(10, -6)*pow(kte, -1.5)*exp(-eni))/x;
	/*3body recombination rate in cm^6/s - this needs to change if not working with H*/
	bfcol[ 1 ] = 4.14133233*pow(10, -16)*pow(te,-3.0/2)*nn*nn*exp(y/nn/nn)*bfcol[ 0 ];
	/*bfcol[ 1 ] = (3.17*pow(10, -27)*pow(kte, -3)*nn*nn)/x;*/


}


//Calulates collision rates for angular momentum changing collisions by protons or electrons*/
//changes needed in dnl if charge of colliding particle != +/-1*/
//Uses formulism of 1964MNRAS.127..165P with dnl defined by 2016MNRAS.459.3498G
void colmomcoef( const int n, double cl[], double cu[], const int te, const double me, const double mp, const double anl[], const double mcoll, const int z)
{
	int l;
	double  redmass, dnll,dnlu, dnl, pc, pc1, lt;



	if(NE==0){
		for( l = 0; l <= n-1; l++ ){
			cu[l] = 0;
			cl[l] = 0;
		}
	}else{

	/*mcoll = mass of colliding particle*/
	redmass = (mcoll*(ATOMN*mp + me) / (mcoll + ATOMN*mp + me));   /*reduced mass of colliding system in g*/

	for( l = 0; l <= n-1; l++ ){
		dnl = 1.0*z*z/(ATOMZ*ATOMZ)*6*n*n*(n*n - l*l - l - 1);
		dnlu = 1.0*z*z/(ATOMZ*ATOMZ)*6.0*n*n*(l+1.0)*(n*n - l*l - 2*l - 1)/(2.0*l+1);

		dnll = 1.0*z*z/(ATOMZ*ATOMZ)*6.0*n*n*l*(n*n - l*l)/(2.0*l+1);
		pc = 1.678 + log10((float)te/NE); //error in Brock71
		if(NE == 0){pc=0;}
		if(anl[l] != 0){    /*lt = 0 for 1s and 2s levels*/
			lt = 1.0/anl[l];
			pc1 = 10.95 + log10(me*te*lt*lt/redmass);
		if(pc1 < pc || NE == 0){
				pc = pc1;
			}
		}
		//c[l]=9.933*pow(10, -6)*sqrt(redmass/me)*dnl/sqrt(te)/*NE*/*(11.538 + log10(te*me/(dnl*redmass)) + pc);
		cu[ l ] = 9.933*pow(10, -6)*sqrt(redmass/me)*dnlu/sqrt(te)/*NE*/*(11.538 + log10(te*me/(dnl*redmass)) + pc);
		//cl[ l ] = 9.933*pow(10, -6)*sqrt(redmass/me)*dnll/sqrt(te)/*NE*/*(11.538 + log10(te*me/(dnl*redmass)) + pc);


	}

	cl[0] = 0;

	for( l = 1; l < n; l++){
		cl[l] = (2.0*(l - 1) + 1.0)/(2*l + 1)*cu[l-1];

	}

	cu[n-1] = 0;

	}
}

/*Calulates collision rates for angular momentum changing collisions by protons or electrons*/
/*changes needed in dnl if charge of colliding particle != +/-1*/
/*Uses method described by 1971MNRAS.153..471B*/
/*1964MNRAS.127..165P claims electron collisions can be neglected wrt protons*/
/*but storey and hummer includes them*/
void colmomcoefbrock( const int n, double cl[], double cu[], const int te, const double me, const double mp, const double anl[], const int nlmaxx, const double mcoll, const int z)
{
	int l;
	double   redmass, dnl, pc, pc1, lt;
	double* c;
	double* cl1;
	double* cu1;


	c = calloc(nlmaxx, sizeof(double));
	cl1 = calloc(nlmaxx, sizeof(double));
	cu1 = calloc(nlmaxx, sizeof(double));
	if(NE==0){
		for( l = 0; l <= n-1; l++ ){
			cu[l] = 0;
			cl[l] = 0;
		}
	}else{

	/*mcoll = mass of colliding particle*/
	redmass = (mcoll*(ATOMN*mp + me) / (mcoll + ATOMN*mp + me));   /*reduced mass of colliding system in g*/

	for( l = 0; l <= n-1; l++ ){
		dnl = 1.0*z*z/(ATOMZ*ATOMZ)*6*n*n*(n*n - l*l - l - 1);
		pc = 1.678 + log10((float)te/NE); //error in Brock71

		if((n != 1) && (n != 2 || l != 0)){    /*lt = 0 for 1s and 2s levels*/
			lt = 1.0/anl[l];
			pc1 = 10.95 + log10(me*te*lt*lt/redmass);
			if(pc1 < pc){
				pc = pc1;
			}
		}
		c[ l ] = 9.933*pow(10, -6)*sqrt(redmass/me)*dnl/sqrt(te)/*NE*/*(11.538 + log10(te*me/(dnl*redmass)) + pc);

	}

	cu[0] = c[0];
	cl[0] = 0;

	for( l = 1; l < n; l++){
		cl[l] = (2.0*(l - 1) + 1.0)/(2*l + 1)*cu[l-1];
		cu[l] = c[l] - cl[l];
	}

	cl1[n-1] = c[n-1];
	cu1[n-1] = 0;

	for( l = n-2; l >= 0; l--){
		cu1[l] = (2.0*(l + 1) + 1.0)/(2*l + 1)*cl1[l+1];
		cl1[l] = c[l] - cu1[l];

	}

	cl1[0] = 0;
	cu[n-1] = 0;
	for( l = 0; l < n; l++){
		cu[l] = (cu[l] + cu1[l]) /2;
		/*cl[l] = (cl[l] + cl1[l]) /2;*/
		cl[l+1] = (2.0*l + 1.0)/(2.0*l + 3.0)*cu[l];
	}
	}
	free(c);
	free(cu1);
	free(cl1);
}



/*calculates the natural log of the bound-free matrix element using method of
Burgess 1965MmRAS..69....1B
output is Theta(n, l; k l') of Burgess*/
void gullc ( const int nn, double gu[], double gl[], double k )
{
	int l, j;
	double pi, g0, fn, n2k2, f1, f2, cu, cl,test;

	k = sqrt(k);
	pi = 2 * acos( 0 );
	f1 = fact( 2*nn - 1 );
	f2 = f1 - log( 2*nn - 1);

	for(l = 0; l < nn; l++){
			gu[l] = 0;
			gl[l] = 0;
	}

	n2k2 = 1 + nn*nn*k*k;

	if ( k == 0 ){
		g0 = 0.5*log( pi / 2 ) + nn*log( 4*nn ) - 2*nn;
		/*gu[l = n-1] = G(n, n-1; 0, n)*/
		gu[ nn-1 ] = g0 - f1 + log( 8*nn );
	}
	else{
		g0 = 0.5*log( pi / 2 ) - 0.5*log( 1 - exp( (-2*pi) / (k) ) ) + nn*log( 4*nn ) -
		(nn + 1)*log( n2k2 ) - ( 2*atan( nn*k ) ) / (k);
		/*gu[l = n-1] = G(n, n-1; k, n)*/
		gu[ nn-1 ] = g0 + log( 8*nn ) - f1 - log( n2k2 );
	}
	if (nn == 1) {
		gl[0] = 0;
	}
	else{
		/*gl[l = n-1] = G(n, n-1; k, n -2)*/
		gl[ nn-1 ] = g0 + log( 4 ) - f1;
		if (nn >= 2){
			/*gu[l = n-2] = G(n, n-2; k, n-1)*/
			gu[ nn-2 ] = g0 + log( 8*nn*nn ) - f2;
			/*gl[l = n-2] = G(n, n-2; k, n-3)*/
				if ( k == 0){
				gl[ nn-2 ] = g0 + log( 4*(nn + 3 )) - f2;
			}
			else{
				gl[ nn-2 ] = g0 + log( 4 + ( nn -1 ) * ( n2k2 ) ) - f2 + log(4);
			}
		}
		if (nn > 2) {
			for (l = nn-3; l >= 1; --l){
				fn = gu[ l + 2 ] - gu[ l + 1 ];
				/*eq (36) in Burgess with l -> l+2*/
				gu[ l ] = gu[ l + 1 ] + log( (4*nn*nn - 4*(l + 2)*(l + 2) + (l + 2)*(2*l + 3)*n2k2)
						- 4.0*nn*nn*(nn*nn - (l + 2)*(l + 2))*(1 + (l + 3)*(l + 3)*k*k)*exp(fn));
					fn = gl[ l + 2 ] - gl[ l + 1 ];
				/*eq (37) in Burgess with l -> l+1*/
				gl[ l ] = gl[ l + 1 ] + log((4*nn*nn - 4*(l + 1)*(l + 1) + (l + 1)*(2*l + 3)*n2k2)
						- 4.0*nn*nn*(nn*nn - (l + 2)*(l + 2))*(1 + (l + 1)*(l + 1)*k*k)*exp(fn));
			}

			fn = gu[ 2 ] - gu[ 1 ];
			gu[ 0 ] = gu[ 1 ] + log( ( 4*nn*nn - 16 + 6*n2k2 )
					- 4.0*nn*nn*( nn*nn -4 )*( 1 + 9*k*k )*exp( fn ) );
			gl[ 0 ] = 0;
		}
	}
	/*convert G -> g*/
	cu = 0.5*(fact( nn ) - fact( nn - 1 ) +log( 1 + k*k)) - nn*log( 2*nn );
	gu[ 0 ] =  cu + gu[ 0 ];

	if ( nn > 1 ){
		cu = cu + 0.5*(log( nn + 1 ) + log( nn - 1 ) + log( 1 + 4*k*k)) + log( 2*nn );
		gu[ 1 ] = cu + gu[ 1 ];
		cl = 0.5*(fact( nn + 1 ) - fact( nn - 2 )) + (1 - nn)*log( 2*nn );
		gl[ 1 ] = cl + gl[ 1 ];

		if (nn > 2){
			for ( l = 2 ; l <= nn - 1; l++ ){
				cu = cu + 0.5*(log( nn + l ) + log( nn - l ) + log( 1 + ( l + 1 )*( l + 1)*k*k)) + log( 2*nn );
				gu[ l ] = cu + gu[ l ];
				cl = cl + 0.5*(log( nn + l ) + log( nn - l ) + log( 1 + ( l - 1 )*( l - 1)*k*k)) + log( 2*nn );
				gl[ l ] = cl + gl[ l ];
			}
		}
	}

	/*computes Theta*/
	for (l = 0; l <= nn - 1; l++){
		gu[ l ] = n2k2*ATOMZ*ATOMZ*(exp( gu[ l ] ))*(exp( gu[ l ] ));
		gl[ l ] = n2k2*ATOMZ*ATOMZ*(exp( gl[ l ] ))*(exp( gl[ l ] ));
	}

	gl[ 0 ] = 0;

}





/*calculates radiative recombination coefficient in cgs
See Eq 9 - 12 of 1965MmRAS..69....1B*/
/*uses Gaussian integration*/
void radrecomb( const int n, const double y, double thu[], double thl[], double rad[] )
{
	int l, i, j,v;
	double h, a1, a2, p, q, r, n2k2, fl, fu, rad0, rad1, ratio;
	double x[5] = {-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664};
	double w[5] = {0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189};

	for(l = 0; l < n; l++)
		rad[ l ] = 0;


	h = 0.0001/n;
	a1 = 0;
	a2 = h;
	p = (a2 - a1)/2;
	q = (a1 + a2)/2;
	r = p*x[0] + q;	/*r = kappa^2 (Burgess)*/

	gullc( n, thu, thl, r);

	n2k2 = (1 + n*n*r)*(1 + n*n*r);
	fu = y*n2k2*thu[0]*exp(-r*y); /*Burgess (10) with l = 0*/
	rad0 = p*w[0]*fu*sqrt(y)*ATOMZ;
	a2 = 0;

	do{
		for(j = 1; j <= 4; j++){  /*change size of h after 4 iterations*/
			a1 = a2;
			a2 = a2 + h;
			p = (a2 - a1)/2;
			q = (a1 + a2)/2;

			for(i = 0; i <=4; i++){
				r = p*x[i] + q;
				gullc( n, thu, thl, r);
				n2k2 = (1 + n*n*r)*(1 + n*n*r);

				for(l = n-1; l >= 0; l--){
					fl = (l)*y*n2k2*thl[l]*exp(-r*y);
					fu = (l+1)*y*n2k2*thu[l]*exp(-r*y);
					rad1 = p*w[i]*(fl + fu)*sqrt(y)*ATOMZ;
					rad[l] = rad[l] + rad1;
				}
			}
		}
		ratio = rad1/rad0;
		h = h*2;
	} while (ratio > 0.000001);

	for(l = 0; l < n; l++){
		rad[ l ] = 5.62597709*pow(10, -15)*rad[l]/(n*n);
	}
}

/*calculates photoionization coefficient in cgs
See Eq 1 - 3 of 1965MmRAS..69....1B*/
/*uses Gaussian integration*/
void photoion( const int n, double thu[], double thl[], double phinl[], double const ryd, const double dw)
{
	int l, i, j;
	double h, a1, a2, p, q, r, fl, fu, ph0, ph1, ratio, mintensity;
	double x[5] = {-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664};
	double w[5] = {0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189};

	for(l = 0; l < n; l++)
		phinl[ l ] = 0;

	if(dw != 0){
	h = 0.0001/n;
	a1 = 0;
	a2 = h;
	p = (a2 - a1)/2;
	q = (a1 + a2)/2;
	r = p*x[0] + q;	/*r = kappa^2 (Burgess)*/
	gullc( n, thu, thl, r);
	mintensity = mintens(n, r, ryd, dw);
	fu = mintensity*thu[ 0 ];
	ph0 = p*w[0]*fu;
	a2 = 0;

	do{
		for(j = 1; j <= 4; j++){  /*change size of h after 4 iterations*/
			a1 = a2;
			a2 = a2 + h;
			p = (a2 - a1)/2;
			q = (a1 + a2)/2;

			for(i = 0; i <=4; i++){
				r = p*x[i] + q;
				gullc( n, thu, thl, r);
				mintensity = mintens(n, r, ryd, dw);
				for(l = n-1; l >= 0; l--){
					fl = l*mintensity*thl[ l ];
					fu = (l+1)*mintensity*thu[ l ];
					ph1 = p*w[i]*(fl + fu);
					phinl[l] = phinl[l] + ph1;
				}
			}
		}
		ratio = ph1/ph0;
		h = h*2;
	} while (ratio > 0.000001);

	for(l = 0; l < n; l++){
		phinl[ l ] = 2.56611968*pow(10, -8)*n*n*ryd*phinl[ l ]/(2*l + 1);
	}
	}
}

void stimrecomb( const int n, const double y, const int te,  double thu[], double thl[], double stimnl[], double const ryd, const double dw )
{
	int l, i, j;
	double h, a1, a2, p, q, r, fl, fu, st0, st1, ratio, mintensity;
	double x[5] = {-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664};
	double w[5] = {0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189};

	for(l = 0; l < n; l++)
		stimnl[ l ] = 0;

	if(dw != 0){

	h = 0.0001/n;
	a1 = 0;
	a2 = h;
	p = (a2 - a1)/2;
	q = (a1 + a2)/2;
	r = p*x[0] + q;	/*r = kappa^2 (Burgess)*/

	gullc( n, thu, thl, r);
	mintensity = mintens(n, r, ryd, dw);
	fu = mintensity*exp(-r*y)*thu[ 0 ];
	st0 = p*w[0]*fu;
	a2 = 0;

	do{
		for(j = 1; j <= 4; j++){  /*change size of h after 4 iterations*/
			a1 = a2;
			a2 = a2 + h;
			p = (a2 - a1)/2;
			q = (a1 + a2)/2;

			for(i = 0; i <=4; i++){
				r = p*x[i] + q;
				gullc( n, thu, thl, r);
				mintensity = mintens(n, r, ryd, dw);
				for(l = n-1; l >= 0; l--){
					fl = l*mintensity*exp(-r*y)*thl[ l ];
					fu = (l+1)*mintensity*exp(-r*y)*thu[ l ];
					st1 = p*w[i]*(fl + fu);
					stimnl[l] = stimnl[l] + st1;
				}
			}
		}
		ratio = st1/st0;
		h = h*2;
	} while (ratio > 0.000001);

	for(l = 0; l < n; l++){
		stimnl[ l ] = 1.06271544*pow(10, -23)*ryd*sqrt( pow(te, -3) )*n*n*stimnl[ l ];
	}
	}
}


void matinv(double a[][SIZE], double x[][SIZE], int depop)
/* Function to invert matrix a[][] with the inverse stored
   in x[][] in the output.  Copyright (c) Tao Pang 2001. */
{
  int i,j,k;
  int indx[SIZE];
  double b[SIZE][SIZE], c[SIZE][SIZE];
  printf("enter matinv\n");
  if(b == NULL || c == NULL/*  || nlmat == NULL  */){
		printf("OUT OF MEMORY!\n");
		exit(1);
	}

	for(i=0+depop;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			c[i][j] = a[i][j];
		}
	}

  for (i = 0+depop; i < SIZE; ++i)
  {
    for (j = 0+depop; j < SIZE; ++j)
    {
      b[i][j] = 0;
    }
  }
  for (i = 0+depop; i < SIZE; ++i)
  {
    b[i][i] = 1;
  }

  elgs (c, indx, depop);

  for (i = 0+depop; i < SIZE-1; ++i)
  {
    for (j = i+1; j < SIZE; ++j)
    {
      for (k = 0+depop; k < SIZE; ++k)
      {
        b[indx[j]][k] = b[indx[j]][k]-c[indx[j]][i]*b[indx[i]][k];
      }
    }
  }

for (i = 0+depop; i < SIZE; ++i)
  {
    x[SIZE-1][i] = b[indx[SIZE-1]][i]/c[indx[SIZE-1]][SIZE-1];
    for (j = SIZE-2; j >= 0+depop; j = j-1)
    {
      x[j][i] = b[indx[j]][i];
      for (k = j+1; k < SIZE; ++k)
      {
        x[j][i] = x[j][i]-c[indx[j]][k]*x[k][i];
      }
	   x[j][i] = x[j][i]/c[indx[j]][j];
    }
  }

}


void elgs (double a[][SIZE], int indx[], int depop)
/* Function to perform the partial-pivoting Gaussian elimination.
   a[][] is the original matrix in the input and transformed
   matrix plus the pivoting element ratios below the diagonal
   in the output.  indx[] records the pivoting order.
   Copyright (c) Tao Pang 2001. */

{
  int i, j, k, itmp;
  double c1, pi, pi1, pj;
  double c[SIZE];

/* Initialize the index */
	printf("enter elgs\n");
 for (i = 0+depop; i < SIZE; ++i)
  {
    indx[i] = i;
  }

/* Find the rescaling factors, one from each row */

 for (i = 0+depop; i < SIZE; ++i)
  {
    c1 = 0;
    for (j = 0+depop; j < SIZE; ++j)
    {
      if (fabs(a[i][j]) > c1) c1 = fabs(a[i][j]);
    }
    c[i] = c1;
  }

/* Search the pivoting (largest) element from each column */

 for (j = 0+depop; j < SIZE-1; ++j)
  {
    pi1 = 0;
    for (i = j; i < SIZE; ++i)
    {
      pi = fabs(a[indx[i]][j])/c[indx[i]];
      if (pi > pi1)
      {
        pi1 = pi;
        k = i;
      }
    }

/* Interchange the rows via indx[] to record pivoting order */

   itmp = indx[j];
    indx[j] = indx[k];
    indx[k] = itmp;
    for (i = j+1; i < SIZE; ++i)
    {
      pj = a[indx[i]][j]/a[indx[j]][j];

/* Record pivoting ratios below the diagonal */

      a[indx[i]][j] = pj;

/* Modify other elements accordingly */

      for (k = j+1; k < SIZE; ++k)
      {
        a[indx[i]][k] = a[indx[i]][k]-pj*a[indx[j]][k];
      }
    }
  }
printf("exit elgs\n");
}

/*multiplies two SIZAExSIZE square matrices with each other*/
void matmultn(double first[][SIZE],double second[][SIZE], double multiply[][SIZE], int depop)
{
	int c,d,k;
	double sum;

    for (c = 0+depop; c < SIZE; c++) {
      for (d = 0+depop; d < SIZE; d++) {
        for (k = 0+depop; k < SIZE; k++) {
          sum = sum + first[c][k]*second[k][d];
        }

        multiply[c][d] = sum;

        sum = 0;
      }

    }


}


/*Multiplies n row vector vec with nxn matrix mat
output is row vector out*/
void matmult(double mat[][SIZE], double vec[], double out[])
{
	int i, j, k;

	for(j=0; j < SIZE; j++){
		out[j] = 0;
	}

	for(j=0; j < SIZE; j++){
		for(k=0; k < SIZE; k++){
			out[j] = out[j] + vec[k]*mat[k][j];
		}
	}
}


/*performs gaussian elimination with partial pivotting of cmat*cbn = veca. The coefficient matrix should be in column form, i.e. veca is a row vector*/
void gausselim(double tmat[][SIZE], double vect[], double cbn[], int del)
{
	int   p,q,r,s,rowx;
	double xfac,temp2,temp1,amax/*, tmat[SIZE][SIZE], vect[SIZE]*/;

	/*reduce tmat and interchange entries in vect so that cmat and veca remain unchanged*/
	/* for(r=0;r<SIZE;r++){
		for(p=0;p<SIZE;p++){
			tmat[r][p] = cmat[r][p];
		}
		vect[r] = veca[r];
	} */

	/* forward reduction*/

	rowx = 0;   /* Keep count of the row interchanges */
	for (r=0+del; r<SIZE-1; ++r){
		amax = (double) fabs(tmat[r][r]) ;
		s = r;
		for (p=r+1; p<SIZE; p++){   /* Find the row with largest pivot */
			xfac = (double) fabs(tmat[r][p]);
			if(xfac > amax) {amax = xfac; s=p;}
		}
		if(s != r){  /* Row interchanges */
			rowx = rowx+1;
			temp1 = vect[r];
			vect[r]  = vect[s];
			vect[s]  = temp1;
			for(q=r; q<SIZE; q++){
				temp2 = tmat[q][r];
				tmat[q][r] = tmat[q][s];
				tmat[q][s] = temp2;
			}
		}
		for (p=r+1; p<SIZE; ++p) {
			xfac = tmat[r][p]/tmat[r][r];
			for (q=r+1; q<SIZE; ++q) {
				tmat[q][p] = tmat[q][p]-xfac*tmat[q][r];
			}
			vect[p] = vect[p]-xfac*vect[r];
		}
	}

	/* back substitution*/

	for (q=1; q<=SIZE-del; ++q) {
		r=SIZE-q;
		cbn[r] = vect[r];
		for(p=r+1; p<SIZE; ++p) {
			cbn[r] = cbn[r]-tmat[p][r]*cbn[p];
		}
		cbn[r] = cbn[r]/tmat[r][r];
	}

}

/*performs iterative refinement procedure on system matrix*cbn = vec where cbn is the first solution and vec and cbn are rowvectors*/
void iterativerefinement(double matrix[][SIZE], double vec[], double cbn[], int del)
{
	int kk,ii,jj,pp, maxout = 1, maxit = 100000;
	double r[SIZE]={0}, y[SIZE], xx[SIZE], jsum, cond,maxx, tol = 0.000000001;
	FILE *irPtr;

	kk=1;
	while(kk<=maxit){
		for(ii=0+del;ii<SIZE;ii++){
			jsum = 0;
			for(jj=0+del;jj<SIZE;jj++){
				jsum = jsum + matrix[jj][ii]*cbn[jj];
			}
			r[ii] = vec[ii] - jsum;
		}

		gausselim(matrix, r, y, del);
		for(ii=0+del;ii<SIZE;ii++){
			xx[ii] = cbn[ii] + y[ii];
		}
		if(kk == 1){
			cond = infnormvec(y,SIZE,del)/infnormvec(cbn,SIZE,del)*pow(10,8);
			printf("condition number = %e\n",cond);
		}
		maxx=0;
		for(pp=0+del;pp<SIZE;pp++){
			if( (double) fabs((cbn[pp]-xx[pp])/xx[pp]) > maxx){
				maxx = (double) fabs(cbn[pp]-xx[pp]);
			}
		}
		if(maxx != 0){
			for(ii=0+del;ii<SIZE;ii++){
				cbn[ii] = xx[ii];
			}
		}
		/*if(maxx < tol){
			if(maxx != 0){printf("number of iterations %d, last max %e\n",kk,maxx);
			}else{printf(" number of iterations %d, last max %e\n",kk-1, maxx);}
			kk = maxit;
			maxout = 0;
		}*/
		kk= kk + 1;
	}
	if(maxout != 0){
		printf("maximum iterations reached %d, last max %e\n",kk-1, maxx);
	}

		irPtr = fopen("bnvalrefined.dat","w"); /* open for writing */

	for ( ii = NMIN; ii <= SIZE; ii++){

		fprintf(irPtr,"%.10e\n", cbn[ii-1]);
	}
	fclose(irPtr);


}



/*Finds the 1-norm of a matrix*/
double matnorm(  double matrix[][SIZE])
{
	int i,j;
	double maxx = 0, colsum;

	for(i=0;i<SIZE;i++){
		colsum = 0;
		for(j=0;j<SIZE;j++){
			colsum = colsum + (double) fabs(matrix[j][i]);
		}
		if(colsum > maxx){
			maxx = colsum;
		}
	}

	return maxx;
}


/*calculates Euclidian norm of vector*/
double vecnorm(const double vec[], const int size)
{
	int i;
	double norm=0;

	for(i=0;i<size;i++){
		norm = norm + vec[i]*vec[i];
	}

	return sqrt(norm);
}

int idx(const int n, const int l)
{
	return n*(n-1)/2.0 + l;
}

int idx2(const int n, const int l,const int min)
{
	return n*(n-1)/2.0 + l + 1 - min;
}



/*calculates the 1-norm of a vector of size "size" */
double normvec1(double vec[], int size, int del)
{
	int i;
	double sum = 0;

	for(i=0+del;i<size;i++){
		sum = sum + (double) fabs(vec[i]);
	}

	return sum;
}


/*calculates the largest absolute value component of a vector of size "s"*/
double infnormvec(double vec[], int s, int del)
{
	int i;
	double max = 0;

	for(i=0+del;i<s;i++){
		if((double) fabs(vec[i]) > max){
			max = (double) fabs(vec[i]);
		}
	}

	return max;
}

void stepgenold(int step[])
{
	int i = 0, m, n, l = 0;

	step[0] = 0;

	/* for(n = 1; n <= SIZE; n ++){
		step[n] = n;
	}  */

	for(n = 1; n <= 9; n += 2){
		for(m = 1; m <= 20; m ++){
			i = i+1;
			l = l + n;
			step[i] = l;
		}
		l = l-1;
	}
}

//creates full path to file if folder path and file name is sent.
char *makepath(const char path[], const char filename[])
{	char *fullpath = malloc(100);

	strcpy(fullpath, path);
	strcat(fullpath,filename);

	return fullpath;
}
