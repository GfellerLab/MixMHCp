/*
############
# Written by David Gfeller
#
# For any question, please contact david.gfeller@unil.ch
#
# To cite MixMHCp, please refer to Bassani-Sternberg M and Gfeller D*, J. Immunol. (2016)
#
# MixMHCp can be used freely by academic groups for non-commercial purposes (see license).
# The product is provided free of charge, and, therefore, on an "as is"
# basis, without warranty of any kind.
#
# FOR-PROFIT USERS
# If you plan to use MixMHCp (version 1.0) in any for-profit
# application, you are required to obtain a separate  license.
# To do so, please contact eauffarth@licr.org or lfoit@licr.org at the Ludwig Institute for  Cancer Research Ltd.
#
# Copyright (2016) David Gfeller
############
*/

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

using namespace std;

//Functions

void init_parameters(int a, int b, int c);

void import_alignment(char * alignment_dir);

double* prob(int k, int *p1);

void normalize(double *pr1, int le);
void normalize_pseudo(double *pr1);

void initialize_comp();

void random_comp(int ncl, int t);

double absv(double a);

void EMsteps();
void initialize_EM(double ***EM_pwm, double *wcl, int ncl);
void expectation(double ***EM_pwm, double *wcl, double **resp, int ncl);
void maximization(double **eprior, double **resp, double ***EM_pwm, double *wcl, double min_error, int ncl);
double loglikelihood(double **eprior, double ***EM_pwm, double *wcl, double *LL, double *sLL, int ncl);
void compute_eprior(double **eprior, int ncl);

void print_EM_pwm(int ncl, double ***EM_pwm, double *wcl);

double test_KLD(double **resp, int ncl);

int position(char s);
void remove_columns(int **tpeptide, int tnaa);

void print_responsibility(int ncl, double **resp);

void import_bias(char * out_dir, int bs);

void make_cluster_pwm(int si, int ncl, double ***m, double **resp);

void best_ncl(int ncl, double *KLD, double ****EM_pwm, double **wcl);
    
//Global variables

int N;    //Alphabet size
int naa; //number of residues (column) in the alignment for each peptide
int **peptide; // peptide composition [domain][position][peptide] (after filtering some columns)
int tnaa;
int **tpeptide; //initial peptides
int kp;   //number of peptides associated with the input sample
char *letter;
int *best_comp;
int sbest_comp;  //number of PWMs chosen by the user
int kpmax;    //size of the largest set of phage peptides
int ***comp_pep;         //peptides in each cluster
int *cl_size;           //size of each cluster
double fsm;
int ncl_max;  //Maximum number of PWMs
int naa_max;
fstream afile;
int decide_comp_number;   //1: the number of PWMs are fixed by the user. 0: The number of PWMs is to be found by the algorithm
char * out_dir;
double *bias;
double pseudo_count;
double pseudo_count_prior;
char * alphabet;

/*
  Run with:
  ./MixMHCp.x 0 0 5 -d ./output/ -b U
*/

int main(int argc, char ** argv)
{
    if (argc < 5) {
	cout << "Invalid arguments. Run through MixMHCp." << endl;
	exit(2);
    }
    char * alignment_dir = new char[4096];
    out_dir = new char[4096];
    alphabet = new char[4096];
    
    int bs;
     
    for (int i=4; i<argc; i+=2) {
	
	if (strcmp(argv[i], "-d") == 0) {
	    //loads the param -d, out_dir
	    strcpy(out_dir, argv[i+1]);
	    sprintf(alignment_dir, "%s/data", out_dir);
	}
 	else if (strcmp(argv[i], "-b") == 0) {
	    //loads the param -b, background frequencies
	    bs=atoi(argv[i+1]);
	}
 	else if (strcmp(argv[i], "-a") == 0) {
	    //loads the param -b, background frequencies
	    strcpy(alphabet, argv[i+1]);
	}
    }
 
    cout<<"MixMHCp input command line:\n";
    for(int i=0; i<argc; i++)
	cout<<argv[i]<<" ";
    cout<<endl<<endl;

    init_parameters(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    import_alignment(alignment_dir);
    import_bias(out_dir, bs);
    
    initialize_comp();
    
    EMsteps();

    return(0);

}



void init_parameters(int a, int b, int c)
{

    //Peptide
    if( a==0 ){
	N=strlen(alphabet);
	letter=new char[N+1];
	strcpy(letter, alphabet);
    }
    
    if(b==0){
	decide_comp_number=0; //1: the number of PWMs are fixed by the user. 0: The number of PWMs is to be found by the algorithm
	if(c==0){
	    ncl_max=10; //Maximal number of clusters
	}
	if(c>0){
	    ncl_max=c;
	}
    }
    if(b>0){
	decide_comp_number=b;
	ncl_max=b; 
    }

    //This is the pseudo_count used in computing KLD (same values as in GibbsCluster).
    pseudo_count=200;
    
    //This parameter determines the strength of the prior. It is a multiplicative factor for the random count added to the frequency
    //Very small random counts (0.1) seem to perform better (note that to compare with pseudo_count, there is a factor 20 that should be used).
    //We could also load it as input... (next version).
    pseudo_count_prior=0.1;
    
}

void EMsteps()
{

    //initialize the PWMs and the mixing coefficients
    double ***EM_pwm;
    double *wcl;
    double **eprior;

    //Final values for the Multiple PWMs
    EM_pwm=new double**[ncl_max];
    wcl=new double[ncl_max];
    eprior=new double*[ncl_max];
    for(int n=0; n<ncl_max; n++) {
	wcl[n]=0;
	EM_pwm[n]=new double*[naa];
	eprior[n]=new double[naa];
	for(int s=0; s<naa; s++) {
	    EM_pwm[n][s]=new double[N];
	    eprior[n][s]=0;
	    for(int j=0; j<N; j++)
		EM_pwm[n][s][j]=0;
	}
    }
    
    double ****full_EM_pwm;
    double **full_wcl;
    
    full_wcl=new double*[ncl_max];
    full_EM_pwm=new double***[ncl_max];
    
    for(int tn=0; tn<ncl_max; tn++) {
	full_EM_pwm[tn]=new double**[tn+1];
	full_wcl[tn]=new double[tn+1];
	for(int n=0; n<tn+1; n++) {
	    full_wcl[tn][n]=0;
	    full_EM_pwm[tn][n]=new double*[naa];
	    for(int s=0; s<naa; s++) {
		full_EM_pwm[tn][n][s]=new double[N];
		for(int j=0; j<N; j++)
		    full_EM_pwm[tn][n][s][j]=0;
	    }
	}
    }

    //temporary values for multiple PWMs (different optimization runs)
    double ***tEM_pwm;      //[component][position][aa]
    double *twcl;
    tEM_pwm=new double**[ncl_max];
    twcl=new double[ncl_max];
    for(int n=0; n<ncl_max; n++) {
	twcl[n]=0;
	tEM_pwm[n]=new double*[naa_max];
    }

    //hidden (latent) variables, useful for the EM
    double ***resp;
    resp=new double**[ncl_max];
    for(int ncl=0; ncl<ncl_max; ncl++){
	resp[ncl]=new double*[kpmax];
	for(int j=0; j<kpmax; j++) {
	    resp[ncl][j]=new double[ncl_max];
	}
    }
    double **tresp;
    tresp=new double*[kpmax];
    for(int j=0; j<kpmax; j++) {
	tresp[j]=new double[ncl_max];
    }
    
    double *LL;
    LL=new double[ncl_max];
    double *sLL;
    sLL=new double[ncl_max];

    int rp=5;  //Number of optimization runs (starting from different initial configuration)

    int st;
    double error=0;
    double min_error=0.005;
    double LLerror;
    int ct;
    int comp_kp;
    int cond=1;
    double dis;
    int step;
    long double max_LL=-100000000000.0;
    long double max_sLL=-100000000000.0;


    int tt=0;
    int corr;

    double MI;
    double pMI;
    int size_corr=50;  //Maximal alignment size to test for MI p-values
    int ttncl;
    int LLct;
    double *KLD;
    int pncl=0;
    double mncl=0;
    int ncl=0;
    KLD=new double[ncl_max+1];
    for(int n=0; n<ncl_max+1; n++){
	KLD[n]=-10000;
    }
  
    char buffer [4096];
    FILE *F;

    int *sz;
    sz=new int[ncl];
    
    
    sprintf(buffer, "%s/project.txt", out_dir);
    afile.open(buffer, ios::out);
    afile<<"#ProjectFile\n";
    afile.close();
    sprintf(buffer, "%s/EM_project.txt", out_dir);
    afile.open(buffer, ios::out);
    afile<<"#ProjectFile\n";
    afile.close();
    	
    cout<<"Number of peptides: "<<kp<<endl;
    
    //Run the EM optimization, 
    
    //Test different number of components
    //Initialize EM_pwm with the values in pwm
    //Here, this differs a bit from the gibbsclustering, but we have to be careful and consistent with the random counts...
    
    for(int j=0; j<kp; j++) {
	resp[0][j][0]=1;
    }
    print_responsibility(1, resp[0]);
	    
    srand(1);
    random_comp( 1, 0);
    initialize_EM(EM_pwm, wcl, 1);
    compute_eprior(eprior, 1);
    loglikelihood(eprior, EM_pwm, wcl, LL, sLL, 1);
     
    print_EM_pwm(1, EM_pwm, wcl);
        
    KLD[1]=test_KLD(resp[0], 1);
        
    for(int ncl=2; ncl<=ncl_max; ncl++) {

	srand(ncl);	    
	max_LL=-10000000000.0;
	
	//Run 'rp' times the optimization
	for(int t=0; t<rp; t++) {
	    
	    random_comp(ncl, t);
	    initialize_EM(tEM_pwm, twcl, ncl);
	    compute_eprior(eprior, ncl);
	    
	    LLerror=1;
	    LL[ncl]=0;
	    
	    loglikelihood(eprior, tEM_pwm, twcl, LL, sLL, ncl);
	    LLct=0;

	    while(LLerror>min_error) {
		LLct++;
		//Compute the responsibilities
		expectation(tEM_pwm, twcl, tresp, ncl);
		
		//Compute the new model coefficients and the mixing coefficients using the responsibilities
		maximization(eprior, tresp, tEM_pwm, twcl, min_error, ncl);
		
		//Compute the new loglikelihood
		LLerror=loglikelihood(eprior, tEM_pwm, twcl, LL, sLL, ncl);
	    }

	    if(LL[ncl]>max_LL) {
		//If the logikelyhood is larger, then keep these values
		max_LL=LL[ncl];
		max_sLL=sLL[ncl];
		for(int n=0; n<ncl; n++) {
		    wcl[n]=twcl[n];
		    for(int s=0; s<naa; s++) {
			for(int j=0; j<N; j++)
			    EM_pwm[n][s][j]=tEM_pwm[n][s][j];
		    }
		}
		for(int j=0; j<kp; j++) {
		    for(int n=0; n<ncl; n++) {
			resp[ncl-1][j][n]=tresp[j][n];
		    }
		}
	    }
	}

	//Keep track of the PWMs (this is to compute the best number of motifs).
	for(int n=0; n<ncl; n++){
	    full_wcl[ncl-1][n]=wcl[n];
	    for(int s=0; s<naa; s++) {
		for(int j=0; j<N; j++)
		    full_EM_pwm[ncl-1][n][s][j]=EM_pwm[n][s][j];
	    }
	}
	
	KLD[ncl]=test_KLD(resp[ncl-1], ncl);
	print_responsibility(ncl, resp[ncl-1]);
	print_EM_pwm(ncl, EM_pwm, wcl);
	    
    }
	
    sprintf(buffer, "%s/KLD/KLD.txt", out_dir);
    F=fopen(buffer, "w");
    fprintf(F, "1\t%.6f\n", KLD[1]);
    for(int n=2; n<=ncl_max; n++){
	fprintf(F, "%d\t%.6f\n", n, KLD[n]);
    }
    fclose(F);

    best_ncl(ncl_max, KLD, full_EM_pwm, full_wcl);
	
}

void best_ncl(int ncl_max, double *KLD, double ****full_EM_pwm, double **full_wcl){

    int ncl_final;

    double thresh1=100/(1.0*N*N*naa);
    double thresh2=200/(1.0*N*N*naa);
    
    //Find the max_KLD
    double Max_KLD=KLD[0];
    int Max_KLD_pos=0;
    
    for(int n=1; n<=ncl_max; n++){
	if(KLD[n]>Max_KLD){
	    Max_KLD=KLD[n];
	    Max_KLD_pos=n;
	}
    }
    ncl_final=Max_KLD_pos;

    double min;
    double mm;
    double d;
    double *lst_b;
    lst_b=new double[ncl_max];
    int *lst_bp;
    lst_bp=new int[ncl_max];
    int cond;
    int ct;
    int pmm;
    double lst_b_max;
    int lst_b_max_pos;
    double lst_b_min;
    int lst_b_min_pos;
    cond=1;
    int minp;
    
    for(int ncl=Max_KLD_pos+1; ncl<=ncl_max && cond==1; ncl++){

	for(int n2=0; n2<ncl; n2++){
	    min=10000;
	    //Find the most similar logo 
	    for(int n1=0; n1<ncl-1; n1++){
		d=0;
		for(int s=0; s<naa; s++) {
		    for(int j=0; j<N; j++){
			d=d+(full_EM_pwm[ncl-2][n1][s][j]-full_EM_pwm[ncl-1][n2][s][j])*(full_EM_pwm[ncl-2][n1][s][j]-full_EM_pwm[ncl-1][n2][s][j]);
		    }
		}
		d = d/(1.0*naa);
		if (d<min) {
		    min=d;
		    minp=n1;
		}
	    }

	    lst_b[n2]=min;
	    lst_bp[n2]=minp;
	}

	//Find the maximum on lst_b
	lst_b_max=-100000;
	lst_b_min=100000;
	for(int n=0; n<ncl; n++){
	    if(lst_b_max < lst_b[n]) {
		lst_b_max=lst_b[n];
		lst_b_max_pos=n;
	    } if(lst_b_min > lst_b[n]) {
		lst_b_min=lst_b[n];
		lst_b_min_pos=n;
	    }
	}
	
	//This is comparing among the new clusters

	min=10000;
	for(int n2=0; n2<ncl; n2++){
	    for(int n1=0; n1<ncl; n1++){ //This is comparing among clusters
		if(n1 != n2){
		    d=0;
		    for(int s=0; s<naa; s++) {
			for(int j=0; j<N; j++){
			    d=d+(full_EM_pwm[ncl-1][n1][s][j]-full_EM_pwm[ncl-1][n2][s][j])*(full_EM_pwm[ncl-1][n1][s][j]-full_EM_pwm[ncl-1][n2][s][j]);
			}
		    }
		    d = d/(1.0*naa);
		    if (d<min) {
			min=d;
		    }
		}
	    }
	}
       	
	//Make sure there is no redundancy (check that the smallest distance is larger than Thresh1
	cond=1;
	if(min < thresh2){
	    cond=0;
	    //cout<<"Redundant motifs: "<<min<<endl;
	}

	//Make sure at least one new motif is present
	//This is an issue, for instance if a cluster is split in two equal parts, because then the two 'new' motifs are still quite similar to the one before
	if(ncl>1){
	    if(lst_b_max < thresh2){
		cond=0;
		//cout<<"No new motifs: "<<lst_b_max<<endl;
	    }
	}
	
	//Make sure the other motifs have been seen before
	for(int n2=0; n2<ncl; n2++){
	    if(lst_b[n2] > thresh1 && n2 != lst_b_max_pos){
		cond=0;
		//cout<<"Too many new motifs: "<<n2<<" "<<lst_b[n2]<<" "<<thresh1<<" "<<thresh2<<endl;
	    }
	}
	//Make sure the new motif is supported by enough data
	if(full_wcl[ncl-1][lst_b_max_pos] < 0.03 || full_wcl[ncl-1][lst_b_max_pos]*kp < 85){
	    cond=0;
	}
	
	if(cond==1){
	    ncl_final=ncl;
	}
    }

    char buffer [4096];
    FILE *F;
    sprintf(buffer, "%s/KLD/best_ncl.txt", out_dir);
    F=fopen(buffer, "w");
    fprintf(F, "Final number of motifs:\t%d\n", ncl_final);
    fclose(F);
    
}

//Compute the Kullback-Leibler distance
//So far we take the responsibilities to decide for a hard clustering
//We could think of another variant of the KLD that only compares with a random model (but similar to maximum likelihood, and likely to favor many small clusters...)
double test_KLD(double **resp, int ncl){

    double KLD=0;
    double max_KLD=0;
    int cl;
    double max_resp;
    double *Sscore;
    Sscore=new double[ncl];

    double lambda=0.8;
    double sigma;

    sigma=10;

    int *sz;
    sz=new int[ncl];
    for(int n=0; n<ncl; n++) {
	sz[n]=0;
    }
    int *gr;
    gr=new int[kp];
        
    //Compute the size of each cluster
    for(int j=0; j<kp; j++) {
	max_resp=0;
	cl=-1;
	for(int n=0; n<ncl; n++) {
	    if(resp[j][n]>max_resp){
		max_resp=resp[j][n];
		cl=n;
	    }
	}
	sz[cl]++;
	gr[j]=cl;
    }

    //Build the global PWMs for each cluster without normalization
    double ***cluster_pwm;
    cluster_pwm=new double**[ncl];
    for(int n=0; n<ncl; n++) {
	cluster_pwm[n]=new double*[naa];
	for(int s=0; s<naa; s++){
	    cluster_pwm[n][s]=new double[N];
	}
    }
    make_cluster_pwm(-1, ncl, cluster_pwm, resp);
 
    double ***cl_pwm;
    cl_pwm=new double**[ncl];
    for(int n=0; n<ncl; n++) {
	cl_pwm[n]=new double*[naa];
	for(int s=0; s<naa; s++){
	    cl_pwm[n][s]=new double[N];
	}
    }
    
    for(int j=0; j<kp; j++) {
	//Find which cluster it belongs to
	max_resp=0;
	cl=gr[j];
	
	//We actually recompute the PWMs only for the cluster which the pepitde is assigned to (excluding the sequence of the peptide itself) and then compute the score.
	//This is to match better the definition in Andreatta et al., although large differences are not expected.

	for(int n=0; n<ncl; n++) {
	    for(int s=0; s<naa; s++){
		for(int j=0; j<N; j++){
		    cl_pwm[n][s][j]=cluster_pwm[n][s][j];
		}
	    }
	}
	
	for(int s=0; s<naa; s++){
	    if(peptide[s][j] != N){
		cl_pwm[cl][s][peptide[s][j]]--;
	    } else {
		for(int p=0; p<N; p++){
		    cl_pwm[cl][s][p]=cl_pwm[cl][s][p]-1.0/N;
		}
	    }
	}
	
	//Include the normalized counts
	for(int n=0; n<ncl; n++) {
	    //If the cluster is empty, give a flat PWM so that the score is 0.
	    if((n==cl && sz[n]==1) || (n != cl && sz[n]==0)){
		for(int s=0; s<naa; s++){
		    for(int p=0; p<N; p++)
			cl_pwm[n][s][p]=bias[p];
		}
	    } else {
		for(int s=0; s<naa; s++){
		    normalize_pseudo(cl_pwm[n][s]);
		}
	    }
	}

	//Compute the scores with each PWM
	for(int n=0; n<ncl; n++) {
	    Sscore[n]=0;
	    for(int s=0; s<naa; s++){
		if(peptide[s][j] != N){
		    Sscore[n]=Sscore[n]+2*log(cl_pwm[n][s][peptide[s][j]]/bias[peptide[s][j]])/log(2);		  
		}
	    }
	    if(n != cl){
		Sscore[n]=Sscore[n]*sz[n]/(sz[n]+sigma);
	    } else {
		Sscore[n]=Sscore[n]*(sz[n]-1)/(sz[n]-1+sigma);
	    }
	}

	//Compute the KLD score for all other clusters
	max_KLD=0;
	for(int n=0; n<ncl; n++) {
	    if(n != cl && max_KLD<Sscore[n]){
		max_KLD=Sscore[n];
	    }
	}
	KLD=KLD+Sscore[cl]-lambda*max_KLD;
    }
    KLD=1.0*KLD/kp; 
    return(KLD);
}

void normalize_pseudo(double *pr1)
{

    double  tot=0;
    for(int i=0; i<N; i++) {
	tot=tot+pr1[i];
    }
    
    //Include the pseudocounts as in Andreatta et al.
    double *g; g=new double[N];
    
    for(int i=0; i<N; i++){
	g[i]=0;
	for(int ii=0; ii<N; ii++){
	    g[i]=g[i]+1.0/N*pr1[ii]/tot;
	}
    }
    for(int i=0; i<N; i++){
	pr1[i]=(pr1[i]+pseudo_count*g[i])/(tot+pseudo_count);
    }
    
    
}

void make_cluster_pwm(int si, int ncl, double ***m, double **resp){

    //Build the pwms for each cluster
    //Do not include the node si itself
    //This is also working if there is only one cluster (all best_comp should be equal to 1).
    //cout<<i<<" "<<si<<" "<<ncl<<" "<<endl;

    for(int n=0; n<ncl; n++){
	for(int s=0; s<naa; s++){
	    for(int p=0; p<N; p++){
		m[n][s][p]=0;
	    }
	}
    }
    double max_resp;
    int cl;

    for(int j=0; j<kp; j++) {
	
	if(j != si){
	    //Find the group of peptide j based on the resp
	    max_resp=0;
	    cl=-1;
	    for(int n=0; n<ncl; n++) {
		if(resp[j][n]>max_resp){
		    max_resp=resp[j][n];
		    cl=n;
		}
	    }
	
	    for(int s=0; s<naa; s++){
		if(peptide[s][j] != N){
		    m[cl][s][peptide[s][j]]++;
		} else {
		    for(int p=0; p<N; p++){
			m[cl][s][p]=m[cl][s][p]+1.0/N;
		    }
		}
	    }
	}
    }
   
}


//Randomly attribute the peptide to different components (initialization of the EM-algorithm).
void initialize_comp()
{

    int comp_kp;

    
    cl_size=new int[ncl_max];
    best_comp=new int[kp];

    comp_pep=new int**[ncl_max];
    for(int n=0; n<ncl_max; n++) {
	comp_pep[n]=new int*[naa];   //[position][cluster_label]
	for(int s=0; s<naa; s++)
	    comp_pep[n][s]=new int[kp];
    }
}


void random_comp(int ncl, int t)
{

    int comp_kp;

    //This part includes random numbers, which changes depending on the platform...
    /*for(int j=0; j<kp; j++) {
      best_comp[j]=rand()%ncl;
      }
    */

    //This part does no depend on random numbers
    if(ncl<kp){
	for(int j=0; j<ncl; j++) {
	    best_comp[j]=j%ncl;
	}
	
	for(int j=ncl; j<kp; j++) {
	    best_comp[j]=int(j/(t+1))%ncl;
	}
    }
    else{
	for(int j=0; j<kp; j++) {
	    best_comp[j]=j%ncl;
	}
    }

    
    for(int n=0; n<ncl; n++)
	cl_size[n]=0;

    for(int j=0; j<kp; j++) {
	if(best_comp[j]>=0) {
	    cl_size[best_comp[j]]++;
	}
    }
    for(int n=0; n<ncl; n++) {
	comp_kp=0;
	for(int j=0; j<kp; j++) {
	    if(best_comp[j]==n) {
		for(int s=0; s<naa; s++) {
		    comp_pep[n][s][comp_kp]=peptide[s][j];
		}
		comp_kp++;
	    }
	}
    }
}


void normalize(double *pr1, int le)
{

    double  tot=0;
    for(int i=0; i<le; i++) {
	tot=tot+pr1[i];
    }
    if(tot>0) {
	for(int i=0; i<le; i++) {
	    pr1[i]=1.0*pr1[i]/tot;
	}
    }
}


//Compute the column of the pwm
double* prob(int k, int *p)
{
    double *f;
    f=new double[N];
    for(int i=0; i<N; i++) {
	f[i]=pseudo_count_prior;
    }

    //Compute the number of occurences of each amino acid in p[]
    for(int i=0; i<k; i++) {
	if(p[i]<N) {
	    f[p[i]]++;
	} else {
	    for(int j=0; j<N; j++) {
		f[j]=f[j]+1.0/N;
	    }
	}
    }

    normalize(f, N);

    return(f);
}

//return the absolute value of a number
double absv(double a)
{
    if(a<0)
	a=-a;
    return(a);
}


void import_bias(char * out_dir, int bs){
    
    char file [4096];
    int t;
    bias=new double[N];
    
    if(bs==0){
	for (int i = 0; i < N; i++) {
	    bias[i] = 1;
	}
	cout << "Flat background" << endl;
	
    } else if(bs==1){
	bias[0]=0.074;
	bias[1]=0.025;
	bias[2]=0.054;
	bias[3]=0.054;
	bias[4]=0.047;
	bias[5]=0.074;
	bias[6]=0.026;
	bias[7]=0.068;
	bias[8]=0.058;
	bias[9]=0.099;
	bias[10]=0.025;
	bias[11]=0.045;
	bias[12]=0.039;
	bias[13]=0.034;
	bias[14]=0.052;
	bias[15]=0.057;
	bias[16]=0.051;
	bias[17]=0.073;
	bias[18]=0.013;
	bias[19]=0.032;
	cout << "Uniprot background" << endl;
    } else if(bs==2){
	
	sprintf(file, "%s/bias.txt", out_dir);
	afile.open(file, ios::in);

	if (afile.is_open()) {
	    cout << "Importing residue biases from " << file << endl;
	    for (int i=0; i<N; i++) afile>>bias[i];      
	} else {
	    cout << "No bias file: Use flat background" << endl;
	    for (int i = 0; i < N; i++) bias[i] = 1;
	}
    }
    
    //Normalize the bias
    double T=0;
    for (int i = 0; i < N; i++)
	T=T+bias[i];
    for (int i = 0; i < N; i++)
	bias[i]=bias[i]/T;

    afile.close();
}

void import_alignment(char * alignment_dir)
{
    int ct=0;
    string line;
    string nada="";
    fstream afile1;

    cout << "Alignment folder:\n" << alignment_dir << endl;
  

    char * pch;
    int tot=0;
    char str [4096]; //Note: We increased this from 40 to accommodate DNA seqs
    int ln;
    char file [4096];
 
   
    kpmax=0;
    naa_max=0;

    sprintf(file, "%s/peptides.fa", alignment_dir);
	
 
    ifstream myfile (file);
    if (myfile.is_open()) {
	ct=0;
	while (! myfile.eof() ) {
	    getline (myfile,line);
	    if(line.compare(nada) != 0) {
		ct++;
	    }
	    //get the peptide length
	    if(ct==2) {
		tnaa=line.length();
		//cout << "peptide length tnaa = " << tnaa << endl;
	    }
	}
	myfile.close();
    } else {
	cout << "Unable to open file: " << file << endl;
	exit(2);
    }
    //Import the peptides
    kp=ct/2;
    
    tot=tot+kp;
    if(kp>kpmax)
	kpmax=kp;

    tpeptide=new int*[tnaa];
    for(int j=0; j<tnaa; j++)
	tpeptide[j]=new int[kp];

    afile1.open(file, ios::in);
    for(int j=0; j<kp; j++) {
	afile1>>str;
	afile1>>str; //Why is this duped?
	for(int s=0; s<tnaa; s++) {
	    tpeptide[s][j]=position(str[s]); //Use a nummerical representation for amino acids
	}
    }
    afile1.close();


    //Remove the positions with very low specificity on both sides of the motif
    remove_columns(tpeptide, tnaa);
  
}

void remove_columns(int **tpeptide, int tnaa)
{

    double *f;
    f=new double[N];
    double ent;
    int *p;
    p=new int[500];
    double ent_cutoff=0.99;  //This is the cut-off for the column to be removed

    int s1, s2, st;
    int stop;

    for(int s=0; s<tnaa; s++) {
	
	p[s]=1;
	for(int t=0; t<N; t++)
	    f[t]=0;
	
	for(int j=0; j<kp; j++) {
	    if(tpeptide[s][j] != N) {
		f[tpeptide[s][j]]++;
	    }
	    if(tpeptide[s][j] == N) {
		for(int t=0; t<N; t++)
		    f[t]=f[t]+1.0/N;
	    }
	}
	ent=0;
	for(int t=0; t<N; t++) {
	    if(f[t]>0) {
		f[t]=1.0*f[t]/kp;
		ent=ent-f[t]*log(f[t])/log(N);
	    }
	}
	if(ent>ent_cutoff) {
	    p[s]=0;
	}
	//p[s]=1;  //DG
    }
    stop=1;
    s1=0;
    //Start removing positions from N-terminal until a position with low enough entropy is found
    for(int s=0; s<tnaa && stop==1; s++) {
	if(p[s] == 0)
	    s1++;
	if(p[s] != 0)
	    stop=0;
    }
    stop=1;
    s2=tnaa-1;
    //Start removing positions from C-terminal until a position with low enough entropy is found
    for(int s=tnaa-1; s>0 && stop==1; s--) {
	if(p[s] == 0)
	    s2--;
	if(p[s] != 0)
	    stop=0;
    }
    s2++;
    st=0;
    naa=s2-s1;

	
    if(naa>naa_max)
	naa_max=naa;

    //Build the new list of peptides
    //cout << "This is the problem: " << naa << " " << i << endl;

    peptide=new int *[naa];
    
    for(int s=s1; s<s2; s++) {
	peptide[st]=new int[kp];
	for(int j=0; j<kp; j++) {
	    peptide[st][j]=tpeptide[s][j];
	}
	st++;
    }
 

}

void print_responsibility(int ncl, double **resp){
  
    char buffer [4096];
    FILE *F;
    int *cluster;
    cluster=new int[kp];
    int ct, pmax;
    double max;
    char file [4096];
     int *size_cluster;
    size_cluster=new int[ncl]; for(int n=0; n<ncl; n++){ size_cluster[n]=0;}
    
  
    sprintf(buffer, "%s/responsibility/resp_%d.txt", out_dir, ncl);
    F=fopen(buffer, "w");
    fprintf(F, "Peptide\t");
    for(int n=0; n<ncl; n++){
	fprintf(F, "%d\t", n+1);
    }
    fprintf(F, "\n");
    for(int j=0; j<kp; j++){
	for(int s=0; s<tnaa; s++){
	    if(tpeptide[s][j] != N)  //Do not print gaps (actually it's better to keep them).
		fprintf(F, "%c", letter[tpeptide[s][j]]);
	    else
		fprintf(F, "-");
	}
	fprintf(F, "\t");
	max=-1;
	for(int n=0; n<ncl; n++){
	    if(resp[j][n]>0.000001)
		fprintf(F, "%.6f\t", resp[j][n]);
	    else
		fprintf(F, "%.6f\t", 0.000001);
	    
	    if(resp[j][n]>max){
		max=resp[j][n];
		cluster[j]=n;
	    }
	}
	size_cluster[cluster[j]]++;
	fprintf(F, "\n");
    }
    fclose(F);
    
    
    //Print the sequences in different clusters in LoLa format

    // WARNING:
    // There is a problem with the X, since they are not treated properly in LoLa...
    // We need to have fake sequence files...
    
    int A=1000;
    int r;
    double **tfr;
    int **int_tfr;
    tfr=new double*[naa];
    int_tfr=new int*[naa];
    for(int s=0; s<naa; s++) {
	tfr[s]=new double[N];
	int_tfr[s]=new int[N];
	for(int j=0; j<N; j++){
	    tfr[s][j]=0;
	    int_tfr[s][j]=0;
	}
    }
    int *index;
    index=new int[naa_max];
    
     
    if(ncl==1){
	sprintf(buffer, "%s/project.txt", out_dir);
	afile.open(buffer, std::ios_base::app);
    } else if(ncl>1){
	sprintf(buffer, "%s/EM_project.txt", out_dir);
	afile.open(buffer, std::ios_base::app);
    }
    
    for(int n=0; n<ncl; n++){
	ct=0;
	if(ncl>1){
	    sprintf(file, "%s/LoLa/LoLa_%d_%d.txt", out_dir,  ncl, n+1);
	} else if(ncl==1){
	    sprintf(file, "%s/LoLa/LoLa_%d.txt", out_dir, n+1);
	}
	F=fopen(file, "w");
	
	if(ncl>1){
	    fprintf(F, "Gene Name	LoLa_%d_%d\nAccession	Refseq:1\nOrganism	H\nNCBITaxonomyID	1\nDomain Number	%d\nDomain Type	HLA\nInterpro ID	1\nTechnique	1\nDomain sequence	A\nDomain Range	1-1\nComment\t\nPeptideName	Peptide	CloneFrequency	QuantData	ExternalIdentifier\n", ncl, n+1, size_cluster[n]);
	} else if(ncl==1){
	    fprintf(F, "Gene Name	LoLa_%d\nAccession	Refseq:1\nOrganism	H\nNCBITaxonomyID	1\nDomain Number	%d\nDomain Type	HLA\nInterpro ID	1\nTechnique	1\nDomain sequence	A\nDomain Range	1-1\nComment\t\nPeptideName	Peptide	CloneFrequency	QuantData	ExternalIdentifier\n", ncl, size_cluster[n]);
	}
	if(size_cluster[n]>0){
	    
	    //Compute the frequency
	    for(int s=0; s<naa; s++) {
		for(int j=0; j<N; j++){
		    tfr[s][j]=0;
		    int_tfr[s][j]=0;
		}
	    }
	    ct=0;
	    for(int j=0; j<kp; j++){
		if(cluster[j]==n){
		    ct++;
		    for(int s=0; s<naa; s++) {
			if(peptide[s][j]<N){
			    tfr[s][peptide[s][j]]++;
			} else {
			    for(int p=0; p<N; p++){
				tfr[s][p]=tfr[s][p]+1.0/N;
			    }
			}
		    }
		}
	    }
	    //normalize
	    for(int s=0; s<naa; s++) {
		index[s]=0;
		for(int p=0; p<N; p++){
		    tfr[s][p]=1.0*tfr[s][p]/ct;
		    int_tfr[s][p]=int(A*tfr[s][p]);
		}
	    }
	    //Check that the PWM columns truly sum up to A
	    for(int s=0; s<naa; s++) {
		ct=0;
		for(int j=0; j<N; j++) {
		    ct=ct+int_tfr[s][j];
		}
		if(ct<A) {
		    for(int t=0; t<A-ct; t++) {
			r=rand()%N;
			int_tfr[s][r]++;
		    }
		}
		//Use the cumulative sum
		for(int j=1; j<N; j++) {
		    int_tfr[s][j]=int_tfr[s][j]+int_tfr[s][j-1];
		}
	    }
		    
	    //Build the fake sequence list
	    for(int t=0; t<A; t++) {
		fprintf(F, "%d\t", t+1);
		for(int s=0; s<naa; s++) {
		    while(t>=int_tfr[s][index[s]]) {
			index[s]++;
		    }
		    fprintf(F, "%c", letter[index[s]]);
		}
		fprintf(F, "\t1\n");
	    }



	    
	   
	} else if(size_cluster[n]==0){
	    fprintf(F, "%d\tEMPTY", 1);
	    for(int s=6; s<naa; s++) {
		fprintf(F, "X");
	    }
	    fprintf(F, "\t1\n");
	}
	fclose(F);
	sprintf(buffer, "LoLa/");
	
	if(ncl>1){
	    afile<<buffer<<"LoLa"<<"_"<<ncl<<"_"<<n+1<<".txt\n";
	} else if (ncl==1){
	    afile<<buffer<<"LoLa"<<"_"<<ncl<<".txt\n";
	}
    }
    afile.close();
}

void print_EM_pwm(int ncl, double ***EM_pwm, double *wcl)
{

    char buffer [4096];
    int ct;
 
    char file [4096];
    FILE *F;
    FILE *F2;
    
    for(int n=0; n<ncl; n++) {
	sprintf(file, "%s/Multiple_PWMs/PWM_%d_%d.txt",out_dir, ncl, n+1);
	F=fopen(file, "w");
	fprintf(F, "PWM_%d_%d\t%.6f\n", ncl, n+1, wcl[n]);
	for(int j=0; j<N; j++) {
	    fprintf(F, "%c\t", letter[j]);
	    for(int s=0; s<naa; s++) {
		if(1.0*EM_pwm[n][s][j]*N>0.000001)
		    fprintf(F, "%.6f\t", 1.0*EM_pwm[n][s][j]*N);
		else
		    fprintf(F, "%.6f\t", 0.000001);
	    }
	    fprintf(F, "\n");
	}
	fclose(F);
    }
    
  
}


//hash to convert the amino acid letters into numbers
int position(char s)
{

    int pp=0;
    for(int i=0; i<strlen(letter); i++){
	if(s == letter[i]){
	    pp=i;
	}
    }  
    return(pp);
    
}



//*********************************
//This part is the initialization and implementation of the EM algorithm

//Compute the Log Likelihood
double loglikelihood(double **eprior, double ***EM_pwm, double *wcl, double *LL, double *sLL, int ncl)
{

    double tLL;
    double sp, p;

    //Compute the loglikelihood. The normalization on the probabilities are not included since they simply correspond to a constant factor
    double LLer=0;
    tLL=0;
    for(int j=0; j<kp; j++) {
	if(best_comp[j]>=0) {
	    p=0;
	    for(int n=0; n<ncl; n++) {
		sp=1;
		for(int s=0; s<naa; s++) {
		    if(peptide[s][j]<N)
			sp=sp*EM_pwm[n][s][peptide[s][j]];
		    else
			sp=sp*1.0/N;
		}
		p=p+wcl[n]*sp;
	    }
	    tLL=tLL+log(p);
	}
    }

    sLL[ncl]=tLL;
    //add the prior
    for(int n=0; n<ncl; n++) {
	for(int s=0; s<naa; s++) {
	    for(int tj=0; tj<N; tj++) {
		if(EM_pwm[n][s][tj]>0)
		    tLL=tLL+eprior[n][s]*log(EM_pwm[n][s][tj]);
		if(EM_pwm[n][s][tj]==0 && eprior[n][s]!=0)
		    cout<<"log(0)\n";
	    }
	}
    }

    LLer=absv(LL[ncl]-tLL);
    LL[ncl]=tLL;

    if(LLer<0.000001) {
	// cout<<1.0*sLL/kp<<endl;
    }
    return(LLer);
}

//compute the responsibilities (E-step)
void expectation(double ***EM_pwm, double *wcl, double **resp, int ncl)
{

    //Compute the prior for the actual EM_pwm
    //***** WARNING: If we want to include the prior, we need to compute the correct normalization
    double tp;
    double gm;

    for(int j=0; j<kp; j++) {
	if(best_comp[j]>=0) {
	    for(int n=0; n<ncl; n++) {
		resp[j][n]=1*wcl[n];
		for(int s=0; s<naa; s++) {
		    if(peptide[s][j]<N)
			resp[j][n]=resp[j][n]*EM_pwm[n][s][peptide[s][j]];
		    else {
			resp[j][n]=resp[j][n]*1.0/N;
		    }
		}
	    }
	    normalize(resp[j], ncl);
	}
    }
}

//Maximize the Log Likelihood (M-step), taking the responsibilities computed at the E-step.
//The maximization can be done analytically using the lagrange Multipliers and is implemented as such
void maximization(double **eprior, double **resp, double ***EM_pwm, double *wcl, double min_error, int ncl)
{

    //Update the mixing coefficients
    for(int n=0; n<ncl; n++) {
	wcl[n]=0;
	for(int j=0; j<kp; j++) {
	    wcl[n]=wcl[n]+resp[j][n];
	}
	wcl[n]=1.0*wcl[n]/kp;
    }

    //Update the model coefficient, i.e. the PWM entries
    for(int n=0; n<ncl; n++) {
	for(int s=0; s<naa; s++) {
	    for(int p=0; p<N; p++) {
		EM_pwm[n][s][p]=eprior[n][s];
		//sum over all peptides that have letter p at position s, excluding the singletons
		for(int j=0; j<kp; j++) {
		    if(peptide[s][j]==p) {
			EM_pwm[n][s][p]=EM_pwm[n][s][p]+resp[j][n];
		    }
		}
	    }
	    //normalize
	    normalize(EM_pwm[n][s], N);
	}
    }
}

//Compute the exponent of the Dirichlet priors. For one component model (single PWM), the exponent of the Dirichlet priors are the random counts
//Use exponents proportional to the entropy of each column (to make it consistent with the random count choice)
void compute_eprior(double **eprior, int ncl)
{

    double *f;
    f=new double[N];
    double tp;
    double xpos;
    double rc;


    for(int n=0; n<ncl; n++) {
	for(int s=0; s<naa; s++) {
	    for(int tj=0; tj<N; tj++) {
		f[tj]=0;
	    }
	    xpos=0;
	    for(int j=0; j<cl_size[n]; j++) {
		//Compute the frequency at each position
		if(comp_pep[n][s][j]!=N) {
		    f[comp_pep[n][s][j]]++;
		} else {
		    for(int tj=0; tj<N; tj++) {
			f[tj]=f[tj]+1.0/N;
		    }
		    xpos++;
		}
	    }
	    eprior[n][s]=(pseudo_count_prior+1.0*xpos/N);
	}
    }

}


void initialize_EM(double ***EM_pwm, double *wcl, int ncl)
{

    int ct;

    if(ncl>0) {
	for(int n=0; n<ncl; n++) {
	    for(int s=0; s<naa; s++) {
		EM_pwm[n][s]=prob(cl_size[n], comp_pep[n][s]);
	    }
	    wcl[n]=1.0*cl_size[n]/kp;
	}
    }
}

