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
# If you plan to use MixMHCp (version 2.0) in any for-profit
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

void init_parameters(int a);

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
void maximization(double **eprior, double **resp, double ***EM_pwm, double *wcl, int ncl);
double loglikelihood(double **eprior, double ***EM_pwm, double *wcl, double *LL, double *sLL, int ncl);
void compute_eprior(double **eprior, int ncl);

void print_EM_pwm(int ncl, double ***EM_pwm, double *wcl);

double test_KLD(double **resp, int ncl);

int position(char s);

void import_bias(char * out_dir, int bs);

void make_cluster_pwm(int si, int ncl, double ***m, double **resp);

void best_ncl(int ncl, double *KLD, double ****EM_pwm, double **wcl);

void predict_other_lengths(int ncl, double ***EM_pwm, double *wcl);

void expectation_all(double ***EM_pwm, double **wcl, double **resp, int ncl, int s1, int **Npos, int **Cpos);
void maximization_all(double ***EM_pwm, double **wcl, double **resp, int ncl, int s1);
double loglikelihood_all(double ***EM_pwm, double **wcl, int ncl, int s1, int **Npos, int **Cpos);


//Global variables

int N;    //Alphabet size
int *naa_all; //number of residues (column) in the alignment for each peptide
int **peptide; // peptide composition [domain][position][peptide] (after filtering some columns)
int **peptide_all; //initial peptides
int kp;   //number of peptides associated with the input sample with length naa_core
int kp_all;   //number of peptides associated with the input sample
char *letter;
int *best_comp;
int sbest_comp;  //number of PWMs chosen by the user
int kpmax;    //size of the largest set of phage peptides
int ***comp_pep;         //peptides in each cluster
int *cl_size;           //size of each cluster
double fsm;
int ncl_max;  //Maximum number of PWMs
int naa_max;
int naa_min;
fstream afile;
char * out_dir;
double * bias;
double pseudo_count;
double pseudo_count_prior;
char * alphabet;
int naa_core;
int trash;
int Nterm;
int Cterm;

double Cterm_pen;
double Nterm_pen;


/*
  Run with:
  ./MixMHCp.x 0 0 5 -d ./output/ -b U
*/

int main(int argc, char ** argv)
{
    if (argc < 3) {
	cout << "Invalid arguments. Run through MixMHCp." << endl;
	exit(2);
    }
    char * alignment_dir = new char[4096];
    out_dir = new char[4096];
    alphabet = new char[4096];
    
    int bs;
     
    for (int i=2; i<argc; i+=2) {
	
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
	    strcpy(alphabet, argv[i+1]);
	}
  	else if (strcmp(argv[i], "-tr") == 0) {
	    trash=atoi(argv[i+1]);
	}
  	else if (strcmp(argv[i], "-lc") == 0) {
	    naa_core=atoi(argv[i+1]);
	}
    }

    cout<<"Trash: "<<trash<<endl;
    cout<<"Core length: "<<naa_core<<endl;
    

    init_parameters(atoi(argv[1]));
    import_alignment(alignment_dir);
    import_bias(out_dir, bs);
    
    initialize_comp();
    
    EMsteps();

    return(0);

}



void init_parameters(int a)
{

    //naa_core=9; //This should be passed as an argument
    //trash=1;

    Nterm=3;
    Cterm=2;
    Cterm_pen=0.2;
    Nterm_pen=0.05;
    
    N=strlen(alphabet);
    letter=new char[N+1];
    strcpy(letter, alphabet);
        
    ncl_max=a;
	
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

    
    int tncl_max=ncl_max+trash;

    //Final values for the Multiple PWMs
    EM_pwm=new double**[tncl_max];
    wcl=new double[tncl_max];
    eprior=new double*[tncl_max];
    for(int n=0; n<tncl_max; n++) {
	wcl[n]=0;
	EM_pwm[n]=new double*[naa_core];
	eprior[n]=new double[naa_core];
	for(int s=0; s<naa_core; s++) {
	    EM_pwm[n][s]=new double[N];
	    eprior[n][s]=0;
	    for(int j=0; j<N; j++)
		EM_pwm[n][s][j]=0;
	}
    }
    
    double ****full_EM_pwm;
    double **full_wcl;
    
    full_wcl=new double*[tncl_max];
    full_EM_pwm=new double***[tncl_max];
    
    for(int tn=0; tn<tncl_max; tn++) {
	full_EM_pwm[tn]=new double**[tn+1+trash];
	full_wcl[tn]=new double[tn+1+trash];
	for(int n=0; n<tn+1+trash; n++) {
	    full_wcl[tn][n]=0;
	    full_EM_pwm[tn][n]=new double*[naa_core];
	    for(int s=0; s<naa_core; s++) {
		full_EM_pwm[tn][n][s]=new double[N];
		for(int j=0; j<N; j++)
		    full_EM_pwm[tn][n][s][j]=0;
	    }
	}
    }

    //temporary values for multiple PWMs (different optimization runs)
    double ***tEM_pwm;      //[component][position][aa]
    double *twcl;
    tEM_pwm=new double**[tncl_max];
    twcl=new double[tncl_max];
    for(int n=0; n<tncl_max; n++) {
	twcl[n]=0;
	tEM_pwm[n]=new double*[naa_core];
	for(int s=0; s<naa_core; s++) {
	    tEM_pwm[n][s]=new double[N];
	}
    }

    //hidden (latent) variables, useful for the EM
    double ***resp;
    resp=new double**[tncl_max];
    for(int ncl=0; ncl<tncl_max; ncl++){
	resp[ncl]=new double*[kpmax];
	for(int j=0; j<kpmax; j++) {
	    resp[ncl][j]=new double[tncl_max];
	}
    }
    double **tresp;
    tresp=new double*[kpmax];
    for(int j=0; j<kpmax; j++) {
	tresp[j]=new double[tncl_max];
    }
    
    double *LL;
    LL=new double[tncl_max];
    double *sLL;
    sLL=new double[tncl_max];

    int rp=5;  //Number of optimization runs (starting from different initial configuration)

    int st;
    double error=0;
    double min_error=0.001;
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
    
    KLD=new double[tncl_max+1];
    for(int n=0; n<tncl_max+1; n++){
	KLD[n]=-10000;
    }
  
    char buffer [4096];
    FILE *F;

    int *sz;
    sz=new int[tncl_max];
    
    
    sprintf(buffer, "%s/project.txt", out_dir);
    afile.open(buffer, ios::out);
    afile<<"#ProjectFile\n";
    afile.close();
    sprintf(buffer, "%s/EM_project.txt", out_dir);
    afile.open(buffer, ios::out);
    afile<<"#ProjectFile\n";
    afile.close();
    	
    cout<<"Number of peptides: "<<kp_all<<endl;
    cout<<"Number of "<<naa_core<<"-mer: "<<kp<<endl;
     
    //Run the EM optimization, 
    
    //Test different number of components
    //Initialize EM_pwm with the values in pwm
    //Here, this differs a bit from the gibbsclustering, but we have to be careful and consistent with the random counts...

    srand(1);
    
    for(int ncl=1; ncl<=ncl_max; ncl++) {

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
		maximization(eprior, tresp, tEM_pwm, twcl, ncl);
		
		//Compute the new loglikelihood
		LLerror=loglikelihood(eprior, tEM_pwm, twcl, LL, sLL, ncl);
		
	    }
	    if(LL[ncl]>max_LL) {
		//If the logikelyhood is larger, then keep these values
		max_LL=LL[ncl];
		max_sLL=sLL[ncl];
		for(int n=0; n<ncl+trash; n++) {
		    wcl[n]=twcl[n];
		    for(int s=0; s<naa_core; s++) {
			for(int j=0; j<N; j++)
			    EM_pwm[n][s][j]=tEM_pwm[n][s][j];
		    }
		}
		
		for(int j=0; j<kp; j++) {
		    for(int n=0; n<ncl+trash; n++) {
			resp[ncl-1][j][n]=tresp[j][n];
		    }
		}
	    }
	}

	//Keep track of the PWMs (this is to compute the best number of motifs).
	for(int n=0; n<ncl+trash; n++){
	    full_wcl[ncl-1][n]=wcl[n];
	    for(int s=0; s<naa_core; s++) {
		for(int j=0; j<N; j++)
		    full_EM_pwm[ncl-1][n][s][j]=EM_pwm[n][s][j];
	    }
	}
	
	KLD[ncl]=test_KLD(resp[ncl-1], ncl);
	print_EM_pwm(ncl, EM_pwm, wcl);
	predict_other_lengths(ncl, EM_pwm, wcl);
	   
    }
	
   

    best_ncl(ncl_max, KLD, full_EM_pwm, full_wcl);
	
}


void best_ncl(int ncl_max, double *KLD, double ****full_EM_pwm, double **full_wcl){

    int ncl_final;

    double thresh1=100/(1.0*N*N*naa_core);
    double thresh2=200/(1.0*N*N*naa_core);
    
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
		for(int s=0; s<naa_core; s++) {
		    for(int j=0; j<N; j++){
			d=d+(full_EM_pwm[ncl-2][n1][s][j]-full_EM_pwm[ncl-1][n2][s][j])*(full_EM_pwm[ncl-2][n1][s][j]-full_EM_pwm[ncl-1][n2][s][j]);
		    }
		}
		d = d/(1.0*naa_core);
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
		    for(int s=0; s<naa_core; s++) {
			for(int j=0; j<N; j++){
			    d=d+(full_EM_pwm[ncl-1][n1][s][j]-full_EM_pwm[ncl-1][n2][s][j])*(full_EM_pwm[ncl-1][n1][s][j]-full_EM_pwm[ncl-1][n2][s][j]);
			}
		    }
		    d = d/(1.0*naa_core);
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
    fprintf(F, "1\t%.6f\n", KLD[1]);
    for(int n=2; n<=ncl_max; n++){
	fprintf(F, "%d\t%.6f\n", n, KLD[n]);
    }
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
    Sscore=new double[ncl+trash];

    double lambda=0.8;
    double sigma;

    sigma=10;

    int *sz;
    sz=new int[ncl+trash];
    for(int n=0; n<ncl+trash; n++) {
	sz[n]=0;
    }
    int *gr;
    gr=new int[kp];
        
    //Compute the size of each cluster
    for(int j=0; j<kp; j++) {
	max_resp=0;
	cl=-1;
	for(int n=0; n<ncl+trash; n++) {
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
    cluster_pwm=new double**[ncl+trash];
    for(int n=0; n<ncl+trash; n++) {
	cluster_pwm[n]=new double*[naa_core];
	for(int s=0; s<naa_core; s++){
	    cluster_pwm[n][s]=new double[N];
	}
    }
    make_cluster_pwm(-1, ncl+trash, cluster_pwm, resp);
    
    
    double ***cl_pwm;
    cl_pwm=new double**[ncl+trash];
    for(int n=0; n<ncl+trash; n++) {
	cl_pwm[n]=new double*[naa_core];
	for(int s=0; s<naa_core; s++){
	    cl_pwm[n][s]=new double[N];
	}
    }
    
    for(int j=0; j<kp; j++) {
	//Find which cluster it belongs to
	max_resp=0;
	cl=gr[j];
	
	//We actually recompute the PWMs only for the cluster which the pepitde is assigned to (excluding the sequence of the peptide itself) and then compute the score.
	//This is to match better the definition in Andreatta et al., although large differences are not expected.

	for(int n=0; n<ncl+trash; n++) {
	    for(int s=0; s<naa_core; s++){
		for(int j=0; j<N; j++){
		    cl_pwm[n][s][j]=cluster_pwm[n][s][j];
		}
	    }
	}
	
	for(int s=0; s<naa_core; s++){
	    if(peptide[j][s] != N){
		cl_pwm[cl][s][peptide[j][s]]--;
	    } else {
		for(int p=0; p<N; p++){
		    cl_pwm[cl][s][p]=cl_pwm[cl][s][p]-1.0/N;
		}
	    }
	}
	
	//Include the normalized counts
	for(int n=0; n<ncl+trash; n++) {
	    //If the cluster is empty, give a flat PWM so that the score is 0.
	    if((n==cl && sz[n]==1) || (n != cl && sz[n]==0)){
		for(int s=0; s<naa_core; s++){
		    for(int p=0; p<N; p++)
			cl_pwm[n][s][p]=bias[p];
		}
	    } else {
		for(int s=0; s<naa_core; s++){
		    normalize_pseudo(cl_pwm[n][s]);
		}
	    }
	}
	
	//Keep flat PWMs for the trash cluster
	if(trash == 1){
	    for(int s=0; s<naa_core; s++){
		for(int p=0; p<N; p++){
		    cl_pwm[ncl][s][p]=bias[p];
		}
	    }
	}
	
	//Compute the scores with each PWM
	for(int n=0; n<ncl+trash; n++) {
	    Sscore[n]=0;
	    for(int s=0; s<naa_core; s++){
		if(peptide[j][s] != N){
		    Sscore[n]=Sscore[n]+2*log(cl_pwm[n][s][peptide[j][s]]/bias[peptide[j][s]])/log(2);		  
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
	for(int n=0; n<ncl+trash; n++) {
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
	for(int s=0; s<naa_core; s++){
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
	    
	    for(int s=0; s<naa_core; s++){
		if(peptide[j][s] != N){
		    m[cl][s][peptide[j][s]]++;
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

    
    cl_size=new int[ncl_max+trash];
    best_comp=new int[kp];

    comp_pep=new int**[ncl_max+trash];
    for(int n=0; n<ncl_max+trash; n++) {
	comp_pep[n]=new int*[naa_core];   //[position][cluster_label]
	for(int s=0; s<naa_core; s++)
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
	for(int j=0; j<ncl+trash; j++) { //This is to make sure that at least one peptide is found in each initial cluster
	    best_comp[j]=j%(ncl+trash);
	}
	
	for(int j=ncl; j<kp; j++) {
	    best_comp[j]=int(j/(t+1))%(ncl+trash);
	}
    }
    else{
	for(int j=0; j<kp; j++) {
	    best_comp[j]=j%(ncl+trash);
	}
    }

    
    for(int n=0; n<ncl+trash; n++)
	cl_size[n]=0;

    for(int j=0; j<kp; j++) {
	cl_size[best_comp[j]]++;
    }
    for(int n=0; n<ncl+trash; n++) {
	comp_kp=0;
	for(int j=0; j<kp; j++) {
	    if(best_comp[j]==n) {
		for(int s=0; s<naa_core; s++) {
		    comp_pep[n][s][comp_kp]=peptide[j][s];
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
    
    if(bs==1){
	bias[0]=0.0702; // A
	bias[1]=0.0230; // C
	bias[2]=0.0473; // D
	bias[3]=0.0710; // E
	bias[4]=0.0365; // F
	bias[5]=0.0657; // G
	bias[6]=0.0263; // H
	bias[7]=0.0433; // I
	bias[8]=0.0572; // K
	bias[9]=0.0996; // L
	bias[10]=0.0213; // M
	bias[11]=0.0359; // N
	bias[12]=0.0631; // P
	bias[13]=0.0477; // Q
	bias[14]=0.0564; // R
	bias[15]=0.0833; // S
	bias[16]=0.0536; // T
	bias[17]=0.0597; // V
	bias[18]=0.0122; // W
	bias[19]=0.0267; // Y
	cout << "Background: Uniprot" << endl;
    } else if(bs==2){
	
	sprintf(file, "%s/bias.txt", out_dir);
	afile.open(file, ios::in);

	char str;
	int p;
	if (afile.is_open()) {
	    cout << "Background: " << file << endl;
	    for (int i=0; i<N; i++) {
		afile>>str;
		p=position(str);
		afile>>bias[p];
	    }
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

    cout << "Input folder:\n" << alignment_dir << endl;
  

    char * pch;
    char str [4096]; //Note: We increased this from 40 to accommodate DNA seqs
    int ln;
    char file [4096];
 
   
    kpmax=0;
    naa_max=0;
    naa_min=100;

    sprintf(file, "%s/peptides.fa", alignment_dir);
	
 
    ifstream myfile (file);
    if (myfile.is_open()) {
	ct=0;
	while (! myfile.eof() ) {
	    getline (myfile,line);
	    if(line.compare(nada) != 0) {
		ct++;
	    }
	}
	myfile.close();
    } else {
	cout << "Unable to open file: " << file << endl;
	exit(2);
    }
    //Import the peptides
    kp_all=ct/2;
    
    if(kp_all>kpmax)
	kpmax=kp_all;

    peptide_all=new int*[kp_all];
    naa_all=new int[kp_all];
    
    std::string delimiter = " ";
    std::string st;
    std::string st2;
    int tmp;
    kp=0;
    afile1.open(file, ios::in);
    for(int j=0; j<kp_all; j++) {
	afile1>>str;
	afile1>>naa_all[j];
	peptide_all[j]=new int[naa_all[j]];

	if(naa_all[j]>naa_max){
	    naa_max=naa_all[j];
	}
	if(naa_all[j]<naa_min){
	    naa_min=naa_all[j];
	}
	
	afile1>>str;
	for(int s=0; s<naa_all[j]; s++) {
	    peptide_all[j][s]=position(str[s]); //Use a nummerical representation for amino acids
	}
	if(naa_all[j]==naa_core){
	    kp++;
	}

    }
    afile1.close();

    peptide=new int*[kp];

    ct=0;
    for(int j=0; j<kp_all; j++) {
	if(naa_all[j]==naa_core){
	    peptide[ct]=new int[naa_core];
	    for(int s=0; s<naa_core; s++) {
		peptide[ct][s]=peptide_all[j][s];
	    }
	    ct++;
	}
    }




    

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
	    for(int s=0; s<naa_core; s++) {
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
	p=0;
	for(int n=0; n<ncl+trash; n++) {
	    sp=1;
	    for(int s=0; s<naa_core; s++) {
		if(peptide[j][s]<N)
		    sp=sp*EM_pwm[n][s][peptide[j][s]];
		else
		    sp=sp*1.0/N;
	    }
	    p=p+wcl[n]*sp;
	}
	tLL=tLL+log(p);
    }

    sLL[ncl]=tLL;
    //add the prior
    for(int n=0; n<ncl+trash; n++) {
	for(int s=0; s<naa_core; s++) {
	    for(int tj=0; tj<N; tj++) {
		if(EM_pwm[n][s][tj]>0)
		    tLL=tLL+eprior[n][s]*log(EM_pwm[n][s][tj]);
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
	for(int n=0; n<ncl+trash; n++) {
	    resp[j][n]=1*wcl[n];
	    for(int s=0; s<naa_core; s++) {
		if(peptide[j][s]<N)
		    resp[j][n]=resp[j][n]*EM_pwm[n][s][peptide[j][s]];
		else {
		    resp[j][n]=resp[j][n]*1.0/N;
		}
	    }
	}
	normalize(resp[j], ncl+trash);
    }
}

//Maximize the Log Likelihood (M-step), taking the responsibilities computed at the E-step.
//The maximization can be done analytically using the lagrange Multipliers and is implemented as such
void maximization(double **eprior, double **resp, double ***EM_pwm, double *wcl, int ncl)
{

    //Update the mixing coefficients
    for(int n=0; n<ncl+trash; n++) {
	wcl[n]=0;
	for(int j=0; j<kp; j++) {
	    wcl[n]=wcl[n]+resp[j][n];
	}
	wcl[n]=1.0*wcl[n]/kp;
    }

    //Update the model coefficient, i.e. the PWM entries
    for(int n=0; n<ncl; n++) {
	for(int s=0; s<naa_core; s++) {
	    for(int p=0; p<N; p++) {
		EM_pwm[n][s][p]=eprior[n][s];
		//sum over all peptides that have letter p at position s, excluding the singletons
		for(int j=0; j<kp; j++) {
		    if(peptide[j][s]==p) {
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


    for(int n=0; n<ncl+trash; n++) {
	for(int s=0; s<naa_core; s++) {
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
	    for(int s=0; s<naa_core; s++) {
		EM_pwm[n][s]=prob(cl_size[n], comp_pep[n][s]);
	    }
	    wcl[n]=1.0*cl_size[n]/kp;
	}
	if(trash==1){
	    
	    wcl[ncl]=1.0*cl_size[ncl]/kp;
	    for(int i=0; i<N; i++){
		for(int s=0; s<naa_core; s++) {
		    EM_pwm[ncl][s][i]=bias[i];
		}
	    }
	}
    }
}



//Compute the Log Likelihood
double loglikelihood_all(double ***EM_pwm, double **wcl, int ncl, int s1, int **Npos, int **Cpos)
{

    double tLL;
    double sp, p;

    //Compute the loglikelihood. The normalization on the probabilities are not included since they simply correspond to a constant factor

    tLL=0;

    for(int j=0; j<kp_all; j++) {
	p=0;
	if(naa_all[j]==s1){
	    
	    for(int n=0; n<ncl+trash; n++) {
		sp=1;
		
		for(int s=Npos[j][n]; s<Npos[j][n]+Nterm; s++) {
		    if(peptide_all[j][s]<N){
			sp=sp*EM_pwm[n][s-Npos[j][n]][peptide_all[j][s]];
		    } else
			sp=sp*1.0/N;
		}
		for(int s=Cpos[j][n]; s<Cpos[j][n]+Cterm; s++) {
		    if(peptide_all[j][s]<N){
			sp=sp*EM_pwm[n][s-Cpos[j][n]-Cterm+naa_core][peptide_all[j][s]];
		    }else
			sp=sp*1.0/N;
		}
		p=p+wcl[s1][n]*sp;
	    }
	    
	    tLL=tLL+log(p);
	}
    }  
    return(tLL);
}

//compute the responsibilities (E-step)
void expectation_all(double ***EM_pwm, double **wcl, double **resp_all, int ncl, int s1, int **Npos, int **Cpos)
{

    //Compute the prior for the actual EM_pwm
    //***** WARNING: If we want to include the prior, we need to compute the correct normalization
    double tp;
    double gm;
    int ps;
    double tresp1;
    double tresp2;
    int tsC;
    int tsN;
    double tr;

    int check=-1;
    
    for(int j=0; j<kp_all; j++) {
	
	if(naa_all[j]==s1){
	    
	    if(naa_all[j]<naa_core){
		
		for(int n=0; n<ncl+trash; n++) {
		    resp_all[j][n]=1*wcl[s1][n];
		    for(int s=0; s<Nterm; s++) {
			if(peptide_all[j][s]<N){
			    resp_all[j][n]=resp_all[j][n]*EM_pwm[n][s][peptide_all[j][s]];
			}
			else {
			    resp_all[j][n]=resp_all[j][n]*1.0/N;
			}
		    }
		    for(int s=naa_core-Cterm; s<naa_core; s++) {
			ps=s-(naa_core-naa_all[j]);
			if(peptide_all[j][ps]<N){
			    resp_all[j][n]=resp_all[j][n]*EM_pwm[n][s][peptide_all[j][ps]];
			}
			else {
			    resp_all[j][n]=resp_all[j][n]*1.0/N;
			}
		    }
		    Cpos[j][n]=naa_all[j]-Cterm;
		    Npos[j][n]=0;
		}   

	    } else if(naa_all[j]>=naa_core){

		for(int n=0; n<ncl+trash; n++) {
		    resp_all[j][n]=0;
		    //scan all positions for the N-terminal motif
		    for(int sN=0; sN<=naa_all[j]-naa_core; sN++){
			
			tresp1=wcl[s1][n];
			for(int s=0; s<Nterm; s++) {
			    if(peptide_all[j][s+sN]<N){
				tresp1=tresp1*EM_pwm[n][s][peptide_all[j][s+sN]];
			    }
			    else {
				tresp1=tresp1*1.0/N;
			    }
			}
			for(int s=0; s<sN; s++){
			    tresp1=tresp1*Nterm_pen;
			}
			//Check all possible C-terminal motifs
			tresp2=0;
			tsC=0;
			for(int sC=sN+naa_core-Cterm; sC<=naa_all[j]-Cterm; sC++){
			    
			    tr=1;
			    for(int s=0; s<Cterm; s++) {
				ps=sC+s;
				
				if(peptide_all[j][ps]<N){
				    tr=tr*EM_pwm[n][naa_core-Cterm+s][peptide_all[j][ps]];
				}
				else {
				    tr=tr*1.0/N;
				}
			    }		   
			   
			    for(int s=sC; s<naa_all[j]-Cterm; s++){
				tr=tr*Cterm_pen;
			    }
			    if(tr>tresp2){
				tresp2=tr;
				tsC=sC;
			    }
			    
			}
			tresp1=tresp1*tresp2;
			if(tresp1>resp_all[j][n]){
			    resp_all[j][n]=tresp1;
			    Npos[j][n]=sN;
			    Cpos[j][n]=tsC;
			}
		    }
		}
	    }
	    normalize(resp_all[j], ncl+trash);
	}
    }
}

//Maximize the Log Likelihood (M-step), taking the responsibilities computed at the E-step.
//The maximization can be done analytically using the lagrange Multipliers and is implemented as such
void maximization_all(double ***EM_pwm, double **wcl, double **resp_all, int ncl, int s1)
{
    int ct;
    
    //Update the mixing coefficients
    for(int n=0; n<ncl+trash; n++) {
	wcl[s1][n]=0;
	ct=0;
	for(int j=0; j<kp_all; j++) {
	    if(naa_all[j]==s1){
		wcl[s1][n]=wcl[s1][n]+resp_all[j][n];
		ct++;
	    }
	}
	wcl[s1][n]=1.0*wcl[s1][n]/ct;
    }
    
}


void predict_other_lengths(int ncl,  double ***EM_pwm, double *wcl){

    ///////////
    //Still need to implement a down-weighting for long C- or N-terminal extensions
    //Still need to change the weight on the noise cluster for longer peptides. Ideally we should try to learn the weights for longer peptides, especially the weight of the trash cluster.
    ///////////
    
    double **resp_all;
    resp_all=new double*[kp_all];
    int ps;
    int ts2;
    double tresp1, tresp2, tr;
    int **Cpos; Cpos=new int*[kp_all];
    int **Npos; Npos=new int*[kp_all];

    double tr_wcl=0.1;
    double **wcl_all; wcl_all=new double*[naa_max+1];
    for(int s1=naa_min; s1<naa_max+1; s1++){
	wcl_all[s1]=new double[ncl+trash];
	if(s1 != naa_core){
	    for(int n=0; n<ncl+trash; n++){
		wcl_all[s1][n]=1.0/(ncl+trash);
	    }
	} else {
	    for(int n=0; n<ncl+trash; n++){
		wcl_all[s1][n]=wcl[n];
	    }
	}
    }

    for(int j=0; j<kp_all; j++){
	Cpos[j]=new int[ncl+trash];
	Npos[j]=new int[ncl+trash];
	resp_all[j]=new double[ncl+1];
	for(int n=0; n<ncl+trash; n++){
	    Cpos[j][n]=naa_all[j]-Cterm;
	    Npos[j][n]=0;
	}
	for(int n=0; n<ncl; n++){
	    resp_all[j][n]=1.0/(ncl+trash);
	}
	if(trash==1){
	    resp_all[j][ncl]=1.0/(ncl+trash);
	} else {
	    resp_all[j][ncl]=0;
	}
    }
    
    double LLerror=1000;
    double min_error=0.00001;

    double LL;
    double sLL;

    int *size;
    size=new int[naa_max+1];
    for(int s1=naa_min; s1<naa_max+1; s1++){
	size[s1]=0;
    }
    for(int i=0; i<kp_all; i++){
	size[naa_all[i]]++;
    }
    
    //////////////
    //Compute the weights by maximizing the log-likelihood
    //////////////
    
    for(int s1=naa_min; s1<naa_max+1; s1++){

	if(s1 != naa_core){

	    if(size[s1]>0){
		LLerror=1000;
		sLL=loglikelihood_all(EM_pwm, wcl_all, ncl, s1, Npos, Cpos);
		
		while(LLerror>min_error) {
		    
		    expectation_all(EM_pwm, wcl_all, resp_all, ncl, s1, Npos, Cpos);
		    maximization_all(EM_pwm, wcl_all, resp_all,ncl,  s1);
		    LL=loglikelihood_all(EM_pwm, wcl_all,ncl, s1, Npos, Cpos);
		    LLerror=absv(LL-sLL);
		    
		    sLL=LL;
		    
		}
	    }
	    
	} else {

	    ////////
	    //This is the special case of the 9-mer where we take the whole PWMs
	    ////////
	    
	    for(int j=0; j<kp_all; j++){
		if(naa_all[j]==naa_core){
		    for(int n=0; n<ncl+trash; n++) {
			resp_all[j][n]=1*wcl[n];
			for(int s=0; s<naa_core; s++) {
			    if(peptide_all[j][s]<N)
				resp_all[j][n]=resp_all[j][n]*EM_pwm[n][s][peptide_all[j][s]];
			    else {
				resp_all[j][n]=resp_all[j][n]*1.0/N;
			    }
			}
			Cpos[j][n]=naa_core-Cterm;
			Npos[j][n]=0;	
		    }
		    
		    normalize(resp_all[j], ncl+trash);
		}
	    }
	    
	}

    }

    
    
    
    char buffer [4096];
    FILE *F;
    int ct, pmax;
    double max;
    char file [4096];

   
    int *cluster;
    cluster=new int[kp_all];
    int t;
    int **size_cluster;
    
    size_cluster=new int*[naa_max+1];
    for(int s=0; s<naa_max+1; s++){
	size_cluster[s]=new int[ncl+trash];
	for(int n=0; n<ncl+trash; n++){
	    size_cluster[s][n]=0;
	}
    }
    
    sprintf(buffer, "%s/responsibility/resp_%d.txt", out_dir, ncl);
    F=fopen(buffer, "w");
    fprintf(F, "Peptide\t");
    for(int n=0; n<ncl; n++){
	fprintf(F, "%d\t", n+1);
    }
    fprintf(F, "Trash\tLength");
   
    for(int n=0; n<ncl; n++){
	fprintf(F, "\tStart_%d\tEnd_%d", n+1, n+1);
    }
    if(trash==1){
	fprintf(F, "\tStart_Trash\tEnd_Trash", n+1, n+1);
    }
    fprintf(F, "\n");

    
    for(int j=0; j<kp_all; j++){

	
	
	for(int s=0; s<naa_all[j]; s++){
	    if(peptide_all[j][s] != N)  //Do not print gaps (actually it's better to keep them).
		fprintf(F, "%c", letter[peptide_all[j][s]]);
	    else
		fprintf(F, "-");
	}
	
	for(int n=0; n<ncl+trash; n++){
	    if(resp_all[j][n]>0.000001)
		fprintf(F, "\t%.6f", resp_all[j][n]);
	    else
		fprintf(F, "\t%.6f", 0.000001);
	}
	if(trash==0){
	    fprintf(F, "\t%.6f", 0.000000);
	}
	max=-1;
	for(int n=0; n<ncl+trash; n++){
	    if(resp_all[j][n]>max){
		max=resp_all[j][n];
		cluster[j]=n;
	    }
	}
	size_cluster[naa_all[j]][cluster[j]]++;
	fprintf(F, "\t%d", naa_all[j]);
	for(int n=0; n<ncl+trash; n++){
	    fprintf(F, "\t%d\t%d", Npos[j][n]+1, Cpos[j][n]+Cterm);
	}
	fprintf(F, "\n");
    }
    fclose(F);


    //Print the weights
    sprintf(buffer, "%s/weights/weights_%d.txt", out_dir, ncl);
    F=fopen(buffer, "w");
    fprintf(F, "Length\tN\t");
    for(int n=0; n<ncl; n++){
	fprintf(F, "%d_weight\t", n+1);
    }
    fprintf(F, "Trash_weight\t");
    for(int n=0; n<ncl; n++){
	fprintf(F, "%d_distr\t", n+1);
    }
    fprintf(F, "Trash_distr\n");

    int *size_cl;
    size_cl=new int[ncl+trash];
    for(int n=0; n<ncl+trash; n++){
	size_cl[n]=0;
	for(int s1=naa_min; s1<naa_max+1; s1++){
	    size_cl[n]=size_cl[n]+size_cluster[s1][n];
	}
    }
    double *w_cl;
    w_cl=new double[ncl+trash];
    for(int n=0; n<ncl+trash; n++){
	w_cl[n]=0;
	for(int s1=naa_min; s1<naa_max+1; s1++){
	    w_cl[n]=w_cl[n]+size[s1]*wcl_all[s1][n];
	}
    }
    
    for(int s1=naa_min; s1<naa_max+1; s1++){
	fprintf(F, "%d\t%d", s1, size[s1]);
	for(int n=0; n<ncl+trash; n++){
	    if(size[s1]>0)
		fprintf(F, "\t%.6f", wcl_all[s1][n]);
	    else
		fprintf(F, "\t%.6f", 0.0);
	}
	if(trash==0){
	    fprintf(F, "\t%.6f", 0.0);
	}

	for(int n=0; n<ncl+trash; n++){
	    if(size[s1]>0)
		//fprintf(F, "\t%.6f", 1.0*wcl_all[s1][n]*size[s1]/w_cl[n]);
		fprintf(F, "\t%d", size_cluster[s1][n]);
	    else
		fprintf(F, "\t%.6f", 0.0);
	}
	if(trash==0){
	    fprintf(F, "\t%.6f", 0.0);
	}
	fprintf(F, "\n");
    }
    fclose(F);
    
    

    //**********
    //Now create the logos for all peptide length
    //**********

    for(int s1=naa_min; s1<=naa_max; s1++){

	if(ncl==1){
	    sprintf(buffer, "%s/project.txt", out_dir);
	    afile.open(buffer, std::ios_base::app);
	} else if(ncl>1){
	    sprintf(buffer, "%s/EM_project.txt", out_dir);
	    afile.open(buffer, std::ios_base::app);
	}
	
	for(int n=0; n<ncl+trash; n++){
	    ct=0;
	    if(n<ncl){
		sprintf(file, "%s/LoLa/LoLa_L%d_%d_%d.txt", out_dir, s1, ncl, n+1);
	    } else if(n==ncl){
		sprintf(file, "%s/LoLa/LoLa_L%d_%d_Trash.txt", out_dir, s1, ncl);
	    }
	    
	    F=fopen(file, "w");
	
	    
	    if(n==ncl){
		fprintf(F, "Gene Name	LoLa_L%d_%d_Trash\nAccession	Refseq:1\nOrganism	H\nNCBITaxonomyID	1\nDomain Number	%d\nDomain Type	HLA\nInterpro ID	1\nTechnique	1\nDomain sequence	A\nDomain Range	1-1\nComment\t\nPeptideName	Peptide	CloneFrequency	QuantData	ExternalIdentifier\n", s1, ncl, size_cluster[s1][n]);
	    } else if(n<ncl){
		fprintf(F, "Gene Name	LoLa_L%d_%d_%d\nAccession	Refseq:1\nOrganism	H\nNCBITaxonomyID	1\nDomain Number	%d\nDomain Type	HLA\nInterpro ID	1\nTechnique	1\nDomain sequence	A\nDomain Range	1-1\nComment\t\nPeptideName	Peptide	CloneFrequency	QuantData	ExternalIdentifier\n", s1, ncl, n+1, size_cluster[s1][n]);
	    }
	    
	    if(size_cluster[s1][n]>0){
	    
		
		//Build the  sequence list
		t=0;
		for(int j=0; j<kp_all; j++) {
		    if(naa_all[j]==s1 && cluster[j]==n){
			fprintf(F, "%d\t", t+1);
			for(int s=0; s<naa_all[j]; s++) {
			    fprintf(F, "%c", letter[peptide_all[j][s]]);
			}
			fprintf(F, "\t1\n");
			t++;
		    }	    
		}
	    } else if(size_cluster[s1][n]==0){
		fprintf(F, "%d\tEMPTY", 1);
		for(int s=6; s<naa_core; s++) {
		    fprintf(F, "X");
		}
		fprintf(F, "\t1\n");
	    }
	    fclose(F);
	    sprintf(buffer, "LoLa/");
	    
	    
	    if(n<ncl){
		afile<<buffer<<"LoLa_L"<<s1<<"_"<<ncl<<"_"<<n+1<<".txt\n";
	    } else if(n==ncl){
		afile<<buffer<<"LoLa_L"<<s1<<"_"<<ncl<<"_Trash.txt\n";
	    }
	   
	}
	afile.close();
    
    }
    
}
