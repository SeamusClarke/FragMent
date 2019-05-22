#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>


// to compile gcc -shared -Wl,-soname,FragMent_C -o FragMent_C.so -fPIC FragMent_C.c -g

// Global definitions 
typedef struct treenode treenode;

struct treenode
{
	//boundary box and centre
	double xmin, xmax;
	double centre;

	//the start and end element of the array which contains particles in this node
	int start;
	int end;

	// pointers to the child nodes
	treenode *leftchild;
	treenode *rightchild;

};

double answer;
int num_direct;



// Internal functions that are not called by Python
// These are minimally commented

// Swap two entries in an array
void swap(float *array, int Index1, int Index2){

	float dummy;

	dummy = array[Index2];
	array[Index2] = array[Index1];
	array[Index1] = dummy;

}

// Split an array around a pivot point, all entries with values greater than the pivot move to the right
int partition(float *array, int lo, int hi){
	
	int pivotIndex = rand() % (hi-lo) + lo ; 	
	float pivotValue = array[pivotIndex];

	swap(array, pivotIndex, hi);

	int storeIndex = lo;
	int ii;

	for(ii=lo;ii<hi;ii++){
		if(array[ii] <= pivotValue){		
			swap(array, ii, storeIndex);
			storeIndex = storeIndex + 1;			
		}
	}

	swap(array, storeIndex, hi);

	return storeIndex;

}

// A Quicksort function
void QuickSort(float *array,int lo, int hi){

	int p;

	if(lo<hi){
	
		p = partition(array,lo, hi);
		QuickSort(array,lo,p-1);
		QuickSort(array,p+1,hi);

	}

}

// Construct a kd tree structure
treenode* kdtree(float *array, int lo, int hi, int maxnumber){

	treenode *node = (treenode*)malloc(sizeof(treenode));

	node->start = lo;
	node->end = hi;

	int medianIndex;

	if(node->end - node->start > maxnumber) {

		medianIndex = (int) ceil(0.5*(node->end-node->start)) + node->start;

		node->leftchild = kdtree(array,node->start,medianIndex,maxnumber);
	 	node->rightchild = kdtree(array,medianIndex,node->end,maxnumber);
	}

	else{
		node->leftchild = NULL;
		node->rightchild = NULL;
	}

	return node;

}

// Define the boundaries of the tree nodes
void BuildBoundaries(float *array, treenode *node){

	int ii;
	double maxx=-1e6,minx=1e6;
	
	treenode *left, *right;

	if((node->leftchild==NULL && node->rightchild!=NULL) || (node->leftchild!=NULL && node->rightchild==NULL)) printf("something is wrong\n");

	else if(node->leftchild!=NULL && node->rightchild!=NULL){

		BuildBoundaries(array,node->leftchild);
		BuildBoundaries(array,node->rightchild);

		left = node->leftchild;
		right = node->rightchild;

		if(left->xmax > right->xmax) node->xmax = left->xmax;
		else node->xmax = right->xmax;

		if(left->xmin < right->xmin) node->xmin = left->xmin;
		else node->xmin = right->xmin;

		node->centre = (left->centre*(left->end - left->start) + right->centre*(right->end - right->start)) / (left->end - left->start + right->end - right->start);

	}
	else{

		node->centre = 0;

		for(ii=node->start;ii<node->end;ii++){

			if(array[ii] > maxx) maxx = array[ii];
			if(array[ii] < minx) minx = array[ii];

			node->centre = node->centre + array[ii];

		}

		node->xmax = maxx;
		node->xmin = minx;

		node->centre = node->centre / (node->end - node->start);


	}


}

// Build a tree
treenode* BuildTree(float *array, int hi, int maxnumber){

	treenode *rootnode = (treenode*)malloc(sizeof(treenode));

	QuickSort(array,0,hi-1);
	rootnode = kdtree(array,0,hi,maxnumber);
	BuildBoundaries(array, rootnode);

	return rootnode;

}

// Delete the tree structure
void DeleteTree(treenode *node){

	if(node->leftchild != NULL && node->rightchild != NULL){

		DeleteTree(node->leftchild);
		DeleteTree(node->rightchild);

		free(node->leftchild);
		free(node->rightchild);

	}

}

// Walk the tree to construct the KDE
void WalkTree(treenode *node, float *array, double position, double std, double error){

	int ii;

	double A,D;

	double opening_crit_dist;

	A = 1./sqrt(2*M_PI*std*std);

	D = 0.5/(std*std);

	opening_crit_dist = error*std*(node->end - node->start);

	if(node->leftchild!=NULL && node->rightchild!=NULL){			

			if(position > node->xmin && position < node->xmax){

				WalkTree(node->leftchild,array,position,std,error);
				WalkTree(node->rightchild,array,position,std,error);
			}

			else if (fabs(position - node->centre) < opening_crit_dist){

				WalkTree(node->leftchild,array,position,std,error);
				WalkTree(node->rightchild,array,position,std,error);

			}

			else{

				answer = answer + (node->end - node->start)*A*exp(-D*(position-node->centre)*(position-node->centre));

			}

	}

	else{
		
		for(ii=node->start;ii<node->end;ii++){

			answer = answer + A*exp(-D*(position-array[ii])*(position-array[ii]));
			num_direct = num_direct + 1;

		}
	
	}


}

// Construct a Gaussian
void Gaussian(double *x, int x_size, const float mean, double std, double *out){

	int ii=0;

	double A,D;

	A = 1./sqrt(2*M_PI*std*std);

	D = 0.5/(std*std);

	for(ii=0;ii<x_size;ii++){

		out[ii] = out[ii] + A*exp(-D*(x[ii]-mean)*(x[ii]-mean));

	}

}

// Determine the bandwidth for the KDE given Scott's method
double Bandwidth_calc(float *array, int a_size){

	int ii=0;

	double std,factor;
	double var=0;
	double mean=0;

	factor = pow((double)a_size,-0.2);

	for(ii=0;ii<a_size;ii++) mean = mean + array[ii];
	mean = mean/a_size;
	for(ii=0;ii<a_size;ii++) var = var + (array[ii]-mean)*(array[ii]-mean);
	var = var/a_size;

	std = factor*1.06*sqrt(var)/1.34;

	return std;

}

// Construct an exact KDE
void KDE(double *x, int x_size, const float *array, int a_size, double bandwidth, double *output){

	int ii=0;

	for(ii=0;ii<a_size;ii++){

		Gaussian(x,x_size,array[ii],bandwidth,output);
	
	}

}

// Construct an approximate KDE
void approxKDE(double *x, int x_size, float *array, int a_size, double bandwidth, double error, double *output){

	treenode *rootnode = (treenode*)malloc(sizeof(treenode));
	rootnode = BuildTree(array, a_size, 8);

	for(int ii=0;ii<x_size;ii++){

		answer = 0;
		num_direct = 0;

		WalkTree(rootnode,array,x[ii],bandwidth,error); 
		output[ii] = answer;

	}

	DeleteTree(rootnode);

}






// These are the external functions which python calls

// Construct the Two-point correlation function by using an exact KDE
void TwoPoint(const double *pos_x, const double *pos_y, int size, int nruns, const double *bounds, int x_size, double lower_lim, double *DD, double *DR, double *RR){

	//struct timeval t1, t2;
	//double elapsedTime;
	//gettimeofday(&t1, NULL);


	// Unpack the boundary array
	double xmin,xmax,ymin,ymax;
	xmin = bounds[0];
	xmax = bounds[1];	
	ymin = bounds[2];
	ymax = bounds[3];

	double max_dist;
	max_dist = sqrt( (xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin) );

	// Initialise counters
	int iDD=0;
	int iDR=0;
	int iRR=0;

	// Create empty arrays to place the random core positions in
	double rpos_x[size];
	double rpos_y[size];

	// Seed for random number generator 
	srand(1);	

	// Determine the size for the DD, DR and RR arrays and then allocate this
	int DD_size=size*(size-1);
	int DR_size=nruns*size*size;
	int RR_size=nruns*size*(size-1);

	float *tot_DD = (float*)malloc(DD_size *sizeof(float));
	float *tot_RR = (float*)malloc(RR_size *sizeof(float));
	float *tot_DR = (float*)malloc(DR_size *sizeof(float));

	// Construct the evaluation point array
	double x[x_size];
	for(int aa=0;aa<x_size;aa++) x[aa] = aa*max_dist/x_size;

	// Loop over all random realisations
	for(int ii=0;ii<nruns;ii++){

		// Randomly place cores and store their locations
		for(int aa=0; aa<size;aa++){
			rpos_x[aa] = (xmax-xmin)*(float)rand()/(float)(RAND_MAX) + xmin;
			rpos_y[aa] = (ymax-ymin)*(float)rand()/(float)(RAND_MAX) + ymin;

		}

		// Loop over all cores to determine separations 
		for(int jj=0;jj<size;jj++){
			for(int kk=0; kk<size; kk++){

				// For DD and RR we do not want to determine the separation between a core and itself
				if(jj!=kk){

					// Only determine the DD array once
					if(ii==0){
						tot_DD[iDD] = sqrt( (pos_x[jj]-pos_x[kk])*(pos_x[jj]-pos_x[kk]) +  (pos_y[jj]-pos_y[kk])*(pos_y[jj]-pos_y[kk]));
						iDD = iDD + 1;

					}

					// Calculate the RR array for each random realisation
					tot_RR[iRR] = sqrt( (rpos_x[jj]-rpos_x[kk])*(rpos_x[jj]-rpos_x[kk]) +  (rpos_y[jj]-rpos_y[kk])*(rpos_y[jj]-rpos_y[kk]));
					iRR = iRR + 1;

				}

				// Calculate the DR array for every random-data point pair
				tot_DR[iDR] = sqrt( (pos_x[jj]-rpos_x[kk])*(pos_x[jj]-rpos_x[kk]) +  (pos_y[jj]-rpos_y[kk])*(pos_y[jj]-rpos_y[kk]));
				iDR = iDR + 1;

			}

		}

	}

	// Determine the bandwidth of the KDE using the DR array
	double bd;
	bd = Bandwidth_calc(tot_DR, DR_size);
	if(bd<lower_lim) bd = lower_lim;
	
	// Construct exact KDEs for DD, DR and RR
	KDE(x,x_size,tot_DD,DD_size,bd,DD);
	KDE(x,x_size,tot_DR,DR_size,bd,DR);
	KDE(x,x_size,tot_RR,RR_size,bd,RR);

	// Normalise the DD, DR and RR arrays before returning them to Python
	for(int bb=0; bb<x_size;bb++){

		DD[bb] = DD[bb] / DD_size;
		DR[bb] = DR[bb] / DR_size;
		RR[bb] = RR[bb] / RR_size;

	}

	// Free the memory for the arrays
	free(tot_DD);
	free(tot_DR);
	free(tot_RR);

	//gettimeofday(&t2, NULL);
	//elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
	//elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
	//printf("The elapsed time in ms for the exact two-point is %lf\n", elapsedTime);

}

// Construct a Two point correlation function using approximate KDEs
void ApproxTwoPoint(const double *pos_x, const double *pos_y, int size, int nruns, const double *bounds, int x_size, double lower_lim, double error, double *DD, double *DR, double *RR){


	//struct timeval t1, t2;
	//double elapsedTime;
	//gettimeofday(&t1, NULL);

	// Unpack the boundary array
	double xmin,xmax,ymin,ymax;

	xmin = bounds[0];
	xmax = bounds[1];	
	ymin = bounds[2];
	ymax = bounds[3];

	double max_dist;
	max_dist = sqrt( (xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin) );

	// Initialise counters
	int iDD=0;
	int iDR=0;
	int iRR=0;

	// Create empty arrays to place the random core positions in
	double rpos_x[size];
	double rpos_y[size];

	// Seed for random number generator 
	srand(1);	

	// Determine the size for the DD, DR and RR arrays and then allocate this
	int DD_size=size*(size-1);
	int DR_size=nruns*size*size;
	int RR_size=nruns*size*(size-1);

	float *tot_DD = (float*)malloc(DD_size *sizeof(float));
	float *tot_RR = (float*)malloc(RR_size *sizeof(float));
	float *tot_DR = (float*)malloc(DR_size *sizeof(float));

	// Construct the evaluation point array
	double x[x_size];
	for(int aa=0;aa<x_size;aa++) x[aa] = aa*max_dist/x_size;

	// Loop over all random realisations
	for(int ii=0;ii<nruns;ii++){

		// Randomly place cores and store their locations
		for(int aa=0; aa<size;aa++){
			rpos_x[aa] = (xmax-xmin)*(float)rand()/(float)(RAND_MAX) + xmin;
			rpos_y[aa] = (ymax-ymin)*(float)rand()/(float)(RAND_MAX) + ymin;

		}

		// Loop over all cores to determine separations 
		for(int jj=0;jj<size;jj++){
			for(int kk=0; kk<size; kk++){

				// For DD and RR we do not want to determine the separation between a core and itself
				if(jj!=kk){

					// Only determine the DD array once
					if(ii==0){
						tot_DD[iDD] = sqrt( (pos_x[jj]-pos_x[kk])*(pos_x[jj]-pos_x[kk]) +  (pos_y[jj]-pos_y[kk])*(pos_y[jj]-pos_y[kk]));
						iDD = iDD + 1;
					}

					// Calculate the RR array for each random realisation
					tot_RR[iRR] = sqrt( (rpos_x[jj]-rpos_x[kk])*(rpos_x[jj]-rpos_x[kk]) +  (rpos_y[jj]-rpos_y[kk])*(rpos_y[jj]-rpos_y[kk]));
					iRR = iRR + 1;

				}

				// Calculate the DR array for every random-data point pair
				tot_DR[iDR] = sqrt( (pos_x[jj]-rpos_x[kk])*(pos_x[jj]-rpos_x[kk]) +  (pos_y[jj]-rpos_y[kk])*(pos_y[jj]-rpos_y[kk]));
				iDR = iDR + 1;

			}

		}

	}


	// Determine the bandwidth of the KDE using the DR array
	double bd;
	bd = Bandwidth_calc(tot_DR, DR_size);
	if(bd<lower_lim) bd = lower_lim;

	// Construct an exact KDE for DD and approximate KDEs for DR and RR arrays
	KDE(x,x_size,tot_DD,DD_size,bd,DD);
	approxKDE(x,x_size,tot_DR,DR_size,bd,error,DR);
	approxKDE(x,x_size,tot_RR,RR_size,bd,error,RR);

	// Normalise the DD, DR and RR arrays before returning them to Python
	for(int bb=0; bb<x_size;bb++){
		DD[bb] = DD[bb] / DD_size;
		DR[bb] = DR[bb] / DR_size;
		RR[bb] = RR[bb] / RR_size;
	}

	// Free the memory for the arrays
	free(tot_DD);
	free(tot_DR);
	free(tot_RR);

	//gettimeofday(&t2, NULL);
	//elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
	//elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
	//printf("The elapsed time in ms for the approximate two-point is %lf\n", elapsedTime);

}



/// A function to return a nrun sets of  npoints random points within a boundary box with separations greater than some lower limit

void ReturnRandomPoints(int nruns, int npoints, const double *bounds, double lower_lim, double *xpos, double *ypos){

	double xmin,xmax,ymin,ymax;

	xmin = bounds[0];
	xmax = bounds[1];	
	ymin = bounds[2];
	ymax = bounds[3];

	int counter = 0;

	int ii,jj;

	double rpos_x[npoints];
	double rpos_y[npoints];

	double dist;
	int flag;

	// Seed for random number generator 
	srand(1);

	while(counter<nruns){

		flag = 0;

		for(ii=0;ii<npoints;ii++){

			rpos_x[ii] = (xmax-xmin)*(float)rand()/(float)(RAND_MAX) + xmin;
			rpos_y[ii] = (ymax-ymin)*(float)rand()/(float)(RAND_MAX) + ymin;

		}

		for(ii=0;ii<npoints-1;ii++){

			for(jj=ii+1;jj<npoints;jj++){

				dist = sqrt( (rpos_x[ii]-rpos_x[jj])*(rpos_x[ii]-rpos_x[jj]) + (rpos_y[ii]-rpos_y[jj])*(rpos_y[ii]-rpos_y[jj]) ); 

				//printf("%lf \n", dist);

				if(dist<lower_lim){

					flag = 1;
					break;
				}

			}

		}

		if(flag==0){

			for(ii=0;ii<npoints;ii++){

				xpos[counter*npoints + ii] = rpos_x[ii];
				ypos[counter*npoints + ii] = rpos_y[ii];

			}

			counter = counter + 1;

		}


	}

}



