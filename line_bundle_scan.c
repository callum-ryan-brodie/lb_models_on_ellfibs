
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// choice for implementing zero slope necessary check
int maxcoeff = 3;
// to be set by input
int 		h11, numkahlergens, numdualeffdivgens, totnumlbs, chirasym, entrymax;
int *		c2X;
int **		kahlergens;
int **		dualeffdivgens;
int ***		inter;
int		mapord; 		// for checking bundles descend to the quotient
int **		mapmat; 		// for checking bundles descend to the quotient
int		largestentryseen; 	// for recording if we see entries larger than entrymax
// to be computed
int ***				transfinter;
float *				c2V;
int					numgoodvecs;
int					nummodels = 0;

int 				numgoodvecguessincrement = 100; 	// for malloc of goodveclist
int 				numgoodvecguesscurrent; 		// for malloc of goodveclist
struct 				goodvec {
						int *	vec;
						float 	index;
						int **	transfmat; 
					};
struct goodvec *	goodveclist;
struct goodvec		finalvecstruct;
// mathematica input and output files
char 	mathout[100];
char 	mathin[100];

// functions for line bundle checks
void 	setc2V(float *c2V, int numbofvecs, int h11, int veclist[numbofvecs][h11]);
float 	getindex(int h11, int vec[h11]);
int 	inmoricheck(int h11, float vec[h11]);
void 	buildtransfinter(int ***inter);
void 	settransfmat(int h11, int numkahlergens, int vec[h11], int **transfmat);
int 	buildgoodveclist(void);
int 	vecinvcheck(int h11, int vec[h11], int mapord, int **mapmat);
int 	subsetaddstozero(int numvecs, int veclength, struct goodvec * passedaddresslist[totnumlbs]);
int 	slopezerocheck1(int numkahlergens, int h11, int **transfmat, int *vec);
int 	slopezerocheck(int numkahlergens, struct goodvec * goodvecaddresslist[totnumlbs], int currentnumlbs);
void 	scannextlevel(struct goodvec * currentaddresslist[totnumlbs], int currentnumlbs, int currenttoplbpos);
void 	scanlasttwolbs(struct goodvec * currentaddresslist[totnumlbs], int currentnumlbs, int currenttoplbpos);
// functions for writing to file
void 	putvect(int veclength, int v[veclength], FILE *pf);
void 	putvectlist(int veclength, int numbofvects, int listofvects[numbofvects][veclength], FILE *pf);
void 	putvectlistlist(int veclength, int listlength, int numboflists, int listoflistofvects[numboflists][listlength][veclength], FILE *pf);
// basic functions
int 	vecdotint(int *vec1, int *vec2, int veclength);
float 	vecdotfloat(float *vec1, float *vec2, int veclength);
int 	matrixcheckint(int nrows, int ncols, int mat[nrows][ncols]);
int 	dotveccheck(int h11, int vec[h11]);
void 	setv1tov2int(int veclength, int *v1, int v2[]);
void 	setm1tom2int(int rowlength, int *m1[rowlength], int m2[][rowlength]);
// functions to malloc space
int 	*create1arrayint(int n);
float 	*create1arrayfloat(int n);
int 	**create2arrayint(int nrows, int ncols);
int 	***create3arrayint(int n1, int n2, int n3);

// for recording progress
int	totnummodelstoscan;
float	progresspercentage = 0.0;
float	percentageleft = 100.0;
float 	tempfactor;

int main()
{
	printf("\n\n");

	int i, j, k; 	// these run up to h11
	int a, b; 	// these run up to numkahlergens

	FILE *fp;

	// read input from mathematica
	// input is (h11,numkahlergens,numdualeffdivgens,c2X,kahlergens,dualeffdivgens,inter,totnumlbs,chirasym,entrymax)
	strcat(strcpy(mathout,getenv("HOME")),"/.mathout");
	strcat(strcpy(mathin,getenv("HOME")),"/.mathin");
	
	
	fp=fopen(mathout,"r");
	if (fp == NULL) {
    		printf("error: couldn't find mathematica input file\n");
		exit(1);
    	return 1;
	}
	
	fscanf(fp,"%d",&h11); fscanf(fp,"%d",&numkahlergens); fscanf(fp,"%d",&numdualeffdivgens); // set by base choice

	c2X =			create1arrayint(h11);
	kahlergens =		create2arrayint(numkahlergens,h11);
	dualeffdivgens =	create2arrayint(numdualeffdivgens,h11);
	inter = 		create3arrayint(h11,h11,h11);
	
	for(i=0;i<h11;i++)
		fscanf(fp,"%d ",&(c2X[i]));
	
	for(a=0;a<numkahlergens;a++)
	for(i=0;i<h11;i++)
		fscanf(fp,"%d ",&(kahlergens[a][i]));

	for(a=0;a<numdualeffdivgens;a++)
	for(i=0;i<h11;i++)
		fscanf(fp,"%d ",&(dualeffdivgens[a][i]));
	
	for(i=0;i<h11;i++)
	for(j=0;j<h11;j++)
	for(k=0;k<h11;k++)
		fscanf(fp,"%d ",&(inter[i][j][k]));
	
	fscanf(fp,"%d",&totnumlbs); fscanf(fp,"%d",&chirasym); fscanf(fp,"%d",&entrymax); // choices for line bundles

	largestentryseen = entrymax;

	// now read (optional) information that represents a map on the space (that we want to quotient the space by)
	mapmat = 		create2arrayint(h11,h11);
	fscanf(fp,"%d",&mapord);
	for(i=0;i<h11;i++)
	for(j=0;j<h11;j++)
		fscanf(fp,"%d ",&(mapmat[i][j]));

	fclose(fp);
	
	// allocate space for other arrays
	c2V = 		create1arrayfloat(h11);
	transfinter = 	create3arrayint(numkahlergens,numkahlergens,h11);
	
	// build transformed intersection numbers
	buildtransfinter(inter);
	
	printf("\n");
	buildgoodveclist();
	printf("numgoodvecs = %d\n",numgoodvecs);
	printf("numgoodvecguesscurrent = %d\n",numgoodvecguesscurrent);
	
	FILE *goodvecoutputfile = fopen("goodvectestoutput.txt", "w");
	int veclist[numgoodvecs][h11];
	for(a = 0; a < numgoodvecs; a++)
		for(i = 0; i < h11; i++)
			veclist[a][i] = goodveclist[a].vec[i];
	putvectlist(h11, numgoodvecs, veclist, goodvecoutputfile);
	fclose(goodvecoutputfile);

	// will use the same struct every time for the final line bundle, and will just overwrite it, so we malloc it here once
	finalvecstruct.vec = create1arrayint(h11);
	finalvecstruct.transfmat = create2arrayint(numkahlergens,numkahlergens);
	
	struct goodvec * blankaddresslist[totnumlbs];
	
	for(i = 0; i < totnumlbs; i++)
		blankaddresslist[i] = NULL;
	
	printf("Starting scan\n");
	
	fp = fopen(mathin,"w");
	fprintf(fp,"{");
	fclose(fp);
	scannextlevel(blankaddresslist, 0, 0); // this is the scan
	fp = fopen(mathin,"a");
	fprintf(fp,"}");
	fclose(fp);
	
	printf("Finished scan, number of models found = %d\n",nummodels);
	
	// free malloc space now that scan is done
	for(i=0;i<numgoodvecs;i++) {
		free(goodveclist[i].vec);
		free(goodveclist[i].transfmat);
	}
	
	free(c2X);free(kahlergens);free(dualeffdivgens);free(inter);
	free(c2V);free(transfinter);
	free(goodveclist);

}



// stores the second Chern class of a set of line bundles in the external array c2V
void setc2V(float c2V[], int numbofvecs, int h11, int veclist[numbofvecs][h11])
{
	int i, j, k, a;
	
	for(k = 0; k < h11; ++k) {
		c2V[k] = 0.;
		for(i = 0; i < h11; ++i)
			for(j = 0; j < h11; ++j)
				for(a = 0; a < numbofvecs; ++a) {
					c2V[k] -= (1.0/2.0)*inter[i][j][k]*veclist[a][i]*veclist[a][j];
					}
	}
} 

// finds the index of a line bundle
float getindex(int h11, int vec[h11])
{
	int temp = 0;
	float index = 0;
	int i, j, k;
	
	for(i = 0; i < h11; ++i)
			temp += 1*vec[i]*c2X[i];
	// index made of two sums, kept separate here just for clarity
	for(i = 0; i < h11; ++i)
		for(j = 0; j < h11; ++j)
			for(k = 0; k < h11; ++k)
				temp += 2*inter[i][j][k]*vec[i]*vec[j]*vec[k];
	index = temp / 12.0; // only divide at the end, to keep the float as accurate as possible
	
	return index;
}

// returns True if a matrix has at least one positive and one negative entry or is all zeros, and False otherwise
// here we assume the matrix is really a symmetric square matrix, since it will be in all our cases
int matrixcheckint(int nrows, int ncols, int mat[nrows][ncols])
{
	int i, j;
	int nneg = 0, npos = 0;
	
	for(i = 0; i < nrows; ++i)
		for(j = i; j < ncols; ++j) {
			if(mat[i][j] < 0)
				++nneg;
			else if(mat[i][j] > 0)
				++npos;
		if(nneg > 0 && npos > 0)
			return 1;
	}
	
	if(nneg == 0 && npos == 0)
		return 1;
	else return 0;
}

// returns True if neither the vector nor its negative are in the cone of effective divisors, and False otherwise
int dotveccheck(int h11, int vec[h11])
{
	int a, i;
	int nneg = 0, npos = 0;
	int dotproduct;
	int gen[h11];
	
	for(a = 0; a < numdualeffdivgens; a++) {
		for(i = 0; i < h11; i++)
			gen[i] = dualeffdivgens[a][i];
		dotproduct = vecdotint(vec,gen,h11);
		if(dotproduct < 0)
			++nneg;
		else if(dotproduct > 0)
			++npos;
		if(nneg > 0 && npos > 0)
			return 1;
	}

	if(nneg == 0 && npos == 0)
		return 1;
	else return 0;
}

// prints a vector to file, in {x1,x2,..} notation
void putvect(int veclength, int v[veclength], FILE *pf)
{
	int i;

	putc('{',pf);
	for(i = 0; i < veclength-1; i++) {
		fprintf(pf,"%d,",v[i]);
	}
	fprintf(pf,"%d",v[veclength-1]);
	putc('}',pf);
}

// prints a list of vectors to file, in {{x1,x2,..},{y1,y2,...},...} notation
void putvectlist(int veclength, int numbofvects, int listofvects[numbofvects][veclength], FILE *pf)
{
	int i;
	
	putc('{',pf);
	for(i = 0; i < numbofvects-1; i++) {
		putvect(veclength,listofvects[i],pf);
		putc(',',pf);
	}
	putvect(veclength,listofvects[numbofvects-1],pf);
	putc('}',pf);
}

// prints a list of a list of vectors to file, in {{{x1,x2,..},...},{{y1,y2,...},...},...} notation
void putvectlistlist(int veclength, int listlength, int numboflists, int listoflistofvects[numboflists][listlength][veclength], FILE *pf)
{
	int i;
	
	putc('{',pf);
	for(i = 0; i < numboflists-1; i++) {
		putvectlist(veclength,listlength,listoflistofvects[i],pf);
		putc(',',pf);
	}
	putvectlist(veclength,listlength,listoflistofvects[numboflists-1],pf);
	putc('}',pf);
}

// checks if vec is in the interior of the Mori cone, by dotting it into the Kahler generators
int inmoricheck(int h11, float vec[h11])
{
	int i, j;
	float currentgenfloat[h11];
	for(i = 0; i < numkahlergens; i++) {
		// have to pass kahlergens[i] as float, so need to convert it
		for(j = 0; j < h11; j++)
			currentgenfloat[j] = (float) kahlergens[i][j];
		if(vecdotfloat(vec,currentgenfloat,h11) < 0) {
			return 0;
		}
	}
	return 1;
}

// dots two integer vectors together, returns the sum
int vecdotint(int *vec1, int *vec2, int veclength)
{	
	int i, sum = 0;
	for(i = 0; i < veclength; i++)
		sum += vec1[i]*vec2[i];
		
	return sum;
}

// dots two floating point vectors together, returns the sum
float vecdotfloat(float *vec1, float *vec2, int veclength)
{	
	int i;
	float sum = 0;
	for(i = 0; i < veclength; i++)
		sum += vec1[i]*vec2[i];
		
	return sum;
}

// sets array v1 equal to array v2, both integer
void setv1tov2int(int veclength, int *v1, int v2[])
{
	int i;
	for(i = 0; i < veclength; i++)
		*(v1+i) = v2[i];
}

// sets matrix m1 equal to array m2, both integer
void setm1tom2int(int rowlength, int *m1[rowlength], int m2[][rowlength])
{
	int i, j;
	for(i = 0; i < rowlength; i++)
		for(j = 0; j < rowlength; j++)
			*(m1[i]+j) = m2[i][j];
}

// transform intersection numbers, store in external variable f, whose pointers have been passed
// we are transforming inter_{i,j,k} to transfinter_{a,b,k}, by contracting with kahler generators
// these are for use in the necessary condition for a common zero slope locus
void buildtransfinter(int ***inter)
{
	int a, b, i, j, k;

	for(a = 0; a < numkahlergens; a++) // set f[a,b,k] for all a,b,k
	for(b = 0; b < numkahlergens; b++)
	for(k = 0; k < h11; k++) {
		transfinter[a][b][k] = 0;
		for(i = 0; i < h11; i++)
		for(j = 0; j < h11; j++)
			transfinter[a][b][k] += kahlergens[a][i]*kahlergens[b][j]*inter[i][j][k];
	}
}

// necessary condition check for existence of zero slope locus for a single vector
// transfmat[a][b] of vec is sum over k of transfinter[a][b][k]*vec[k]
int slopezerocheck1(int numkahlergens, int h11, int **transfmat, int *vec)
{
	int i, j;

	int matrix[numkahlergens][numkahlergens];
	
	for(i = 0; i < numkahlergens; i++)
	for(j = 0; j < numkahlergens; j++)
		matrix[i][j] = transfmat[i][j];

	if(dotveccheck(h11,vec))
		return 1;
	else
		return 0;
}

// checks if a line bundle is invariant under the action of a map on the space
// the map action is specified by a matrix, 'mapmat', that acts on a vector representing a line bundle
int vecinvcheck(int h11, int vec[h11], int mapord, int **mapmat)
{
	int i, j;
	int actedonvec[h11];

	if(mapord>1) {
		for(i = 0; i < h11; i++) {
			actedonvec[i] = 0;
			for(j = 0; j < h11; j++)
				actedonvec[i] += mapmat[i][j]*vec[j];
			if(actedonvec[i] != vec[i])
				return 0;
		}	
	}
	
	return 1;
}

// checks if the index of a line bundle is divisible by the order of a map on the space
int inddivcheck(int index, int mapord)
{
	if(mapord>1) {
		if(fmodf(index,mapord) == 0.0)
			return 1;
		else
			return 0;
	}
	else
		return 1;
}


// cycles through possible line bundle vectors, looks for those with good index and that pass slope condition
int buildgoodveclist(void)
{
	int i, j, k, a, b;

	float 	currentindex;
	int * 	currentvec		= create1arrayint(h11);
	int ** 	currenttransfmat	= create2arrayint(numkahlergens,numkahlergens);
	
	int ** tempmatrix;
	
	printf("Building good vector list\n");
	numgoodvecs = 0;
	numgoodvecguesscurrent = numgoodvecguessincrement;
	goodveclist = malloc(numgoodvecguesscurrent * sizeof(struct goodvec));
	
	if(goodveclist == NULL) {
		printf("error: malloc failure for goodveclist\n");
		exit(1);
	}
	
	// start with 'lowest' possible vector (line bundle)
	for(i = 0; i < h11; i++)
		currentvec[i] = -entrymax;
	
	// now step through possible vectors (line bundles)
	for(i = 0; i < pow((2*entrymax+1),h11); i++) {
	
		// don't want the (0,...,0) vector, note this occurs one less than halfway though
		j=0;
		while(j<h11 && currentvec[j]==0)
			j++;
		if(j!=h11) {
		
			// here is where to do stuff with the vectors being stepped through
			currentindex = getindex(h11, currentvec);
			if(currentindex <= 0 && currentindex >= chirasym && inddivcheck(currentindex,mapord)) { // this assumes that chirasym is negative
				for(a = 0; a < numkahlergens; a++)
				for(b = 0; b < numkahlergens; b++) {
					currenttransfmat[a][b] = 0;
					for(k = 0; k < h11; k++) {
						currenttransfmat[a][b] += transfinter[a][b][k]*currentvec[k];
					}
				}
				if(slopezerocheck1(numkahlergens, h11, currenttransfmat, currentvec) && vecinvcheck(h11, currentvec, mapord, mapmat)) {
					goodveclist[numgoodvecs].vec = create1arrayint(h11);
					for(j = 0; j < h11; j++)
						goodveclist[numgoodvecs].vec[j] = currentvec[j];
					goodveclist[numgoodvecs].index = currentindex;
					goodveclist[numgoodvecs].transfmat = create2arrayint(numkahlergens,numkahlergens);
					for(a = 0; a < numkahlergens; a++)
						for(b = 0; b < numkahlergens; b++)
							goodveclist[numgoodvecs].transfmat[a][b] = currenttransfmat[a][b];
					numgoodvecs++;
				
					if(numgoodvecs == numgoodvecguesscurrent-2) {
						numgoodvecguesscurrent = numgoodvecguesscurrent + numgoodvecguessincrement;
						goodveclist = realloc(goodveclist, numgoodvecguesscurrent * sizeof(struct goodvec)); 
						if(goodveclist == NULL) {
							printf("error: realloc failure for goodveclist\n");
							exit(1);
						}
					}
				}
			}
		
		}

		// handle stepping to the next vector
		j = h11-1; // start at the end
		while(currentvec[j] == entrymax && j >= 0) // work backwards to find entry < entrymax
			j--;
		if(j < 0) // if all are entrymax, then we're done
			break;
		else {
			currentvec[j]++;	   // otherwise increment current slot
			for(k = j+1; k < h11; k++) // and reset all those to the right
				currentvec[k] = -entrymax;
		}
	}
	// finished stepping through possible vectors
	
	free(currentvec);free(currenttransfmat);
	
	return 1;
}



// index check used when we have multiple vectors that are all themselves 'good'
// gets passed currentnumlbs vector (line bundles) addresses
int indexcheck(struct goodvec * goodvecaddresslist[totnumlbs], int currentnumlbs)
{
	float 	twicefinalvecindex, twiceothervecindex, twicerestofpairindex, totindex = 0;
	int 	finalvec[h11], othervec[h11];
	
	// check that total index is not too negative
	// (index of line bundle sum is sum of indices, and each line bundle has non-positive index, so if we've already overshot then we're done for)
	int a;
	for(a = 0; a < currentnumlbs; a++)
		totindex += goodvecaddresslist[a]->index;
	if(totindex < chirasym)
		return 0;

	// want to check that all tensor products of pairs of line bundles have non-positive index
	// we will have checked at an earlier stage all pairs except pairs involving the final (new) vector, so only look at pairs involving the final vector
	int i, j, k;
	for(i = 0; i < h11; i++)
		finalvec[i] = goodvecaddresslist[currentnumlbs-1]->vec[i];
	twicefinalvecindex = 2*(goodvecaddresslist[currentnumlbs-1]->index); 
	for(a = 0; a < currentnumlbs-1; a++) {
		for(i = 0; i < h11; i++)
			othervec[i] = goodvecaddresslist[a]->vec[i];
		twiceothervecindex = 2*(goodvecaddresslist[a]->index);
		twicerestofpairindex = 0;
		for(i = 0; i < h11; i++)
			for(j = 0; j < h11; j++)
				for(k = 0; k < h11; k++)
					twicerestofpairindex += inter[i][j][k]*(finalvec[i]*othervec[j]*othervec[k]+othervec[i]*finalvec[j]*finalvec[k]);
		if(twiceothervecindex + twicerestofpairindex + twicefinalvecindex > 0)
			return 0;
	}

	// if we get here, then the sum has passed all the index checks
	return 1;
}



// checks necessary condition for existence of a common zero slope locus
// gets passed currentnumlbs vector (line bundles) addresses
int slopezerocheck(int numkahlergens, struct goodvec * goodvecaddresslist[totnumlbs], int currentnumlbs)
{
	int coeff[currentnumlbs];
	int vecsum[h11];
	int matrixsum[numkahlergens][numkahlergens];
	
	int i, j, k, x, y, a;
	
	// sets up initial coefficients. these will be incremented as we go
	for(i = 0; i < currentnumlbs; i++)
		coeff[i] = -maxcoeff;
	// here we loop over possible coefficients in the matrix sum
	int teststatus = 1 ;
	for(i = 0; i < pow((2*maxcoeff+1),currentnumlbs); i++) {

		// construct vector sum using coefficients
		for(j = 0; j < h11; j++) {
			vecsum[j] = 0;
			for(a = 0; a < currentnumlbs; a++)
				vecsum[j] += coeff[a]*(goodvecaddresslist[a]->vec[j]);
		}

		if(!dotveccheck(h11,vecsum)) {
			teststatus = 0;
			break;
		};

		// set up the coefficients for the next sum
		// finds rightmost element not yet at max value, increments it, and resets those to its right
		j = currentnumlbs-1;
		while(coeff[j] == maxcoeff && j >= 0) 
			j--;
		if(j < 0)
			break;
		else { 
			coeff[j]++;
			for(k = j+1; k < currentnumlbs; k++)
				coeff[k] = -maxcoeff;
		}
	}
	
	return teststatus;
}


// when have 'good' set of n line bundles, call this to look for extensions to 'good' set of n-1 line bundles
// function calls itself to move up levels, so to scan we just call it once, with initial conditions
void scannextlevel(struct goodvec * passedaddresslist[totnumlbs], int currentnumlbs, int currenttoplbpos)
{
	int i, j;
	
	struct goodvec * currentaddresslist[totnumlbs];
	for(i = 0; i < currentnumlbs; i++) {
		currentaddresslist[i] = passedaddresslist[i];
	}
	
	currentnumlbs++;
	for(i = currenttoplbpos; i < numgoodvecs; i++) {

		if(currentnumlbs == 1) {
			tempfactor = (float)(numgoodvecs-i)/(numgoodvecs+(totnumlbs-1)-i);
			percentageleft = percentageleft*tempfactor;
			if(floor(100.0-percentageleft) > floor(progresspercentage)) {
				progresspercentage = 100.0-percentageleft;
				printf("Progress: %.0f%%, number of models found = %d\n",floor(progresspercentage),nummodels);
			}
		}
	
		currentaddresslist[currentnumlbs-1] =  &goodveclist[i];
		if(currentnumlbs > 1) {
			if(indexcheck(currentaddresslist, currentnumlbs) && slopezerocheck(numkahlergens, currentaddresslist, currentnumlbs) && !subsetaddstozero(currentnumlbs, h11, currentaddresslist)) {
				if(currentnumlbs < (totnumlbs-2)) {
					scannextlevel(currentaddresslist, currentnumlbs, i);
				}
				else
					scanlasttwolbs(currentaddresslist, currentnumlbs, i);
			}
		}
		else {
			if(totnumlbs == 3)
				scanlasttwolbs(currentaddresslist, currentnumlbs, i);
			else
				scannextlevel(currentaddresslist, currentnumlbs, i);
		}
	}
}

// putting in final vector is different, as it gets built from the others to impose c1=0; this function sorts this
void scanlasttwolbs(struct goodvec * passedaddresslist[totnumlbs], int currentnumlbs, int currenttoplbpos)
{
	int i, j, k, a, b;
	
	float index, finalindex;
	float c2Xminusc2V[h11];
	int finalvec[h11], numzerosinfinalvec;
	int veclist[totnumlbs][h11];
	int finaltransfmat[numkahlergens][numkahlergens];
	int goodmodel[totnumlbs][h11];
	FILE * fp;

	struct goodvec * currentaddresslist[totnumlbs];
	for(i = 0; i < currentnumlbs; i++) {
		currentaddresslist[i] = passedaddresslist[i];
	}
	
	int teststatus = 0;
	
	currentnumlbs++;
	for(i = currenttoplbpos; i < numgoodvecs; i++) {
		currentaddresslist[currentnumlbs-1] = &goodveclist[i];
		
		// note: we could get the same model more than once, by chance
		// construct final vector
		numzerosinfinalvec = 0;
		for(k = 0; k < h11; k++) {
			finalvec[k] = 0;
			for(j = 0; j < currentnumlbs; j++)
				finalvec[k] -= currentaddresslist[j]->vec[k];
			if(finalvec[k] == 0)
				numzerosinfinalvec++;
		}
		
		// don't want an all-zeros vector
		// this check not necessary if check as go along for subsets of vects that add to zero
		
		finalindex = getindex(h11, finalvec);
		// check index of final vector first, as may be rubbish
		if(finalindex > 0 || finalindex < chirasym || !inddivcheck(finalindex,mapord))
			continue;

		// must send all vectors to slope check, not just totnumlbs-1, as we use nec. cond. not suf. cond. (else would be guaranteed for totnumlbs if held for numtotlbs-1)
		// need final transformed matrices
		for(a = 0; a < numkahlergens; a++)
			for(b = 0; b < numkahlergens; b++) {
				finaltransfmat[a][b] = 0;
				for(j = 0; j < h11; j++) {
					finaltransfmat[a][b] += transfinter[a][b][j]*finalvec[j];
				}
			}

		// write the struct for the final vector
		for(j = 0; j < h11; j++)
			finalvecstruct.vec[j] = finalvec[j];
		finalvecstruct.index = finalindex;
		for(a = 0; a < numkahlergens; a++)
			for(b = 0; b < numkahlergens; b++)
				finalvecstruct.transfmat[a][b] = finaltransfmat[a][b];
		currentaddresslist[totnumlbs-1] = &finalvecstruct;

		index = 0;
		for(j = 0; j < totnumlbs; j++)
			index += currentaddresslist[j]->index;
			
		// doing the two full indexchecks is a bit redundant, but we need them to make sure tensor pairs involving the final vector don't have non-negative index
		if(index == chirasym && slopezerocheck(numkahlergens, currentaddresslist, totnumlbs) && indexcheck(currentaddresslist, totnumlbs-1) && indexcheck(currentaddresslist, totnumlbs) && !subsetaddstozero(currentnumlbs, h11, currentaddresslist)) {
			for(j = 0; j < totnumlbs; j++)
				for(k = 0; k < h11; k++)
					veclist[j][k] = currentaddresslist[j]->vec[k];

			setc2V(c2V, totnumlbs, h11, veclist);				
			for(j = 0; j < h11; j++) {
				c2Xminusc2V[j] = c2X[j]-c2V[j];
			}
			
			if(inmoricheck(h11, c2Xminusc2V)) {
				for(j = 0; j < totnumlbs; j++)
					for(k = 0; k < h11; k++)
						goodmodel[j][k] = veclist[j][k];
				
				// write model to file
				fp = fopen(mathin,"a");
				if(nummodels > 0)
					fprintf(fp, ",");
				putvectlist(h11, totnumlbs, goodmodel, fp);
				fclose(fp);
				nummodels++;

				for(j = 0; j < h11; j++)
					if(abs(finalvec[j]) > largestentryseen) {
						largestentryseen = abs(finalvec[j]);
						printf("New largest (absolute) line bundle entry seen: %d\n",largestentryseen);
					}
			}
		}
	}
}


// check if any subset of goodvecs in a list of goodvecs add to the zero vector
// if there is such a subset, then this list is not interesting, so we throw it away
int subsetaddstozero(int numvecs, int veclength, struct goodvec * passedaddresslist[totnumlbs])
{
	int coeff[numvecs], maxcoeff = 1, mincoeff = 0;
	int veclist[numvecs][veclength];

	int numzeros, currentelement;	
	int i, j, k, x, a;
	
	for(a = 0; a < numvecs; a++)
		for(x = 0; x < veclength; x++)
			veclist[a][x] = passedaddresslist[a]->vec[x];
	
	for(i = 0; i < numvecs; i++)
		coeff[i] = mincoeff;
	coeff[numvecs-1] = mincoeff+1; // don't want all the coefficients to be zero
	
	// these loops go over possible coefficients
	for(i = 0; i < pow(2,numvecs)-1; i++) { // need the -1 as we skip the coefficients-all-zero case
		
		numzeros = 0;
		for(x = 0; x < veclength; x++) {
			currentelement = 0;
			for(a = 0; a < numvecs; a++)
					currentelement += coeff[a]*veclist[a][x];
			if(currentelement == 0)
				numzeros++;
		}
		if(numzeros == veclength)
			return 1;
				
		// set up the coefficients for the next sum
		// finds rightmost element not yet at max value, increments it, and resets those to its right
		j = numvecs-1;
		while(coeff[j] == maxcoeff && j >= 0) 
			j--;
		if(j < 0)
			break;
		else {
			coeff[j]++;
			for(k = j+1; k < numvecs; k++)
				coeff[k] = mincoeff;
		}
	}

	return 0;
}

int *create1arrayint(int n)
{
	int *p = malloc(n * sizeof(int));
	
	if(p == NULL) {
		printf("error: malloc failure\n");
		exit(1);
	}
		
	return p;
}

float *create1arrayfloat(int n)
{
	float *p = malloc(n * sizeof(float));
	
	if(p == NULL) {
		printf("error: malloc failure\n");
		exit(1);
	}
		
	return p;
}

int **create2arrayint(int nrows, int ncols)
{
	int i;
	int *p = malloc(nrows * ncols * sizeof(int));
	int **pp;
	
	pp = malloc(nrows * sizeof(int *));
	
	if(p == NULL || pp == NULL) {
		printf("error: malloc failure\n");
		exit(1);
	}
	
	for (i = 0; i < nrows; i++)
		pp[i] = &p[i * ncols];
	
	return pp;
}

int ***create3arrayint(int n1, int n2, int n3)
{
	int i;
	int *p = malloc(n1 * n2 * n3 * sizeof(int));
	int **pp;
	int ***ppp;
	
	pp = malloc(n1 * n2 * sizeof(int *));
	
	if(p == NULL || pp == NULL) {
		printf("error: malloc failure\n");
		exit(1);
	}
	
	for (i = 0; i < n1 * n2; i++)
		pp[i] = &p[i * n3];
	
	ppp = malloc(n1 * sizeof(int **));
	
	if(ppp == NULL) {
		printf("error: malloc failure\n");
		exit(1);
	}
	
	for (i = 0; i < n1; i++)
		ppp[i] = &pp[i * n2];
	
	return ppp;
}
