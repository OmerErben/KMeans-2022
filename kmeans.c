#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
typedef struct node
{
    double* vec;
    struct node* next;
    struct node* prev;
} node;

typedef struct cluster{
    double* clustVec;
    node* vecLst;
} cluster;

typedef struct modify{
    node* nodeVec;
    int oldClustIndex;
    int newClustIndex;
} modify;

#define DEFAULT_MAX_ITER 200
#define EPSILON 0.001

/* function decleration*/
double** initMat(FILE* fp, int dim, int numLines);
cluster* initClusters(double** mat, int K, int dim);
int argmin(double* vec, cluster* clusters, int K, int dim);
void addVecToClust(node* nodeVec, cluster* clustVec);
void removeVecFromClust(node* nodeVec, cluster* clustVec);
void swap(node* nodeVec, int oldClust, int newClust, cluster* clusters);
double euqlideNorm(double* oldClust, double* newClust, int dim);
double updateClusters(cluster* clusters, int K, int dim);
int write(cluster* clusters, int K, FILE* fpOut, int dim);
int matLen(FILE *fp);
int dimCalc(FILE *fp);
int initClustForVec( int numLines, int K, int dim, cluster* clusters, double** mat);
modify* updModify(cluster* clusters, modify* modArray, int dim, int K );
void impChanges(cluster* clusters, int numLines, modify* modArray);
int KMean(int K, int maxIter, FILE* fpIn, FILE* fpOut);
void freeMemory(double** mat, cluster* clusters, int matDim, int k);
int submitArgs(int argc, char **argv, FILE** fpIn, FILE** fpOut, double* k, double* maxIter);



/* create an array of pointers to vectors called "mat"
 output: read file and update the array of vertex , init all clusters to null*/
double** initMat(FILE* fp, int dim, int numLines){
     double** mat = (double**)calloc(numLines , sizeof(double*));
     int lineNumber, i;
     double* vec;
     char* line;
     char* currNumStart;
     char* currComma;
     rewind(fp);

        line = (char*)calloc(1024,sizeof(char));
        if (line == NULL)
        {
            printf("An Error Has Occurred\n");
            return NULL;
        }
        
        for(lineNumber =0 ; lineNumber < numLines; lineNumber++)
        {
            currNumStart = line;
            fscanf(fp, "%s",line);
            vec = (double*) calloc(dim,sizeof(double));
            if (vec == NULL)
            {
                printf("An Error Has Occurred\n");
                return NULL;
            }
            for ( i =0; i<dim -1; i++){
                currComma = strchr((char*)currNumStart,','); /* point to the first "," in line */
                *currComma = '\0'; 
                vec[i] = atof(currNumStart); /* insert word to the vector */
                currNumStart = ++currComma; /* update line to tne next char after , */
            }
            vec[dim-1] = atof(currNumStart);
            mat[lineNumber]= vec; 
       
        }
        free(line);
        fclose(fp);
        return mat;
}
/* free memory at the end of KMeans */
void freeMemory(double** mat, cluster* clusters, int matDim, int k){
    node* curr;
    node* step;
    int i;
    for (i = 0; i < matDim; i++)
    {
        free(mat[i]);
    }
    free(mat);
    for (i = 0; i < k; i++)
    {
        free(clusters[i].clustVec);
        curr = clusters[i].vecLst;
        while(curr){
            step = curr->next;
            free(curr);
            curr = step;
        }

    }
    free((void*)clusters);
}

/* create mu_array from first K xi from X
 output: return mu_array - the firt K xi from X */
cluster* initClusters(double** mat, int K, int dim){
    cluster* clusters = (cluster*)calloc(K,sizeof(cluster));
    int i, word;
    for(i =0; i< K; i++){
        clusters[i].clustVec= (double*)calloc(dim,sizeof(double));
        if(clusters[i].clustVec == NULL){
            printf("An Error Has Occurred\n");
            return NULL;
        }
        clusters[i].vecLst = NULL;
        for(word =0; word < dim; word++){
            clusters[i].clustVec[word] = mat[i][word];
        }
    }
    return clusters;
}

/* outpot: the index of its clostest cluster of vec*/
int argmin(double* vec, cluster* clusters, int K, int dim){
   double minSum, sum;
   int minIndex=0, clustIndex, curDim;

    for(clustIndex =0; clustIndex < K; clustIndex++ ){
       sum =0;
       for(curDim = 0; curDim< dim; curDim++){
           sum += pow((vec[curDim] - clusters[clustIndex].clustVec[curDim]),2);
       }
       if (clustIndex ==0) minSum = sum;
        else if(sum<minSum) {
            minSum = sum;
            minIndex = clustIndex;
        } 
    }
    return minIndex;
}

/* output: add node of vec to vecLst of clustVec*/
void addVecToClust(node* nodeVec, cluster* clustVec){
    nodeVec->prev = NULL;
    if(clustVec->vecLst == NULL){ /* vecLst empty*/
        clustVec->vecLst = nodeVec; 
        nodeVec->next = NULL;
    }else{
    clustVec->vecLst->prev = nodeVec;
    nodeVec->next = clustVec->vecLst;
    clustVec->vecLst = nodeVec;
    }
}

/* output: remove vec from vecLst of clustVec*/
void removeVecFromClust(node* nodeVec, cluster* clustVec){
    if(nodeVec->prev == NULL){ /* first node in list*/
        if(nodeVec->next == NULL){ /* nodeVec is the only node in list*/
            clustVec->vecLst = NULL;
        }else{
            nodeVec->next->prev = NULL;
            clustVec->vecLst= nodeVec->next;
        }
    } 
    else if (nodeVec->next ==NULL) /* last node in list*/
    {
        nodeVec->prev->next=  NULL;
    }
    else{ 
        nodeVec->prev->next = nodeVec->next;
        nodeVec->next->prev = nodeVec->prev;
    }
}

/* output: remote nodeVec from oldCluster list and add it to the start of newCluster list*/
void swap(node* nodeVec, int oldClust, int newClust, cluster* clusters){
    removeVecFromClust(nodeVec, &clusters[oldClust]);
    addVecToClust(nodeVec,&clusters[newClust]);
}

/*euqlideNorm
output: for every dim (sum of (oldClust[curDim]**2 -newClust[curDim])**2 )**0.5*/
double euqlideNorm(double* oldClust, double* newClust, int dim){
    double sum =0, power;
    int i;
    for(i =0; i< dim ; i++){
        power = pow(oldClust[i] - newClust[i],2);
        sum += power;
    }
    sum = pow(sum,0.5);
    return sum;
}

/* output: compute and update newClustArray and returns the maximum delta */
double updateClusters(cluster* clusters, int K, int dim){
    double deltaMax =0, delta;
    int clusterIndex, vecLstLen, curDim;
    double* newClust;
    double* oldClust;
    node* currVec;
    
    for(clusterIndex =0; clusterIndex< K; clusterIndex++){
        vecLstLen =0;
        newClust = (double*)calloc(dim, sizeof(double));
        if(newClust == NULL){
            printf("An Error Has Occurred\n");
            return -1;
        }
        currVec = clusters[clusterIndex].vecLst;
        while (currVec != NULL)
        {
            for(curDim =0; curDim<dim; curDim++){
                newClust[curDim] += currVec->vec[curDim];
            }
            vecLstLen++;
            currVec = currVec->next;
        }
        /* we assume vecLstLen != 0 becuse clustVec dosent have an empty list */
        for(curDim =0; curDim< dim; curDim++){
            newClust[curDim] = newClust[curDim]/vecLstLen;
        }
        oldClust =  clusters[clusterIndex].clustVec;
        delta = euqlideNorm(oldClust, newClust, dim);
        if (delta > deltaMax) deltaMax= delta;
        free(oldClust);
        clusters[clusterIndex].clustVec = newClust;
    }
    return deltaMax;
}

/* output: create output.txt and write to it*/
int write(cluster* clusters, int K, FILE* fpOut, int dim ){
    int cluster, curDim;
    for(cluster =0; cluster <K; cluster++){
        for(curDim =0; curDim< dim; curDim++){
             fprintf(fpOut,"%.4f",clusters[cluster].clustVec[curDim]);
            if (curDim == dim -1 && cluster != K-1 )
            {
                fprintf(fpOut,"%s","\n");
            }else if (curDim != dim-1)
            {
                fprintf(fpOut,"%s",",");
            }     
        }
    }
    fprintf(fpOut,"%s","\n");
    fclose(fpOut);
    return 0; 
}

/*output: returns number of lines*/
int matLen(FILE *fp ){
    char ch;
    int vecCount =0;
    rewind(fp); /* goes to the beginning of the file*/
     do{
        ch = fgetc(fp);
        if(ch =='\n') vecCount ++;
    } while (ch != EOF);
    return vecCount ;
}
/* output: return the dimension*/
int dimCalc(FILE *fp){
    char ch = '0';
    int dimCount = 1;
    rewind(fp); /* goes to the beginning of th file*/
    while (ch!='\n')
    {
        ch = fgetc(fp);
        if(ch == ',') dimCount++;
    }
  return dimCount;
}
/*create node for each vec and insert it to it's closest cluster */
int initClustForVec( int numLines, int K, int dim, cluster* clusters, double** mat){
    int vec, newClust;
    node* newNodeVec;
    for(vec = 0; vec< numLines; vec++){
        newClust = argmin(mat[vec],clusters,K,dim);
        newNodeVec = (node*)malloc(sizeof(node));
        if(newNodeVec == NULL ) {
            printf("An Error Has Occurred\n");
            return 0;
        }
        newNodeVec->vec = mat[vec];
        addVecToClust(newNodeVec, &clusters[newClust]);
    }
    return 1;
}

/* for every vec save its old cluster & new cluster */ 
modify* updModify( cluster* clusters, modify* modArray, int dim, int K ){
    int newClust, clustIndex, vecCount =0;
    node* curVec;
    for(clustIndex =0; clustIndex< K; clustIndex++){
            curVec = clusters[clustIndex].vecLst;
            while (curVec != NULL)
            {
                newClust = argmin(curVec->vec,clusters,K,dim);
                modArray[vecCount].newClustIndex = newClust;
                modArray[vecCount].oldClustIndex=clustIndex;
                modArray[vecCount].nodeVec = curVec;
                vecCount++;
                curVec = curVec->next;
            }
    }
    return modArray;
}

/* update cluster's linked list according modArray */ 
void impChanges(cluster* clusters, int numLines, modify* modArray){
    int modNum,newClust,oldClust;
    for(modNum = 0; modNum< numLines; modNum++){
            newClust = modArray[modNum].newClustIndex;
            oldClust = modArray[modNum].oldClustIndex;
            if(newClust != oldClust ){
                swap(modArray[modNum].nodeVec, oldClust,newClust, clusters);
            }
        }
}

int KMean(int K, int maxIter, FILE* fpIn, FILE* fpOut){ 
    double maxDelta = EPSILON;
    int  iter =0, dim, numLines, memAloc;
    modify* modArray;
    double** mat; 
    cluster* clusters;
    dim = dimCalc(fpIn);
    numLines = matLen(fpIn);

    if(K > numLines){
        printf("Invalid Input!\n");
        return 1;
    }

    mat =initMat(fpIn,dim,numLines);
    if(mat == NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }
    clusters = initClusters(mat,K,dim);
    if(clusters == NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }
    memAloc= initClustForVec(numLines, K, dim, clusters,mat);
    if(memAloc == 0) {
        printf("An Error Has Occurred\n");
        return 1;
    }
    while (iter < maxIter && maxDelta >= EPSILON)
    {
        modArray =(modify*) calloc(numLines,sizeof(modify));
        if(modArray == NULL ) {
            printf("An Error Has Occurred\n");
            return 1;
        }
        modArray = updModify(clusters,modArray,dim,K); /* update modArray according new clusters */
        impChanges(clusters, numLines,modArray); /* link evrey vec to its new cluster according to modArray */ 
        maxDelta = updateClusters(clusters,K,dim); /* compute new max deltas */ 
        if(maxDelta == -1) {
            printf("An Error Has Occurred\n");
            return 1;
        }
        iter++;
        
        /*free modArray*/
        free(modArray);
        
        
    }
    write(clusters,K,fpOut,dim);
    freeMemory(mat,clusters,numLines,K);
    return 0;
}

/* submit args to vars, return 1 if successed else 0 */
int submitArgs(int argc, char **argv, FILE** fpIn, FILE** fpOut, double* k, double* maxIter){
    char* inputFile;
    char* outputFile;
    char* eptr;
    if (argc != 4 && argc != 5)
    {
        printf("Invalid Input!\n");
        return 1;
    }

    *k = strtod(argv[1], &eptr);

    /* if maxIter is not given */
    if (argc == 4)
    {
        *maxIter = DEFAULT_MAX_ITER;
        inputFile = argv[2];
        outputFile = argv[3];
    }

    /* if maxIter is given */
    if (argc == 5)
    {
        *maxIter = strtod(argv[2], &eptr);
        inputFile = argv[3];
        outputFile = argv[4];
    }

    *fpIn  = fopen(inputFile,"r");
    *fpOut = fopen(outputFile,"w");
    /* one of the file didnt open */
    if(fpIn == NULL || fpOut == NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }


    /* input check */
    if (*k <= 0|| *k != (int)*k || *maxIter <= 0|| *maxIter != (int)*maxIter || *fpIn == NULL || fpOut == NULL)
    {
        printf("Invalid Input!\n");
        return 1;
    }
    return 1;
}

int main(int argc, char **argv){ 
    double KDouble, maxIterDouble ;
    int K, maxIter;
    FILE* fpIn;
    FILE* fpOut;
    if(submitArgs(argc,argv,&fpIn,&fpOut,&KDouble,&maxIterDouble) == 0) return 1;
    K = (int) KDouble;
    maxIter = (int)maxIterDouble;
    KMean(K,maxIter,fpIn, fpOut);
    return 0;
}