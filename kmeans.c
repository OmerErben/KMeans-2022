#define PY_SSIZE_T_CLEAN
#include <Python.h>
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

/* function decleration*/
double** readDataFromFile(FILE* fp, int dim, int numLines);
double* readLineIntoVec(int dim, char* currNumberStart, double* vec );
cluster* initClusters(double** mat, int K, int dim);
int argmin(double* vec, cluster* clusters, int K, int dim);
void addVecToClust(node* nodeVec, cluster* clustVec);
void removeVecFromClust(node* nodeVec, cluster* clustVec);
void swap(node* nodeVec, int oldClust, int newClust, cluster* clusters);
double euqlideNorm(double* oldClust, double* newClust, int dim);
double updateClusters(cluster* clusters, int K, int dim);
int writeClust(cluster* clusters, int K, FILE* fpOut, int dim );
int matLen(FILE *fp);
int dimCalc(FILE *fp);
int initClustForVec(int numLines, int K, int dim, cluster* clusters, double** mat);
modify* updModify(cluster* clusters, modify* modArray, int dim, int K );
void impChanges(cluster* clusters, int numLines, modify* modArray);
int KMean(int K, int maxIter, double epsilon, char* dataFilename, char* clustFilename);
void freeMemory(double** mat, double** mat2, cluster* clusters, int matDim, int K);
int submitArgs(int argc, char **argv, FILE** fpIn, FILE** fpOut, double* K, double* maxIter);
static PyObject* fit(PyObject *self, PyObject *args);


static PyObject* fit(PyObject *self, PyObject *args){
    int K,maxIter,epsilon;
    char *dataFilename, *clustFilename;
    if(!PyArg_ParseTuple(args, "iidss", &K,&maxIter,&epsilon,&dataFilename, &clustFilename)) {
        printf("An Error Has Occurred\n");
        return NULL; 
    }
    return Py_BuildValue("i", KMean(K,maxIter,epsilon,dataFilename,clustFilename)); 
}

static PyMethodDef capiMethods[] = {
        {"k_means",                   /* the Python method name that will be used */
                (PyCFunction) fit, /* the C-function that implements the Python function and returns static PyObject*  */
                METH_VARARGS,           /* flags indicating parameters accepted for this function */
                        PyDoc_STR("kmeans ++")}, /*  The docstring for the function */
        {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

/* This initiates the module using the above definitions. */
static PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp", /* name of module */
        NULL, /* module documentation, may be NULL */
        -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};

/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_mykmeanssp(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}

double** readDataFromFile(FILE* fp, int dim, int numLines){
    int lineNumber;
    char *currNumberStart, *line;
    double* vec;
    double** mat = (double**)calloc(numLines , sizeof(double*));
    rewind(fp);
    line = (char*)calloc(1024,sizeof(char));
    if (line == NULL || mat == NULL){
        printf("An Error Has Occurred\n");
        return NULL;
    }
    for(lineNumber =0 ; lineNumber < numLines; lineNumber++){
        currNumberStart = line;
        if(fscanf(fp, "%s",line)== EOF){
            printf("An Error Has Occurred\n");
            return NULL;
        }
        vec = (double*) calloc(dim,sizeof(double));
        if (vec == NULL){
            printf("An Error Has Occurred\n");
            return NULL;
        }
        mat[lineNumber]= readLineIntoVec(dim,currNumberStart,vec);
    }
    free(line);
    return mat;
}


/* read line from file and write it to vec */
double* readLineIntoVec(int dim, char* currNumberStart, double* vec ){
    char* currComma;
    int wordNum;
    for ( wordNum =0; wordNum < dim -1; wordNum++){
        currComma = strchr((char*)currNumberStart,','); /* point to the first "comma" in line */
        *currComma = '\0'; 
        vec[wordNum] = atof(currNumberStart); /* insert word to vec */
        currNumberStart = ++currComma; /* update line to tne next char after , */
    }
    vec[dim-1] = atof(currNumberStart);
    return vec;
}

/* free memory at the end of KMeans */
void freeMemory(double** mat, double** mat2, cluster* clusters, int matDim, int K){
    node* curr;
    node* step;
    int i;
    for (i = 0; i < matDim; i++)
    {
        free(mat[i]);
    }
    free(mat);
    for (i = 0; i < K; i++)
    {
        free(mat2[i]);
    }
    free(mat2);
    for (i = 0; i < K; i++)
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
       for(curDim = 0; curDim < dim; curDim++){
           sum += pow((vec[curDim] - clusters[clustIndex].clustVec[curDim]),2);
       }
       if (clustIndex == 0) minSum = sum;
        else if(sum < minSum) {
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
    double deltaMax = 0, delta;
    int clusterIndex, vecLstLen, curDim;
    double* newClust;
    double* oldClust;
    node* currVec;
    
    for(clusterIndex = 0; clusterIndex< K; clusterIndex++){
        vecLstLen = 0;
        newClust = (double*)calloc(dim, sizeof(double));
        if(newClust == NULL){
            printf("An Error Has Occurred\n");
            return -1;
        }
        currVec = clusters[clusterIndex].vecLst;
        while (currVec != NULL)
        {
            for(curDim = 0; curDim < dim; curDim++){
                newClust[curDim] += currVec->vec[curDim];
            }
            vecLstLen++;
            currVec = currVec->next;
        }
        /* we assume vecLstLen != 0 becuse clustVec dosent have an empty list */
        for(curDim = 0; curDim < dim; curDim++){
            newClust[curDim] = newClust[curDim] / vecLstLen;
        }
        oldClust = clusters[clusterIndex].clustVec;
        delta = euqlideNorm(oldClust, newClust, dim);
        if (delta > deltaMax) deltaMax = delta;
        free(oldClust);
        clusters[clusterIndex].clustVec = newClust;
    }
    return deltaMax;
}

/* output: create output.txt and write to it*/
int writeClust(cluster* clusters, int K, FILE* fpOut, int dim ){
    int cluster, curDim;
    for(cluster = 0; cluster < K; cluster++){
        for(curDim = 0; curDim < dim; curDim++){
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
    return 0; 
}

/*output: returns number of lines*/
int matLen(FILE *fp ){
    char ch;
    int vecCount = 0;
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
    while (ch != '\n')
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
    for(vec = 0; vec < numLines; vec++){
        newClust = argmin(mat[vec],clusters,K,dim);
        newNodeVec = (node*)malloc(sizeof(node));
        if(newNodeVec == NULL) {
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
    int newClust, clustIndex, vecCount = 0;
    node* curVec;
    for(clustIndex = 0; clustIndex < K; clustIndex++){
            curVec = clusters[clustIndex].vecLst;
            while (curVec != NULL)
            {
                newClust = argmin(curVec->vec,clusters,K,dim);
                modArray[vecCount].newClustIndex = newClust;
                modArray[vecCount].oldClustIndex = clustIndex;
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
    for(modNum = 0; modNum < numLines; modNum++){
            newClust = modArray[modNum].newClustIndex;
            oldClust = modArray[modNum].oldClustIndex;
            if(newClust != oldClust ){
                swap(modArray[modNum].nodeVec, oldClust, newClust, clusters);
            }
        }
}

int KMean(int K, int maxIter, double epsilon, char* dataFilename, char* clustFilename){ 
    double maxDelta = epsilon;
    int  iter = 0, dim, numLines, memAloc;
    modify* modArray;
    double** mat; 
    double** mat2;
    cluster* clusters;
    FILE *dataFile, *clustFile;
    dataFile = fopen(dataFilename, "r");
    clustFile = fopen(clustFilename, "r");
    dim = dimCalc(dataFile);
    numLines = matLen(dataFile);
    
    if(K > numLines){
        printf("Invalid Input!\n");
        return 1;
    }
    mat =readDataFromFile(dataFile, dim, numLines);
    if(mat == NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }
    mat2 = readDataFromFile(clustFile, dim, K);
    clusters = initClusters(mat2,K,dim);
    if(clusters == NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }
    memAloc= initClustForVec(numLines, K, dim, clusters, mat);
    if(memAloc == 0) {
        printf("An Error Has Occurred\n");
        return 1;
    }
    fclose(clustFile);
    fclose(dataFile);
    
    while (iter < maxIter && maxDelta >= epsilon)
    {
        modArray =(modify*) calloc(numLines,sizeof(modify));
        if(modArray == NULL ) {
            printf("An Error Has Occurred\n");
            return 1;
        }
        modArray = updModify(clusters, modArray, dim, K); /* update modArray according new clusters */
        impChanges(clusters, numLines, modArray); /* link evrey vec to its new cluster according to modArray */ 
        maxDelta = updateClusters(clusters, K, dim); /* compute new max deltas */ 
        if(maxDelta == -1) {
            printf("An Error Has Occurred\n");
            return 1;
        }
        iter++;
        
        /*free modArray*/
        free(modArray);
        
        
    }
    clustFile = fopen(clustFilename, "w");
    writeClust(clusters, K, clustFile, dim);
    fclose(clustFile);
    freeMemory(mat, mat2, clusters, numLines, K);
    return 0;
}

/* submit args to vars, return 1 if successed else 0 */
int submitArgs(int argc, char **argv, FILE** fpIn, FILE** fpOut, double* K, double* maxIter){
    char* inputFile = NULL;
    char* outputFile = NULL;
    char* eptr;
    if (argc != 4 && argc != 5)
    {
        printf("Invalid Input!\n");
        return 1;
    }

    *K = strtod(argv[1], &eptr);

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
    if (*K <= 0 || *K != (int)*K || *maxIter <= 0 || *maxIter != (int)*maxIter || *fpIn == NULL || fpOut == NULL)
    {
        printf("Invalid Input!\n");
        return 1;
    }
    return 1;
}