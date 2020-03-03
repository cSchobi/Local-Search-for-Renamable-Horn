#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <string.h>

#ifdef DEBUG
#define DEBUG_PRINT(x) printf x
#else
#define DEBUG_PRINT(x) do {} while (0)
#endif

#ifdef DEBUG
#define DEBUG_PRINT_iARRAY(x, n) {int i=0; for(i=0; i < n; i++) printf("%i ", x[i]);} printf("\n")
#else
#define DEBUG_PRINT_iARRAY(x, n) do {} while (0)
#endif


#ifdef DEBUG
#define DEBUG_PRINT_dARRAY(x, n) {int i=0; for(i=0; i < n; i++) printf("%f ", x[i]);} printf("\n")
#else
#define DEBUG_PRINT_dARRAY(x, n) do {} while (0)
#endif


/*
 * Parameters
 * */

int skipNegLit = 0;
int nPoints = 10;
int print_remaining_clauses_flag = 0;
int print_optimum_flag = 0;
int print_status_flag = 0;
double pWalkSATSKC = 0.567;
double cbProbSATexp = 2.5;
double cbProbSATpoly = 2.3;
const double epsilon = 1;
int mode = -1; 

#define MAX_FLIPS 200
#define START_SIZE 16 // start size of occurrence lists of literals
#define TRIES 10 // number of trials for experiments
/*
 * ***Data structures***
 * */
typedef struct literal{
	int nPosClauses, nNegClauses;
	int *posClauses, *negClauses;
	int posSize, negSize;
} Literal;

typedef struct solver {
	int nClauses, nVars;
	int *memory, *flipped;
	int **clauses; // array with pointers into memory at the start and end of the clauses
	int *posLiterals; // counter for positive literals in each clause
	double *cumProbs; // buffer for cumulative probabilities in probSAT variant
	int *bufferVars; // buffer used in choosing literal
	int *bestFlipped, bestNoHorn;
	Literal *literals;
} Solver;

// Parsing
int parse(Solver* solver, char *fileName);
int skip_comment(FILE* file);

//Memory
void init(Solver* solver);
void allocate_memory(Solver* solver);
void free_solver(Solver* solver);

//Setup of solver
void add_clause(Solver* solver, int* buffer, int size, int nClause);
void addClauseToLiteral(Solver* solver, int nLit, int clause);
void setupSolver(Solver* solver);
void reset(Solver* solver);

//Printing
void print(Solver* solver);
void print_status(Solver* solver);
void print_remaining_clauses(Solver* solver, int i);
void print_usage();

// Search
void solve(Solver* solver, int (* getLiteral)(Solver* solver, int nClause), int maxFlips);
int getNonHornClause(Solver* solver);
int getBreakCount(Solver* solver, int nLit, int minBreak, int earlyBreak);
int getMakeCount(Solver* solver, int nLit);
int countPosLit(Solver* solver, int i);
int isHorn(Solver* solver, int nClause);
int isNonHorn(Solver* solver, int nClause);
int countHornClauses(Solver *solver);
int isPosLit(Solver* solver, int lit);
int isPosInClause(Solver* solver, int nClause, int lit);
void flipLiteral(Solver* solver, int nLit);
int areAllClausesHorn(Solver* solver);

// update solver
void updateOptimum(Solver* solver);
void updateSolver(Solver* solver, int lit);

// Literal choice
int getLitWalkSATSKC(Solver* solver, int nClause);
int getLitProbSAT(Solver* solver, int nClause, double (*f)(double breakCount));
int getLitProbSATexp(Solver* solver, int nClause);
int getLitProbSATpoly(Solver* solver, int nClause);
double probSATexp(double breakCount);
double probSATpoly(double breakCount);

//Experiments
void trendExperiment(Solver* solver, int maxFlips);
void trendSkipExperiment(Solver* solver, int maxFlips);
void walkSATExperiment(Solver* solver, int maxFlips, int nPoints);
void probSATexpExperiment(Solver* solver, int maxFlips, int nPoints);
void probSATpolyExperiment(Solver* solver, int maxFlips, int nPoints);
void parse_parameters(int argc, char *argv[]);

/* 
 * ************************************************************************
 * */
int main(int argc, char **argv){
	Solver solver;
	srand(time(NULL));
	init(&solver);
	parse_parameters(argc, argv);
	if(parse(&solver, argv[1]) == 1){
		setupSolver(&solver);
		DEBUG_PRINT(("skip neg literal: %d\n", skipNegLit));
		if(mode == 0){
			print_status_flag=1;
			print_remaining_clauses_flag=0;
			solve(&solver, getLitWalkSATSKC, 20);
		}else if(mode == 1){
			trendExperiment(&solver, MAX_FLIPS);
		}else if(mode == 2){
			walkSATExperiment(&solver, MAX_FLIPS, nPoints);
		}else if(mode == 3){
			nPoints = 20;
			probSATexpExperiment(&solver, MAX_FLIPS, nPoints);
		}else if(mode == 4){
			nPoints = 20;
			probSATpolyExperiment(&solver, MAX_FLIPS, nPoints);
		}else if(mode == 5){
			trendSkipExperiment(&solver, MAX_FLIPS);
		}else{
			printf("invalid mode\n");
			exit(0); 
		}
	}else
		printf("Error while parsing file\n");
	free_solver(&solver);
}

void parse_parameters(int argc, char *argv[]){
	int i;
	for(i = 2; i < argc; i++){ // first argument is filename
		if(strcmp("-skip", argv[i]) == 0){
			i++;
			if(strcmp("0", argv[i]) == 0){
				skipNegLit = 0;
			}else if(strcmp("1", argv[i]) == 0){
				skipNegLit = 1;
			}else{
				printf("Illegal argument for -skip\n");
				print_usage();
				exit(0);
			}
		}else if(strcmp("-mode", argv[i]) == 0){
			i++;
			mode = atoi(argv[i]);
		}else if(strcmp("-h", argv[i]) == 0){
			exit(0); 
		}
	}
}

void print_usage(){
	printf("Usage: fileName -mode x [optional parameters]\n");
	printf("Parameters:\n");
	printf("-mode = 0 testing, 1 trend experiment, 2 walkSATexperiment, \n");
	printf("\t3 probSatexponential experiment, 4  probSATpoly experiment\n");	
	printf("Optional parameters:\n");
	printf("-skip = 0 do not skip negative literals for literal selection, 1 skip negative literals\n");
}

int parse(Solver* solver, char *fileName){
	FILE* file;
	int c, nClause = 0, lit, size = 0;
	int* buffer;
	file = fopen(fileName, "r");
	if(!file){
		printf("Could not open file %s\n", fileName);
		return -1;
	}
	c = fgetc(file);
	while(c == 'c'){
		if(skip_comment(file) != 1){
			fclose(file);
			return -1;
		}
		c = fgetc(file);
	}
	ungetc(c, file);
	if(fscanf(file, "p cnf %i %i", &solver->nVars, &solver->nClauses) != 2)
		return -1;
	allocate_memory(solver);
	buffer = (int *) calloc(solver->nVars, sizeof(int));
	while(nClause < solver->nClauses){
		c = fgetc(file);
		if( c == 'c') {
			if(skip_comment(file) != 1) {
				free(buffer); fclose(file); return -1;
			} 
			continue;
		}
		else if (c == ' ' || c == '\n') continue;
		else {
			ungetc(c, file);
			if(fscanf(file, "%i", &lit) != 1) {free(buffer); fclose(file); return -1;}
			if(lit == 0){
				add_clause(solver, buffer, size, nClause);
				nClause++;
				size = 0;
			}else{
				buffer[size++] = lit;
			}
		}
	}	
	free(buffer);
	fclose(file);
	return 1;
}

void init(Solver* solver){
	solver->memory = NULL;
	solver->clauses = NULL;
	solver->flipped = NULL;
	solver->posLiterals = NULL;
	solver->cumProbs = NULL;
	solver->bestFlipped = NULL;
	solver->literals = NULL;
	solver->bufferVars = NULL;
	solver->nVars = 0;
	solver->nClauses = 0;
}

void allocate_memory(Solver* solver){
	int i;
	solver->memory = (int *) calloc(solver->nVars * solver->nClauses, sizeof(int)); // memory for all clauses
	solver->clauses = (int **) calloc(solver->nClauses + 1, sizeof(int *)); // pointers into memory at start and end of each clause (last pointer points to end of last clause)
	solver->clauses[0] = solver->memory;
	solver->flipped = (int *) calloc(solver->nVars , sizeof(int)); // flip status of literal (-1: literal is flipped, 1: literal is not flipped)
	// number of positive literals in each clause
	solver->posLiterals = (int *) calloc(solver->nClauses, sizeof(int));
	solver->cumProbs = (double *) calloc(solver->nVars, sizeof(double));
	solver->bestFlipped = (int *) calloc(solver->nVars, sizeof(int));
	solver->literals = (Literal *) calloc(solver->nVars, sizeof(Literal));
	solver->bufferVars = (int *) calloc(solver->nVars, sizeof(int));

	for(i = 0; i < solver->nVars; i++){
		solver->flipped[i] = 1;
		solver->bestFlipped[i] = 1;
		solver->literals[i].nPosClauses = 0;
		solver->literals[i].nNegClauses = 0;
		solver->literals[i].posSize = START_SIZE;
		solver->literals[i].negSize = START_SIZE;
		solver->literals[i].posClauses = calloc(START_SIZE, sizeof(int));
		solver->literals[i].negClauses = calloc(START_SIZE, sizeof(int));

	}
}

void add_clause(Solver* solver, int* buffer, int size, int nClause){
	int* clause = solver->clauses[nClause];
	int i;
	for(i = 0; i < size; i++){
		assert(-solver->nVars <= buffer[i] && buffer[i] <= solver->nVars);
		clause[i] = buffer[i];
	}	
	solver->clauses[nClause+1] = solver->clauses[nClause] + size;
	for(i = 0; i < size; i++){
		addClauseToLiteral(solver, clause[i], nClause);
	}
}

/*
 * nLit is in the range -nVars <= nLit <= -1 && 1 <= nLit <= nVars;
 * */
void addClauseToLiteral(Solver* solver, int nLit, int clause){
	Literal* lit = &solver->literals[abs(nLit)-1];
	if(nLit > 0){
		if(lit->posSize == lit->nPosClauses){
			lit->posSize *= 2;
			lit->posClauses = realloc(lit->posClauses , sizeof(int) * lit->posSize);
		}
		lit->posClauses[lit->nPosClauses] = clause;
		lit->nPosClauses++;
	}else{
		if(lit->negSize == lit->nNegClauses){
			lit->negSize *= 2;
			lit->negClauses = realloc(lit->negClauses , sizeof(int) * lit->negSize);
		}
		lit->negClauses[lit->nNegClauses] = clause;
		lit->nNegClauses++;
	}
}

int skip_comment(FILE* file){
	int ch;
	while((ch = getc(file)) != '\n');
	if(ch == EOF) return -1;
	else return 1;
}

void free_solver(Solver* solver){
	int i;
	for(i = 0; i < solver->nVars; i++){
		free(solver->literals[i].posClauses);
		free(solver->literals[i].negClauses);
	}
	free(solver->literals); solver->literals = NULL;
	free(solver->clauses); solver->clauses = NULL;
	free(solver->memory); solver->memory = NULL;
	solver->nClauses = solver->nVars = 0;
	free(solver->flipped); solver->flipped = NULL;
	free(solver->posLiterals); solver->posLiterals = NULL;
	free(solver->cumProbs); solver->cumProbs = NULL;
	free(solver->bestFlipped); solver->bestFlipped = NULL;
	free(solver->bufferVars); solver->bufferVars = NULL;
}

/*
 * initialize number of positive Literals and optimum number of horn clauses
 * */
void setupSolver(Solver* solver){
	int i,j, size;
	int* clause;
	solver->bestNoHorn = 0;
	for(i = 0; i < solver->nClauses; i++){
		solver->posLiterals[i] = 0;
		clause = solver->clauses[i];
		size = solver->clauses[i+1] - solver->clauses[i];
		for(j = 0; j < size; j++){
			solver->posLiterals[i] += clause[j] > 0 ? 1 : 0;
		}	
		solver->bestNoHorn += solver->posLiterals[i] < 2 ? 1 : 0;
	}	
}

void print(Solver* solver){
	int size, i, j;
	int* clause;
	if(!print_status_flag) return;	
	printf("p cnf %i %i\n", solver->nVars, solver->nClauses);
	for(i = 0; i < solver->nClauses; i++){
		printf("%3d ", i);
		clause = solver->clauses[i];
		size = solver->clauses[i+1] - solver->clauses[i];
		for(j = 0; j < size; j++){
			printf("%3d ", clause[j] * solver->flipped[abs(clause[j])-1]);
		}
		printf("|%2d\n", solver->posLiterals[i]);
	}
}

void solve(Solver* solver, int (* getLiteral)(Solver* solver, int nClause), int maxFlips){
	int i, lit, clause;
	print_status(solver);
	for(i = 0; i < maxFlips; i++){
		print_remaining_clauses(solver, i);
		DEBUG_PRINT(("Trie: %i\n", i));
		clause = getNonHornClause(solver);
		if(clause == -1){
		   DEBUG_PRINT(("Done\n"));
	   	   break;
		}
		DEBUG_PRINT(("chosen clause: %d\n", clause));
		lit = (*getLiteral)(solver, clause);
		updateSolver(solver, lit);
		flipLiteral(solver, lit);
		updateOptimum(solver);
		print_status(solver);
	}
	if(print_optimum_flag) printf(",%d", solver->nClauses - solver->bestNoHorn); // print optimum at the end
	DEBUG_PRINT(("Optimal number of horn clauses found: %i\n", solver->bestNoHorn));
	DEBUG_PRINT_iARRAY(solver->bestFlipped, solver->nVars);
}

void print_remaining_clauses(Solver* solver, int i){
		if(!print_remaining_clauses_flag) return;
		if(i != 0) printf(",");
		printf("%d", solver->nClauses - countHornClauses(solver));
}

int areAllClausesHorn(Solver* solver){
	int i;
	for(i = 0; i < solver->nClauses; i++){
		if(solver->posLiterals[i] > 1) return 0;
	}
	return 1;
}

/*
 * update solver such that counter for positive literals, break and create count
 * are still up to date when lit is flipped
 * */
void updateSolver(Solver* solver, int nLit){
	int i;
	Literal* lit = &solver->literals[abs(nLit)-1];
	for(i = 0; i < lit->nPosClauses; i++){
		if(solver->flipped[abs(nLit)-1] == 1)
			solver->posLiterals[lit->posClauses[i]]--;
		else
			solver->posLiterals[lit->posClauses[i]]++;
	}
	for(i = 0; i < lit->nNegClauses; i++){
		if(solver->flipped[abs(nLit)-1] == 1)
			solver->posLiterals[lit->negClauses[i]]++;
		else
			solver->posLiterals[lit->negClauses[i]]--;
	}
}

void updateOptimum(Solver* solver){
	int cntHorn = countHornClauses(solver), i;
	if(cntHorn > solver->bestNoHorn){
		solver->bestNoHorn = cntHorn;
		for(i = 0; i < solver->nVars; i++){
			solver->bestFlipped[i] = solver->flipped[i];
		}
	}
}

/*
int getMaxCreateValue(Solver* solver){
	int max = 0, i;
	for(i = 0; i < solver->nVars; i++){
		if(solver->createCount[i] > max) max = solver->createCount[i];
	}
	return max;
}
*/


void flipLiteral(Solver* solver, int nLit){
	DEBUG_PRINT(("Flipping literal: %i\n", abs(nLit)));
	solver->flipped[abs(nLit)-1] *= -1;
}

/*
 * iterate over all clauses and return the index of a random
 * nonHorn clause
 * return -1 if no nonHorn clause exists
 * */
int getNonHornClause(Solver* solver){
	int i, nNonHorn = 0;
	int isNonHornClauses[solver->nClauses];
	for(i = 0; i < solver->nClauses; i++){
		if(isNonHorn(solver, i)){
			isNonHornClauses[nNonHorn++] = i;
		}
	}
	DEBUG_PRINT(("Number of NonHorn clauses: %i\n", nNonHorn));
	if(nNonHorn > 0) return isNonHornClauses[rand() % nNonHorn];
	else return -1;
	
}

// return 1 if clause  is not a horn clause, 0 otherwise
int isNonHorn(Solver* solver, int nClause){
	return countPosLit(solver, nClause) > 1;
}

int isHorn(Solver* solver, int nClause){
	return 1 - isNonHorn(solver, nClause);
}

/*
 * count the number of positive literals in the clause
 * */
int countPosLit(Solver* solver, int nClause){
	int* clause = solver->clauses[nClause];
	int size = solver->clauses[nClause+1] - solver->clauses[nClause], i, literal, posLit = 0;
	for(i = 0; i < size; i++){
		literal = clause[i];
		if(solver->flipped[abs(literal)-1] * literal > 0) posLit++;
	}
	return posLit;
}

/*
 * return true if given literal is positive
 * */
int isPosLit(Solver* solver, int lit){
	return solver->flipped[abs(lit)-1] * lit > 0;
}

/*
 * returns true if given literal is positive in given clause
 * */
int isPosInClause(Solver* solver, int nClause, int lit){
	int	i = 0;
	int *clause = solver->clauses[nClause];
	int size = solver->clauses[nClause+1] - solver->clauses[nClause];
	while(abs(clause[i]) != abs(lit) && i < size) i++; // use size to avoid infinite loop
	if(i < size && abs(clause[i]) == abs(lit))
		return clause[i] * solver->flipped[abs(lit)-1] > 0;
	else // literal not found
		return 0;
}

/*
 * count the break value between previous flipping state and new flipping state;
 * posLiterals are the number of positive literals in previous flipping state;
 * posLiterals2 are for the new flipping state
 * return the number of horn clauses in posLiterals that become nonHorn in posLiterals
 * */
int getBreakValue(Solver* solver, int* posLiterals, int* posLiterals2){
	int i, breakValue = 0;
	for(i = 0; i < solver->nClauses; i++){
		if(posLiterals[i] <= 1 && posLiterals2[i] > 1)
			breakValue++;	
	}
	return breakValue;
}

void print_status(Solver* solver){
	if(!print_status_flag) return;
	int i;
	print(solver);
	printf("Horn clauses: %i\n", countHornClauses(solver));
	for(i = 0; i < solver->nVars; i++)
		printf("--");
	printf("\n");
	for(i = 0; i < solver->nVars; i++)
		printf("%3d ", i+1);
	printf("\n");
	for(i = 0; i < solver->nVars; i++)
		printf("%3d ", solver->flipped[i]);
	printf("| flip Status\n");
	for(i = 1; i <= solver->nVars; i++)
		printf("%3d ", getBreakCount(solver, i, 0, 0));
	printf("| break count\n");
	for(i = 1; i <= solver->nVars; i++)
		printf("%3d ", getMakeCount(solver, i));
	printf("| make count\n");
}

int countHornClauses(Solver *solver){
	int i, nHorn = 0;
	for(i = 0; i < solver->nClauses; i++){
		nHorn += isHorn(solver, i);
	}
	return nHorn;
}

void trendExperiment(Solver* solver, int maxFlips){
	int i,n;
	int (*heuristics[])(Solver* solver, int nClause) = {getLitWalkSATSKC, getLitProbSATexp, getLitProbSATpoly};
	char* heuristic_names[] = {"WalkSATSKC", "probSATexp", "probSATpoly"};
	n = sizeof(heuristics) / sizeof(heuristics[0]);
	print_remaining_clauses_flag = 1;
	print_optimum_flag = 1;
	for(i = 0; i < n; i++){
		printf("%s,", heuristic_names[i]);
		solve(solver, heuristics[i], maxFlips);
		printf("\n");
		reset(solver);
	}
}

void trendSkipExperiment(Solver* solver, int maxFlips){
	int i,skip,n;
	int (*heuristics[])(Solver* solver, int nClause) = {getLitWalkSATSKC, getLitProbSATexp, getLitProbSATpoly};
	char* heuristic_names[] = {"WalkSATSKC", "probSATexp", "probSATpoly"};
	n = sizeof(heuristics) / sizeof(heuristics[0]);
	print_remaining_clauses_flag = 1;
	print_optimum_flag = 1;
	for(i = 0; i < n; i++){
		for(skip = 0; skip <= 1; skip++){
			printf("%s[skip=%d],", heuristic_names[i], skip);
			skipNegLit = skip;
			solve(solver, heuristics[i], maxFlips);
			printf("\n");
			reset(solver);
		}
	}
}

void walkSATExperiment(Solver* solver, int maxFlips, int nPoints){
	pWalkSATSKC = 0.0;
	int i, j, sum;
	print_optimum_flag = 0;
	for(i = 0; i < nPoints; i++){
		sum = 0;
		for(j = 0; j < TRIES; j++){
			solve(solver, getLitWalkSATSKC, maxFlips);
			sum += solver->nClauses - solver->bestNoHorn;
			reset(solver);
		}
		printf(",%d", sum / TRIES);
		pWalkSATSKC += (double) 1 / (double) (nPoints-1);
	}	
	printf("\n");
}

void probSATexpExperiment(Solver* solver, int maxFlips, int nPoints){
	double startValue = 1.0, endValue = 4.0;
	int i, j, sum;
	print_optimum_flag = 0;
	cbProbSATexp = startValue;
	for(i = 0; i < nPoints; i++){
		sum = 0;
		for(j = 0; j < TRIES; j++){
			solve(solver, getLitProbSATexp, maxFlips);
			sum += solver->nClauses - solver->bestNoHorn;
			reset(solver);
		}
		printf(",%d", sum / TRIES);
		cbProbSATexp += (double) (endValue - startValue) / (double) (nPoints-1);
	}	
	printf("\n");
}

void probSATpolyExperiment(Solver* solver, int maxFlips, int nPoints){
	double startValue = 1.0, endValue = 4.0;
	int i, j, sum;
	print_optimum_flag = 0;
	cbProbSATpoly = startValue;
	for(i = 0; i < nPoints; i++){
		sum = 0;
		for(j = 0; j < TRIES; j++){
			solve(solver, getLitProbSATpoly, maxFlips);
			sum += solver->nClauses - solver->bestNoHorn;
			reset(solver);
		}
		printf(",%d", sum / TRIES);
		cbProbSATpoly += (double) (endValue - startValue) / (double) (nPoints-1);
	}	
	printf("\n");
}

void reset(Solver* solver){
	int i;
	for(i = 0; i < solver->nVars; i++){
		solver->flipped[i] = 1;
		solver->bestFlipped[i] = 1;
	}
	setupSolver(solver);
}

int getBreakCount(Solver* solver, int nLit, int minBreak, int earlyBreak){
	int breakCount = 0, i;
	Literal* lit = &solver->literals[abs(nLit)-1];
	if(solver->flipped[abs(nLit)-1] == 1){
		for(i = 0; i < lit->nNegClauses; i++){
			if(solver->posLiterals[lit->negClauses[i]] == 1) breakCount++;
			if(earlyBreak && breakCount > minBreak) return breakCount;
		}
	}else{
		for(i = 0; i < lit->nPosClauses; i++){
			if(solver->posLiterals[lit->posClauses[i]] == 1) breakCount++;
			if(earlyBreak && breakCount > minBreak) return breakCount;
		}
	}
	return breakCount;
}

/*
 * as an alternatie to break count, the number of non Horn clauses
 * that become horn might be a good heuristic
 * */
int getMakeCount(Solver* solver, int nLit){
	int makeCount = 0, i;
	Literal* lit = &solver->literals[abs(nLit)-1];
	if(solver->flipped[abs(nLit)-1] == -1){
		for(i = 0; i < lit->nNegClauses; i++){
			if(solver->posLiterals[lit->negClauses[i]] == 2) makeCount++;
		}
	}else{
		for(i = 0; i < lit->nPosClauses; i++){
			if(solver->posLiterals[lit->posClauses[i]] == 2) makeCount++;
		}
	}
	return makeCount;
}

int getLitWalkSATSKC(Solver* solver, int nClause){
	int* clause = solver->clauses[nClause];
	int size, i, min, lit, cnt = 0, breakCount;
	min = solver->nClauses + 1;
	size = solver->clauses[nClause+1] - solver->clauses[nClause];
	for(i = 0; i < size; i++){
		lit = clause[i];
		if(skipNegLit && !isPosLit(solver, lit)) continue; // should only negative literals be considered?
		breakCount = getBreakCount(solver, lit, min, 1);
		if(breakCount < min){
			min = breakCount;
			cnt = 0;
		}
		if(breakCount == min){
			solver->bufferVars[cnt] = lit;
			cnt++;
		}
	}
	DEBUG_PRINT(("Candidate literals: "));
	DEBUG_PRINT_iARRAY(solver->bufferVars, cnt);
	assert(cnt > 0);
	if(min == 0 || ((float) rand() / (float) RAND_MAX > pWalkSATSKC))
		return solver->bufferVars[rand() % cnt];
	else // if minimum is not zero take a random step with propability pWalkSATSKC
	{
		DEBUG_PRINT(("Taking random step\n"));
		return clause[rand() % size];				
	}
}

double probSATexp(double breakCount){
	return pow(cbProbSATexp, -breakCount);
}

int getLitProbSATexp(Solver* solver, int nClause){
	return getLitProbSAT(solver, nClause, probSATexp);
}

double probSATpoly(double breakCount){
	return pow(epsilon + breakCount, -cbProbSATpoly);
}

int getLitProbSATpoly(Solver* solver, int nClause){
	return getLitProbSAT(solver, nClause, probSATpoly);
}

int getLitProbSAT(Solver* solver, int nClause, double (*f)(double breakCount)){
	int* clause = solver->clauses[nClause];
	int size, i, lit, breakCount;
	double prob, r;
	size = solver->clauses[nClause+1] - solver->clauses[nClause];
	lit = clause[0];
	breakCount = getBreakCount(solver, lit, 0, 0);
	prob = f(breakCount);
	solver->cumProbs[0] = prob;
	for(i = 1; i < size; i++){
		lit = clause[i];
		breakCount = getBreakCount(solver, lit, 0, 0);
		prob = f(breakCount);
		solver->cumProbs[i] = solver->cumProbs[i-1] + prob;
	}
	r = ((double) rand() / (float) RAND_MAX) * solver->cumProbs[size-1];
	i = 0;
	while(solver->cumProbs[i] < r) i++;
	DEBUG_PRINT(("probabilities: "));
	DEBUG_PRINT_dARRAY(solver->cumProbs, size);
	DEBUG_PRINT(("random double: %f\n", r));
	return clause[i];
}
