#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <string.h>

/*
 * Parameters
 * */

int skipNegLit = 0;
int nPoints = 20;
int print_remaining_clauses_flag = 0;
int print_optimum_flag = 0;
int print_status_flag = 0;
double pWalkSATSKC = 0.0; //0.567; //best values for SAT according to literature
double pWalkSATWithMake = 1.0;
double cmProbSATexp = 8.0;
double cbProbSATexp = 4.0; //2.5;
double cmProbSATpoly = 5.0;
double cbProbSATpoly = 6.0; //2.3;
const double epsilon = 1;
int mode = -1; 
const int trackBestSolution_flag = 0;

#define MAX_FLIPS 100000
#define MAX_CLAUSE_SIZE 6
#define START_SIZE 16 // start size of occurrence lists of literals
#define TRIES 10 // number of trials for experiments


/* ***********************
 * *** Data structures ***
 * ***********************
 * */
typedef struct literal{
	int nPosClauses, nNegClauses; // number of clauses the literal appears in (positively and negatively)
	int *posClauses, *negClauses; // arrays with numbers of clauses where literal appears positively/negatively
	int posSize, negSize; // size of respective array
} Literal;

typedef struct solver {
	int nClauses, nVars;
	int *memory, *flipped;
	int mem_size;
	int **clauses; // array with pointers into memory at the start and end of the clauses
	int *posLiterals; // counter for positive literals in each clause
	int *nonHorn, *whereNonHorn; // array with clauses that are nonHorn and array which contains position of each clause in nonHorn
	int numNonHorn;
	double *cumProbs; // buffer for cumulative probabilities in probSAT variant
	int *bufferVars; // buffer used in choosing literal
	int *bestFlipped, bestNumHorn; // best number of horn clauses
	Literal *literals;
} Solver;

// Parsing
void parse_parameters(int argc, char *argv[]);
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
int countHornClauses(Solver *solver);
int isPosLit(Solver* solver, int lit);
void flipLiteral(Solver* solver, int nLit);

// update solver
void updateOptimum(Solver* solver);
void updateSolver(Solver* solver, int lit);
void addNonHornClause(Solver* solver, int nClause);
void removeNonHornClause(Solver* solver, int nClause);

// Literal choice
int getLitWalkSATSKC(Solver* solver, int nClause);
int getLitWalkSATWithMake(Solver* solver, int nClause);
int getLitProbSAT(Solver* solver, int nClause, double (*f)(int makeCount, int breakCount));
int getLitProbSATexp(Solver* solver, int nClause);
int getLitProbSATpoly(Solver* solver, int nClause);
double probSATexp(int makeCount, int breakCount);
double probSATpoly(int makeCount, int breakCount);

//Experiments
void trendExperiment(Solver* solver, int maxFlips);
void boxPlotExperiment(Solver* solver, int maxFlips);
void walkSATExperiment(Solver* solver, int maxFlips, int nPoints);
void walkSATWithMakeExperiment(Solver* solver, int maxFlips, int nPoints);
void probSATexpExperiment(Solver* solver, int maxFlips, int nPoints);
void probSATpolyExperiment(Solver* solver, int maxFlips, int nPoints);

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
		if(mode == 0){
			printf("Testing mode\n");
			print_status_flag=1;
			print_remaining_clauses_flag=1;
			print_optimum_flag = 1;
			solve(&solver, getLitProbSATexp, 200);
		}else if(mode == 1){
			trendExperiment(&solver, 3 * solver.nVars);
		}else if(mode == 2){
			walkSATExperiment(&solver, 3 * solver.nVars, nPoints);
		}else if(mode == 3){
			probSATexpExperiment(&solver, 3 * solver.nVars, nPoints);
		}else if(mode == 4){
			probSATpolyExperiment(&solver, 3 * solver.nVars, nPoints);
		}else if(mode == 5){
			boxPlotExperiment(&solver, 3 * solver.nVars);
		}else if(mode == 6){
			walkSATWithMakeExperiment(&solver, 3 * solver.nVars, nPoints);
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
			print_usage();
			exit(0); 
		}
	}
}

void print_usage(){
	printf("Usage: fileName -mode x [optional parameters]\n");
	printf("Parameters:\n");
	printf("-mode = 0 testing, 1 trend experiment, 2 walkSATexperiment, \n");
	printf("\t3 probSatexponential experiment, 4  probSATpoly experiment\n");	
	printf("\t5 poxPLot experiment, 6 walkSAT with make count\n");
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
	solver->nonHorn = NULL;
	solver->whereNonHorn = NULL;
	solver->cumProbs = NULL;
	solver->bestFlipped = NULL;
	solver->literals = NULL;
	solver->bufferVars = NULL;
	solver->nVars = 0;
	solver->nClauses = 0;
	solver->numNonHorn = 0;
}

void allocate_memory(Solver* solver){
	int i;
	solver->mem_size = MAX_CLAUSE_SIZE * solver->nClauses;
	solver->memory = (int *) calloc(solver->mem_size, sizeof(int)); // memory for all literals
	solver->clauses = (int **) calloc(solver->nClauses + 1, sizeof(int *)); // pointers into memory at start and end of each clause (last pointer points to end of last clause)
	solver->clauses[0] = solver->memory;
	solver->flipped = (int *) calloc(solver->nVars , sizeof(int)); // flip status of literal (-1: literal is flipped, 1: literal is not flipped)
	solver->posLiterals = (int *) calloc(solver->nClauses, sizeof(int)); // number of positive literals in each clause
	solver->nonHorn = (int *) calloc(solver->nClauses, sizeof(int)); // number of positive literals in each clause
	solver->whereNonHorn = (int *) calloc(solver->nClauses, sizeof(int)); // number of positive literals in each clause
	solver->cumProbs = (double *) calloc(solver->nVars, sizeof(double));
	solver->bestFlipped = (int *) calloc(solver->nVars, sizeof(int));
	solver->literals = (Literal *) calloc(solver->nVars, sizeof(Literal));
	solver->bufferVars = (int *) calloc(solver->nVars, sizeof(int));

	for(i = 0; i < solver->nVars; i++){ // init literals
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

/*
 * add new clause contained in buffer to the solver
 * */
void add_clause(Solver* solver, int* buffer, int size, int nClause){
	int* clause = solver->clauses[nClause];
	assert(clause - solver->clauses[0] + size < solver->mem_size);
	int i;
	for(i = 0; i < size; i++){
		assert((-solver->nVars <= buffer[i] && buffer[i] <= 1) || (1 <= buffer[i] && buffer[i] <= solver->nVars));
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
	solver->mem_size = solver->nClauses = solver->nVars = 0;
	free(solver->flipped); solver->flipped = NULL;
	free(solver->posLiterals); solver->posLiterals = NULL;
	free(solver->nonHorn); solver->nonHorn = NULL;
	free(solver->whereNonHorn); solver->whereNonHorn = NULL;
	free(solver->cumProbs); solver->cumProbs = NULL;
	free(solver->bestFlipped); solver->bestFlipped = NULL;
	free(solver->bufferVars); solver->bufferVars = NULL;
}

/*
 * initialize number of positive Literals, nonHornClauses and optimum number of horn clauses
 * */
void setupSolver(Solver* solver){
	int i,j, size;
	int* clause;
	solver->bestNumHorn = 0;
	solver->numNonHorn = 0;
	for(i = 0; i < solver->nClauses; i++){
		solver->posLiterals[i] = 0;
		clause = solver->clauses[i];
		size = solver->clauses[i+1] - solver->clauses[i];
		for(j = 0; j < size; j++){
			solver->posLiterals[i] += clause[j] > 0 ? 1 : 0;
		}	
		if(solver->posLiterals[i] > 1)
			addNonHornClause(solver, i);
		solver->bestNumHorn += solver->posLiterals[i] < 2 ? 1 : 0;
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
		clause = getNonHornClause(solver);
		if(clause == -1){
	   	   break;
		}
		lit = (*getLiteral)(solver, clause);
		updateSolver(solver, lit);
		flipLiteral(solver, lit);
		updateOptimum(solver);
		print_status(solver);
	}
	if(print_optimum_flag) printf(",%d", solver->nClauses - solver->bestNumHorn); // print optimum at the end
}

void print_remaining_clauses(Solver* solver, int i){
		if(!print_remaining_clauses_flag) return;
		if(i != 0) printf(",");
		printf("%d", solver->numNonHorn);
}

/*
 * update solver such that counters for positive literals
 * are still up to date when lit is flipped
 * */
void updateSolver(Solver* solver, int nLit){
	int i, nClause;
	Literal* lit = &solver->literals[abs(nLit)-1];

	// update number of positive literals and the array with non horn clauses
	for(i = 0; i < lit->nPosClauses; i++){
		nClause = lit->posClauses[i];
		if(solver->flipped[abs(nLit)-1] == 1){
			solver->posLiterals[nClause]--;
			if(solver->posLiterals[nClause] == 1)
				removeNonHornClause(solver, nClause);
		}else{
			solver->posLiterals[nClause]++;
			if(solver->posLiterals[nClause] == 2)
				addNonHornClause(solver, nClause);
		}
	}
	for(i = 0; i < lit->nNegClauses; i++){
		nClause = lit->negClauses[i];
		if(solver->flipped[abs(nLit)-1] == 1){
			solver->posLiterals[nClause]++;
			if(solver->posLiterals[nClause] == 2)
				addNonHornClause(solver, nClause);
		}else{
			solver->posLiterals[nClause]--;
			if(solver->posLiterals[nClause] == 1)
				removeNonHornClause(solver, nClause);
		}
	}
}

/*
 * add the clause to the array which contains all nonHorn clauses;
 * also keep track where the clause is stored in the array
 * */
void addNonHornClause(Solver* solver, int nClause){
	solver->nonHorn[solver->numNonHorn] = nClause;
	solver->whereNonHorn[nClause] = solver->numNonHorn;
	solver->numNonHorn++;
}

void removeNonHornClause(Solver* solver, int nClause){
	solver->numNonHorn--;
	solver->nonHorn[solver->whereNonHorn[nClause]] = solver->nonHorn[solver->numNonHorn];
	solver->whereNonHorn[solver->nonHorn[solver->numNonHorn]] = solver->whereNonHorn[nClause];
}

void updateOptimum(Solver* solver){
	int cntHorn = solver->nClauses - solver->numNonHorn, i;
	if(cntHorn > solver->bestNumHorn){
		solver->bestNumHorn = cntHorn;
		if(trackBestSolution_flag){
			for(i = 0; i < solver->nVars; i++){
				solver->bestFlipped[i] = solver->flipped[i];
			}
		}
	}
}

void flipLiteral(Solver* solver, int nLit){
	solver->flipped[abs(nLit)-1] *= -1;
}

/*
 * iterate over all clauses and return the index of a random
 * nonHorn clause
 * return -1 if no nonHorn clause exists
 * */
int getNonHornClause(Solver* solver){
	if(solver->numNonHorn == 0) return -1;
	else return solver->nonHorn[rand() % solver->numNonHorn];
}

/*
 * return true if given literal is positive
 * */
int isPosLit(Solver* solver, int lit){
	return solver->flipped[abs(lit)-1] * lit > 0;
}

void print_status(Solver* solver){
	if(!print_status_flag) return;
	int i;
	print(solver);
	printf("Horn clauses: %i\n", countHornClauses(solver));
	for(i = 0; i < solver->nVars; i++)
		printf("--");
	printf("\n    ");
	for(i = 0; i < solver->nVars; i++)
		printf("%3d ", i+1);
	printf("\n    ");
	for(i = 0; i < solver->nVars; i++)
		printf("%3d ", solver->flipped[i]);
	printf("| flip Status\n    ");
	for(i = 1; i <= solver->nVars; i++)
		printf("%3d ", getBreakCount(solver, i, 0, 0));
	printf("| break count\n    ");
	for(i = 1; i <= solver->nVars; i++)
		printf("%3d ", getMakeCount(solver, i));
	printf("| make count\n");
	printf("numNonHorn: %d\n", solver->numNonHorn);
	printf("Non Horn:\n");
	for(i = 0; i < solver->nClauses; i++)
		printf("%3d ", solver->nonHorn[i]);
	printf("\nWhere non Horn\n");
	for(i = 0; i < solver->nClauses; i++)
		printf("%3d ", solver->whereNonHorn[i]);
	printf("\n");
}

int countHornClauses(Solver *solver){
	return solver->nClauses - solver->numNonHorn;
}

/* **************************************
 * *******    Experiments ***************
 * **************************************
 * */

void trendExperiment(Solver* solver, int maxFlips){
	int i,n;
	int (*heuristics[])(Solver* solver, int nClause) = 
		{getLitWalkSATSKC, getLitWalkSATWithMake, getLitProbSATexp, getLitProbSATpoly};
	char* heuristic_names[] = {"WalkSATSKC", "WalkSATWithMake", "probSATexp", "probSATpoly"};
	n = sizeof(heuristics) / sizeof(heuristics[0]);
	print_remaining_clauses_flag = 1;
	print_optimum_flag = 0;
	for(i = 0; i < n; i++){
		printf("%s,", heuristic_names[i]);
		solve(solver, heuristics[i], maxFlips);
		printf("\n");
		reset(solver);
	}
}

void boxPlotExperiment(Solver* solver, int maxFlips){
	int i,j,n;
	int (*heuristics[])(Solver* solver, int nClause) = 
		{getLitWalkSATSKC, getLitWalkSATWithMake, getLitProbSATexp, getLitProbSATpoly};
	char* heuristic_names[] = {"WalkSATSKC", "WalkSATWithMake", "probSATexp", "probSATpoly"};
	n = sizeof(heuristics) / sizeof(heuristics[0]);
	print_remaining_clauses_flag = 0;
	print_optimum_flag = 1;
	for(i = 0; i < n; i++){
		for(skipNegLit = 0; skipNegLit < 2; skipNegLit++){
			printf("%s[skip=%d]", heuristic_names[i], skipNegLit);
			for(j = 0; j < TRIES; j++){
				solve(solver, heuristics[i], maxFlips);
				reset(solver);
			}
			printf("\n");
		}
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
			sum += solver->nClauses - solver->bestNumHorn;
			reset(solver);
		}
		printf(",%d", sum / TRIES);
		pWalkSATSKC += (double) 1 / (double) (nPoints-1);
	}	
	printf("\n");
}

void walkSATWithMakeExperiment(Solver* solver, int maxFlips, int nPoints){
	pWalkSATWithMake = 0.0;
	int i, j, sum;
	print_optimum_flag = 0;
	for(i = 0; i < nPoints; i++){
		sum = 0;
		for(j = 0; j < TRIES; j++){
			solve(solver, getLitWalkSATWithMake, maxFlips);
			sum += solver->nClauses - solver->bestNumHorn;
			reset(solver);
		}
		printf(",%d", sum / TRIES);
		pWalkSATWithMake += (double) 1 / (double) (nPoints-1);
	}	
	printf("\n");
}

void probSATexpExperiment(Solver* solver, int maxFlips, int nPoints){
	double cbStart = 1.0, cbEnd = 20.0;
	double cmStart = 0.0, cmEnd = 20.0;
	int i, j, k, sum;
	print_optimum_flag = 0;
	cbProbSATexp = cbStart;
	for(k = 0; k < nPoints; k++){
		cmProbSATexp = cmStart;
		for(i = 0; i < nPoints; i++){
			sum = 0;
			for(j = 0; j < TRIES; j++){
				solve(solver, getLitProbSATexp, maxFlips);
				sum += solver->nClauses - solver->bestNumHorn;
				reset(solver);
			}
			if(i != 0) printf(",");
		   	printf("%d", sum / TRIES);
			cmProbSATexp += (cmEnd - cmStart) / (double) (nPoints - 1);
		}	
		cbProbSATexp += (cbEnd - cbStart ) / (double) (nPoints-1);
		printf("\n");
	}
}

void probSATpolyExperiment(Solver* solver, int maxFlips, int nPoints){
	double cbStart = -2.0, cbEnd = 10.0;
	double cmStart = -2.0, cmEnd = 10.0;
	int i, j, k, sum;
	print_optimum_flag = 0;
	cbProbSATpoly = cbStart;
	for(k = 0; k < nPoints; k++){
		cmProbSATpoly = cmStart;
		for(i = 0; i < nPoints; i++){
			sum = 0;
			for(j = 0; j < TRIES; j++){
				solve(solver, getLitProbSATpoly, maxFlips);
				sum += solver->nClauses - solver->bestNumHorn;
				reset(solver);
			}
			if(i != 0) printf(",");
			printf("%d", sum / TRIES);
			cmProbSATpoly += (cmEnd - cmStart) / (double) (nPoints - 1);
		}	
		cbProbSATpoly += (cbEnd - cbStart) / (double) (nPoints-1);
		printf("\n");
	}
}

void reset(Solver* solver){
	int i;
	for(i = 0; i < solver->nVars; i++){
		solver->flipped[i] = 1;
		solver->bestFlipped[i] = 1;
	}
	setupSolver(solver);
}

/*
 * get break count of literal
 * if earlyBreak is TRUE then function will return as soon as the break count 
 * exceeds the minimum break count
 * */
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
 * as an alternative to break count, the number of non Horn clauses
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

/* ******************************************
 * **** heuristics for literal selection ****
 * ******************************************
 * */
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
	assert(cnt > 0);
	if(min == 0 || (((double) rand() / (double) RAND_MAX) > pWalkSATSKC))
		return solver->bufferVars[rand() % cnt];
	else // if minimum is not zero take a random step with propability pWalkSATSKC
	{
		return clause[rand() % size];				
	}
}

int getLitWalkSATWithMake(Solver* solver, int nClause){
	int* clause = solver->clauses[nClause];
	int size, i, max , lit, cnt = 0, value;
	max = -solver->nClauses;
	size = solver->clauses[nClause+1] - solver->clauses[nClause];
	for(i = 0; i < size; i++){
		lit = clause[i];
		if(skipNegLit && !isPosLit(solver, lit)) continue; // should only negative literals be considered?
		value = getMakeCount(solver, lit) - getBreakCount(solver, lit, 0, 0);
		if(value > max){
			max = value;
			cnt = 0;
		}
		if(value == max){
			solver->bufferVars[cnt] = lit;
			cnt++;
		}
	}
	assert(cnt > 0);
	if(max >= 0 || ((double) rand() / (double) RAND_MAX > pWalkSATWithMake)) // README: maybe try different cutoff values for max
		return solver->bufferVars[rand() % cnt];
	else // if minimum is not zero take a random step with propability pWalkSATSKC
	{
		return clause[rand() % size];				
	}
	
}

int getLitProbSAT(Solver* solver, int nClause, double (*f)(int makeCount, int breakCount)){
	int* clause = solver->clauses[nClause];
	int size, i, lit, breakCount, makeCount, cnt;
	double prob, r;
	size = solver->clauses[nClause+1] - solver->clauses[nClause];
	for(cnt = 0, i = 0; i < size; i++){
		lit = clause[i];
		if(skipNegLit && !isPosLit(solver, lit)) continue;
		makeCount = getMakeCount(solver, lit);
		breakCount = getBreakCount(solver, lit, 0, 0);
		prob = f(makeCount, breakCount);
		if(cnt == 0)
			solver->cumProbs[cnt] = prob;
		else
			solver->cumProbs[cnt] = solver->cumProbs[cnt-1] + prob;
		cnt++;
	}
	r = ((double) rand() / (double) RAND_MAX) * solver->cumProbs[cnt-1];
	i = 0;
	while(solver->cumProbs[i] < r) i++;
	if(skipNegLit){
		cnt = 0;
		lit = 0;
		while(!isPosLit(solver, clause[lit])) lit++;
		while(cnt < i){
			lit++;
			while(!isPosLit(solver, clause[lit])) lit++;
			cnt++;
		}
		return clause[lit];
	}
	else
		return clause[i];
}

double probSATexp(int makeCount, int breakCount){
	return pow(cmProbSATexp, makeCount) / pow(cbProbSATexp, breakCount);
}

int getLitProbSATexp(Solver* solver, int nClause){
	return getLitProbSAT(solver, nClause, probSATexp);
}

double probSATpoly(int makeCount, int breakCount){
	return pow(makeCount, cmProbSATpoly) / pow(epsilon + breakCount, cbProbSATpoly);
}

int getLitProbSATpoly(Solver* solver, int nClause){
	return getLitProbSAT(solver, nClause, probSATpoly);
}

