#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

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

typedef struct literal{
	int nClauses;
	int *clauses;
	int size;
} Literal;

typedef struct solver {
	int nClauses, nVars;
	int *memory, *flipped;
	int **clauses; // array with pointers into memory at the start and end of the clauses
	int *breakCount;
	int *posLiterals;
	int *cumProbs; // buffer for cumulative probabilities in probSAT variant
	int *optFlip;
	int optNHorn;
	int *createCount;
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
void updateSolver(Solver* solver, int lit);

//Printing
void print(Solver* solver);
void print_status(Solver* solver);

// Search
void solve_incrementally(Solver* solver, int maxTries);
int getNonHornClause(Solver* solver);
int getLitInClause(Solver* solver, int nClause);
int getLit(Solver* solver);
int getLitProb(Solver* solver);
int chooseLiteral(int* clause, int* breakValues, int size, int min);
int getBreakValue(Solver* solver, int* posLiterals, int* posLiterals2);
int* countPosLits(Solver* solver);
int countPosLit(Solver* solver, int i);
int isHorn(Solver* solver, int nClause);
int isNonHorn(Solver* solver, int nClause);
int countHornClauses(Solver *solver);
int isPosLit(Solver* solver, int lit);
int isPosInClause(Solver* solver, int nClause, int lit);
void flipLiteral(Solver* solver, int nLit);
int getMinBreakCount(Solver* solver);
void updateClause(Solver* solver, int nClause, int lit);
void updateBreakCreateCount(Solver* solver, int nClause, int lit);
void updateOptimum(Solver* solver);
int allClausesAreHorn(Solver* solver);

//Experiments
void run_experiments(Solver* solver, int maxFlips);
void solve(Solver* solver, int (* heuristic)(Solver* solver), int maxFlips);

//Heuristics
int getLitByBreak(Solver* solver);
int getLitByBreakProb(Solver* solver);
int getLitByCreate(Solver* solver);
int getLitByCreateProb(Solver* solver);

/* 
 * ************************************************************************
 * */
int main(int argc, char **argv){
	Solver solver;
	if(argc != 2){ printf("Usage: %s fileName\n", argv[0]); exit(0); }
	init(&solver);
	if(parse(&solver, argv[1]) == 1){
		setupSolver(&solver);
		//run_experiments(&solver, 100);
		solve_incrementally(&solver, 500);
		
	}else
		printf("Error while parsing file\n");
	free_solver(&solver);
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
	solver->breakCount = NULL;
	solver->posLiterals = NULL;
	solver->cumProbs = NULL;
	solver->optFlip = NULL;
	solver->createCount = NULL;
	solver->literals = NULL;
	solver->nVars = 0;
	solver->nClauses = 0;
}

void allocate_memory(Solver* solver){
	int i;
	solver->memory = (int *) calloc(solver->nVars * solver->nClauses, sizeof(int)); // memory for all clauses
	solver->clauses = (int **) calloc(solver->nClauses + 1, sizeof(int *)); // pointers into memory at start and end of each clause (last pointer points to end of last clause)
	solver->clauses[0] = solver->memory;
	solver->flipped = (int *) calloc(solver->nVars , sizeof(int)); // flip status of literal (-1: literal is flipped, 1: literal is not flipped)
	// breakCount for each literal 
	solver->breakCount = (int *) calloc(solver->nVars, sizeof(int));
	// number of positive literals in each clause
	solver->posLiterals = (int *) calloc(solver->nClauses, sizeof(int));
	solver->cumProbs = (int *) calloc(solver->nVars, sizeof(int));
	solver->optFlip = (int *) calloc(solver->nVars, sizeof(int));
	solver->createCount = (int *) calloc(solver->nVars, sizeof(int));
	solver->literals = (Literal *) calloc(solver->nVars, sizeof(Literal));

	for(i = 0; i < solver->nVars; i++){
		solver->flipped[i] = 1;
		solver->optFlip[i] = 1;
		solver->breakCount[i] = 0;
		solver->createCount[i] = 0;
		solver->literals[i].nClauses = 0;
		solver->literals[i].size = 1;
		solver->literals[i].clauses = calloc(1, sizeof(int));

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
	if(lit->size == lit->nClauses){
		lit->size *= 2;
		lit->clauses = realloc(lit->clauses, sizeof(int) * lit->size);
	}
	lit->clauses[lit->nClauses] = clause;
	lit->nClauses++;
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
		free(solver->literals[i].clauses);
	}
	free(solver->literals); solver->literals = NULL;
	free(solver->clauses); solver->clauses = NULL;
	free(solver->memory); solver->memory = NULL;
	solver->nClauses = solver->nVars = 0;
	free(solver->flipped); solver->flipped = NULL;
	free(solver->breakCount); solver->breakCount = NULL;
	free(solver->posLiterals); solver->posLiterals = NULL;
	free(solver->cumProbs); solver->cumProbs = NULL;
	free(solver->optFlip); solver->optFlip = NULL;
	free(solver->createCount); solver->createCount = NULL;
}

/*
 * initialize breakCount and number of positive Literals
 * */
void setupSolver(Solver* solver){
	int i,j, size, lit;
	int* clause;
	solver->optNHorn = 0;
	for(i = 0; i < solver->nClauses; i++){
		solver->posLiterals[i] = 0;
		clause = solver->clauses[i];
		size = solver->clauses[i+1] - solver->clauses[i];
		for(j = 0; j < size; j++){
			solver->posLiterals[i] += clause[j] > 0 ? 1 : 0;
		}	
		solver->optNHorn += solver->posLiterals[i] < 2 ? 1 : 0;
	}	
	for(i = 0; i < solver->nClauses; i++){
		clause = solver->clauses[i];
		size = solver->clauses[i+1] - solver->clauses[i];
		for(j = 0; j < size; j++){
			lit = clause[j];
			solver->breakCount[abs(lit)-1] += solver->posLiterals[i] == 1 && lit < 0;
			solver->createCount[abs(lit)-1] += solver->posLiterals[i] == 2 && lit  > 0;

		}
	}
}

void print(Solver* solver){
	int nVar, nClause, size;
	int* clause;
	
	printf("p cnf %i %i\n", solver->nVars, solver->nClauses);
	for(nClause = 0; nClause < solver->nClauses; nClause++){
		clause = solver->clauses[nClause];
		size = solver->clauses[nClause+1] - solver->clauses[nClause];
		for(nVar = 0; nVar < size; nVar++){
			printf("%i ", clause[nVar]);
		}
		printf("\n");
	}

}

void solve_incrementally(Solver* solver, int maxFlips){
	int i, lit, clause;
	//print_status(solver);
	for(i = 0; i < maxFlips; i++){
		DEBUG_PRINT(("Trie: %i\n", i));
		clause = getNonHornClause(solver);
		if(clause == -1){
		   DEBUG_PRINT(("Done\n"));
	   	   break;
		}
		DEBUG_PRINT(("Chosen clause: %d\n", clause));
		lit = getLitInClause(solver, clause);
		DEBUG_PRINT(("lit: %3d\n", lit));
		printf("%d ", lit);
		updateSolver(solver, lit);
		flipLiteral(solver, lit);
		updateOptimum(solver);
		//printf("%d ", solver->nClauses - countHornClauses(solver));
		//print_status(solver);
	}
	DEBUG_PRINT(("Optimal number of horn clauses found: %i\n", solver->optNHorn));
	DEBUG_PRINT_iARRAY(solver->optFlip, solver->nVars);
}

int allClausesAreHorn(Solver* solver){
	int i;
	for(i = 0; i < solver->nClauses; i++){
		if(solver->posLiterals[i] > 1) return 0;
	}
	return 1;
}

/*
 * get a random literal with minimum break count
 * */
int getLit(Solver* solver){
	int min, i = 0, litIndex, nLiterals = 1;
	min = getMinBreakCount(solver);
	while( i < solver->nVars && solver->breakCount[i] != min) i++; // get first literal with minimum break count
	litIndex = i;
	for(; i < solver->nVars; i++){
		if(solver->breakCount[i] == min && rand() < RAND_MAX/(nLiterals+1)){ // randomly keep current literal or replace with new one
		   	litIndex = i;
			nLiterals++;
		}
	}
	return litIndex+1;
}

/*
 * choose a literal using the breakCount and the method used in probSAT
 * */
int getLitProb(Solver* solver){
	int i, prob, randInt;
	prob = solver->nClauses + 1 - solver->breakCount[0];
	solver->cumProbs[0] = prob;
	for(i = 1; i < solver->nVars; i++){
		prob = solver->nClauses + 1 - solver->breakCount[i];
		solver->cumProbs[i] = solver->cumProbs[i-1] + prob;
	}
	randInt = rand() % solver->cumProbs[solver->nVars-1];
	i = 0;
	while(solver->cumProbs[i] <= randInt) i++;
	DEBUG_PRINT(("Rand Int: %i\n", randInt));
	DEBUG_PRINT(("cumulative probabilities:\n"));
	DEBUG_PRINT_iARRAY(solver->cumProbs, solver->nVars);
	return i+1; // "i" is index of the literal "i+1"
	
}

void updateSolver(Solver* solver, int lit){
	int i;
	for(i = 0; i < solver->literals[abs(lit)-1].nClauses; i++){
		updateClause(solver, solver->literals[abs(lit)-1].clauses[i], lit);
	}
}

/*
 * update number of positive literals and breakcount for all variables in clause
 * if necessary
 * */
void updateClause(Solver* solver, int nClause, int lit){
	updateBreakCreateCount(solver, nClause, lit);
	if(isPosInClause(solver, nClause, lit)) 
		solver->posLiterals[nClause]--;
	else 
		solver->posLiterals[nClause]++;
}

void updateBreakCreateCount(Solver* solver, int nClause, int flipLit){
	int i, nPosLiterals, size, lit;
	int* clause;
	nPosLiterals = solver->posLiterals[nClause];
	clause = solver->clauses[nClause];
	size = solver->clauses[nClause+1] - solver->clauses[nClause];
	if(isPosInClause(solver, nClause, flipLit)){ // flipLit is changed from positive to negative
		if(nPosLiterals != 1 && nPosLiterals != 2 && nPosLiterals != 3)
			return;
		for(i = 0; i < size; i++){ 
			lit = clause[i];
			if(abs(lit) != abs(flipLit)){
				if(!isPosLit(solver, lit)){
					if(nPosLiterals == 1) // clause cannot be broken by one flip anymore because there will be 0 positive literals in the clause
						solver->breakCount[abs(lit)-1] -= 1;
					else if(nPosLiterals == 2)// clause becomes horn because number of positive literals decreases from 2 to 1
						solver->breakCount[abs(lit)-1] += 1;
				}else{
					if(nPosLiterals == 2) // clause becomes horn -> literal i cannot make clause horn
						solver->createCount[abs(lit)-1] -= 1;
					else if(nPosLiterals == 3) // clause stays nonhorn but can become horn with one flip of literal 1
						solver->createCount[abs(lit)-1] += 1;
				}
			}
			else{
				if(nPosLiterals == 2){ // number of literals decreases from 2 to 1 and lit will be negated in the clause -> clause will be able to be broken
					solver->breakCount[abs(lit)-1] += 1;
					solver->createCount[abs(lit)-1] -= 1;
				}
			}
		}
	}else{ // lit is changed from negated to positive
		if(nPosLiterals != 0 && nPosLiterals != 1 && nPosLiterals != 2)
			return;
		for(i = 0; i < size; i++){
			lit = clause[i];
			if(abs(lit) != abs(flipLit)){
				if(!isPosLit(solver, lit)){
					if(nPosLiterals == 0) // number of positive literals increases from 0 to 1 -> clause can be broken
						solver->breakCount[abs(lit)-1] += 1;
					else if(nPosLiterals == 1) // number of positive literals increases from 1 to 2 -> clause cannot be broken anymore
						solver->breakCount[abs(lit)-1] -= 1;
				}else{
					if(nPosLiterals == 1) // there will be 2 pos. literals in clause -> by flipping literal i this clause will become horn again
						solver->createCount[abs(lit)-1] += 1;
					else if(nPosLiterals == 2) // there will be 3 pos. literals in clause -> flipping literal i cannot make this clause horn anymore
						solver->createCount[abs(lit)-1] -= 1;
				}
			}else{ // the flip broke the clause -> decrease breakCount
				if(nPosLiterals == 1){
					solver->breakCount[abs(lit)-1] -= 1;
					solver->createCount[abs(lit)-1] += 1;
				}
			}
		}
	}
}

void updateOptimum(Solver* solver){
	int cntHorn = countHornClauses(solver), i;
	if(cntHorn > solver->optNHorn){
		solver->optNHorn = cntHorn;
		for(i = 0; i < solver->nVars; i++){
			solver->optFlip[i] = solver->flipped[i];
		}
	}
}

/*
 * returns minimal break value of all variables
 * */
int getMinBreakCount(Solver* solver){
	int min = solver->nClauses, i;
	for(i = 0; i < solver->nVars; i++){
		if(solver->breakCount[i] < min) {
			min = solver->breakCount[i];
		}
	}
	return min;
}

int getMaxCreateValue(Solver* solver){
	int max = 0, i;
	for(i = 0; i < solver->nVars; i++){
		if(solver->createCount[i] > max) max = solver->createCount[i];
	}
	return max;
}


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
 * Returns a literal in the given clause that should be flipped
 * candidate literals are those with the lowest break count
 * if multiple of those exist a random one is chosen
 * */
int getLitInClause(Solver* solver, int nClause){
	int* clause = solver->clauses[nClause];
	int size = solver->clauses[nClause+1] - solver->clauses[nClause], i, min = solver->nClauses +1;
	int lit, cnt, chosenLit;
	for(i = 0; i < size; i++){
		lit = clause[i];
		//if(!isPosLit(solver, lit)) continue;
		if(solver->breakCount[abs(lit)-1] < min){
			cnt = 1;
			min = solver->breakCount[abs(lit)-1];
			chosenLit = lit;
		}else if(solver->breakCount[abs(lit)-1] == min){ // literal with minimum break count found
			if((float) rand() < (float) RAND_MAX/(cnt+1)){ // choose literal according to random uniform distribution
				chosenLit = lit;
			}
			cnt++;	
		}
	}
	return chosenLit;
}

/*
int isPosLit(Solver* solver, int* clause, int i){
	return clause[i] * solver->flipped[abs(clause[i])-1] > 0;

}
*/

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

/*
 * choose a random literal among those with the minimum break value
 * */
int chooseLiteral(int* clause, int* breakValues, int size, int min){
	int i = 0, litIndex;
	while(i < size && breakValues[i] != min) i++; // get first literal candidate
	litIndex = i;
	for(; i < size; i++){ // choose one random literal candidate
		if(breakValues[i] == min)
			if(rand() <  RAND_MAX / 2) litIndex = i;
	}
	return clause[litIndex]; 
}

void print_status(Solver* solver){
	int i,j, size;
	int *clause;
	printf("Horn clauses: %i\n", countHornClauses(solver));
	for(i = 0; i < solver->nClauses; i++){
		clause = solver->clauses[i];
		size = solver->clauses[i+1] - solver->clauses[i];
		for(j = 0; j < size; j++){
			printf("%3d ", clause[j] * solver->flipped[abs(clause[j])-1]);
		}
		printf("|%2d\n", solver->posLiterals[i]);
	}
	for(i = 0; i < solver->nVars; i++)
		printf("--");
	printf("\n");
	for(i = 0; i < solver->nVars; i++)
		printf("%3d ", i+1);
	printf("\n");
	for(i = 0; i < solver->nVars; i++)
		printf("%3d ", solver->flipped[i]);
	printf("| flip Status\n");
	for(i = 0; i < solver->nVars; i++)
		printf("%3d ", solver->breakCount[i]);
	printf("| break count\n");
	for(i = 0; i < solver->nVars; i++)
		printf("%3d ", solver->createCount[i]);
	printf("| create count\n");

}

int countHornClauses(Solver *solver){
	int i, nHorn = 0;
	for(i = 0; i < solver->nClauses; i++){
		nHorn += isHorn(solver, i);
	}
	return nHorn;
}


void run_experiments(Solver* solver, int maxFlips){
	int i,n;
	int (*heuristics[])(Solver* solver) = {getLitByBreak, getLitByBreakProb, getLitByCreate, getLitByCreateProb};
	char* heuristic_names[] = {"break_count", "break_count_prob", "create_count", "create_count_prob"};
	n = sizeof(heuristics) / sizeof(heuristics[0]);
	for(i = 0; i < n; i++){
		printf("%s ", heuristic_names[i]);
		solve(solver, heuristics[i], maxFlips);
		printf("\n");
		reset(solver);
	}
}


void solve(Solver* solver, int (* heuristic)(Solver* solver), int maxFlips){
	int i, lit;
	for(i = 0; i < maxFlips; i++){
		if(allClausesAreHorn(solver)){
	   	   break;
		}
		lit = (*heuristic)(solver);
		updateSolver(solver, lit);
		flipLiteral(solver, lit);
		updateOptimum(solver);
		printf("%i ", countHornClauses(solver));
	}
}


void reset(Solver* solver){
	int i;
	for(i = 0; i < solver->nVars; i++){
		solver->flipped[i] = 1;
		solver->optFlip[i] = 1;
	}
	setupSolver(solver);
}

int getLitByBreak(Solver* solver){
	int min, i = 0, litIndex, nLiterals = 1;
	min = getMinBreakCount(solver);
	while( i < solver->nVars && solver->breakCount[i] != min) i++; // get first literal with minimum break count
	litIndex = i;
	for(; i < solver->nVars; i++){
		if(solver->breakCount[i] == min && rand() < RAND_MAX/(nLiterals+1)){ // randomly keep current literal or replace with new one
		   	litIndex = i;
			nLiterals++;
		}
	}
	return litIndex+1;
}

int getLitByBreakProb(Solver* solver){
	int i, prob, randInt;
	prob = solver->nClauses + 1 - solver->breakCount[0];
	solver->cumProbs[0] = prob;
	for(i = 1; i < solver->nVars; i++){
		prob = solver->nClauses + 1 - solver->breakCount[i];
		solver->cumProbs[i] = solver->cumProbs[i-1] + prob;
	}
	randInt = rand() % solver->cumProbs[solver->nVars-1];
	i = 0;
	while(solver->cumProbs[i] <= randInt) i++;
	return i+1; // "i" is index of the literal "i+1"
}

int getLitByCreate(Solver* solver){
	int max, i = 0, litIndex, nLiterals = 1;
	max = getMaxCreateValue(solver);
	while( i < solver->nVars && solver->createCount[i] != max) i++; // get first literal with minimum break count
	litIndex = i;
	for(; i < solver->nVars; i++){
		if(solver->createCount[i] == max && rand() < RAND_MAX/(nLiterals+1)){ // randomly keep current literal or replace with new one
		   	litIndex = i;
			nLiterals++;
		}
	}
	return litIndex+1;

}

int getLitByCreateProb(Solver* solver){
	int i, prob, randInt;
	prob = solver->createCount[0];
	solver->cumProbs[0] = prob;
	for(i = 1; i < solver->nVars; i++){
		prob = solver->createCount[i];
		solver->cumProbs[i] = solver->cumProbs[i-1] + prob;
	}
	randInt = rand() % solver->cumProbs[solver->nVars-1];
	i = 0;
	while(solver->cumProbs[i] <= randInt) i++;
	return i+1; // "i" is index of the literal "i+1"
}
