#include <stdio.h>
#include <stdlib.h>

#ifdef DEBUG
#define DEBUG_PRINT(x) printf x
#else
#define DEBUG_PRINT(x) do {} while (0)
#endif

typedef struct solver {
	int nClauses, nVars;
	int *memory, *flipped;
	int **clauses; // array with pointers into memory at the start and end of the clauses
} Solver;

// Parsing
int parse(Solver* solver, char *fileName);
int skip_comment(FILE* file);
//Memory
void allocate_memory(Solver* solver);
void free_solver(Solver* solver);
void add_clause(Solver* solver, int* buffer, int size, int nClause);
void print(Solver* solver);
void print_status(Solver* solver);

// Search
void solve_globally(Solver* solver, int maxTries);
int findNonHornClause(Solver* solver);
int getLit(Solver* solver, int nClause);
int chooseLiteral(int* clause, int* breakValues, int size, int min);
int getBreakValue(Solver* solver, int* posLiterals, int* posLiterals2);
int* countPosLits(Solver* solver);
int countPosLit(Solver* solver, int i);
int isHorn(Solver* solver, int nClause);
int isNonHorn(Solver* solver, int nClause);
int countHornClauses(Solver *solver);
int isPosLit(Solver* solver, int* clause, int i);
void flipLiteral(Solver* solver, int nLit);

int main(int argc, char **argv){
	Solver solver;
	if(argc != 2){ printf("Usage: %s fileName\n", argv[0]); exit(0); }
	if(parse(&solver, argv[1]) == 1){
		print(&solver);
		printf("\n\n");
		solve_globally(&solver, 10);
	}else
		printf("Error while parsing file\n");
	free_solver(&solver);
}

int parse(Solver* solver, char *fileName){
	FILE* file;
	int c, i = 0, lit, size = 0;
	int* buffer;
	file = fopen(fileName, "r");
	if(!file){
		printf("Could no open file %s", fileName);
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
	while(i < solver->nClauses){
		c = fgetc(file);
		if( c == 'c') {
			if(skip_comment(file) != 1) {
				free(buffer); fclose(file); return -1;
			} 
			continue;}
		else if (c == ' ' || c == '\n') continue;
		else {
			ungetc(c, file);
			if(fscanf(file, "%i", &lit) != 1) {free(buffer); fclose(file); return -1;}
			if(lit == 0){
				add_clause(solver, buffer, size, i);
				i++;
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

void allocate_memory(Solver* solver){
	int i;
	solver->memory = (int *) calloc(solver->nVars * solver->nClauses, sizeof(int));
	solver->clauses = (int **) calloc(solver->nClauses + 1, sizeof(int *));
	solver->clauses[0] = solver->memory;
	// store nVars + 1 such that the absolute literal value can be used as index (index 0 will not be used)
	solver->flipped = (int *) calloc(solver->nVars + 1, sizeof(int));
	for(i = 0; i < solver->nVars + 1; i++)
		solver->flipped[i] = 1;
}

void add_clause(Solver* solver, int* buffer, int size, int nClause){
	int* clause = solver->clauses[nClause];
	int i;
	for(i = 0; i < size; i++){
		clause[i] = buffer[i];
	}	
	solver->clauses[nClause+1] = solver->clauses[nClause] + size;
}

int skip_comment(FILE* file){
	int ch;
	while((ch = getc(file)) != '\n');
	if(ch == EOF) return -1;
	else return 1;
}

void free_solver(Solver* solver){
	free(solver->clauses);
	solver->clauses = NULL;
	free(solver->memory);
	solver->memory = NULL;
	solver->nClauses = solver->nVars = 0;
	free(solver->flipped);
	solver->flipped = NULL;
}


void print(Solver* solver){
	int nClause, nVar, size;
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

void solve_globally(Solver* solver, int maxFlips){
	int i, nClause, nLit;
	print_status(solver);
	for(i = 0; i < maxFlips; i++){
		DEBUG_PRINT(("Trie %i\n", i));
		nClause = findNonHornClause(solver);
		DEBUG_PRINT(("Chosen Non Horn Clause: %i\n", nClause));
		if(nClause < 0) // no clause can be flipped -> process is done
			return;
		nLit = getLit(solver, nClause);
		flipLiteral(solver, nLit);
		print_status(solver);
	}
}

void flipLiteral(Solver* solver, int nLit){
	DEBUG_PRINT(("Flipping literal: %i\n", abs(nLit)));
	solver->flipped[abs(nLit)] *= -1;
}

/*
 * iterate over all clauses and return the index of a random
 * nonHorn clause
 * return -1 if no nonHorn clause exists
 * */
int findNonHornClause(Solver* solver){
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

// return 1 if clause  is isNonHorn, 0 otherwise
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
		if(solver->flipped[abs(literal)] * literal > 0) posLit++;
	}
	return posLit;
}

/*
 * Returns a literal in the given clause that should be flipped
 * candidate literals are thos with the lowest break count
 * if multiple of those exist a random one is chosen
 * */
int getLit(Solver* solver, int nClause){
	int* clause = solver->clauses[nClause], *posLiterals, *posLiterals2;
	int size = solver->clauses[nClause+1] - solver->clauses[nClause], i, min = nClause +1;
	int breakValues[size];
	int lit;
	posLiterals = countPosLits(solver);
	printf("Pos Literals: ");
	for(i=0; i < solver->nClauses; i++){
		printf("%i ", posLiterals[i]);
	} printf("\n");
	for(i = 0; i < size; i++){
		lit = clause[i];
		if(isPosLit(solver, clause, i)){ 
			solver->flipped[abs(lit)] *= -1;
			posLiterals2 = countPosLits(solver);	
			breakValues[i] = getBreakValue(solver, posLiterals, posLiterals2);
			if(breakValues[i] < min)
				min = breakValues[i];
			free(posLiterals2); posLiterals2 = NULL;
			solver->flipped[abs(lit)] *= -1;
		}else{ // negative literal cannot be flipped
			breakValues[i] = solver->nClauses + 1; // nClauses+1 clauses cannot be broken -> literal will not be considered afterwards
		}
	}
	printf("Break values: ");
	for(i = 0; i < size; i++){
		printf("%i ", breakValues[i]);
	}
	printf("\n");
    lit = chooseLiteral(clause, breakValues, size, min);		
	free(posLiterals);	
	return lit;
}

int isPosLit(Solver* solver, int* clause, int i){
	return clause[i] * solver->flipped[abs(clause[i])] > 0;
}

/*
 * count positive literals for all clauses
 * returns dynamically allocated array with number of positive
 * literals for each clause
 * */
int* countPosLits(Solver* solver){
	int* clauses = (int *) calloc(solver->nClauses, sizeof(int));
	int i;
	for(i = 0; i < solver->nClauses; i++){
		clauses[i] = countPosLit(solver, i);
	}
	return clauses;
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
	int i;
	printf("Horn clauses: %i\n", countHornClauses(solver));
	printf("Flip status: ");
	for(i = 1; i <= solver->nVars; i++){
		printf("%i ", solver->flipped[i]);
	}
	printf("\n\n");

}

int countHornClauses(Solver *solver){
	int i, nHorn = 0;
	for(i = 0; i < solver->nClauses; i++){
		nHorn += isHorn(solver, i);
	}
	return nHorn;
}
