/** @file 
 * QAOA
 *
 * @author Tom O'Leary
 */

# include <stdio.h>
# include <math.h>

# include "QuEST.h" 


// C should be a collection of Pauli strings e.g. z1z2z3, x2z4, ...
// Think I will store my own in a file in this way then convert to the quEST format

/* This will use some optimisation algorithm to
 * find the qaoa angles which
 * minimise (or maximise) the cost function
 */
int analyticOptimalAngles(int n, int p, C, init);

typedef struct {
    double costOperatorAngle;
    double mixerOperatorAngle;
} anglePair;

int analyticOptimalAngles(int n, int p, C, init){
    //   
}



int main (int narg, char** varg) {


    /* 	
     * PREPARE QuEST
     */

    // model parameters
    int numQubits = 9;
    
    // prepare QuEST
    QuESTEnv env = createQuESTEnv();

    // create qureg; let zeroth qubit be ancilla
    Qureg qureg = createQureg(numQubits, env);
    initZeroState(qureg);


   /*
     * FREE MEMORY
     */

    destroyQureg(qureg, env); 
    destroyQuESTEnv(env);
    return 0;
}
