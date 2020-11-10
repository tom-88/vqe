/** @file 
 * A demo of QuEST
 *
 * @author Ania Brown
 * @author Tyson Jones
 */

#include <stdio.h>
#include "QuEST.h"

int main (int narg, char *varg[]) {

    

    /*
     * PREPARE QuEST environment
     * (Required only once per program)
     */

    QuESTEnv env = createQuESTEnv();

    printf("-------------------------------------------------------\n");
    printf("WOOOOOOOO its vqe:\n\t baby.\n");
    printf("-------------------------------------------------------\n");



    /*
     * PREPARE QUBIT SYSTEM
     */

    Qureg qubits = createQureg(1, env);
    initZeroState(qubits);



    /*
     * REPORT SYSTEM AND ENVIRONMENT
     */
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    /*
     * HAMILTONIAN
     */
    int z[2][2] = {{1,0},{0,-1}};

    /*
     * ANSATZ PARAMS 
     *
     * Single qubit variational form : U3(theta, thi, lambda) = [[cos(theta/2), -e^(i*lambda)sin(theta/2)],[e^(i*phi)sin(theta/2),e^(i*lambda + i*phi)cos(theta/2)]] 
     */
    qreal theta = M_PI / 2;
    qreal phi = 0;;
    qreal lambda = 0;

    qreal cosTheta = cos(theta / 2);
    qreal sinTheta = sin(theta / 2);
 
    Complex expIPhi;
    expIPhi.real = cos(phi);
    expIPhi.imag = sin(phi);

    Complex expILambda;
    expILambda.real = cos(lambda);
    expILambda.imag = sin(lambda);

    ComplexMatrix2 u3 = {{cosTheta, -1*expIPhi*sinTheta},{expIPhi*sinTheta, expIPhi*expILambda*cosTheta}} 

    /*
     * INITIAL STATE
     */

    /*
     * APPLY ANSATZ GATES
     */ 
    Vector v = {.x=0, .y=1, .z=0};
    rotateAroundAxis(qubits, 0, theta, v); 
    /*
     * MEASUREMENTS : MEASURE EACH SUB HAMILTONIAN 
     */

    /*
     * PAULI Z MEASUREMENTS : JUST MEASURE IN COMPUTATIONAL BASIS
     */
    int outcome = measure(qubits, 0);
    printf("Qubit 0 was measured in state %d\n", outcome);

    /*
     * PAULI X MEASUREMENTS : ROTATE BASIS ABOUT Y AXIS BY -PI/2, MEASURE IN COMPUTATIONAL BASIS
     */

    /*
     * PAULI Y MEASUREMENTS : ROTATE BASIS ABOUT Y AXIS BY PI/2, MEASURE IN COMPUTATIONAL BASIS
     */
     
    /*
     * MULTIPLY / SUM PAULI MEASUREMENT OUTCOMES TO GET TERM EXPECTATION VALUE
     */

    /*
     * SUM EXPECTATION VALUES OF SUB HAMILTONIAN TO GIVE ENERGY EXPECTATION
     */

     /*
      * PASS ENERGY EXPECTATION TO CLASSICAL OPTIMIZER
      */



   



    /*
     * FREE MEMORY
     */

    destroyQureg(qubits, env); 



    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    destroyQuESTEnv(env);
    return 0;
}
