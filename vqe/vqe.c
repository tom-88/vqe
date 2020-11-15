/** @file 
 * A demo of QuEST
 *
 * @author Ania Brown
 * @author Tyson Jones
 */

#include <stdio.h>
#include "QuEST.h"

/*
 *  1 -> 1
 *  0 -> -1
 *
 */

int measurementToEvalue(outcome){
    return 2*outcome - 1;
}

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
    
    Complex u300;
    u300.real = cosTheta;
    u300.imag = 0;
    Complex u301;
    u301.real = -1*expILambda.real*sinTheta; 
    u301.imag = -1*expILambda.imag*sinTheta;
    Complex u310;
    u310.real = expIPhi.real*sinTheta;
    u310.imag = expIPhi.imag*sinTheta;
    Complex u311;
    u311.real = (expILambda.real*expIPhi.real - expIPhi.imag*expILambda.imag)*cosTheta;
    u311.imag = (expILambda.imag*expIPhi.real + expIPhi.imag*expILambda.real)*cosTheta;


    ComplexMatrix2 u3 = {
                           .real={{u300.real,u301.real},{u310.real,u311.real}},
                           .imag={{u300.imag,u301.imag},{u310.imag,u311.imag}}
                        };
    
    /*
     * INITIAL STATE
     */

    /*
     * APPLY ANSATZ GATES
     */ 
    unitary(qubits,0,u3);
    Qureg qubitsX = createQureg(qubits, env); //copy circuit
    Qureg qubitsY = createQureg(qubits, env); //copy circuit
    /*
     * MEASUREMENTS : MEASURE EACH SUB HAMILTONIAN 
     */

    /*
     * PAULI Z MEASUREMENTS : JUST MEASURE IN COMPUTATIONAL BASIS
     */
    int outcome = measure(qubits, 0);
    printf("Qubit 0 was measured in state %d\n", outcome);
    printf("Qubit 0 output measurement as pauli z eigenvalue: %d\n", measurementToEvalue(outcome));

    /*
     * PAULI X MEASUREMENTS : ROTATE BASIS ABOUT Y AXIS BY -PI/2, MEASURE IN COMPUTATIONAL BASIS
     */
    rotateY(qubitsX,0,-1*M_PI/2);  
    int outcome = measure(qubitsX, 0);
    printf("Qubit 0 was measured in state %d\n", outcome);
    printf("Qubit 0 output measurement as pauli x eigenvalue: %d\n", measurementToEvalue(outcome));


    /*
     * PAULI Y MEASUREMENTS : ROTATE BASIS ABOUT Y AXIS BY PI/2, MEASURE IN COMPUTATIONAL BASIS
     */
    rotateY(qubitsY,0,M_PI/2);  
    int outcome = measure(qubitsY, 0);
    printf("Qubit 0 was measured in state %d\n", outcome);
    printf("Qubit 0 output measurement as pauli y eigenvalue: %d\n", measurementToEvalue(outcome));


    
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
    destroyQureg(qubitsX, env); 
    destroyQureg(qubitsZ, env); 



    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    destroyQuESTEnv(env);
    return 0;
}
