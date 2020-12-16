
// TODO: assuming Hamiltonian is reduced 1 term pauli entries with consts, need a way of representing this so that measurements can be decided


#include <stdio.h>
#include <stdlib.h>
#include "QuEST.h"
#include <QuEST_complex.h>

int measurementToEvalue(int outcome);
ComplexMatrix2 singleQubitVariationalForm(qreal theta, qreal phi, qreal lambda);
int energyExpectation(qreal theta,qreal phi,qreal lambda);
int gradientDescent();
int symmetricDiffQuotient(int(*cost)(qreal, qreal, qreal),qreal x1,qreal x2, qreal x3, int target_param, double step);




/*
* PREPARE QuEST environment
* (Required only once per program)
*/

QuESTEnv env;

/*
 *  1 -> 1
 *  0 -> -1
 *
 */

int measurementToEvalue(outcome){
    return 2*outcome - 1;
}

/* Returns a pointer to a universal gate for a single qubit
 *
 * Single qubit variational form : U3(theta, thi, lambda) = [[cos(theta/2), -e^(i*lambda)sin(theta/2)],[e^(i*phi)sin(theta/2),e^(i*lambda + i*phi)cos(theta/2)]] 
 *
 */

ComplexMatrix2 singleQubitVariationalForm(qreal theta, qreal phi, qreal lambda){
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

    return u3;
}
 

int energyExpectation(qreal theta,qreal phi,qreal lambda){
   /*
     * ANSATZ PARAMS 
     *
     */
   
    /*
     * INITIAL STATE
     */

    /*
     * APPLY ANSATZ GATES
     */ 
    ComplexMatrix2 u3 = singleQubitVariationalForm(theta, phi, lambda);
 
     /*
     * PREPARE QUBIT SYSTEM
     */

    Qureg qubitsZ = createQureg(1, env);
    
    initZeroState(qubitsZ);
    unitary(qubitsZ,0,u3);
   
    Qureg qubitsX = createQureg(1, env);
    cloneQureg(qubitsX,qubitsZ); //copy circuit
    Qureg qubitsY = createQureg(1, env);
    cloneQureg(qubitsY, qubitsZ); //copy circuit
 
     /*
     * REPORT SYSTEM AND ENVIRONMENT
     */
    reportQuregParams(qubitsZ);
    reportQuregParams(qubitsX);
    reportQuregParams(qubitsY);
    reportQuESTEnv(env);

    /*
     * MEASUREMENTS : MEASURE EACH SUB HAMILTONIAN 
     */

    /*
     * PAULI Z MEASUREMENTS : JUST MEASURE IN COMPUTATIONAL BASIS
     */
    int outcomeZ = measure(qubitsZ, 0);
    printf("Qubit 0 was measured in state %d\n", outcomeZ);
    printf("Qubit 0 output measurement as pauli z eigenvalue: %d\n", measurementToEvalue(outcomeZ));

    /*
     * PAULI X MEASUREMENTS : ROTATE BASIS ABOUT Y AXIS BY -PI/2, MEASURE IN COMPUTATIONAL BASIS
     */
    rotateY(qubitsX,0,-1*M_PI/2);  
    int outcomeX = measure(qubitsX, 0);
    printf("Qubit 0 was measured in state %d\n", outcomeX);
    printf("Qubit 0 output measurement as pauli x eigenvalue: %d\n", measurementToEvalue(outcomeX));


    /*
     * PAULI Y MEASUREMENTS : ROTATE BASIS ABOUT Y AXIS BY PI/2, MEASURE IN COMPUTATIONAL BASIS
     */
    rotateY(qubitsY,0,M_PI/2);  
    int outcomeY = measure(qubitsY, 0);
    printf("Qubit 0 was measured in state %d\n", outcomeZ);
    printf("Qubit 0 output measurement as pauli y eigenvalue: %d\n", measurementToEvalue(outcomeZ));


    
    /*
     * MULTIPLY / SUM PAULI MEASUREMENT OUTCOMES TO GET TERM EXPECTATION VALUE
     */

    /*
     * SUM EXPECTATION VALUES OF SUB HAMILTONIAN TO GIVE ENERGY EXPECTATION
     */
    


   



    /*
     * FREE MEMORY
     */

    destroyQureg(qubitsY, env); 
    destroyQureg(qubitsX, env); 
    destroyQureg(qubitsZ, env); 

     /*
      * PASS ENERGY EXPECTATION TO CLASSICAL OPTIMIZER
      */
    
    return measurementToEvalue(outcomeX) + measurementToEvalue(outcomeY) + measurementToEvalue(outcomeZ);


}

/* Performs gradient descent with params for 1 qubit variational form
 * 
 * TODO: If possible, make this modular & generic for any function and a list of params
 */

int gradientDescent(){
    qreal theta = M_PI / 2;
    qreal phi = 0;
    qreal lambda = 0;

    double step_size = 0.1;
    double diff_step_size = 0.1;
    int prev_cost = energyExpectation(theta, phi, lambda);
    int current_cost = 0;

    while (abs(current_cost - prev_cost) > 0.1){
        theta = theta - step_size*symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 1, diff_step_size);
        phi = phi - step_size*symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 2, diff_step_size);
        lambda = lambda - step_size*symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 3, diff_step_size);
        prev_cost = current_cost;
        current_cost = energyExpectation(theta, phi, lambda);
    }
    
    
    
    return current_cost;
}

int symmetricDiffQuotient(int(*cost)(qreal, qreal, qreal),qreal x1,qreal x2, qreal x3, int target_param, double step){
    if (target_param == 1){
        return (cost(x1+step,x2,x3) - cost(x1-step,x2,x3)) / (2*step);
    }
    else if (target_param == 2){
        return (cost(x1,x2+step,x3) - cost(x1,x2-step,x3)) / (2*step);
    }
    else{
        return (cost(x1,x2,x3+step) - cost(x1,x2,x3-step)) / (2*step);
    }
}



int main (int narg, char *varg[]) {

    env = createQuESTEnv();



    printf("-------------------------------------------------------\n");
    printf("WOOOOOOOO its vqe:\n\t baby.\n");
    printf("-------------------------------------------------------\n");

    int groundStateEnergyBound = gradientDescent();
    

      
    printf("Ground state energy bound %d \n", groundStateEnergyBound);
    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    destroyQuESTEnv(env);
    return 0;
}
