
// TODO: assuming Hamiltonian is reduced 1 term pauli entries with consts, need a way of representing this so that measurements can be decided

#include <stdio.h>
#include <stdlib.h>
#include "QuEST.h"
#include <QuEST_complex.h>

int measurementToEvalue(int outcome);
ComplexMatrix2 singleQubitVariationalForm(qreal theta, qreal phi, qreal lambda);
double energyExpectation(qreal theta,qreal phi,qreal lambda);
double gradientDescent();
double symmetricDiffQuotient(double(*cost)(qreal, qreal, qreal),qreal x1,qreal x2, qreal x3, int target_param, qreal step);
void testSymmetricDiffQuotient();
double addUs(qreal x, qreal y, qreal z);
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
 

double energyExpectation(qreal theta,qreal phi,qreal lambda){
    int measurementRepeats = 1000;
  
    /*
     * INITIAL STATE
     */

    /*
     * APPLY ANSATZ GATES
     */ 
    ComplexMatrix2 u3 = singleQubitVariationalForm(theta, phi, lambda);
 
    int outcomeZTotal = 0;
    int outcomeXTotal = 0;
    int outcomeYTotal = 0;
    
    for (int i = 0; i < measurementRepeats; i++){
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
         * PAULI Z MEASUREMENTS : JUST MEASURE IN COMPUTATIONAL BASIS
         */
        int outcomeZ = measure(qubitsZ, 0);
 //       printf("Qubit 0 was measured in state %d\n", outcomeZ);
 //       printf("Qubit 0 output measurement as pauli z eigenvalue: %d\n", measurementToEvalue(outcomeZ));
        outcomeZTotal = outcomeZTotal + measurementToEvalue(outcomeZ);
        /*
         * PAULI X MEASUREMENTS : ROTATE BASIS ABOUT Y AXIS BY -PI/2, MEASURE IN COMPUTATIONAL BASIS
         */
        rotateY(qubitsX,0,-1*M_PI/2);  
        int outcomeX = measure(qubitsX, 0);
 //       printf("Qubit 0 was measured in state %d\n", outcomeX);
 //       printf("Qubit 0 output measurement as pauli x eigenvalue: %d\n", measurementToEvalue(outcomeX));
        outcomeXTotal = outcomeXTotal + measurementToEvalue(outcomeX);

        /*
         * PAULI Y MEASUREMENTS : ROTATE BASIS ABOUT X AXIS BY PI/2, MEASURE IN COMPUTATIONAL BASIS
         */
        rotateX(qubitsY,0,M_PI/2);  
        int outcomeY = measure(qubitsY, 0);
 //       printf("Qubit 0 was measured in state %d\n", outcomeZ);
 //       printf("Qubit 0 output measurement as pauli y eigenvalue: %d\n", measurementToEvalue(outcomeZ));
        outcomeYTotal = outcomeYTotal + measurementToEvalue(outcomeY);

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
    }
         /*
      * PASS ENERGY EXPECTATION TO CLASSICAL OPTIMIZER
      */
    
    return (outcomeZTotal + outcomeXTotal + outcomeYTotal) / (double)measurementRepeats;
}

/* Performs gradient descent with params for 1 qubit variational form
 * 
 * Exit condition : after the cost function stops changing by some amount
 *
 * TODO: If possible, make this modular & generic for any function and a list of params
 *
 * could I derive the analytic derivative for the cost function of smaller Hamiltonians?
 */

double gradientDescentExitCondCostTolerance(){
    qreal theta = 0;
    qreal phi = 0;
    qreal lambda = 0;

    double step_size = 0.001;
    double diff_step_size = 0.1; // how to set this, tends to derivative when this tends to 0, so as small as possible?
    double tolerance = 0.0001;
    double prev_cost = 1000;
    double current_cost = 0;
    
    while(fabs(prev_cost-current_cost) > tolerance){
        printf("cost : %f \n",current_cost);
        double deriv = symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 0, diff_step_size);
        printf("deriv : %f \n",deriv);
        theta = theta - step_size*symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 0, diff_step_size);
        phi = phi - step_size*symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 1, diff_step_size);
        lambda = lambda - step_size*symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 2, diff_step_size);
        prev_cost = current_cost;
        current_cost = energyExpectation(theta, phi, lambda);
    }   
    
    return current_cost;
}

/* Performs gradient descent with params for 1 qubit variational form
 * 
 * Exit condition : after the gradient cost function stops changing by some amount
 *
 */

double gradientDescentExitCondGradientTolerance(){
    qreal theta = 0;
    qreal phi = 0;
    qreal lambda = 0;

    double step_size = 0.001;
    double diff_step_size = 0.1; 
    double tolerance = 0.0001;
    double prev_cost_gradient = 1000;
    double current_cost_gradient = 0;
    
    while(fabs(prev_cost_gradient-current_cost_gradient) > tolerance){
        double deriv_param_0 = symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 0, diff_step_size);
        double deriv_param_1 = symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 1, diff_step_size);
        double deriv_param_2 = symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 2, diff_step_size);
        
        theta = theta - step_size*deriv_param_0;
        phi = phi - step_size*deriv_param_1;
        lambda = lambda - step_size*deriv_param_2;
        
        prev_cost_gradient = current_cost_gradient;
        current_cost_gradient = deriv_param_0 + deriv_param_1 + deriv_param_2;   
    }

    return energyExpectation(theta, phi, lambda);
}

/* Performs gradient descent with params for 1 qubit variational form
 * exit condition : after some fixed number of iterations
 */
double gradientDescentExitCondIterations(){
    qreal theta = 0;
    qreal phi = 0;
    qreal lambda = 0;

    double step_size = 0.001;
    double diff_step_size = 0.1; // how to set this, tends to derivative when this tends to 0, so as small as possible?
    int max_it = 10000;
    int i = 0;
    double prev_cost = 1000;
    double current_cost = 0;
    
    while(i < max_it){
        printf("%i, ",i);
        double deriv = symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 0, diff_step_size);
        theta = theta - step_size*symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 0, diff_step_size);
        phi = phi - step_size*symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 1, diff_step_size);
        lambda = lambda - step_size*symmetricDiffQuotient(energyExpectation, theta, phi, lambda, 2, diff_step_size);
        prev_cost = current_cost;
        current_cost = energyExpectation(theta, phi, lambda);
        i++;
    }   
    
    return current_cost;
}

double symmetricDiffQuotient(double(*cost)(qreal, qreal, qreal),qreal x1,qreal x2, qreal x3, int target_param, qreal step){
    if (target_param == 0){
        return (cost(x1+step,x2,x3) - cost(x1-step,x2,x3)) / (2*step);
    }
    else if (target_param == 1){
        return (cost(x1,x2+step,x3) - cost(x1,x2-step,x3)) / (2*step);
    }
    else{
        return (cost(x1,x2,x3+step) - cost(x1,x2,x3-step)) / (2*step);
    }
}

void testSymmetricDiffQuotient(){
    qreal x = 2;
    qreal y = 3;
    qreal z = 4 ;

    double analyticResultDerivFByX = 4;
    double numericalResultDerivFByX = symmetricDiffQuotient(addUs,x,y,z,0,10.0);

    printf("Analytic result = %f \n",analyticResultDerivFByX);
    printf("Numerical result = %f \n",numericalResultDerivFByX);
}

// deriv of this is 1 + the other two variables
double addUs(qreal x, qreal y, qreal z){
    return (x*x)+y+z;
}


int main (int narg, char *varg[]) {
    testSymmetricDiffQuotient();
    
    env = createQuESTEnv();

    printf("-------------------------------------------------------\n");
    printf("WOOOOOOOO its vqe:\n\t baby.\n");
    printf("-------------------------------------------------------\n");

    double groundStateEnergyBound = gradientDescentExitCondIterations();
    
    printf("Expected ground state energy : -1.73205...");
    printf("Ground state energy bound %f \n", groundStateEnergyBound);
    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    destroyQuESTEnv(env);
    return 0;
}
