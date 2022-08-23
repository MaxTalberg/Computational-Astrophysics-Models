#include <stdio.h> // pre-compiler statement
#include <stdlib.h> // free
#include <math.h> // pre-compiler statement, remember to LINK EXPLICITLY using "-lm"
#include <time.h> // clock

/// introduce parameters
#define n_cells 200              // number of cells // Questions 2: n_cells 1000 //
#define n_x (n_cells + 2)          // number of cells + ghost cells
#define x_initial 0          // length of tube
#define x_final 1

/// introduce variables
int i, j;
double t_final;     // final time changes with shocktube
double dx = (double) (x_final - x_initial) / n_cells;
double current_time, dt;    // timestep and current (total) time in loop
double X[n_x], S[n_x][3];   // array for x count from 0 to 1 and L/R wavespeed values
double F[n_x][3], Q[n_x][3], Qn[n_x][3]; // state cell vector arrays
double F_plus[n_x-1][3], F_minus[n_x-1][3]; // boundary vector arrays
                                // defining dynamic variables
double *r0, *p0, *v0;        // r0 = rho, p0 = pressure, v0 = velocity
double *E, *c_s, *a;        // E = energy density, c_s = speed of sound, a = max speed
double r0L, p0L, v0L, EL, c_sL; // defining variables for the HLL method
double r0R, p0R, v0R, ER, c_sR;
double sL, sR, p0_star;

double gama = 1.4;  // gamma value for shocktube A and B, changes for C
double CFL = 1;      // to reduce stepsize in case of supersonic shock
int A = 0;          // variables and constants to allow smooth operation
int B = 1;
int C = 2;
int LF = 3;
int HLL = 4;
int tube;


/// generate the maximum value of any array of any size int size
double maximum(double const array[], int size)
{
    double max;
    max = array[0];
    for(i=1 ;i < size - 1; i++)
    {
        if(array[i] > max)
            max = array[i];
    }
    return max;
}

/// function to allocate dynamic memory
void Allocate_Memory() {
    size_t memory_size = n_x * sizeof(double);
    r0 = malloc(memory_size);
    p0 = (double *) malloc(memory_size);
    v0 = (double *) malloc(memory_size);
    E = (double *) malloc(memory_size);
    c_s = (double *) malloc(memory_size);
    a = (double *) malloc(memory_size);
}

/// generate the initial conditions of the shock tube grid
/// dynamic input for shocktube A, B or C
void Initial_Conditions(int shocktube) {
    /// initialise cell arrays
    for(j=0; j < 3; j++) {
        for (i = 0; i < n_x ; i++)
        {
            Q[i][j] = 0;
            Qn[i][j] = 0;
            F[i][j] = 0;
            S[i][j] = 0;
            r0[i] = 0;
            p0[i] = 0;
            v0[i] = 0;
            E[i] = 0;
            c_s[i] = 0;
            a[i] = 0;
        }
        /// initialise boundary arrays
        for (i = 0; i < n_x-1; i++)
        {
            F_minus[i][j] = 0;
            F_plus[i][j] = 0;
        }
    }

    /// initial conditions of Shocktube A
    if(shocktube == A) {
        tube = A;
        t_final = 0.2;
        for (i = 0; i < n_x; i++) {
            X[i] = i * dx;
            if (X[i] < 0.3) {
                r0[i] = 1;
                p0[i] = 1;
                v0[i] = 0.75;
            } else {
                r0[i] = 0.125;
                p0[i] = 0.1;
                v0[i] = 0;
            }
        }
    }
    /// initial conditions of Shocktube B
    else if(shocktube == B) {
        tube = B;
        t_final = 0.012;
        for (i = 0; i < n_x; i++) {
            X[i] = i * dx;
            if (X[i] < 0.8) {
                r0[i] = 1;
                p0[i] = 1000;
                v0[i] = -19.59745;
            } else {
                r0[i] = 1;
                p0[i] = 0.01;
                v0[i] = -19.59745;
            }
        }
    }
    /// initial conditions of Shocktube C
    else if(shocktube == C) {
        tube = C;
        gama = 1.667;       // gamma value is 5/3
        t_final = 0.1;     // max t value calculated in Update_Timestep function
                        // t_final stops at 0.008799 using the if loop in Update_Timestep
                        // multiple values were used in the report
        for (i = 0; i < n_x; i++) {
            X[i] = (10 * i * dx);   // box size of 10
            if (X[i] < 5) {         // dimensionless values
                r0[i] = 10;         // calculated in report appendix
                p0[i] = 20631;
                v0[i] = 50;
            } else {
                r0[i] = 30;
                p0[i] = 618.93;
                v0[i] = -50;
            }
        }
    }
}

/// set boundary conditions/ ghost cell states
/// constant boundary conditions used for simplicity
void Boundary_Conditions() {
    /// Boundary conditions A
    if(tube == A) {
        /// boundary conditions x = 0
        r0[0] = 1;
        p0[0] = 1;
        v0[0] = 0.75;
        E[0] = p0[0] / ((gama - 1)) + r0[0] * v0[0] * v0[0] * 0.5;
        c_s[0] = sqrt((gama * p0[0]) / r0[0]);

        /// boundary conditions at x = 1
        r0[n_x - 1] = 0.125;
        p0[n_x - 1] = 0.1;
        v0[n_x - 1] = 0;
        E[n_x - 1] = p0[n_x - 1] / ((gama - 1)) + r0[n_x - 1] * v0[n_x - 1] * v0[n_x - 1] * 0.5;
        c_s[n_x - 1] = sqrt((gama * p0[n_x - 1]) / r0[n_x - 1]);
    }
    /// Boundary conditions B
    else if(tube == B) {
        /// boundary conditions x = 0
        r0[0] = 1;
        p0[0] = 1000;
        v0[0] = -19.59745;
        E[0] = p0[0] / ((gama - 1)) + r0[0] * v0[0] * v0[0] * 0.5;
        c_s[0] = sqrt((gama * p0[0]) / r0[0]);

        /// boundary conditions at x = 1
        r0[n_x - 1] = 1;
        p0[n_x - 1] = 0.01;
        v0[n_x - 1] = -19.59745;
        E[n_x - 1] = p0[n_x - 1] / ((gama - 1)) + r0[n_x - 1] * v0[n_x - 1] * v0[n_x - 1] * 0.5;
        c_s[n_x - 1] = sqrt((gama * p0[n_x - 1]) / r0[n_x - 1]);
    }
        /// Boundary conditions C
    else if(tube == C) {
        /// boundary conditions x = 0
        r0[0] = 10;
        p0[0] = 20631;
        v0[0] = 50;
        E[0] = p0[0] / ((gama - 1)) + r0[0] * v0[0] * v0[0] * 0.5;
        c_s[0] = sqrt((gama * p0[0]) / r0[0]);

        /// boundary conditions at x = 10
        r0[n_x - 1] = 30;
        p0[n_x - 1] = 618.93;
        v0[n_x - 1] = -50;
        E[n_x - 1] = p0[n_x - 1] / ((gama - 1)) + r0[n_x - 1] * v0[n_x - 1] * v0[n_x - 1] * 0.5;
        c_s[n_x - 1] = sqrt((gama * p0[n_x - 1]) / r0[n_x - 1]);
    }
}

/// determine the stepsize dt based on speed of sounf and wavespeed
void Update_Timestep(){

    /// speed of sound and max possible wave speed
    for(i=1; i < n_x - 1; i++) {
        c_s[i] = sqrt((gama * p0[i]) / r0[i]);
        a[i] = fabs(v0[i]) + c_s[i];
    }

    /// determine step size using the systems largest eigenvalue
    dt = (double) CFL * dx / maximum(a, n_x);

    /// condition for C to end script when impact parameters
    /// are felt at the end of the zone
    if (tube == C){
        if(round(r0[1]) != r0[0])
        {
            t_final = 0;
        }
    }
}

/// update cell values of state vectors
void Update_Matrix(){
    /// update values for Q, F, E and c_s
    for(i=0; i < n_x; i++)
    {
        /// specific internal energy
        E[i] = p0[i] / ((gama-1) * r0[i]) + v0[i] * v0[i] * 0.5;

        /// update state cell vector Q values
        Q[i][0] = r0[i];
        Q[i][1] = r0[i] * v0[i];
        Q[i][2] = E[i] * r0[i];

        /// update state cell vector F values
        F[i][0] = r0[i] * v0[i];
        F[i][1] = r0[i] * v0[i] * v0[i] + p0[i];
        F[i][2] = (E[i] * r0[i] + p0[i]) * v0[i];
    }
}

/// determine flux between cells using different methods
/// Lax-Friedrich or the HLL method
void Determine_Flux(int method) {
    /// implementation of the HLL method
    if(method == HLL){
        /// calculate initial cell properties Left and Right
        for (i = 0; i < n_x; i++) {
            /// initial conditions
            /// left i-1
            r0L = Q[i - 1][0];
            v0L = Q[i - 1][1] / r0L;
            EL = Q[i - 1][2];
            p0L = (gama - 1) * (EL - (r0L * v0L * v0L * 0.5));
            c_sL = sqrt((gama * p0L) / r0L);

            /// right i
            r0R = Q[i][0];
            v0R = Q[i][1] / r0R;
            ER = Q[i][2];
            p0R = (gama - 1) * (ER - (r0R * v0R * v0R * 0.5));
            c_sR = sqrt((gama * p0R) / r0R);

            /// estimate pressure between waves
            p0_star = fmax(0, 0.5 * (p0L + p0R) - 0.5 * (v0R - v0L) * 0.5 * (r0L + r0R) * 0.5 * (c_sL + c_sR));

            /// calculate the wave speed values L and R
            /// left
            if (p0_star <= p0L) {
                sL = v0L - c_sL;
            } else {
                sL = v0L - c_sL * sqrt(1 + ((gama + 1) / (2 * gama) * ((p0_star / p0L) - 1)));
            }

            /// right
            if (p0_star <= p0R) {
                sR = v0R + c_sR;
            } else {
                sR = v0R + c_sR * sqrt(1 + ((gama + 1) / (2 * gama) * ((p0_star / p0R) - 1)));
            }

            /// store values
            S[i][0] = sL;
            S[i][1] = sR;

        }

        /// 1/2 flux values
        /// calculating flux boundary values
        for (j = 0; j < 3; j++) {
            for (i = 1; i < n_x; i++) {
                sL = S[i][0];
                sR = S[i][1];
                if ((sL >= 0) && (sR >= 0)) {
                    F_minus[i][j] = F[i - 1][j];
                } else if ((sL <= 0) && (sR >= 0)) {
                    F_minus[i][j] = (sR * F[i - 1][j] - sL * F[i][j] +
                                     sR * sL * (Q[i][j] - Q[i - 1][j])) / (sR - sL);
                } else if ((sR <= 0) && (sL <= 0)) {
                    F_minus[i][j] = F[i][j];
                }
            }
        }
        for (j = 0; j < 3; j++) {
            for (i = 1; i < n_x - 1; i++) {
                F_plus[i][j] = F_minus[i + 1][j];
            }
        }
    }
    /// implementation of the Lax-Friedrichs algorithm
    else {
        for(j=0; j < 3; j++) {
            for (i = 1; i < n_x - 1 ; i++)
            {
            /// F + 1/2
            F_plus[i][j] = ((F[i][j] + F[i+1][j])*0.5) + ((dx/dt)*(Q[i][j] - Q[i+1][j])*0.5);
            /// F - 1/2
            F_minus[i][j] = ((F[i-1][j] + F[i][j])*0.5) - ((dx/dt)*(Q[i][j] - Q[i-1][j])*0.5);
            }
        }
    }
}

/// function to update the flux values
/// after one of either the LF or HLL method has been applied
/// to find flux values between cells
void Update_Flux() {
    for (j = 0; j < 3; j++) {
        for (i = 1; i < n_x-1 ; i++) {
            /// Q n+1 the new Q value
            Qn[i][j] = Q[i][j] - ((dt/dx) * (F_plus[i][j] - F_minus[i][j]));
        }
    }
}

/// update initial variables with the updated state vector
void Update_Values(){
    /// assign output to a new array
    for(j=0; j < 3; j++) {
        for (i = 0; i < n_x ; i++) {
            Q[i][j] = Qn[i][j];
        }
    }
    /// update variables
    for(i=1; i< n_x-1; i++)
    {
        r0[i] = Qn[i][0];
        v0[i] = Qn[i][1] / r0[i];
        E[i] = Qn[i][2];
        p0[i] = (gama - 1) * ((E[i]) - r0[i] * v0[i] * v0[i] *0.5);
        c_s[i] = sqrt((gama * p0[i])/r0[i]);
    }
}

/// algorithm undergos a while loop
///  until the final time has been acheived by the system
void Algorithm(method){
    while(current_time < t_final)
    {
        Boundary_Conditions();

        Update_Timestep();

        Update_Matrix();

        Determine_Flux(method);

        Update_Flux();

        Update_Values();

        current_time += dt;     // update the current time
    }
}

/// save final arrays to then be printed in python
void Save_Values(){
    /// store output as a txt file to plot in python
    FILE *fp=NULL;
    fp=fopen("24511PH30110.text", "w");

    for(i=1; i < n_x-1; i++)
    {
        fprintf(fp, "%1.6f\t%1.6f\t%1.6f\t%1.6f\t%1.6f\n", r0[i], p0[i], v0[i], p0[i]/((gama-1)*r0[i]), X[i]);
    }
}

/// free allocated memory
void Free_Memory(){
    free(r0);
    free(p0);
    free(v0);
    free(E);
    free(c_s);
    free(a);
}

/// run programme
int main() {

    Allocate_Memory();

    Initial_Conditions(A);  // assign shocktube - A, B or C

    Algorithm(HLL);     // assign method - LF or HLL

    Save_Values();

    Free_Memory();

    return 0;
}