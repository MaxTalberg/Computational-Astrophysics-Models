#include <stdio.h> // pre-compiler statement
#include <stdlib.h> // free
#include <math.h> // pre-compiler statement, remember to LINK EXPLICITLY using "-lm"
#include <time.h> // clock



/**
 * Random number generator
 * @return
 * a random number in the range 0 to 1
 */
double random_number()
{
    static int64_t seed = 1;
    static int32_t m = 2147483647;
    static int16_t a = 16807;
    static int16_t c;

    seed = (a * seed + c) % m;
    return (double) seed / m;
}

/**
 * Exponential distribution function
 * @param x
 * @return
 * y values obeying this function
 */
double f(double x)
{
    return exp(-x);
}


/**\
 * Cumulative distribtuion function
 * @param y
 * @return
 * x values obeying this cumulaitve function
 */
double f_cum(double y)
{
    return -log(1-y);
}

/**
 * polar vector used to store polar coordinates
 */
typedef struct polar_vector
{
    double r;
    double theta;
    double phi;
    double z;
} polar_vector;

/**
 * data type polar vector
 */
struct polar_vector r1;

/**
 * initialised vector for polar coordinates
 */
double outputPolarVector[4] = {0, 0, 0, 0};

/**
 * function calculating polar coordinates
 */
void generatePolarVector()
{
    /// r is the ratio of cumulative distribution and optical depth
    outputPolarVector[0] = (double) f_cum(random_number())/10;
    /// theta calculated from a non-flat distribution
    outputPolarVector[1] = (double) acos(2*random_number()-1);
    /// phi calculated using a random number generator
    outputPolarVector[2] = (double) random_number()*2*M_PI - M_PI;
    /// z is r multiplied by cos theta
    outputPolarVector[3] = (double) ((outputPolarVector[0])  * cos(outputPolarVector[1]));
}

/// un-normalised Rayleigh function
double rayleigh(double x)
{
    return (3/(16*M_PI))*(1+(cos(x)*cos(x)));
}

/**
 * / initialised vector for cartesian coordinates
 */
double outputCartesianVector[3] = {0, 0, 0};

/**
 * function calculating cartesian coordinates
 * @param r
 * @param theta
 * @param phi
 */
void generateCartesianVector(double r, double theta, double phi)
{
    outputCartesianVector[0] = (double) r*sin(theta)*cos(phi);
    outputCartesianVector[1] = (double) r*sin(theta)*sin(phi);
    outputCartesianVector[2] = (double) r*cos(theta);
}

/**
 * // xyz vector used to store cartesian coordinates
 */
typedef struct xyz_vector
{
    double x;
    double y;
    double z;
} xyz_vector;

/**
 * data type cartesian vector
 */
struct xyz_vector r2;

/**
 * initialised matrix for rotation
 */
double rotationMatrix[4][4];
double inputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};
double outputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};

/**
 * matrix multiplication function
 */
void multiplyMatrix()
{
    /// combine inverse and original rotation to give the relative output
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            outputMatrix[i][j] = 0;
            for(int k = 0; k < 4; k++){
                outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
            }
        }
    }
}

/**
 * rotation matrix function
 * @param angle
 *
 * @param u
 * @param v
 * @param w
 */
void RotationMatrix(double angle, double u, double v, double w)
{
    double L = (u*u + v*v + w*w);
    double u2 = u*u;
    double v2 = v*v;
    double w2 = w*w;
    /// ref - https://stackoverflow.com/questions/45160580/rotation-about-an-arbitrary-axis-in-3-dimensions-using-matrix
    /// rotation about arbitrary axis in 3 dimensions
    rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
    rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][3] = 0.0;
    rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
    rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][3] = 0.0;
    rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
    rotationMatrix[2][3] = 0;
    rotationMatrix[3][0] = 0;
    rotationMatrix[3][1] = 0;
    rotationMatrix[3][2] = 0;
    rotationMatrix[3][3] = 1;
}

/**
 * initialised vector for polar coordinates
 */
double outputRayleighScattering[3] = {0, 0, 0};

/**
 * function generating Rayleigh scattering polar coordinates
 * @param tau
 */
void generateRayleighScattering(float tau)
{

    /// r is the ratio of cumulative distribution and optical depth
    outputRayleighScattering[0] = (double) f_cum(random_number())/tau;

    /// theta calculated from a non-flat distribution using the rejection method
    /// obeying the Rayleigh phase function
    double theta_vals = (double) acos(2*random_number()-1);
    double y_vals = (double) (3/(8*M_PI))*random_number();
    while (y_vals > rayleigh(theta_vals))
    {
        theta_vals = (double) acos(2*random_number()-1);
        y_vals = (double) (3/(8*M_PI))*random_number();
    }
    outputRayleighScattering[1] = (double) theta_vals;

    /// phi calculated using a random number generator
    outputRayleighScattering[2] = (double) random_number()*2*M_PI - M_PI;
}

/**
 * initialised vector for polar coordinates
 */
double Generate_rotation[3] = {0, 0, 0};

/**
 * function rotating Rayleigh coordinates and tracking values
 * @param r
 * @param theta
 * @param phi
 * @param x
 * @param y
 * @param z
 */
void generateRotation(double r, double theta, double phi, double x, double y, double z) {
    /**
     * xyz vector used to store initial and final cartesian co-ordinates
     */
    typedef struct xyz_vector {
        double x1;
        double y1;
        double z1;
        double x2;
        double y2;
        double z2;
    } xyz_vector;

    /**
     * data type cartesian vector
     */
    struct xyz_vector r3;

    /// old xyz points
    r3.x1 = x;
    r3.y1 = y;
    r3.z1 = z;

    /// generate new cartesian co-ordinates from polars
    generateCartesianVector(r, theta, phi);

    /// new xyz points
    r3.x2 = outputCartesianVector[0];
    r3.y2 = outputCartesianVector[1];
    r3.z2 = outputCartesianVector[2];

    /// initial points to transform
    inputMatrix[0][0] = r3.x2;
    inputMatrix[1][0] = r3.y2;
    inputMatrix[2][0] = r3.z2;
    inputMatrix[3][0] = 1;

    /// axis vector
    double u = r3.x2 - r3.x1;
    double v = r3.y2 - r3.y1;
    double w = r3.z2 - r3.z1;

    /// rotation to generate arbitrary axis and combine matrixs to give output
    RotationMatrix(theta, u, v, w);
    multiplyMatrix();

    /// output co-ordiantes
    Generate_rotation[0] = outputMatrix[0][0];
    Generate_rotation[1] = outputMatrix[1][0];
    Generate_rotation[2] = outputMatrix[2][0];

}

int main()
{
    // 1 a

    /// initial variables
    int i, i_bin;
    int Naccept = 0;
    int N_samples = 1000000;
    int fmax = 1;

    /// call malloc to allocate heap space
    double *X = malloc(N_samples*sizeof(double));
    double *Y = malloc(N_samples*sizeof(double));

    // 1 b

    /// initial variables
    int N = 1000000;
    int N_remove = 0;
    int N_bins = 10;
    int binned[N_bins];
    double zmax = 1.0;
    double albedo = 1.0;
    double dx = 1./N_bins;

    // 1 c

    /// initial variables
    double x_init = 0;
    double y_init = 0;
    double z_init = 0;
    double blue_t = 10;
    double other_t = 0.1;
    double u, v, w;

    // #############################################################################
    /////////////////////////////////////  1 a   ///////////////////////////////////
    printf("############## 1 a ##############\n");
    // #############################################################################


    /// handling errors
    if ( (X == NULL) || (Y == NULL) ){
        printf("ERROR: malloc failed! \n");
        exit(0);
    }

    /// Generate 1000 random numbers for X and Y in the range 0-1 and 0-y_max
    for (i=0; i<N_samples; i++)
    {
        X[i] = (double) random_number()*6;
        Y[i] = (double) random_number() * fmax;
    }

    /// clock begins to measure time for the rejection method
    clock_t rbegin = clock();

    /// store output as a txt file to plot in python
    FILE *fp=NULL;
    fp=fopen("rejection_method.text", "w");

    /// Rejection Method
    for (i = 0; i < N_samples; i++) {

        /// rejecting and reassigning x values
        while (Y[i] > f(X[i])){
            X[i] = (double) random_number()*6;
            Y[i] = (double) random_number() * fmax;
        }

        /// accepted x values
        fprintf(fp, "%lf\n", X[i]);

    }

    /// clock ends
    clock_t rend = clock();
    double rejection_time_spent = (double)(rend - rbegin) / CLOCKS_PER_SEC;
    printf("rejection method: %1.8f us\n", rejection_time_spent*1000000);


    /// re-initialise values
    for (i=0; i<N_samples; i++)
    {
        Y[i] = (double) random_number() * fmax;
    }

    /// clock begins to measure time for the cumulative method
    clock_t begin = clock();

    /// store output as a txt file, to then plot
    FILE *fp1=NULL;
    fp1=fopen("cumulative_method.text", "w");

    /// Cumulative method
    for (i = 0; i < N_samples; i++) {
        X[i] = f_cum(Y[i]);

        fprintf(fp1, "%lf\n", X[i]);
    }

    /// clock ends
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("cumulative method: %1.8f us\n", time_spent*1000000);


    /// deallocating memory
    free(X);
    X = NULL;
    free(Y);
    Y = NULL;

    // #############################################################################
    /////////////////////////////////////  1 b   ///////////////////////////////////
    printf("############## 1 b ##############\n");
    // #############################################################################

    /// clear out the bins first and initialise values
    for (i = 0; i < N_bins; i++)
        binned[i] = 0;

    Naccept = 0;


    /// Generate photons
    while (Naccept < N) {
        /// While the total number of photons that have beena accepted is less than N = 1,000,000
        /// photons are generated in polar and z co-ordinates
        generatePolarVector();
        r1.r = outputPolarVector[0];
        r1.theta = outputPolarVector[1];
        r1.phi = outputPolarVector[2];
        r1.z = outputPolarVector[3];
        /// while the photon is between z = z_min = 0 and z = z_max = 1
        while (r1.z <= zmax && r1.z >= 0) {
            /// check for scattering or absorption
            if (albedo >= random_number()) {
                /// if the albedo is larger than a random probability
                /// the photon undergoes scattering and new polar and z values are generated
                generatePolarVector();
                r1.r = outputPolarVector[0];
                r1.theta = outputPolarVector[1];
                r1.phi = outputPolarVector[2];
                /// z value is updated
                r1.z = outputPolarVector[3] + r1.z;
            } else {
                /// if the albedo is less than the random probability
                /// the photon undergoes absorption and is removed
                r1.z = 0;
                break;
            }
        }
        /// photons have escaped the while loop condition
        if (r1.z >= zmax) {
            /// photons in the correct direction
            /// bin the photons respectively between 10 bins
            i_bin = (int) (cos(r1.theta) / dx);
            binned[i_bin]++;
            Naccept++;

        }
            /// photons have gone in the wrong direction or have been absorbed
        else {
            /// counter of removed photons increases
            N_remove++;
        }
    }

    /// update these values for bins
    printf("midpoint, Number of photons binned\n");
    for (i = 0; i < N_bins; i++) {
        printf("%1.3f,  %d\n", ((i + 0.5) * dx), binned[i]);
    }



    // #############################################################################
    /////////////////////////////////////  1 c   ///////////////////////////////////
    printf("############## 1 c ##############\n");
    // #############################################################################

    /// initialize values
    Naccept = 0;
    /// store output as a txt file, to then plot
    FILE *fp5=NULL;
    fp5=fopen("XYZ_other.text", "w");


    while (Naccept < N)
    {
        /// generate initial polar values for scattering
        generateRayleighScattering(other_t);
        r1.r = outputRayleighScattering[0];
        r1.theta = 0;
        r1.phi = outputRayleighScattering[2];
        /// generated rotation matrix following the initial polar conditions
        generateRotation(r1.r, r1.theta, r1.phi, x_init, y_init, z_init);
        /// initial injected xyz co-ordinates
        r2.x = Generate_rotation[0];
        r2.y = Generate_rotation[1];
        r2.z = Generate_rotation[2];

        while ( r2.z <= zmax && r2.z >= 0 ){
            /// check for scattering or absorption
            if (albedo >= random_number()) {
                /// if the albedo is larger than a random probability
                /// the photon undergoes Rayleigh scattering
                generateRayleighScattering(other_t);
                /// polar values re-generated
                r1.r = outputRayleighScattering[0];
                r1.theta = outputRayleighScattering[1];
                r1.phi = outputRayleighScattering[2];
                /// rotation matrix updated according to polar co-ordinates
                generateRotation(r1.r, r1.theta, r1.phi, r2.x, r2.y, r2.z);
                r2.x = Generate_rotation[0] + r2.x;
                r2.y = Generate_rotation[1] + r2.y;
                r2.z = Generate_rotation[2] + r2.z;


            } else {
                /// if the albedo is less than the random probability
                /// the photon undergoes absorption
                r2.z = 0;
                break;
            }
        }
        if (r2.z >= zmax){
            /// Tracking the number of accepted photons and storing cartesian co-ordiantes
            Naccept++;
            fprintf(fp5, "%lf\t%lf\t%lf\n", r1.z, r2.y, r2.z);

        }
            /// photons have gone in the wrong direction or have been absorbed
        else {
            /// counter of removed photons increases
            N_remove++;
        }

        }
    printf("The sky is blue");

    return 0;
}
