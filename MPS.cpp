/*** 
    Program to simulate a solitary wave in a straight channel in 2D using MPS 
     - All the quantites are in SI units.  
***/

#include <bits/stdc++.h>
using namespace std;

#define PI 3.14159265

//Parameters for the problem
const double d = 0.5, l = 7.515, w = 1, r = 0.045;  // depth and length of the water column.
const int Ny_walls = round(w/r); // number of nodes inx & y direction inside the domain.

double rho = 1000, g = 9.81, beta = 0.97, n_0 = 1000.0, mu = 1e-6, re = 2*r, n_air = 1.204; // parameters pre-defined.
double H,L; // wave parameters taken in input from user.
double total_time, dt; // stores total time of sim and time step respectively.

// Data structure of a particle in the simulation.
typedef struct particle 
{
    int ID;    // for fluid particles - 1 and for wall particles - 0
    double x;  // x coordinate of the particle.
    double y;  // y coordinate of the particle.
    double u;  // x velocity component of the particle.
    double v;  // x velocity component of the particle.
    double p;  // pressure of the particle.
    double n;  // particle density of the particle.
}particle;

// Pre-defining the functions required in main() program:
void Domain_discretisation (vector<particle>& Domain);
void Initialise (vector<particle>& Domain);
void Intermediate (vector<particle>& Domain);
void pressure_update (vector<particle>& Domain_old, vector<particle>& Domain_new);
void correction (vector<particle>& Domain_old, vector<particle>& Domain_new);
void Introduce_particles (vector<particle>& Domain_new, double t);
void Impose_BC (vector<particle>& Domain_new);

// Pre-defining the functions required for implementation of programs used in main():
double particle_density (vector<particle>& Domain, int i);
bool ifPresent (vector<particle>& Domain_new, double x, double y);
bool ifOut(particle ele);
double delP(vector<particle>& Domain, int component, int i);
double delsquare (vector<particle>& Domain, int function, int i);
double weight (double r);

int main ()
{
    // Taking input the wave characteristics.
    cout << "Enter the following parameters for the wave - \n";
    cout << "Wavelength(L):\n";
    cin >> L;
    cout << "Wave Height(H):\n";
    cin >> H;

    // Taking input the total time of simulation and time step.
    cout << "Enter the total time -\n";
    cin >> total_time;
    cout << "Enter the time step -\n";
    cin >> dt;

    // starting a output file to store output data.
    FILE* f;
    f = fopen("Output-SW.dat", "w");
    fprintf(f, "Time\tID\tx\t\ty\tVel_X\tVel_y\tPress\tParticle Density\n");

    // Matrices to store the field variables values for surface and domain nodes for previous and current time step.
    // They are stored in the matrix in the t,x,y,Vx,Vy,P columns in this order
    vector<particle> Domain_old, Domain_new;

    // First discretise the domain and free surface of the wave.
    Domain_discretisation(Domain_old);
    Domain_discretisation(Domain_new);

    // Initialising the solution domain.
    Initialise(Domain_old);
    Initialise(Domain_new);

    // Storing the information for the initial stage, i.e, time = 0
    fprintf (f, "\n");
    for (int i=0; i<Domain_old.size(); i++)
        fprintf(f, "%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n", 0.0, Domain_old[i].ID, Domain_old[i].x, Domain_old[i].y, Domain_old[i].u, Domain_old[i].v, Domain_old[i].p, Domain_old[i].n);

    cout << "\nThe discretisation and initialisation is done!\n";
    cout << "Total number of particles initially are - " << Domain_new.size() << endl;

    // Looping through all the time steps.
    for (double t=dt; t<=total_time; t=t+dt)
    {   
        Intermediate (Domain_old);
        
        pressure_update (Domain_old, Domain_new);

        correction (Domain_old, Domain_new);

        Introduce_particles (Domain_new, t);

        Impose_BC (Domain_new);

        // Storing the information for the current time step, i.e, time = t
        fprintf (f, "\n");
        for (int i=0; i<Domain_new.size(); i++)
            fprintf(f, "%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n", t, Domain_new[i].ID, Domain_new[i].x, Domain_new[i].y, Domain_new[i].u, Domain_new[i].v, Domain_new[i].p, Domain_new[i].n);

        // Copying particles in Domain_new to Domain_old vector for next time step.
        Domain_old = Domain_new;

        cout << "\n" << (int)t/dt << "th time step is completed!\n";
        cout << "Total number of particles now are - " << Domain_new.size() << endl;
    }

    fclose(f);

    cout << "The program execution if finished!!\n";

    return 0;
}

// function to discretising the internal domain nodes for the initial condition.
void Domain_discretisation (vector<particle>& Domain)
{
    for (double x=0.0; x<=(l+3*r); x=x+r)
    {
        // For the particles between the side two walls
        if ( (x > r) && (x < (l+2*r) ) )
        {
            double eta = H/pow(cosh((x-2*r)/L), 2);
            for (double y=0.0; y<=d+r+eta; y=y+r)
            {
                // Fluid particles.
                if (y > r)
                {
                    particle p;
                    p.ID = 1;
                    p.x = x;
                    p.y = y;
                    Domain.push_back(p); 
                }
                // floor particles. 
                else 
                {
                    particle p;
                    p.ID = 0;
                    p.x = x;
                    p.y = y;
                    Domain.push_back(p);            
                }
            }
        }
        
        // For the side walls particles.
        else
        {
            for (double y=0.0; y<=w+r; y=y+r)
            {
                particle p;
                p.ID = 0;
                p.x = x;
                p.y = y;
                Domain.push_back(p); 
            }
        }
    }
}

// function to initialise all the nodes.
void Initialise (vector<particle>& Domain)
{
    double c = sqrt(g*d)*(1+H/(2*d));
    for (int i=0; i<Domain.size(); i++)
    {
        // For wall particles.
        if (Domain[i].ID == 0)
        {
            // For particles at the inlet / outlet walls.
            if (Domain[i].x == 0.0 || Domain[i].x == r || Domain[i].x == l+2*r || Domain[i].x == l+3*r)
            {
                // the wave elevation that could have been there if wave extended.
                double eta = H/pow(cosh((Domain[i].x-2*r)/L), 2);
                // if particle is below wave elevation then its considered water particle
                if (Domain[i].y < d+r+eta)
                {
                    Domain[i].u = 0.0;
                    Domain[i].v = 0.0; 
                    Domain[i].p = 101325 + (rho*c*c)/(2*pow(cosh(2*PI*(Domain[i].x-2*r)/L), 2)); // pressure value based on the water particles which could have been there for the wave considered.
                    Domain[i].n = n_0;
                }
                // if particle is above wave elevation then its considered air particle
                else
                {
                    Domain[i].u = 0.0;
                    Domain[i].v = 0.0; 
                    Domain[i].p = 0.0; // pressure value based on the air particles which could have been there for the wave considered.
                    Domain[i].n = n_air;
                }
            }

            // For particles below the water level depth.
            else 
            {

                Domain[i].u = 0.0;
                Domain[i].v = 0.0; 
                Domain[i].p = 101325 + (rho*c*c)/(2*pow(cosh(2*PI*(Domain[i].x-2*r)/L), 2)); // pressure value based on the water particles which could have been there for the wave considered.
                Domain[i].n = n_0;
            }
        }

        // For water particles in the domain.
        else if (Domain[i].ID == 1)
        {
            double x = Domain[i].x - 2*r;
            double y = Domain[i].y - 2*r;
            double eta = H/pow(cosh(x/L), 2);

            Domain[i].u = (eta*sqrt(g*d))/d; // x-velocity based on Korteweg–de Vries (KdV).
            Domain[i].v = 0.0; // y-velocity based on Korteweg–de Vries (KdV).
            
            Domain[i].p = 101325 + (rho*c*c)/(2*pow(cosh(2*PI*(Domain[i].x-2*r)/L), 2)); // pressure value.
            Domain[i].n = n_0;
        }
        
        else
            cout << "Particle with incorrect(random) ID is present!" << endl;
    }
}

// function to calculate intermediate state of all the particles since wall particle also acts as fluid during pressure calculation.
void Intermediate (vector<particle>& Domain)
{
    for (int i=0; i<Domain.size(); i++)
    {
        Domain[i].u = Domain[i].u + dt*(mu*delsquare(Domain, 3, i));
        Domain[i].v = Domain[i].v + dt*(mu*delsquare(Domain, 4, i));
        Domain[i].x = Domain[i].x + dt*Domain[i].u;
        Domain[i].y = Domain[i].y + dt*Domain[i].v;
    }
    for (int i=0; i<Domain.size(); i++)
        Domain[i].n = particle_density(Domain, i);
}

// function to calculate particle density present at ith position in given particles vector.
double particle_density (vector<particle>& Domain, int i)
{
    double sum = 0.0, size = Domain.size();
    for (int j=0; j<size; j++)
    {
        if (j != i)
        {
            double distance = sqrt(pow(Domain[j].x-Domain[i].x, 2) + pow(Domain[j].y-Domain[i].y, 2));
            sum = sum + weight(distance);
        }
    }
}

// function to calculate pressure for next time step in Domain_new vector from values at intermediate time step in Domain_old.
void pressure_update (vector<particle>& Domain_old, vector<particle>& Domain_new)
{
    int size = Domain_old.size();
    double coeff[size][size], p[size], constant[size];
    
    // Calculating the coeff matrix containing the distances values as given in derivation.
    for (int i=0; i<size; i++)
    {
        for (int j=0; j<size; j++)
        {
            // Diagonal elements which is sum of distances with all other particles.
            if (i==j)
            {
                double sum = 0.0;
                for (int k=0; k<size; k++)
                {
                    if (i!=k)
                    {
                        double distance = sqrt(pow(Domain_old[k].x-Domain_old[i].x, 2) + pow(Domain_old[k].y-Domain_old[i].y, 2));
                        sum = sum + weight(distance);
                    }
                }
                coeff[i][j] = sum;
            }
            // Non-diagonal elements which is distance between i-j position element
            else
            {
                double distance = sqrt(pow(Domain_old[i].x-Domain_old[j].x, 2) + pow(Domain_old[i].y-Domain_old[j].y, 2));
                coeff[i][j] = weight(distance);
            }
        }
    }

    // Calculating the RHS matrix based on particle density diff for each ith particle.
    for (int i=0; i<size; i++)
    {
        double sum_n=0.0, sum_d=0.0, lambda;
        for (int j=0; j<size; j++)
        {
            double distance = sqrt(pow(Domain_old[j].x-Domain_old[i].x, 2) + pow(Domain_old[j].y-Domain_old[i].y, 2));
            sum_n = sum_n + pow(distance,2)*weight(distance);
            sum_d = sum_d + weight(distance);
        }
        lambda = sum_n / sum_d; 

        constant[i] = ((lambda*n_0)/4)*(rho/pow(dt,2))*(Domain_old[i].n/n_0 - 1);
    }

    //Using Gauss elimination to solve the linear equation [coeff] x [p] = [constant]
    double A[size][size+1];
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j <= size; j++)
        {
            if (j == size)
                A[i][j] = constant[i];
            else
                A[i][j] = coeff[i][j];
        }
    }

    double l;
    for (int k = 0; k < size; k++)
    {
        for (int i = 0; i <= size; i++)
        {
            l = A[i][k];

            for (int j = 0; j <= size; j++)
            {
                if (i != k)
                    A[i][j] = (A[k][k] * A[i][j]) - (l * A[k][j]);
            }
        }
    }

    for (int m=0; m<size; m++)
        p[m] = A[m][3] / A[m][m];

    // Substituting the calculated pressure values in the Domain_new vector.
    for (int i=0; i<size; i++)
    {
        Domain_new[i].p = p[i];
    }
}

// function to calculate correct position for fluid particles into Domain_new and wall particles have correct x,y,u,v,n already.
void correction (vector<particle>& Domain_old, vector<particle>& Domain_new)
{
    for (int i=0; i<Domain_new.size(); i++)
    {
        if (Domain_new[i].ID == 1)
        {
            Domain_new[i].u = Domain_old[i].u - (dt/rho)*delP(Domain_new, 0, i);
            Domain_new[i].v = Domain_old[i].v - (dt/rho)*delP(Domain_new, 1, i);
            Domain_new[i].x = Domain_old[i].x + dt*Domain_old[i].u;
            Domain_new[i].y = Domain_old[i].y + dt*Domain_old[i].v;
        }
    }
    // updating particle density of all particles once new positions are updated in Domain_new.
    for (int i=0; i<Domain_new.size(); i++)
    {
        if (Domain_new[i].ID == 1)
            Domain_new[i].n = particle_density(Domain_new, i);
    }
}

// function to introduce new particles from left wall inlet and remove particles now outside right wall outlet.
void Introduce_particles (vector<particle>& Domain_new, double t)
{
    // removing particles now outside the domain from all four corners.
    Domain_new.erase(remove_if(Domain_new.begin(), Domain_new.end(), [](particle x) { return ifOut(x); }), Domain_new.end());
    
    // position in Domain_new vector where we start introducing new elements in the vector.
    int position = 2*Ny_walls+2;

    // the elevation of the surface particle at that time step at inlet boundary.
    double c = sqrt(g*d)*(1+H/(2*d));
    double y = H/pow(cosh( ((-1)*c*t)/L ), 2);

    // looping vertically till the elevation reached using particles of size r.
    for (double j = 2*r; j<=y; j=j+r)
    {
        double x = 2*r;
        if (ifPresent(Domain_new, x, j))
            continue;

        // if particle at that calculated position not present.
        else
        {
            particle p;
            p.ID = 1;
            p.x = x;
            p.y = j;

            double eta = y;
            p.u = (eta*sqrt(g*d))/d; // x-velocity based on Korteweg–de Vries (KdV) equation.
            p.v = 0.0; // y-velocity based on Korteweg–de Vries (KdV) equation.
            
            p.p = 101325 + (rho*c*c)/(2*pow(cosh((-1)*2*PI*c*t/L), 2)); // pressure value.

            Domain_new.insert(Domain_new.begin()+position, p); // particles added just after first set of the wall particles.
        }
    }

}

// function to impose free surface and bottom surface BC for fluid particles.
void Impose_BC (vector<particle>& Domain_new)
{
    for (int i=0; i<Domain_new.size(); i++)
    {
        if (Domain_new[i].n < beta*n_0 && Domain_new[i].ID == 1)
            Domain_new[i].p = 0.0;
        else if (Domain_new[i].y == 2*r && Domain_new[i].ID == 1)
            Domain_new[i].v = 0.0;
    }
}

// function to check if given coordinates of a particle present in given Domain vector.
bool ifPresent (vector<particle>& Domain_new, double x, double y)
{
    for (int i=0; i<Domain_new.size(); i++)
    {
        if (x == Domain_new[i].x && y == Domain_new[i].y)
            return true;
        else
            return false;
    }
}

// function to see if given particle is outside the four extreme locations in considered domain lengths.
bool ifOut(particle ele)
{
    if (ele.ID==1 && (ele.x >= l+r || ele.y >= w-r || ele.y < 2*r || ele.x < 2*r ) )
        return true;
    else 
        return false;
}

// function to calculate divergence of pressure of given position particle in Domain vector in component direction (x or y)
double delP(vector<particle>& Domain, int component, int i)
{
    // to calculate the min value of pressure in the region of influence.
    int min=i;
    for (int j=0; j<Domain.size(); j++)
    {
        double distance = sqrt(pow(Domain[i].x - Domain[j].x, 2) + pow(Domain[i].y - Domain[j].y, 2));
        if (distance <= re)
        {
            if (Domain[j].p < Domain[min].p)
                min = j;
        }
    }

    // calcuate the divergence value based on derivation mentioned in the given direction component.
    double sum=0.0;
    for (int j=0; j<Domain.size(); j++)
    {
        double distance = sqrt(pow(Domain[i].x - Domain[j].x, 2) + pow(Domain[i].y - Domain[j].y, 2));
        if (j!=i)
        {
            if (component == 0)
                sum = sum + (2/n_0)*( (Domain[j].p-Domain[i].p)/pow(distance,2) )*weight(distance)*(Domain[j].x-Domain[i].x);
            else if (component == 1)
                sum = sum + (2/n_0)*( (Domain[j].p-Domain[i].p)/pow(distance,2) )*weight(distance)*(Domain[j].y-Domain[i].y);
        }
    }

    return sum;
}

// function to calculate laplacian of u or v at particle present in ith location of Domain vector.
double delsquare (vector<particle>& Domain, int function, int i)
{
    double sum_n=0.0, sum_d=0.0, lambda, value=0.0, size = Domain.size();
    
    // calculate the lambda value based on the mentioned derivation.
    for (int j=0; j<size; j++)
    {
        double distance = sqrt(pow(Domain[j].x-Domain[i].x, 2) + pow(Domain[j].y-Domain[i].y, 2));
        sum_n = sum_n + pow(distance,2)*weight(distance);
        sum_d = sum_d + weight(distance);
    }
    lambda = sum_n / sum_d;

    // calculate the laplacian scalar value.
    for (int j=0; j<size; j++)
    {
        if (function == 3)
        {
            double distance = sqrt(pow(Domain[j].x-Domain[i].x, 2) + pow(Domain[j].y-Domain[i].y, 2));
            value = value + (4/(lambda*n_0))*(Domain[j].u - Domain[i].u)*weight(distance);
        }
        if (function == 4)
        {
            double distance = sqrt(pow(Domain[j].x-Domain[i].x, 2) + pow(Domain[j].y-Domain[i].y, 2));
            value = value + (4/(lambda*n_0))*(Domain[j].v - Domain[i].v)*weight(distance);
        }
    }

    return value;
}

// function to calcuate weight for a given distance r between two particles based on considered radius of influence (re).
double weight (double r)
{
    double weight;

    if (r<=re)
        weight = (re/r) - 1;
    else
        weight = 0.0;

    return weight;
}

