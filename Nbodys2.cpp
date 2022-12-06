/*!
 @mainpage Title
 
 \section intro_sec Introduction
The program (written in C++) analyzes an N-body system. It takes in a CSV 
file with initial conditionsin the following format: 
mass, x, y, z, Vx, Vy, Vz, q.
Where (x, y, z) are the coordinates of the body, (Vx, Vy, Vz) are its velocity 
components, q is its charge, and mass is the mass.The other inputs taken by 
the program are choosing between the gravity or electric simulators(g/e) and 
then inputting the time of analysis(ex. 400 days in the gravity simulator).
Both of these are terminal inputs. Upon completion, the program outputs N+1 
CSV files. N CSV files contain the position and velocity of each body after
every calculation. The last CSV file contains the position of every body 
after each calculation set. This last file, called “data.tmp” can then be 
used to graph the path taken by each body. The program also output the time 
taken to run the complete analysis and the final positions and velocities 
of each body(terminal outputs).This is in the following format: 
(x, y, z) | (Vx, Vy, Vz). 

 @author Samarth Kumar
 @date December 5th, 2022
*/


#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <getopt.h>
#include <stdlib.h>

using std::string;

/*!
 @brief Vector Class

This class is used to denote a vector in 3d. It consists of three doubles,
corresponding to the components of a vector in the x, y, and z directions.
It also contains functions to add and subtract vectors, and also to take 
the dot and cross products of two vectors. There is also a function that 
overloads the operator << to allow direct output of the vector in the format:
(x, y, z).
*/
class Vector
{
private:
    double _x; /*!< x component*/
    double _y; /*!< y component*/
    double _z; /*!< z component*/

public:
    /*!
     @brief Constructor for Vector class.

    Takes in three double values and initialises the vector with them. 

    @param x the x component of the vector
    @param y the y component of the vector
    @param z the z component of the vector
    @return void
    */
    Vector(double x, double y, double z)
    {
        _x=(x);
        _y=(y);
        _z=(z);
    }

    /*!
    @brief Manipluator function for the x component of the vector
    
    Changes the value of the x component to the new value.
    @param x the new x component
    @return void
    */
    void x(double x) 
    {
        // manipulator functions(to change specific values)
        _x = x;
    }

    /*!
    @brief Manipluator function for the y component of the vector
    
    Changes the value of the y component to the new value.
    @param y the new y component
    @return void
    */
    void y(double y)
    {
        _y = y;
    }

    /*!
    @brief Manipluator function for the z component of the vector
    
    Changes the value of the z component to the new value.
    @param z the new z component
    @return void
    */
    void z(double z)
    {
        _z = z;
    }

    /*!
    @brief Accessor function for the x component of the vector
    
    @return value of the x component
    */
    double x(void) const 
    {
        // accessor functions
        return _x;
    }

    /*!
    @brief Accessor function for the y component of the vector
    
    @return value of the y component
    */
    double y(void) const
    {
        return _y;
    }

    /*!
    @brief Accessor function for the z component of the vector
    
    @return value of the z component
    */
    double z(void) const
    {
        return _z;
    }

    /*!
    @brief Addition of Two Vectors
    
    Overloads the '+' operator to add vectors the correct way.
    */
    Vector operator+(const Vector &obj)
    { 
        // overloading the + operator
        Vector temp = Vector(0, 0, 0);
        temp._x = _x + obj._x;
        temp._y = _y + obj._y;
        temp._z = _z + obj._z;
        return temp;
    }
    
    /*!
    @brief Subtraction of Two Vectors
    
    Overloads the '-' operator to subtract vectors the correct way.
    */
    Vector operator-(const Vector &obj)
    { 
        // overloading the - operator
        Vector temp = Vector(0, 0, 0);
        temp._x = _x - obj._x;
        temp._y = _y - obj._y;
        temp._z = _z - obj._z;
        return temp;
    }

    /*!
    @brief Overloads Output Operator
    
    This function overloads the output operator 'std::cout<<'. Vector can be
    outputted directly using std::cout<<. The format will be (x, y, z).
    */
    friend std::ostream &operator<<(std::ostream &os, const Vector &obj)
    { 
        // overloading the << operator
        os << "(" << obj._x << ", " << obj._y << ", " << obj._z << ")";
        return os;
    }

    /*!
    @brief Dot Product Function
    
    Function created to calculate the dot product of two vectors. 
    @param a One of the vectors for dot product
    @param b The other vector for dot product
    @return A number. 
    */
    double dot(Vector a, Vector b);

    /*!
    @brief Cross Product Function
    
    Function created to calculate the cross product of two vectors. 
    @param a One of the vectors for cross product
    @param b The other vector for cross product
    @return Another vector, of the form (x, y, z). 
    */
    Vector cross(Vector a, Vector b);
};


double Vector::dot(Vector a, Vector b) 
{
    // dot product
    double temp = 0;
    temp = a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
    return temp;
}
Vector Vector::cross(Vector a, Vector b) 
{
    // cross product
    Vector temp = Vector(0, 0, 0);
    temp._x = a._y * b._z - a._z * b._y;
    temp._y = a._z * b._x - a._x * b._z;
    temp._z = a._x * b._y - a._y * b._x;
    return temp;
}

/*!
 @brief Particle Class

This class is used to denote a particle(celestial body/atomic body).
It consists of two vectors,position and velocity; and two doubles; charge and mass.
It contains a function that overloads the output operator '<<' to directly 
output the properties of the particle. 
*/
class Particle
{
private:
    double _mass; /*!<Mass of the particle*/
    double _charge; /*!<Charge of the particle*/
    Vector Position = Vector(0, 0, 0); /*!<Position of the particle*/
    Vector Velocity = Vector(0, 0, 0); /*!<Velocity of the particle*/

public:
    /*!
    Constructor for Particle Class

    Takes in 8 doubles and initialises the class. 
    @param mass mass of the particle
    @param x x component of position
    @param y y component of position
    @param z z component of position
    @param vx x component of velocity
    @param vy y component of velocity
    @param vz z component of velocity
    @param charge charge of the particle
    @return void
    */
    Particle(double mass, double x, double y, double z, double vx,
             double vy, double vz, double charge)
    {
        _mass = mass;
        Position.x(x), Position.y(y), Position.z(z);
        Velocity.x(vx), Velocity.y(vy), Velocity.z(vz);
        _charge = charge;
    }

    /*!
    Output Operator Overloader 

    Overloads the output operator to directly output the particle properties. 
    Can use 'std::cout<<' directly for output .
    */
    friend std::ostream &operator<<(std::ostream &os, const Particle &obj)
    { 
        // overloading the << operator
        os << "Mass = " << obj._mass << "\nPosition is " << obj.Position <<
         "\nVelocity is " << obj.Velocity << "Charge = " << obj._charge << std::endl;
        return os;
    }

    /*!
    Manipulator function for Particle Class

    Allows multiple manipulations at the same time.(This funciton is redundent)
    @param mass mass of the particle
    @param x x component of position
    @param y y component of position
    @param z z component of position
    @param vx x component of velocity
    @param vy y component of velocity
    @param vz z component of velocity
    @param charge charge of the particle
    @return void
    */
    void input(double mass, double x, double y, double z, double vx,
               double vy, double vz, double charge)
    {
        _mass = mass;
        Position.x(x), Position.y(y), Position.z(z);
        Velocity.x(vx), Velocity.y(vy), Velocity.z(vz);
        _charge = charge;
    }

    /*!
    Acessor function for position

    @return Ouputs the position vector of the particle.
    */
    Vector P(void) const 
    {
        // accesors
        return Position;
    }
    
    /*!
    Acessor function for velocity

    @return Ouputs the velocity vector of the particle.
    */
    Vector V(void) const
    {
        return Velocity;
    }

    /*!
    Acessor function for mass

    @return Ouputs the mass of the particle.
    */
    double M(void) const
    {
        return _mass;
    }

    /*!
    Acessor function for charge

    @return Ouputs the charge of the particle.
    */
    double Q(void) const
    {
        return _charge;
    }

    /*!
    Manipluator function for the mass of the particle
    
    Changes the value of the mass.
    @param mass the new mass
    @return void
    */
    void M(double mass)
    {
        _mass = mass;
    }

    /*!
    Manipluator function for the position of the particle
    
    Changes the values of the position vector of the particle.
    @param x x component of position
    @param y y component of position
    @param z z component of position
    @return void
    */
    void P(double x, double y, double z)
    {
        Position = Vector(x, y, z);
    }

    /*!
    Manipluator function for the velocity of the particle
    
    Changes the values of the velcoity vector of the particle.
    @param x x component of velocity
    @param y y component of velocity
    @param z z component of velocity
    @return void
    */
    void V(double x, double y, double z)
    {
        Velocity = Vector(x, y, z);
    }

    /*!
    Manipluator function for the charge of the particle
    
    Changes the value of the charge.
    @param charge the new charge
    @return void
    */
    void Q(double charge)
    {
        _charge = charge;
    }
   
};

/*!
 @brief Base Class

This class consits of a vector called Force corresponding with the force on 
the particle. It also contain virutal functions that the inherited classes 
define. This class basically gives an outline of for how the force calculations
on the particles is to be done. 
*/
class Simulation_Engine 
{
    // Force
private:
    Vector Force = Vector(0, 0, 0);/*!<Force on a particle*/

public:

    /*!
    Constructor for Simulation_Engine Class

    Takes in 3 doubles and initialises the class. It basically initialises 
    the force vector.  
    @param fx x component of force
    @param fy y component of force
    @param fz z component of force
    @return void
    */
    Simulation_Engine(double fx, double fy, double fz)
    {
        Force.x(fx), Force.y(fy), Force.z(fz);
    }

    /*!
    Overloads Output Operator
    
    This function overloads the output operator 'std::cout<<'. The force can be
    outputted directly as a vector using std::cout<<. 
    */
    friend std::ostream &operator<<(std::ostream &os,
                                    const Simulation_Engine &obj)
    {
        os << obj.Force << std::endl;
        return os;
    }
    
    /*!
    Acessor function for force

    @return Ouputs the force vector
    */
    Vector F(void)
    {
        return Force;
    }

    /*!
    Manipulator function for force

    Takes in 3 doubles corresponding to the components of force and changes 
    their values. 

    @param fx x component of force
    @param fy y component of force
    @param fz z component of force
    @return void
    */
    void F(double fx, double fy, double fz)
    {
        Force = Vector(fx, fy, fz);
    }

    /*!
    Calculates conditions after time 't'

    This function calcualtes the position and velocities of the particle
    after time t has passed, using Newton's Laws of motion.  It takes in
    the acceleration calculated by the force function. This function along 
    with the force function are responsible for the main analysis. The 
    smaller the t, the more correct the answer.(t=25 for gravity sim) 
    @param a the particle that has to be analysed
    @param F_net the net force on that particle
    @param t the time interval of analysis
    @return return the same particle with new position and velocity components
    */
    virtual Particle run(Particle a, Vector F_net, double t) = 0;

    /*!
    Calculates the force on a particle due to another

    Calculates the force on particle a due to particle b.
    @param a Particle that is being affected. 
    @param b Particle that is causing the effect.
    @return A force vector. 
    */
    virtual Vector force(Particle a, Particle b) = 0;
};

/*!
 @brief Inherted class specifcally for gravitational simulations

Defines the virtual functions stated in the base class using basic laws of 
gravitation.
*/
class Gravity_Sim : public Simulation_Engine
{
    // for particle 'a'
    // run the force function N times
    // for *each* particle, then sum
    //  to get the total.
public:
    Gravity_Sim(double fx, double fy, double fz)
        : Simulation_Engine(fx, fy, fz){};

    Vector force(Particle a, Particle b);

    Particle run(Particle a, Vector F_net, double t);
};


Vector Gravity_Sim::force(Particle a, Particle b)
    { 
        // force on a due to b
        Vector obj = Vector(b.P().x() - a.P().x(), b.P().y() - a.P().y(),
                            b.P().z() - a.P().z());
        // std::cout << obj.x()<< " "<< obj.y() << " " << obj.z() << std::endl;
        double fm = (6.6743e-11) * a.M() * b.M() / pow(obj.dot(obj, obj), 1.5); 
        // replace 1 with 6.67e-11
        // force magnitude
        F(fm * obj.x(),
          fm * obj.y(),
          fm * obj.z()); // force vector
                         // is the force vector (with x, y, z, components) 
                         //on 'a' due to 'b'
        // std::cout << F() << std::endl;
        return F();
    }
Particle Gravity_Sim::run(Particle a, Vector F_net, double t)
{
    Particle old = Particle(a.M(), a.P().x(), a.P().y(), a.P().z(),
                            a.V().x(), a.V().y(), a.V().z(), 0.0);
    // std::cout << old.P() << " | " << old.V() << std::endl;
    a.V(old.V().x() + F_net.x() * t / old.M(),
        old.V().y() + F_net.y() * t / old.M(),
        old.V().z() + F_net.z() * t / old.M());
    a.P(old.P().x() + old.V().x() * t + F_net.x() * (pow(t, 2)) / (2 * old.M()),
        old.P().y() + old.V().y() * t + F_net.y() * (pow(t, 2)) / (2 * old.M()),
        old.P().z() + old.V().z() * t + F_net.z() * (pow(t, 2)) / (2 * old.M()));
    // std::cout << a.P() << " | " << a.V() << std::endl;
    return a;
}

/*!
 @brief Inherted class specifcally for electrostatic simulations

Defines the virtual functions stated in the base class using columbs law of
elcectrostatics.
*/
class Electro_Sim : public Simulation_Engine
{
    // for particle 'a'
    // run the force function N times
    // for *each* particle, then sum
    //  to get the total.
public:
    Electro_Sim(double fx, double fy, double fz)
        : Simulation_Engine(fx, fy, fz){};

    Vector force(Particle a, Particle b);

    Particle run(Particle a, Vector F_net, double t);
};


Vector Electro_Sim::force(Particle a, Particle b)
{ 
        // force on a due to b
        Vector obj = Vector(b.P().x() - a.P().x(), b.P().y() - a.P().y(),
                            b.P().z() - a.P().z());
        // std::cout << obj.x()<< " "<< obj.y() << " " << obj.z() << std::endl;
        double fm = (8.99e+9) * a.Q() * b.Q() / pow(obj.dot(obj, obj), 1.5); 
        // replace 1 with 6.67e-11
        // force magnitude
        F(fm * obj.x(),
          fm * obj.y(),
          fm * obj.z()); // force vector
                         // is the force vector (with x, y, z, components) 
                         //on 'a' due to 'b'
        // std::cout << F() << std::endl;
        return F();
}
Particle Electro_Sim::run(Particle a, Vector F_net, double t)
{
        Particle old = Particle(a.M(), a.P().x(), a.P().y(), a.P().z(),
                                a.V().x(), a.V().y(), a.V().z(), 0.0);
        // std::cout << old.P() << " | " << old.V() << std::endl;
        a.V(old.V().x() + F_net.x()*t / old.M(),
            old.V().y() + F_net.y()*t / old.M(),
            old.V().z() + F_net.z()*t / old.M());
        a.P(old.P().x() + old.V().x()*t + F_net.x()*(pow(t,2)) / (2 * old.M()),
            old.P().y() + old.V().y()*t + F_net.y()*(pow(t,2)) / (2 * old.M()),
            old.P().z() + old.V().z()*t + F_net.z()*(pow(t,2)) / (2 * old.M()));
         //std::cout << a.P() << " | " << a.V() << std::endl;
        return a;
}


/*!
 @brief Execution Time Function

 A varidiac function that takes in any function and its inputs, and outputs 
 results of the function along with the time taken to run it.
 @param A Pointer to a funtion
 @param var1 All the other input variables required for the function
 @return Time taken to run and outputs of funciton
*/
template <typename A, typename... Types>
    void Exec_Time(A(*func)(Types...), Types... var1)
    {
        time_t start, end;
        time(&start);
        std::ios_base::sync_with_stdio(false);
        std::cout << func(var1...) << std::endl;
        time(&end);
        double time_taken = double(end - start);
        std::cout << "Time taken by function is : " << std::fixed
                  << time_taken << std::setprecision(5);
        std::cout << " sec " << std::endl;
    }

/*!
 @brief Execution Function

 Function that has all the code required to execute the program.
 
*/
int start(void)
{
    int N = 0;
    // Vector a = Vector(1, 2, 3);
    // Vector b = Vector(2, 3, 4);
    // Vector c = Vector(0, 0, 0);
    // c = a - b;
    // std::cout << c << std::endl;
    // c = a + b;
    // std::cout << c << std::endl;
    // std::cout << pow(c.dot(a-b, a-b), 1.5) << std::endl;
    //  std::cout << " \n";
    //  Particle A = Particle(5, 1, 1, 1, 0, 0, 0, 0);
    //  Particle B = Particle(5, 2, 2, 2, 0, 0, 0, 0);
    //  Gravity_Sim test = Gravity_Sim(0, 0, 0);
    //  c = test.force(A, B);
    //  std::cout << c.x() << " " << c.y() << " " << c.z() << std::endl;
    //  Particle C = test.run(A, c, 500);
    //  A.input(C.M(), C.P().x(), C.P().y(), C.P().z(), C.V().x(), C.V().y(), C.V().z(), C.Q());
    //  std::cout << A.P() << " | " << A.V() << std::endl;
    //  // << std::endl;

     std::ifstream myFile("settings.txt");
     string line;
     char thing;
     std::cout << "What simulation do you want to run[g/e]: " << std::endl;
     std::cin >> thing;
     while (std::getline(myFile, line))
     {
         N++;
    }
   // std::cout << N << std::endl;
    myFile.clear();
    myFile.seekg(0);

    Particle *Par = (Particle *)malloc(sizeof(Particle) * (N - 1));
    // number of particles
    std::getline(myFile, line);
     //std::cout << N << std::endl;
    int i = 0;
    while (std::getline(myFile, line))
    {
        double mass, x, y, z, vx, vy, vz, charge;
        string temp;
        std::stringstream ss(line);
        std::getline(ss, temp, ',');
        mass = std::stod(temp);
        //std::cout << N << std::endl;
        std::getline(ss, temp, ',');
        x = std::stod(temp);
        //std::cout << N << std::endl;
        std::getline(ss, temp, ',');
        y = std::stod(temp);
        //std::cout << N << std::endl;
        std::getline(ss, temp, ',');
        z = std::stod(temp);
        //std::cout << N << std::endl;
        std::getline(ss, temp, ',');
        vx = std::stod(temp);
        //std::cout << N << std::endl;
        std::getline(ss, temp, ',');
        vy = std::stod(temp);
        //std::cout << N << std::endl;
        std::getline(ss, temp, ',');
        vz = std::stod(temp);
        //std::cout << N << std::endl;
        std::getline(ss, temp, ',');
        charge = std::stod(temp);
        //std::cout << N << std::endl;
        Par[i] = Particle(mass, x, y, z, vx, vy, vz, charge);
       // std::cout << N << std::endl;
        i++;
    }
    //std::cout << N << std::endl;
    // std::cout << Par[0].P() << " | " << Par[0].V() << std::endl;
    // std::cout << Par[1].P() << " | " << Par[1].V() << std::endl;
    //std::cout << i << std::endl;
    // for (int j = 0; j < N - 1; j++)
    // {
    //     std::cout << Par[j].P() << " | " << Par[j].V() << std::endl;
    // }
    myFile.close();
    //std::cout << i << std::endl;
    i = 0;
    // Vector *forceM = (Vector *)malloc(sizeof(Vector) * (pow(N - 1, 2))); 
    // 2d matrix of forces
    // for (int j = 0; j < (N - 1)*(N-1); j++)
    // {
    //     forceM[j] = Vector(0, 0, 0);
    // }
    // for (i, j) force(i*j+j);

    Vector *force_net = (Vector *)malloc(sizeof(Vector) * (N - 1)); 
    // net force of each particle
    for (int j = 0; j < N - 1; j++)
    {
        force_net[j] = Vector(0, 0, 0);
    }

    Particle ph = Particle(0, 0, 0, 0, 0, 0, 0, 0); // place holder

    FILE *fp = NULL;
    FILE *gnupipe = NULL;
    char *GnuCommands[] = { "set title \"N-Body System\"", "plot 'data.tmp'"};
    fp = fopen("data.tmp", "w");
    gnupipe = popen("gnuplot -persitent", "w");
//  std::cout << N << std::endl;
//  std::cout << N << std::endl;

    //std::cout << Par[1].P() << " | " << Par[1].V() << std::endl;
    //float lx = 0, ly = 0, lz = 0, lvx = 0, lvy, lvz = 0;
    if (thing == 103)
    {

        int T = 0;
        double time_step = 0;
        std::cout << "Time of measurement(days): " << std::endl;
        std::cin >> T;
        std::cout << "Time step(seconds): " << std::endl;
        std::cin >> time_step;

        Gravity_Sim grav_force = Gravity_Sim(0, 0, 0);

        std::ofstream outFile;
        char filename[] = "0planet.txt";

        for (int i = 0; i < N - 1; i++)
        {
            fprintf(fp, "%lf %lf %lf %d\n", Par[i].P().x(), Par[i].P().y(), Par[i].P().z(), i);
        }

        // fprintf(fp, "%lf %lf\n", 2e+12, 0);

        // std::cout << Par[0].P() << " | " << Par[0].V() << std::endl;
        // std::cout << Par[1].P() << " | " << Par[1].V() << std::endl;
        // std::cout << "--------------------" << std::endl;

        Vector *forceM = (Vector *)malloc(sizeof(Vector) * (pow(N - 1, 2)));
        for (int j = 0; j < (N - 1) * (N - 1); j++)
        {
            forceM[j] = Vector(0, 0, 0);
        }
        int time = 0;
        while (time < T * 86400 / time_step)
        {
            // T*86400
            // std::cout << "--------------------" << std::endl;
            for (int j = 0; j < (N - 1) * (N - 1); j++)
            {
                forceM[j].x(0);
                forceM[j].y(0);
                forceM[j].z(0);

            }
            for (int i = 0; i < N - 1; i++)
            {
                force_net[i].x(0);
                force_net[i].y(0);
                force_net[i].z(0);
            }

        for (int i = 0; i < N - 1; i++)
        {
          // std::cout << Par[i].P() << " | " << Par[i].V() << std::endl;
            for (int j = 0; j < N - 1; j++)
            {
                //std::cout << i << ", " << j << std::endl;

                if(j>i){
                    forceM[i * (N - 1) + j] = grav_force.force(Par[i], Par[j]);
                    forceM[j * (N - 1) + i].x(-1 * (forceM[i * (N - 1) + j].x()));
                    forceM[j * (N - 1) + i].y(-1 * (forceM[i * (N - 1) + j].y()));
                    forceM[j * (N - 1) + i].z(-1 * (forceM[i * (N - 1) + j].z()));
                }
            }
        }

        // for (int i = 0; i < N - 1; i++)
        // {
        //     for (int j = 0; j < N - 1; j++)
        //     {
        //         std::cout << "( " << forceM[i * (N - 1) + j].x() << ", " << forceM[i * (N - 1) + j].y()
        //                   << ", " << forceM[i * (N - 1) + j].z() << ")"
        //                   << " | ";
        //     }
        //     std::cout << "\n";
        // }

        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
               //if(i!=0){

                    force_net[i] = force_net[i] + forceM[i * (N - 1) + j];
                //}
            }
            // std::cout << "( " << force_net[i].x() << ", " << force_net[i].y()
            //               << ", " << force_net[i].z() << ")" << std::endl;
        }


        for (int i = 0; i < N - 1; i++){
            ph = grav_force.run(Par[i], force_net[i], time_step);
            // std::cout << "( " << force_net[i].x() << ", " << force_net[i].y()
            //               << ", " << force_net[i].z() << ")" << std::endl;
            Par[i].input(ph.M(), ph.P().x(), ph.P().y(), ph.P().z(),
                         ph.V().x(), ph.V().y(), ph.V().z(), ph.Q());
            // std::cout << Par[i].P() << " | " << Par[i].V() << std::endl;

            filename[0] = i + 1;
            outFile.open(filename, std::ios::app);
            outFile << time << " | ";
            outFile << Par[i].M() << ", " << Par[i].P().x() << ", " << Par[i].P().y() << ", ";
            outFile << Par[i].P().z() << ", " << Par[i].V().x() << ", " << Par[i].V().y() << ", ";
            outFile << Par[i].V().z() <<"\n "<<std::endl;
            outFile.close();

            fprintf(fp, "%lf %lf %lf %d\n", Par[i].P().x(), Par[i].P().y(), Par[i].P().z(), i);
        }

        //std::cout << time << std::endl;
        time++;
    }
        for (int i = 0; i < N - 1; i++)
        {
            std::cout << Par[i].P() - Par[0].P() << " | " << Par[i].V() - Par[0].V() << std::endl;
        }
    }

    if( thing == 101){

    int T = 0;
    std::cout << "Time of measurement(10^(-20)s): " << std::endl;
    std::cin >> T;

    Electro_Sim electro_force = Electro_Sim(0, 0, 0);

    std::ofstream outFile;
    char filename[] = "0electro.txt";

    for (int i = 0; i < N - 1; i++){
        fprintf(fp, "%lf %lf %lf %d\n", Par[i].P().x(), Par[i].P().y(), Par[i].P().z(), i);
    }

    //fprintf(fp, "%lf %lf\n", 2e+12, 0);

    // std::cout << Par[0].P() << " | " << Par[0].V() << std::endl;
    // std::cout << Par[1].P() << " | " << Par[1].V() << std::endl;
   // std::cout << "--------------------" << std::endl;
    Vector *forceM = (Vector *)malloc(sizeof(Vector) * (pow(N - 1, 2)));
    for (int j = 0; j < (N - 1) * (N - 1); j++)
            {
                forceM[j] = Vector(0, 0, 0);
            }
    int time = 0;
    while(time<T*10e4){
      //std::cout << "--------------------" << std::endl;
         for (int j = 0; j < (N - 1) * (N - 1); j++)
            {
                forceM[j].x(0);
                forceM[j].y(0);
                forceM[j].z(0);

            }
        for (int i = 0; i < N - 1; i++)
            {
                force_net[i].x(0);
                force_net[i].y(0);
                force_net[i].z(0);
            }

        for (int i = 0; i < N - 1; i++)
        {
          // std::cout << Par[i].P() << " | " << Par[i].V() << std::endl;
            for (int j = 0; j < N - 1; j++)
            {
                //std::cout << i << ", " << j << std::endl;

                if(j>i){
                forceM[i * (N-1) + j] = electro_force.force(Par[i], Par[j]);
                forceM[j * (N - 1) + i].x(-1 * (forceM[i * (N - 1) + j].x()));
                forceM[j * (N - 1) + i].y(-1 * (forceM[i * (N - 1) + j].y()));
                forceM[j * (N - 1) + i].z(-1 * (forceM[i * (N - 1) + j].z()));
                }
            }
        }

        // for (int i = 0; i < N - 1; i++)
        // {
        //     for (int j = 0; j < N - 1; j++)
        //     {
        //         std::cout << "( " << forceM[i * (N - 1) + j].x() << ", " << forceM[i * (N - 1) + j].y()
        //                   << ", " << forceM[i * (N - 1) + j].z() << ")"
        //                   << " | ";
        //     }
        //     std::cout << "\n";
        // }

        for (int i = 0; i < N - 1; i++)
        {
            for (int j = 0; j < N - 1; j++)
            {
               //if(i!=0){

                    force_net[i] = force_net[i] + forceM[i * (N - 1) + j];
                //}
            }
            // std::cout << "( " << force_net[i].x() << ", " << force_net[i].y()
            //               << ", " << force_net[i].z() << ")" << std::endl;
        }


        for (int i = 0; i < N - 1; i++){
            ph = electro_force.run(Par[i], force_net[i], 10e-16);
            
            // std::cout << "( " << force_net[i].x() << ", " << force_net[i].y()
            //               << ", " << force_net[i].z() << ")" << std::endl;
            Par[i].input(ph.M(), ph.P().x(), ph.P().y(), ph.P().z(), ph.V().x(),
                         ph.V().y(), ph.V().z(), ph.Q());
            // std::cout << Par[i].P() << " | " << Par[i].V() << std::endl;

            filename[0] = i + 1;
            outFile.open(filename, std::ios::app);
            outFile << time << " | ";
            outFile << Par[i].M() << ", " << Par[i].P().x() << ", " << Par[i].P().y() << ", ";
            outFile << Par[i].P().z() << ", " << Par[i].V().x() << ", " << Par[i].V().y() << ", ";
            outFile << Par[i].V().z() << "\n " << std::endl;
            outFile.close();

            fprintf(fp, "%lf %lf %lf %d\n", Par[i].P().x(), Par[i].P().y(), Par[i].P().z(), i);

            // std::cout << time << std::endl;
            time++;
    }
    }

    //std::cout << Par[0].P() << " | " << Par[0].V() << std::endl;`
   // std::cout << Par[1].P() << " | " << Par[1].V() << std::endl;

    // for (int i = 0; i < 2; i++){
    //     fprintf(gnupipe, "%s\n", GnuCommands[i]);
    // }

    // fprintf(gnupipe, "set xlabel 'X'\n");
    // fprintf(gnupipe, "set ylabel 'Y'\n");
    // fflush(gnupipe);
    // getchar();
    }

    std::cout << "DONE" << std::endl;
    return 0;
}

int main(){
    Exec_Time(&start);
    return 0;
}
