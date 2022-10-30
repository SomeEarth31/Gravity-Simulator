#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

using std::string;

class Vector
{
private:
    double _x;
    double _y;
    double _z;

public:
    Vector(double x, double y, double z)
    {
        _x=(x);
        _y=(y);
        _z=(z);
    }
    void x(double x) 
    {
        // manipulator functions(to change specific values)
        _x = x;
    }
    void y(double y)
    {
        _y = y;
    }
    void z(double z)
    {
        _z = z;
    }
    double x(void) const 
    {
        // accessor functions
        return _x;
    }
    double y(void) const
    {
        return _y;
    }
    double z(void) const
    {
        return _z;
    }
    Vector operator+(const Vector &obj)
    { 
        // overloading the + operator
        Vector temp = Vector(0, 0, 0);
        temp._x = _x + obj._x;
        temp._y = _y + obj._y;
        temp._z = _z + obj._z;
        return temp;
    }
    Vector operator-(const Vector &obj)
    { 
        // overloading the - operator
        Vector temp = Vector(0, 0, 0);
        temp._x = _x - obj._x;
        temp._y = _y - obj._y;
        temp._z = _z - obj._z;
        return temp;
    }

    friend std::ostream &operator<<(std::ostream &os, const Vector &obj)
    { 
        // overloading the << operator
        os << "(" << obj._x << ", " << obj._y << ", " << obj._z << ")";
        return os;
    }
    double dot(Vector a, Vector b);
    Vector cross(Vector a, Vector b);
   // Vector *Getgrav_force(int a, int b);
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

class Particle
{
private:
    double _mass;
    double _charge;
    Vector Position = Vector(0, 0, 0);
    Vector Velocity = Vector(0, 0, 0);

public:
    Particle(double mass, double x, double y, double z, double vx,
             double vy, double vz, double charge)
    {
        _mass = mass;
        Position.x(x), Position.y(y), Position.z(z);
        Velocity.x(vx), Velocity.y(vy), Velocity.z(vz);
        _charge = charge;
    }
    friend std::ostream &operator<<(std::ostream &os, const Particle &obj)
    { 
        // overloading the << operator
        os << "_mass = " << obj._mass << "\nPosition is " << obj.Position << "\nVelocity is " << obj.Velocity << std::endl;
        return os;
    }
    void input(double mass, double x, double y, double z, double vx,
               double vy, double vz, double charge)
    {
        _mass = mass;
        Position.x(x), Position.y(y), Position.z(z);
        Velocity.x(vx), Velocity.y(vy), Velocity.z(vz);
        _charge = charge;
    }
    Vector P(void) const 
    {
        // accesors
        return Position;
    }
    Vector V(void) const
    {
        return Velocity;
    }
    double M(void) const
    {
        return _mass;
    }
    double Q(void) const
    {
        return _charge;
    }
    void M(double mass)
    {
        _mass = mass;
    }
    void P(double x, double y, double z)
    {
        Position = Vector(x, y, z);
    }
    void V(double x, double y, double z)
    {
        Velocity = Vector(x, y, z);
    }
    void Q(double charge)
    {
        _charge = charge;
    }
    // Particle operator+(const Particle &obj)
    // {
    //     Particle temp = Particle(0, 0, 0, 0, 0, 0, 0, 0);
    //     // temp.M() = temp.M() + obj.M();
    //     temp.V() = temp.V() + obj.V();
    //     temp.P() = temp.P() + obj.P();
    //     return temp;
    // }
    // Particle operator-(const Particle &obj)
    // {
    //     Particle temp = Particle(0, 0, 0, 0, 0, 0, 0, 0);
    //     // temp.M() = temp.M() - obj.M();
    //     temp.V() = temp.V() - obj.V();
    //     temp.P() = temp.P() - obj.P();
    //     return temp;
    // }
    // Particle operator=(const Particle &obj)
    // {
    //     Particle temp = Particle(obj.M(), obj.P().x(), obj.P().y(), obj.P().z(), obj.V().x(), obj.V().y(), obj.V().z(), 0);
    //     return temp;
    // }
};

class Simulation_Engine 
{
    // Force
private:
    Vector Force = Vector(0, 0, 0);

public:
    Simulation_Engine(double fx, double fy, double fz)
    {
        Force.x(fx), Force.y(fy), Force.z(fz);
    }
    friend std::ostream &operator<<(std::ostream &os,
                                    const Simulation_Engine &obj)
    {
        os << obj.Force << std::endl;
        return os;
    }
    Vector F(void)
    {
        return Force;
    }
    void F(double fx, double fy, double fz)
    {
        Force = Vector(fx, fy, fz);
    }
    virtual Particle run(Particle a, Vector F_net, double t) = 0;
    virtual Vector force(Particle a, Particle b) = 0;
};

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


int main(void)
{
    time_t start, end;
    time(&start);
    std::ios_base::sync_with_stdio(false);
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
    if( thing == 103){

    int T = 0;
    std::cout << "Time of measurement(days): " << std::endl;
    std::cin >> T;

    Gravity_Sim grav_force = Gravity_Sim(0, 0, 0);

    std::ofstream outFile;
    char filename[] = "0planet.txt";

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
    while(time<T*86){ 
        //T*86400
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
                forceM[i * (N-1) + j] = grav_force.force(Par[i], Par[j]);
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
            ph = grav_force.run(Par[i], force_net[i], 1004);
                                                   // ^100^
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
    time(&end);
    double time_taken = double(end - start);

    std::cout << "Time taken by simulation : " << std::fixed
              << time_taken << " sec " << std::endl;
    // std::cout << " sec " << std::endl;
    return 0;
}
