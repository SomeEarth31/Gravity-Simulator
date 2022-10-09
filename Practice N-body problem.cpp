#include<iostream>
#include<cmath>
using namespace std;

//given n bodies(mass, position in 3d, and velocity components) this finds their positions after time t. All in SI units.

float G=6.67*pow(10, -11);

struct prop{
    float mass;
    float x[3];
    float V[3];
};

struct acc{
    float ax, ay, az;
};

float distance(float x1, float x2, float y1, float y2, float z1, float z2){
  float d;
  d=sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2)) ;
    return d;
}

float Force_Coeff(float m1, float m2, float d){
  float f;
  f=(G*m1*m2)/pow(d, 2);
    return f; 
}

float Force_s(float x, float y, float z){
    float Fx;
    Fx=x/sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    return Fx;
}

float dist_s(float x, float a, float V, int t ){
  float Sx;
  Sx=x+V*t+(a*t*t)/2;
  return Sx;
}

float vel_s(float u, float a, int t){
    float v;
    v=u+a*t;
    return v;
}


int main(){
 float T;
 int N, t; //Get N from file hopefully
 struct prop values[N], place[N];
cout << "What is the time: ";
cin >> T;
t=T/1;
 for(int i=0; i<N; i++){

    cin>>values[i].mass;
    place[i].mass=values[i].mass;

      for(int j=0; j<3; j++){
        cin>>values[i].x[j];
        place[i].x[j]=values[i].x[j];
      }
      for(int j=0; j<3; j++){
        cin>>values[i].V[j];
        place[i].V[j]=values[i].V[j];
      }
 }
  
float r=T-t;
int i=0, l=0, k=0;
struct acc A[N][N];

 for(i=0; i<t+1; i++){
    for(l=0; l<N; l++){
        for(k=0; k<N; k++) {
         if(k<l){ continue; }
             float d = distance(place[k].x[0], place[k].x[1], place[k].x[2], place[l].x[0], place[l].x[1], place[l].x[2]);
             float Fc=Force_Coeff(place[0].mass, place[1].mass, d);
             A[k][l].ax=Fc*Force_s(place[k].x[0],  place[k].x[1], place[k].x[2])/place[k].mass;
               A[l][k].ax=A[k][l].ax;
             A[k][l].ay=Fc*Force_s(place[k].x[1], place[k].x[0], place[k].x[2])/place[k].mass;
               A[l][k].ax=A[k][l].ay;
             A[k][l].az=Fc*Force_s(place[k].x[2], place[k].x[1], place[k].x[0])/place[k].mass; 
               A[l][k].ax=A[k][l].az; }}

    struct acc An[N];
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            An[i].ax=A[i][j].ax;
            An[i].ay=A[i][j].ay;
            An[i].az=A[i][j].az;
        }
    }
    
    for(int l=0; l<N; l++){
    float x, y, z, Vx, Vy, Vz;
    x=place[l].x[0];
    y=place[l].x[1];
    z=place[l].x[2];
    Vx=place[l].V[0];
    Vy=place[l].V[1];
    Vz=place[l].V[2];
    if(i<t){
                place[l].x[0]=dist_s(x, An[l].ax, Vx, 1);
                place[l].x[1]=dist_s(y, An[l].ay, Vy, 1);
                place[l].x[2]=dist_s(z, An[l].az, Vz, 1);

                place[l].V[0]=vel_s(Vx, An[l].ax, 1);
                place[l].V[1]=vel_s(Vy, An[l].ay, 1);
                place[l].V[2]=vel_s(Vz, An[l].az, 1); }
    else if(i==t){
                place[l].x[0]=dist_s(x, An[l].ax, Vx, r);
                place[l].x[1]=dist_s(y, An[l].ay, Vy, r);
                place[l].x[2]=dist_s(z, An[l].az, Vz, r);

                place[l].V[0]=vel_s(Vx, An[l].ax, r);
                place[l].V[1]=vel_s(Vy, An[l].ay, r);
                place[l].V[2]=vel_s(Vz, An[l].az, r);   }
    }
     cout<<sqrt(pow(place[0].V[0],2)+pow(place[0].V[1],2)+pow(place[0].V[2],2))<<", ("<<place[0].x[0]<<", "<<place[0].x[1]<<", "<<place[0].x[2]<<")\n";
    }
 

return 0;
}