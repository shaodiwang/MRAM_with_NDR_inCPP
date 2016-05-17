#include <iostream>
#include "NDR_Solver.cpp"
using namespace std;

int main(int argc, char* argv[])
{
    if(argc < 5){
        cout <<" no enough input, [VI of NDR] [VI of mos] [resistance of p-MTJ] [vdd]"<<endl;
        return 1;
    }
    double** VIGndr =  Read_voltage_current(argv[1]);
    double** VIGmos = Read_voltage_current(argv[2]);
    double Rmtj = atof(argv[3]);
    double vdd = atof(argv[4]);
    double Vndr = Solve_stable_vndr(VIGndr,VIGmos,0.75*Rmtj,vdd);
    cout<<"Vndr: "<<Vndr<<endl;
    
    //Starts to solve the transient equation
    //Parameters
    double Cload = 17e-15;
    double t_step =1e-12;
    double t_sim = 1e-9;
    //find the ic status
    double Imtj = IG_V(Vndr,VIGndr,1); //current through MTJ and nmos
    double Indr = Imtj; //current through NDR
    double Vmos = V_I(Imtj, VIGmos);
    double d_Rmtj = 0; // delta Rmtj
    double d_Imtj = 0; // delta Imtj
    double Rmos = 0; // dVm/dImtj
    // Solve the equation: Indr = Imtj + Cload*(Imtj*Rmtj + Vmos(Imtj))
    double t = 0;
    while(t < t_sim){
        t += t_step;
        if ( abs(Imtj) > 1e-6) d_Rmtj = - 5e11*t_step;
        d_Imtj = ( t_step*Indr - t_step*Imtj - Cload*Imtj*d_Rmtj) / ( Cload/IG_V(Vmos,VIGmos,2) + Cload*Rmtj);
        Rmtj += d_Rmtj;
        Imtj += d_Imtj;
        Vmos = V_I(Imtj, VIGmos);
        Vndr = vdd - Imtj*Rmtj - Vmos;
        Indr = IG_V(Vndr,VIGndr,1);
        std::cout<<"Vndr Indr Vmos Imtj Rmtj: "<<Vndr<<" "<<Indr<<" "<<Vmos<<" "<<Imtj<<" "<<Rmtj<<std::endl;
    }
    return 0;
}


