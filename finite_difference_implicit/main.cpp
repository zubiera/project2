//
//  main.cpp
//  Finite_difference_implict for the heston model PDE
// the crank nicolson one is saved as Finite_differenc
//  Created by amit arfan on 27/08/2019
//  Copyright © 2019 amit arfan. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
using namespace std;

// A generic lagrange interpolation
double lagrangeInterpolation(const vector<double>& y, const vector<double>& x,double x0, unsigned int n )
{
    if(x.size()<n) return lagrangeInterpolation(y, x, x0, x.size());
    if(n==0) throw;
    int nHalf = n/2;
    int jStar;
    double dx = x[1]-x[0];
    if(n%2 == 0)
        jStar = int((x0-x[0])/dx) - (nHalf-1);
    else
        jStar = int((x0-x[0])/dx+0.5) - (nHalf);
    jStar = std::max(0,jStar);
    jStar = std::min(int(x.size()-n),jStar);
    if (n==1) return y[jStar];
    double temp = 0;
    for(unsigned int i = jStar; i < jStar+n;i++)
    {
        double int_temp;
        int_temp = y[i];
        for(unsigned int j=jStar;j<jStar+n;j++)
        {
            if(j==i){continue;}         //This is part fo lagrange formula for i != j
            int_temp *= (x0-x[j]) / (x[i]-x[j]); // lagrange formula
        }
        temp += int_temp;
    }
    
    return temp;
}


//LU deccomposition function Thomas algorithm
vector<double> tridag(const vector<double>& a, const vector<double>& beta, const vector<double>& c, const vector<double>& d )
{
    int n = a.size();
    vector<double> b(n);
    //move d to rhs
    vector<double> rhs(d);
    //initilise first value of b
    b[0] = beta[0];
    
    for (int j=1; j<n;j++)
    {
        b[j]   = beta[j] - c[j-1]*a[j]/b[j-1];
        rhs[j] = rhs[j]  - rhs[j-1]*a[j]/b[j-1];
    }
    //calc solution
    rhs[n-1]  = rhs[n-1]/b[n-1];
    for(int j=n-2;j>=0;j--)
    {
        rhs[j] = (rhs[j] - c[j]*rhs[j+1])/b[j];
    }
    return rhs;
}

int calculateExpectedUtility(
                             double T,  //final time t
                             double mu, //asset returns
                             double sigma, //asset voaltility
                             double gamma, // tisk aversion
                             double A,   //volume of sales param
                             double k,   // market depth param
                             int qMax,  // number of assets
                             int n,   //grid size in fund
                             vector<double>& S,
                             vector<vector<double>>& deltat,
                             vector<vector<double>>& Ut)
{
    int iMax = T*n;
    int jMax  = n;
    //local parameter setup
    double dt = (T)/iMax;
    double ds = (S[jMax] - S[0])/jMax;
    //S.resize(jMax+1);
    //deltat.clear();
    //Ut.clear();
    //deltat.resize(qMax+1,vector<double>(jMax+1));
    //Ut.resize(qMax+1,vector<double>(jMax+1));
    
    
    //star looping through time levels now (from terminal)
    for (int i=iMax-1;i>=0;i--)
    {
        //store the values at t+dt as Old known values
        vector<vector<double>> UtOld(Ut);
        //lopp through q values
        for(int q=1;q<=qMax;q++)
        {
            //firt we calculate the the optimal delta values using values at t+dt (UtOld)
            //and store it in a vector
            vector<double> SupTerm(jMax+1);
            for (int j=0;j<=jMax;j++)
            {
                //optimal choice of delta when we solved sup[Ae^(-k delta)( e^{-gamma(s+delta)}U(t+dt,q-1,S) - U(t+dt,q,s)) ]
                deltat[q][j] = 1./gamma * log( ( (k+gamma)*UtOld[q-1][j] ) /(k*UtOld[q][j]) ) - S[j];
                //then sub back in to find supTerm value
                SupTerm[j] = A*exp(-k*deltat[q][j])*( exp(-gamma*(S[j]+deltat[q][j])) *UtOld[q-1][j] - UtOld[q][j] );
            }
            //supTerm is calcualted and known so we can put onto rhs ie we include it in the d term used in the thomas algorithm for LU decomp of tridiagonal
            //vecotrs for matrix exuations
            vector<double> a(jMax+1);
            vector<double> b(jMax+1);
            vector<double> c(jMax+1);
            vector<double> d(jMax+1);
            
            //boundary conditions assumes U = e^{-gamma*q*s}f
            //U_s=-gamma*q*U  U_ss = gamma^2q^2U sub these back into the pde to get
            a[0] = 0;
            b[0] = -1./dt + (0.5*sigma*sigma*gamma*gamma*q*q - mu*gamma*q);
            c[0] = 0;
            d[0] = -UtOld[q][0]/dt - SupTerm[0];
            
            for(int j =1;j<=jMax;j++)
            {
                a[j] = 0.5*(sigma*sigma/ds/ds -mu/ds);
                b[j] = -sigma*sigma/ds/ds - 1./dt;
                c[j] = 0.5*(sigma*sigma/ds/ds + mu/ds);
                d[j] = -1./dt*UtOld[q][j] -SupTerm[j];
            }
            
            //same as lower boundary
            a[jMax] = 0;
            b[jMax] = -1./dt+(0.5*sigma*sigma*gamma*gamma*q*q - mu*gamma*q);
            c[jMax] = 0;
            d[jMax] = -UtOld[q][jMax]/dt - SupTerm[jMax];
            
            //solve amtrix equation
            Ut[q] = tridag(a,b,c,d);
        }
    }
    
    return 0;
}

int Solve_stochastic_vol(double T, double mu,
                         double theta, //reversion rate
                         double alpha, //lng term mean of variance
                         double xi, // volatility of vvariance
                         double rho, //correlation
                         double gamma, double A, double k,
                         int qMax, // the most posive we want to be in assets
                         int n, vector<double>& Nu,
                         vector<vector<double>>& deltat_a,  //store the ask side delta
                         vector<vector<double>>& deltat_b, //store the bid side delta
                         vector<vector<double>>& Ut)
{
    int iMax = max(1,int(T*n)/10);  //time steps per time interval
    int jMax = n;
    //local parameter setup
    double dt = (T)/iMax;
    double dnu = (Nu[jMax]-Nu[0])/jMax;
    //    Nu.resize(jMax+1);
    //    deltat_a.clear();
    //    deltat_b.clear();
    //    Ut.clear();
    //    deltat_a.resize(2*qMax+1, vector<double>(jMax+1));
    //    deltat_b.resize(2*qMax+1, vector<double>(jMax+1));
    //    Ut.resize(2*qMax+1,vector<double>(jMax+1));
    
    
    //    //setup and initilise stockprice
    //    for(int j=0;j<=jMax;j++)
    //    {
    //        Nu[j] = nuMin + j*dnu;
    //    }
    //
    //
    //    //initilise terminal conditoon
    //    //here we take into account the absolute inventort not jsut postitive
    //    for (int q=-qMax; q<=qMax;q++)
    //    {
    //        for (int j=0;j<=jMax;j++)
    //        {
    //            Ut[q+qMax][j]  =  1;    // assumng l=0 or -exp(-gamma*(l(|q|))) i.e -exp(-gamma(q*q))
    //        }
    //    }
    
    //start looping through time levels now (from terminal)
    for (int i=iMax-1;i>=0;i--)
    {
        //store the values at t+dt as Old known values
        vector<vector<double>> UtOld(Ut);
        
        //creat a storage matrix for the sup terms instead of a vector.
        //vector<vector<double>> SupTerm_bid(2*qMax+1, vector<double>(n+1));
        //vector<vector<double>> SupTerm_ask(2*qMax+1, vector<double>(n+1));
        
        
        //******for -Q  ******************
        {
            //firt we calculate the the optimal delta values using values at t+dt (UtOld)
            //and store it in a vector
            vector<double> SupTerm_bid(jMax+1);
            for (int j=0;j<=jMax;j++)
            {
                
                //optimal choice of delta^a and delta^b
                deltat_b[0][j] = (1./gamma) * (log(1 + gamma/k) + log(UtOld[1][j]/UtOld[0][j]));
                //then sub back in to find supTerm value
                SupTerm_bid[j] = A*exp(-k*deltat_b[0][j])*( exp(-gamma*deltat_b[0][j]) *UtOld[1][j] - UtOld[0][j] )  ;
            }
            
            
            //supTerm is calcualted and known so we can put onto rhs ie we include it in the d term used in the thomas algorithm for LU decomp of tridiagonal
            //vecotrs for matrix exuations
            vector<double> a(jMax+1);
            vector<double> b(jMax+1);
            vector<double> c(jMax+1);
            vector<double> d(jMax+1);
            
            //boundary conditions assumes at nu=0 (j=0), one sided (forward) implicit nicolson ..........
            //
            a[0] = 0;
            b[0] = -1./dt - 0.5*theta*alpha/dnu;
            c[0] = 0.5*theta*alpha/dnu;
            d[0] =  - UtOld[0][0]/dt - SupTerm_bid[0];
            
            //Implicit nicolson scheme with upwond and downwid takin into accont
            for(int j =1;j<=jMax;j++)
            {
//                double a_coeff = 0.5*(theta*(alpha - Nu[j]) - rho*xi*gamma*(-qMax)*(-qMax)*Nu[j])/(dnu) ;
                double a_coeff = 0.5*(theta*(alpha - Nu[j]) - rho*xi*gamma*(-qMax)*Nu[j])/(dnu) ;
                double b_coeff = 0.5* xi*xi*Nu[j]/dnu/dnu;
                a[j] = b_coeff;
                b[j] = -2.*b_coeff;
                c[j] = b_coeff;
                
                if(a_coeff > b_coeff)
                {
                    //Change to forward difference scheme
                    b[j] = b[j] - 2.*a_coeff;
                    c[j] = c[j] + 2.*a_coeff;
                }
                else if(a_coeff< -b_coeff)
                {
                    //Change to backwards difference scheme
                    a[j] = a[j] - 2.*a_coeff;
                    b[j] = b[j] + 2.*a_coeff;
                }
                else
                {
                    //Keep with central dfference scheme
                    a[j] = a[j] - a_coeff;
                    c[j] = c[j] + a_coeff;
                }
                //add back the time difference portion and the U_i^j term to b[j]
                b[j] = b[j] - 1./dt + Nu[j]*gamma*gamma*(-qMax)*(-qMax)/2.;
                
                d[j] =  UtOld[0][j]*(  -1/dt  ) - SupTerm_bid[j];

            }
            
            // boundary conditions approximating nu->infinity using crank nicolson aprox
            //here we assume that U_nu = 0
            a[jMax] = -1.;
            b[jMax] = 1.;
            c[jMax] = 0.;
            d[jMax] = 0;
            
            //solve Matrix equation
            Ut[0] = tridag(a,b,c,d);
        }
        
        //lopp through q values
        for(int q=-qMax+1; q<=qMax-1; q++)  //cannot includde -Q or Q
        {
            //firt we calculate the the optimal delta values using values at t+dt (UtOld)
            //and store it in a vector
            vector<double> SupTerm_ask(jMax+1);
            vector<double> SupTerm_bid(jMax+1);
            for (int j=0;j<=jMax;j++)
            {
                
                //optimal choice of delta^a and delta^b
                deltat_a[q+qMax][j] = (1./gamma) * (log(1. + gamma/k) + log(UtOld[q-1 +qMax][j]/UtOld[q + qMax][j]));
                deltat_b[q+qMax][j] = (1./gamma) * (log(1. + gamma/k) + log(UtOld[q+1 +qMax][j]/UtOld[q +qMax][j]));
                //then sub back in to find supTerm value
                SupTerm_ask[j] = A*exp(-k*deltat_a[q +qMax][j])*( exp(-gamma*deltat_a[q+qMax][j]) *UtOld[q-1+qMax][j] - UtOld[q+qMax][j] );
                SupTerm_bid[j] = A*exp(-k*deltat_b[q+qMax][j])*( exp(-gamma*deltat_b[q+qMax][j]) *UtOld[q+1+qMax][j] - UtOld[q+qMax][j] )  ;
                
            }
            
            
            //supTerm is calcualted and known so we can put onto rhs ie we include it in the d term used in the thomas algorithm for LU decomp of tridiagonal
            //vecotrs for matrix exuations
            vector<double> a(jMax+1);
            vector<double> b(jMax+1);
            vector<double> c(jMax+1);
            vector<double> d(jMax+1);
            
            //boundary conditions assumes at nu=o, one sided (forward) implicit nicolson ..........
            //
            a[0] = 0;
            b[0] = -1./dt - 0.5*theta*alpha/dnu;
            c[0] = 0.5*theta*alpha/dnu;
            d[0] =  - UtOld[q+qMax][0]/dt  - SupTerm_ask[0] - SupTerm_bid[0];
            
            //Implicit nicolson scheme with upwond and downwid takin into accont
            for(int j =1;j<=jMax;j++)
            {
                
                //double a_coeff = 0.5*(theta*(alpha - Nu[j]) - rho*xi*gamma*q*q*Nu[j])/(dnu) ;
                double a_coeff = 0.5*(theta*(alpha - Nu[j]) - rho*xi*gamma*q*Nu[j])/(dnu) ;
                double b_coeff = 0.5* xi*xi*Nu[j]/dnu/dnu;
                a[j] = b_coeff;
                b[j] = -2.*b_coeff;
                c[j] = b_coeff;
                
                if(a_coeff > b_coeff)
                {
                    //Change to forward difference scheme
                    b[j] = b[j] - 2.*a_coeff;
                    c[j] = c[j] + 2.*a_coeff;
                }
                else if(a_coeff<-b_coeff)
                {
                    //Change to backwards difference scheme
                    a[j] = a[j] - 2.*a_coeff;
                    b[j] = b[j] + 2.*a_coeff;
                }
                else
                {
                    //Keep with central dfference scheme
                    a[j] = a[j] - a_coeff;
                    c[j] = c[j] + a_coeff;
                }
                //add back the time difference portion and the U_i^j term to b[j]
                b[j] = b[j] - 1./dt + Nu[j]*gamma*gamma*q*q/2.;
                
                d[j] =  UtOld[q+qMax][j]*(  -1/dt  )  - SupTerm_ask[j] - SupTerm_bid[j];
                
            }
            
            // boundary conditions approximating nu->infinity using crank nicolson aprox
            //here we assume that U_nu = 0
            a[jMax] = -1.;
            b[jMax] = 1.;
            c[jMax] = 0.;
            d[jMax] = 0.;
            
            //solve Matrix equation
            Ut[q+qMax] = tridag(a,b,c,d);
        }
        //******for Q  ******************
        {
            //firt we calculate the the optimal delta values using values at t+dt (UtOld)
            //and store it in a vector
            vector<double> SupTerm_ask(jMax+1);
            for (int j=0;j<=jMax;j++)
            {
                
                //optimal choice of delta^a and delta^b
                deltat_a[2*qMax][j] = (1./gamma) * (log(1 + gamma/k) + log(UtOld[2*qMax-1][j]/UtOld[2*qMax][j]));
                //then sub back in to find supTerm value
                SupTerm_ask[j] = A*exp(-k*deltat_a[2*qMax][j])*( exp(-gamma*deltat_a[2*qMax][j]) *UtOld[2*qMax-1][j] - UtOld[2*qMax][j] )  ;
            }
            
            
            //supTerm is calcualted and known so we can put onto rhs ie we include it in the d term used in the thomas algorithm for LU decomp of tridiagonal
            //vecotrs for matrix exuations
            vector<double> a(jMax+1);
            vector<double> b(jMax+1);
            vector<double> c(jMax+1);
            vector<double> d(jMax+1);
            
            //boundary conditions assumes at nu=0 (j=0), one sided crank nicolson ..........
            //
            a[0] = 0;
            b[0] = -1./dt - 0.5*theta*alpha/dnu;
            c[0] = 0.5*theta*alpha/dnu;
            d[0] =  - UtOld[2*qMax][0]/dt - SupTerm_ask[0];
            
            //Implicit nicolson scheme with upwond and downwid takin into accont
            for(int j = 1;j<=jMax;j++)
        {
                
//                double a_coeff = 0.5*(theta*(alpha - Nu[j]) - rho*xi*gamma*(qMax)*(qMax)*Nu[j])/(dnu) ;
            
                double a_coeff = 0.5*(theta*(alpha - Nu[j]) - rho*xi*gamma*(qMax)*Nu[j])/(dnu) ;
                double b_coeff = 0.5* xi*xi*Nu[j]/dnu/dnu;
                a[j] = b_coeff;
                b[j] = -2.*b_coeff;
                c[j] = b_coeff;
                
                if(a_coeff > b_coeff)
                {
                    //Change to forward difference scheme
                    b[j] = b[j] - 2.*a_coeff;
                    c[j] = c[j] + 2.*a_coeff;
                }
                else if(a_coeff<-b_coeff)
                {
                    //Change to backwards difference scheme
                    a[j] = a[j] - 2.*a_coeff;
                    b[j] = b[j] + 2.*a_coeff;
                }
                else
                {
                    //Keep with central dfference scheme
                    a[j] = a[j] - a_coeff;
                    c[j] = c[j] + a_coeff;
                }
                //add back the time difference portion and the U_i^j term to b[j]
                b[j] = b[j] - 1./dt + Nu[j]*gamma*gamma*(qMax)*(qMax)/2.;
                
                d[j] =  UtOld[2*qMax][j]*(  -1/dt  )  -SupTerm_ask[j]  - SupTerm_ask[j];
                
            }
            
            // boundary conditions approximating nu->infinity using crank nicolson aprox
            //here we assume that U_nu = 0
            a[jMax] = -1.;
            b[jMax] = 1.;
            c[jMax] = 0.;
            d[jMax] = 0.;
            
            //solve Matrix equation
            Ut[2*qMax] = tridag(a,b,c,d);
        }
        
    }
    
    return 0;
}




int main(int argc, char *argv[]) {

    
    //   //Parameters for the stochastic vol PDE
    double nuMin    =   0.;
    double nuMax    =   std::stod(argv[9]);//12000;//1.;
    double mu       =   0.;
    double theta    =   std::stod(argv[1]); //0.02;
    double alpha    =   std::stod(argv[2]); //0.3;    //lng term mean of variance
    double xi       =   std::stod(argv[3]); //0.01;  // volatility of variance
    double rho      =   std::stod(argv[4]); // 0.7; //correlation
    double gamma    =   std::stod(argv[5]);  //0.05;
    double A        =   std::stod(argv[6]); //.1;
    double k        =   std::stod(argv[7]);  //0.3;
    int qMax        =   10;    // the most posive we want to be in assets
    double T        =   std::stod(argv[8]);//1380; //300.;
    int N           =   500;
    vector<double> Nu(N+1);
    vector<vector<double>> deltat_a(2*qMax+1, vector<double>(N+1));  //store the ask side
    vector<vector<double>> deltat_b(2*qMax+1, vector<double>(N+1)); //bid side delta
    vector<vector<double>> Ut(2*qMax+1, vector<double>(N+1));
    
    
    int TimePeriods = 1380; //300;
    double DeltaT   = T/TimePeriods;
    double dnu      = (nuMax-nuMin)/N;
    
    //create the files here to store delta_a
    std::string path_b = "/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_b_implicit_Qchange.txt";
    std::string path_a = "/Users/amitarfan1/Documents/Phd/3yr/Python_Code/heston_solns/delta_a_implicit_Qchange.txt";
    ofstream output_b(path_b);
    ofstream output_a(path_a);
    
    //setup and initlise the nu's
    for(int j=0; j<=N; j++)
    {
        Nu[j] = nuMin + j*dnu;
    }
    // setup and initlise the final conditions on U
    //here we assumes that U(T,nu,q) = -exp(-l(|q|)) where l=0
    for (int q=-qMax; q<=qMax;q++)
    {
        for (int j=0;j<=N;j++)
        {
            Ut[q+qMax][j]  =  1.;    // assumng l=0 or -exp(-gamma*(l(|q|))) i.e -exp(-gamma(q*q))
            // calculate delta at T
        }
    }
    
    // output terminal stuff here for t=T
    for(int j=0; j<=N; j++)
    {
        deltat_b[0][j] = (1./gamma) * (log(1 + gamma/k) + log(Ut[1][j]/Ut[0][j]));
        output_b    << deltat_b[0][j] << ",";
        output_a    << 0. << ",";
        for(int q=-qMax+1; q<=qMax-1; q++)
        {
            deltat_a[q+qMax][j] = (1./gamma) * (log(1. + gamma/k) + log(Ut[q-1 +qMax][j]/Ut[q + qMax][j]));
            deltat_b[q+qMax][j] = (1./gamma) * (log(1. + gamma/k) + log(Ut[q+1 +qMax][j]/Ut[q +qMax][j]));
            output_b << deltat_b[q+qMax][j] << ",";
            output_a << deltat_a[q+qMax][j] << ",";
        }
        deltat_a[2*qMax][j] = (1./gamma) * (log(1 + gamma/k) + log(Ut[2*qMax-1][j]/Ut[2*qMax][j]));
        output_b << 0. << ",";
        output_a << deltat_a[2*qMax][j] << ",";
    }
    
    // Here we output delta_b and delta_a as cout or to a file.
    for (int T_k =0; T_k<=TimePeriods; T_k++)
    {
        Solve_stochastic_vol(DeltaT, mu, theta, alpha, xi, rho, gamma, A, k, qMax, N, Nu, deltat_a, deltat_b, Ut);
        
        for(int j=0; j<=N; j++)
        {
            //cout << (TimePeriods-T_k)*DeltaT << "," << Nu[j] << ",";
            for(int q=-qMax; q<=qMax; q++)
            {
                //double U_Nu   =   (Ut[q+qMax][j+1] - Ut[q+qMax][j-1])/(2*(1/N));
                //double U_NuNu =   (Ut[q+qMax][j+1] - 2.*Ut[q+qMax][j] + Ut[q+qMax][j-1]) / (1./N)/(1./N) ;
                output_b << deltat_b[q+qMax][j] << ",";
                output_a << deltat_a[q+qMax][j] << ",";
                
                //cout  << deltat_b[q+qMax][j] << ","; //remove spaces here for pyhton
                
            }
            //cout << endl;
        }
    }
    output_b << endl;
    output_a << endl;
    //cout << endl; //Only keep this endl when I want to output a long string.
    
    
    
    return 0;
}

