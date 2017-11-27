//
//  main.cpp
//  exam
//
//  Created by XingSong  Lin on 11/13/17.
//  Copyright © 2017 XingSong  Lin. All rights reserved.
//

#include <iostream>
#include <cmath>
using namespace std;

double cum_norm(double x)
{
    const double root = sqrt(0.5);
    return 0.5*(1.0 + erf(x*root));
}
int BS_call(const double &S,
            const double &K,
            const double &r,
            const double &q,
            const double &sigma,
            const double &T,
            const double &t,
            double &call,
            double &delta)
{
    call = 0;
    delta = 0;
    //const double pi = 4.0*atan2(1.0,1.0);
    double tau = T-t;
    if (tau < 0) return 1;
    if (tau == 0) {
        call = (S > K) ? (S-K) : 0.0;
        delta = (S > K) ? 1.0 : 0.0;
        return 0;
    }
    if (sigma < 0) return 1;
    if (K <= 0) {
        call = S*exp(-q*tau) -K*exp(-r*tau);
        delta = exp(-q*tau);
        return 0;
    }
    if (S <= 0) return 1;
    double psi = sigma*sqrt(tau);
    double d1 = (log(S/K) +(r-q)*tau)/psi +0.5*psi;
    double d2 = d1 - psi;
    double Nd1 = cum_norm(d1);
    double Nd2 = cum_norm(d2);
    call = S*exp(-q*tau)*Nd1 -K*exp(-r*tau)*Nd2;
    delta = exp(-q*tau)*Nd1;
    return 0;
}

void compute_X(double &c1, double &c2, double &c3, double &X){
    double B = c1- 2*c2+c3;
    X =100*B;
}

int main(int argc, const char * argv[]) {
    // insert code here...
    cout<<"5.3: Butterfly Spread"<<endl;
    cout<<"S = 100, q = 1%, sigma = 45%, r = 5%, t=0, T=1"<<endl;
    double S=100;
    double q=0.01;
    double sigma = 0.45;
    double r=0.05;
    double t=0.0;
    double T=1.0;
    double call, delta;
    double c1, c2, c3, K, X;
    
    cout<<"--------------------"<<endl;
    BS_call(S, 49, r, q, sigma, T, t, call, delta);
    c1 = call;
    BS_call(S, 50, r, q, sigma, T, t, call, delta);
    c2 = call;
    BS_call(S, 51, r, q, sigma, T, t, call, delta);
    c3= call;
    compute_X(c1,c2,c3,X);
    cout<<"K=50, X = "<<X<<endl;
    
     cout<<"--------------------"<<endl;
    BS_call(S, 74, r, q, sigma, T, t, call, delta);
    c1 = call;
    BS_call(S, 75, r, q, sigma, T, t, call, delta);
    c2 = call;
    BS_call(S, 76, r, q, sigma, T, t, call, delta);
    c3= call;
    compute_X(c1,c2,c3,X);
    cout<<"K=75, X = "<<X<<endl;
    
    cout<<"--------------------"<<endl;
    BS_call(S, 99, r, q, sigma, T, t, call, delta);
    c1 = call;
    BS_call(S, 100, r, q, sigma, T, t, call, delta);
    c2 = call;
    BS_call(S, 101, r, q, sigma, T, t, call, delta);
    c3= call;
    compute_X(c1,c2,c3,X);
    cout<<"K=100, X = "<<X<<endl;
    
    cout<<"--------------------"<<endl;
    BS_call(S, 124, r, q, sigma, T, t, call, delta);
    c1 = call;
    BS_call(S, 125, r, q, sigma, T, t, call, delta);
    c2 = call;
    BS_call(S, 126, r, q, sigma, T, t, call, delta);
    c3= call;
    compute_X(c1,c2,c3,X);
    cout<<"K=126, X = "<<X<<endl;
    
    cout<<"--------------------"<<endl;
    BS_call(S, 149, r, q, sigma, T, t, call, delta);
    c1 = call;
    BS_call(S, 150, r, q, sigma, T, t, call, delta);
    c2 = call;
    BS_call(S, 151, r, q, sigma, T, t, call, delta);
    c3= call;
    compute_X(c1,c2,c3,X);
    cout<<"K=150, X = "<<X<<endl;
    cout<<"================================"<<endl;
    cout<<"5.4: Implied volatility"<<endl;
    cout<<"S = 100, q = 1%, r = 5%, t=0, T=1, K=150"<<endl;
    K=150;

    BS_call(S, K, r, q, 0.1, T, t, call, delta);
    cout<<"sigma=0.1, Call = "<<call<<endl;
    BS_call(S, K, r, q, 0.2, T, t, call, delta);
    cout<<"sigma=0.2, Call = "<<call<<endl;
    BS_call(S, K, r, q, 0.3, T, t, call, delta);
    cout<<"sigma=0.3, Call = "<<call<<endl;
    BS_call(S, K, r, q, 0.31, T, t, call, delta);
    cout<<"sigma=0.31, Call = "<<call<<endl;
    BS_call(S, K, r, q, 0.4, T, t, call, delta);
    cout<<"sigma=0.4, Call = "<<call<<endl;
    BS_call(S, K, r, q, 0.5, T, t, call, delta);
    cout<<"sigma=0.5, Call = "<<call<<endl;
    cout<<"-----------------------"<<endl;
    cout<<"Finding sigma that 1.99<c<2.01"<<endl;
    BS_call(S, K, r, q, 0.3035, T, t, call, delta);
    cout<<"sigma=0.3035, Call = "<<call<<endl;
    BS_call(S, K, r, q, 0.304, T, t, call, delta);
    cout<<"sigma=0.304, Call = "<<call<<endl;
    BS_call(S, K, r, q, 0.3042, T, t, call, delta);
    cout<<"sigma=0.3042, Call = "<<call<<endl;
    BS_call(S, K, r, q, 0.3043, T, t, call, delta);
    cout<<"sigma=0.3043, Call = "<<call<<endl;
    BS_call(S, K, r, q, 0.305, T, t, call, delta);
    cout<<"sigma=0.305, Call = "<<call<<endl;
    cout<<"---------------------------"<<endl;
    cout<<"Calculate the Delta of a call option using the sigma_0, where sigma_0 = 0.304"<<endl;
    BS_call(S, K, r, q, 0.304, T, t, call, delta);
    cout<<"sigma_0=0.304, Delta = "<<delta<<endl;
    cout<<"Finish."<<endl;
    
    return 0;
}
