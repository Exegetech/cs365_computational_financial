//
//  main.cpp
//  HW8_365
//
//  Created by XingSong  Lin on 11/27/17.
//  Copyright © 2017 XingSong  Lin. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

int binomial_simple(double S, double K, double r, double q, double sigma, double T, double t0, bool call, bool American, int n, double & V){
    //7.2: Validation tests
    if(n<1 || S<=0 || T <= t0 || sigma <= 0.0){
        return 1;
    }
    
    //7.3: Parameters
    double dt = (T-t0)/double(n);
    double df = exp(-r*dt);
    double growth = exp((r-q)*dt);
    double u = exp(sigma*sqrt(dt));
    double d = 1.0/u;
    double p_prob = (growth - d)/(u-d);
    double q_prob = 1.0 - p_prob;
    
    if(p_prob<0.0 || p_prob>1.0){
        return 1;
    }
    
    double** stock_nodes = new double*[n+1];
    double** option_nodes = new double*[n+1];
    double* S_tmp;
    double* V_tmp;
    for (int i = 0; i <= n; ++i) {
        stock_nodes[i] = new double[n+1];
        option_nodes[i] = new double[n+1];
        S_tmp = stock_nodes[i];
        V_tmp = option_nodes[i];
        for (int j = 0; j <= n; ++j) {
            S_tmp[j] = 0;
            V_tmp[j] = 0;
        }
    }
    
    S_tmp = stock_nodes[0];
    S_tmp[0]=S;
    
    for (int i = 1; i <= n; ++i) {
        double * prev = stock_nodes[i-1];
        S_tmp = stock_nodes[i];
        S_tmp[0] = prev[0] * d;
        for (int j = 1; j <= n; ++j) {
            S_tmp[j] = S_tmp[j-1]*u*u;
        }
    }
    
    int i=n;
    S_tmp = stock_nodes[i];
    V_tmp = option_nodes[i];
    for(int j=0; j<=n; ++j){
        double instrinsic =0;
        if(call){
            if(S_tmp[j]>K){
                instrinsic = S_tmp[j]-K;
            }
        }else{
            if(S_tmp[j]<K){
                instrinsic = K-S_tmp[j];
            }
        }
        V_tmp[j]= instrinsic;
    }
    
    for(int i=n-1; i>=0; --i){
        S_tmp = stock_nodes[i];
        V_tmp = option_nodes[i];
        
        double * V_next = option_nodes[i+1];
        for (int j = 0; j <= i; ++j) {
            V_tmp[j] = df*(p_prob*V_next[j+1] + q_prob*V_next[j]);
            
            //early exercise test
            if(American){
                if(!call){
                    double temp = max(V_tmp[j], max(K-S_tmp[j], 0.0));
                    V_tmp[j]=temp;
                }else{
                    double temp = max(V_tmp[j], max(S_tmp[j]-K, 0.0));
                    V_tmp[j]= temp;
                }
            }
        }
    }
    
    int index=0;
    V_tmp = option_nodes[index];
    V=V_tmp[0];
    
    for (int i = 0; i <= n; ++i) {
        delete [] stock_nodes[i];
        delete [] option_nodes[i];
    }
    delete [] stock_nodes;
    delete [] option_nodes;
    
    return 0;
}

class Derivative
{
public:
    virtual ~Derivative() {}
    virtual double TerminalPayoff(double S) const { return 0; }
    virtual int ValuationTests(double S, double & V) const { return 0; }
    // data
    double r;
    double q;
    double sigma;
    double T;
protected:
    Derivative() { r = 0; q = 0; sigma = 0; T = 0; }
};

class Option : public Derivative
{
public:
    Option() { K = 0; isCall = false; isAmerican = false; }
    virtual ~Option() {}
    virtual double TerminalPayoff(double S) const;
    virtual int ValuationTests(double S, double &V) const;
    // data
    double K;
    bool isCall;
    bool isAmerican;
};
double Option::TerminalPayoff(double S) const
{
    // *** RETURN TERMINAL PAYOFF FOR PUT OR CALL OPTION ***
    if(isCall){
        if(S>K){
            return S-K;
        }
    }else{
        if(S<K){
            return K-S;
        }
    }
    return 0;
}
int Option::ValuationTests(double S, double &V) const
{
    // *** TEST IF THE VALUE OF V SHOULD BE UPDATED TO THE INTRINSIC VALUE ***
    if(isAmerican){
        if(isCall){
            double temp =max(V, max(S-K, 0.0));
            V= temp;
        }else{
            double temp = max(V, max(K-S, 0.0));
            V= temp;
        }
    }
    return 0;
}


class BinomialModel
{
public:
    BinomialModel(int n);
    ~BinomialModel();
    int FairValue(int n, Derivative * p_derivative, double S, double t0, double & V);
private:
    // methods
    void Clear();
    int Allocate(int n);
    // data
    int n_tree;
    double **stock_nodes;
    double **value_nodes;
};

BinomialModel::BinomialModel(int n)
{
    n_tree = 0;
    stock_nodes = NULL;
    value_nodes = NULL;
    Allocate(n);
    
}
BinomialModel::~BinomialModel()
{
    Clear();
}
void BinomialModel::Clear()
{
    // *** WRITE THE FUNCTION TO RELEASE ALLOCATED MEMORY ***
    stock_nodes = NULL;
    value_nodes = NULL;
    /*
     for (int i = 0; i <= n; ++i) {
     delete [] stock_nodes[i];
     delete [] value_nodes[i];
     }
     delete [] stock_nodes;
     delete [] value_nodes;
     */
}

int BinomialModel::Allocate(int n)
{
    if (n <= n_tree) return 0;
    // deallocate old tree
    
    Clear();
    // allocate memory
    n_tree = n;
    
    // *** WRITE THE FUNCTION TO ALLOCATE NEW MEMORY ***
    stock_nodes = new double*[n_tree+1];
    value_nodes = new double*[n_tree+1];
    for (int i = 0; i <= n_tree; ++i) {
        stock_nodes[i] = new double[n_tree+1];
        value_nodes[i] = new double[n_tree+1];
    }
    
    return 0;
}

int BinomialModel::FairValue(int n, Derivative * p_derivative, double S, double t0, double & V)
{
    int rc = 0;
    V = 0;
    // validation checks
    if(n<1 || S<=0 || p_derivative == NULL || p_derivative->T <= t0 || p_derivative->sigma <= 0.0 ){
        return 1;
    }
    // declaration of local variables (I use S_tmp and V_tmp)
    // declarated in constoctor.
    double* S_tmp=NULL;
    double* V_tmp=NULL;
    // calculate parameters
    double dt = (p_derivative->T - t0)/double(n);
    double df = exp(-p_derivative->r*dt);
    double growth = exp((p_derivative->r - p_derivative->q)*dt);
    double u = exp(p_derivative->sigma*sqrt(dt));
    double d = 1.0/u;
    double p_prob = (growth - d)/(u-d);
    double q_prob = 1.0 - p_prob;
    
    // more validation checks
    if(p_prob<0.0 || p_prob>1.0){
        return 1;
    }
    
    // allocate memory if required (call Allocate(n))
    Allocate(n);
    for(int i=0; i<=n; ++i){
        S_tmp = stock_nodes[i];
        V_tmp = value_nodes[i];
        for (int j = 0; j <= n; ++j) {
            S_tmp[j] = 0;
            V_tmp[j] = 0;
        }
    }
    // set up stock prices in tree
    S_tmp = stock_nodes[0];
    S_tmp[0]=S;
    
    for (int i = 1; i <= n; ++i) {
        double * prev = stock_nodes[i-1];
        S_tmp = stock_nodes[i];
        S_tmp[0] = prev[0] * d;
        for (int j = 1; j <= n; ++j) {
            S_tmp[j] = S_tmp[j-1]*u*u;
        }
    }
    
    // set terminal payoff (call virtual function in derivative class to calculate payoff)
    int i = n;
    S_tmp = stock_nodes[i];
    V_tmp = value_nodes[i];
    for (int j = 0; j <= n; ++j) {
        V_tmp[j] = p_derivative->TerminalPayoff(S_tmp[j]);
    }
    
    // valuation loop (call virtual function in derivative class for valuation tests)
    for (int i = n-1; i >= 0; --i) {
        S_tmp = stock_nodes[i];
        V_tmp = value_nodes[i];
        double * V_next = value_nodes[i+1];
        for (int j = 0; j <= i; ++j) {
            V_tmp[j] = df*(p_prob*V_next[j+1] + q_prob*V_next[j]);
            p_derivative->ValuationTests(S_tmp[j], V_tmp[j]); // VALUATION TESTS
        }
    }
    
    // derivative fair value
    V_tmp = value_nodes[0];
    V = V_tmp[0];
    
    // deallocate memory (if necessary)
    
    
    return 0;
}


int binomial_test()
{
    int rc = 0;
    
    // output file
    ofstream ofs;
    ofs.open("output.txt");
    
    double S = 100;
    double K = 100;
    double r = 0.05;
    double q = 0.01;
    double sigma = 0.5;
    double T = 1.0;
    double t0 = 0;
    
    Option Eur_put;
    Eur_put.r = r;
    Eur_put.q = q;
    Eur_put.sigma = sigma;
    Eur_put.T = T;
    Eur_put.K = K;
    Eur_put.isCall = false;
    Eur_put.isAmerican = false;
    
    Option Am_put;
    Am_put.r = r;
    Am_put.q = q;
    Am_put.sigma = sigma;
    Am_put.T = T;
    Am_put.K = K;
    Am_put.isCall = false;
    Am_put.isAmerican = true;
    
    Option Eur_call;
    Eur_call.r = r;
    Eur_call.q = q;
    Eur_call.sigma = sigma;
    Eur_call.T = T;
    Eur_call.K = K;
    Eur_call.isCall = true;
    Eur_call.isAmerican = false;
    
    Option Am_call;
    Am_call.r = r;
    Am_call.q = q;
    Am_call.sigma = sigma;
    Am_call.T = T;
    Am_call.K = K;
    Am_call.isCall = true;
    Am_call.isAmerican = true;
    
    double FV_Am_put = 0;
    double FV_Eur_put = 0;
    double FV_Am_call = 0;
    double FV_Eur_call = 0;
    
    int n = 100;
    
    BinomialModel binom(n);
    
    double dS = 0.1;
    int imax = 2000;
    int i;
    for (i = 1; i <= imax; ++i) {
        S = i*dS;
        rc = binom.FairValue(n, &Am_put, S, t0, FV_Am_put);
        rc = binom.FairValue(n, &Eur_put, S, t0, FV_Eur_put);
        rc = binom.FairValue(n, &Am_call, S, t0, FV_Am_call);
        rc = binom.FairValue(n, &Eur_call, S, t0, FV_Eur_call);
        
        ofs << setw(16) << S << " ";
        ofs << setw(16) << FV_Am_put << " ";
        ofs << setw(16) << FV_Eur_put << " ";
        ofs << setw(16) << FV_Am_call << " ";
        ofs << setw(16) << FV_Eur_call << " ";
        ofs << endl;
    }
    binom.~BinomialModel();
    ofs.close();
    return 0;
}

int main(int argc, const char * argv[]) {
    /*
     //double check for code from hw7. result is same
     double S = 100;
     double K = 100;
     double r = 0.1;
     double q = 0.0;
     double sigma = 0.5;
     double T = 0.3;
     double t0 = 0;
     double n = 3;
     double FV_Am_put = 0;
     double FV_Eur_put = 0;
     double FV_Am_call = 0;
     double FV_Eur_call = 0;
     binomial_simple(S, K, r, q, sigma, T, t0, false, true, n, FV_Am_put);
     binomial_simple(S, K, r, q, sigma, T, t0, false, false, n, FV_Eur_put);
     binomial_simple(S, K, r, q, sigma, T, t0, true, true, n, FV_Am_call);
     binomial_simple(S, K, r, q, sigma, T, t0, true, false, n, FV_Eur_call);
     cout << FV_Am_put << " "<<FV_Eur_put << " "<< FV_Am_call << " "<< FV_Eur_call << " "<< endl;
     */
    ofstream outfile;
    outfile.open("result_8.2.txt");
    double S = 0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.01;
    double sigma = 0.5;
    double T = 1.0;
    double t0 = 0.0;
    double FV_Am_put = 0;
    double FV_Eur_put = 0;
    double FV_Am_call = 0;
    double FV_Eur_call = 0;
    int n = 100;
    double dS = 0.1;
    int imax = 2000;
    for (int i = 1; i <= imax; ++i) {
        S = i*dS;
        binomial_simple(S, K, r, q, sigma, T, t0, false, true, n, FV_Am_put);
        binomial_simple(S, K, r, q, sigma, T, t0, false, false, n, FV_Eur_put);
        binomial_simple(S, K, r, q, sigma, T, t0, true, true, n, FV_Am_call);
        binomial_simple(S, K, r, q, sigma, T, t0, true, false, n, FV_Eur_call);
        // print output to file
        outfile << setw(16) << S << " ";
        outfile << setw(16) << FV_Am_put << " ";
        outfile << setw(16) << FV_Eur_put << " ";
        outfile << setw(16) << FV_Am_call << " ";
        outfile << setw(16) << FV_Eur_call << " ";
        outfile << endl;
    }
    outfile.close();
    
    binomial_test();
    
    return 0;
}



