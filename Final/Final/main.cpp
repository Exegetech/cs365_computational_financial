//
//  main.cpp
//  Final
//
//  Created by XingSong  Lin on 12/6/17.
//  Copyright © 2017 XingSong  Lin. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

//Question 1
void price_from_yield(double F, double c[], double y, int n, double & B)
{
    double y_decimal = 0.01*y;
    B=0.0;
    for(int i=1; i<=n; i++){
        if(i!=n){
            B+= (c[i-1]/2)/(pow(1+0.5*y_decimal, i));
        }else{
            B+=(F+c[i-1]/2)/(pow(1+0.5*y_decimal, i));
        }
    }
}

int yield_from_price(double F, double c[], int n, double B_market, double tol, int max_iter, double & y)
{
    double B, B_y_low, B_y_high, y_low, y_high;
    y_low=0;
    price_from_yield(F,c, y_low, n, B_y_low);
    if(abs(B_y_low-B_market) <= tol){
        y=y_low;
        return 0;
    }
    
    y_high=100;
    price_from_yield(F, c, y_high, n, B_y_high);
    if(abs(B_y_high-B_market)<=tol){
        y=y_high;
        return 0;
    }
    
    if(B_y_low< B_market || B_y_high>B_market){
        y=0;
        return 1;
    }
    //step 1-8 end.
    
    for(int i=0; i<max_iter; ++i){
        y=(y_low+y_high)/2.0;
        price_from_yield(F, c, y, n, B);
        
        if(abs(B-B_market)<=tol){
            return 0;
        }else if (B>B_market){
            y_low=y;
        }else{
            y_high=y;
        }
        
        if((y_high-y_low)<=tol){
            return 0;
        }
        
        if(i==(max_iter-1)){
            y=0;
            return 1;
        }
    }
    return 0;
}

void compute_Question1(double &y){
    ofstream outfile;
    outfile.open("Solution_for_Q1.txt");
    outfile<<"Student ID: 23083576 --> array of{2,3,0,8,3,5,7,6}\n\n";
    double c[8] = {2,3,0,8,3,5,7,6};
    double F = 100;
    y = 5;
    int max_iter=100;
    int n = 8;
    double tol = 1.0e-4;
    double B;
    price_from_yield(F, c, y, n, B);
    outfile<<"The theoretical fair value of the bond is "<<fixed<<setprecision(2)<<B<<" if yield of the bond is 5%.\n\n";
    double B_market = 100;
    yield_from_price(F, c, n, B_market, tol, max_iter, y);
    outfile<<"The yield of the bond is "<<fixed<<setprecision(2)<<y<<"% if market price of the bond is 100.\n\n";
    outfile<<"\n---------------------------Separated Line-----------------------------------\n\n";
    outfile.close();
}

//Question 2
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

//Question 3
class BinomialModel
{
public:
    BinomialModel(int n);
    ~BinomialModel();
    int FairValue(int n, Derivative * p_derivative, double S, double t0, double & V);
    
    int ImpliedVolatility(int n, Derivative * p_derivative, double S, double t0, double target, double & implied_vol, int & num_iter);
private:
    // methods
    void Clear();
    int Allocate(int n);
    int ImpliedVolatilityPrivate(int n, Derivative * p_derivative, double S, double t0, double target, double & implied_vol, int & num_iter);
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
    if (n_tree <= 0) return;
    
    for(int i=0; i <= n_tree; i++){
        delete[] stock_nodes[i];
        delete[] value_nodes[i];
    }
    delete[] stock_nodes;
    delete[] value_nodes;
    stock_nodes= NULL;
    value_nodes= NULL;
    n_tree = 0;
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

int BinomialModel::ImpliedVolatility(int n, Derivative * p_derivative, double S, double t0,
                                     double target, double & implied_vol, int & num_iter)
{
    int rc = 0;
    const double saved_vol = p_derivative->sigma;
    rc = ImpliedVolatilityPrivate(n, p_derivative, S, t0, target, implied_vol, num_iter);
    p_derivative->sigma = saved_vol;
    return rc;
}

int BinomialModel::ImpliedVolatilityPrivate(int n, Derivative * p_derivative, double S, double t0, double target,
                                            double & implied_vol, int & num_iter)
{
    const double tol = 1.0e-6;
    const int max_iter = 100;
    
    implied_vol = 0;
    num_iter = 0;
    
    double sigma_low = 0.01; // 1%
    double sigma_high = 3.0; // 300%
    double FV_low = 0;
    double FV_high = 0;
    double FV = 0;
    
    p_derivative->sigma = sigma_low;
    FairValue(n, p_derivative, S, t0, FV_low);  //FV_low is 0 after computed!!!
    double diff_FV_low = FV_low - target;
    if (fabs(diff_FV_low) <= tol) {  //false
        implied_vol = p_derivative->sigma;
        return 0;
    }
    
    p_derivative->sigma = sigma_high;
    FairValue(n, p_derivative, S, t0, FV_high);//FV_high is 0 after computed!!!
    double diff_FV_high = FV_high-target;
    if (fabs(diff_FV_high) <= tol) {  //false
        implied_vol = p_derivative->sigma;
        return 0;
    }
    
    if (diff_FV_low * diff_FV_high > 0) {  //break here
        implied_vol = 0;
        return 1; // fail
    }
    
    for (int i = 0; i < max_iter; ++i) {
        
        p_derivative->sigma = 0.5*(sigma_low + sigma_high);
        FairValue(n, p_derivative, S, t0, FV);
        double diff_FV = FV - target;
        
        if(abs(diff_FV)<=tol){
            implied_vol = p_derivative->sigma;
            num_iter = i;
            return 0;
        }else if(diff_FV_low * diff_FV >0){
            sigma_low =p_derivative->sigma;
        }else{
            sigma_high=p_derivative->sigma;
        }
        
        if(abs(sigma_low - sigma_high) <= tol){
            implied_vol = p_derivative->sigma;
            num_iter = i;
            return 0;
        }
        
        if(i==(max_iter-1)){
            num_iter=max_iter;
            implied_vol=0;
            return 1;
        }
    }
    
    return 0;
}

//Question 4
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

//Question 5
class straddle : public Derivative
{
public:
    straddle(){K = 0; isCall = false; isAmerican = false;}
    virtual ~straddle(){}
    virtual double TerminalPayoff(double S) const;
    virtual int ValuationTest(double S, double &V) const;
    
    double K;
    bool isCall;
    bool isAmerican;
};

double straddle::TerminalPayoff(double S) const
{
    double temp =0.0;
    if(isCall){
        if(S>K){
            temp+=S-K;
        }
    }else{
        if(S<K){
            temp+=K-S;
        }
    }
    return temp;
}

int straddle::ValuationTest(double S, double &V) const
{
    if(isAmerican){
        double temp =0.0;
        if(isCall){
            temp += max(V, max(S-K, 0.0));
        }else{
            temp += max(V, max(K-S, 0.0));
        }
        return temp;
    }
    
    return 0;
}

//Question 6
class BinaryOption : public Derivative
{
public:
    BinaryOption() { K = 0; isCall = false; isAmerican = false; }
    virtual ~BinaryOption() {}
    virtual double TerminalPayoff(double S) const;
    virtual int ValuationTests(double S, double &V) const;
    // data
    double K;
    bool isCall;
    bool isAmerican;
};
double BinaryOption::TerminalPayoff(double S) const
{
    // *** RETURN TERMINAL PAYOFF FOR PUT OR CALL OPTION ***
    if(isCall){
        if(S>=K){
            return 1;
        }else{
            return 0;
        }
    }else{
        if(S<K){
            return 1;
        }else{
            return 0;
        }
    }
}
int BinaryOption::ValuationTests(double S, double &V) const
{
    // *** TEST IF THE VALUE OF V SHOULD BE UPDATED TO THE INTRINSIC VALUE ***
    //only do European, no need for this American early exercise.
    return 0;
}

class UpOutBarrierCallOption: public Derivative
{
public:
    UpOutBarrierCallOption(){K=0; B=0; isCall=false; isAmerican=false;}
    virtual ~UpOutBarrierCallOption(){}
    virtual double TerminalPayoff(double S) const;
    virtual int ValuationTests(double S, double &V) const;
    // data
    double K;
    double B; //barrier
    bool isCall;
    bool isAmerican;
};

double UpOutBarrierCallOption::TerminalPayoff(double S) const
{
    //K<B , and S depends.
    if(K<=S && S<B){
        return S-K;
    }
    if(S>=B){
        return B-K;
    }
    if(S<K){
        return 0;
    }
    return 0;
}

int UpOutBarrierCallOption::ValuationTests(double S, double &V) const
{
    if(isAmerican){
        double V_barrier = B-K;
        if(K<S && S<B && V_barrier<(S-K)){
            return V_barrier;
        }
    }
    return 0;
}

int Question4_test(){
    int rc = 0;
    
    // output file
    ofstream ofs;
    ofs.open("Solution_for_Q4.txt");
    
    //srand((unsigned)time(NULL));
    
    double S = 590;
    double K = 150;
    double r = 0.12;
    double q = 0.04;
    double sigma = 0.8;
    double T = 2.0;
    double t0 = 1.3;
    /*
    double S = rand();
    double K = rand();
    double r = (double)(rand()%101)/101;
    double q = (double)(rand()%101)/101;
    double sigma = (double)(rand()%101)/101;
    double T = (double)(rand()%101)/101;
    double t0 = (double)(rand()%101)/101;
    
    ofs<<S<<" "<<K<<" "<<r<<" "<<q<<" "<<sigma<<" "<<T<<" "<<t0<<endl;
    */
    
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
    for (int i = 1; i <= imax; ++i) {
        S = i*dS;
        rc = binom.FairValue(n, &Am_put, S, t0, FV_Am_put);
        rc = binom.FairValue(n, &Eur_put, S, t0, FV_Eur_put);
        rc = binom.FairValue(n, &Am_call, S, t0, FV_Am_call);
        rc = binom.FairValue(n, &Eur_call, S, t0, FV_Eur_call);
    }
    ofs << setw(16) << "S" << " ";
    ofs << setw(16) << "FV_Am_put" << " ";
    ofs << setw(16) << "FV_Eur_put" << " ";
    ofs << setw(16) << "FV_Am_call" << " ";
    ofs << setw(16) << "FV_Eur_call" << " ";
    ofs << endl;
    ofs << setw(16) << S << " ";
    ofs << setw(16) << FV_Am_put << " ";
    ofs << setw(16) << FV_Eur_put << " ";
    ofs << setw(16) << FV_Am_call << " ";
    ofs << setw(16) << FV_Eur_call << " ";
    ofs << endl;
    ofs <<"---------------------------Separated Line-----------------------------------\n\n";
    
    //test 4.1
    ofs << "Test for Equation 4.1: \n\n";
    double tol = 1.0e-6;
    double V_c_sub_p = FV_Eur_call - FV_Eur_put;
    double V_C_sub_P = FV_Am_call - FV_Am_put;
    double dt = T-t0;
    double e1 = exp(-q*dt);
    double V_Se = S * e1;
    double e2 = exp(-r*dt);
    double V_Ke = K * e2;
    bool check_result=false;
    /*
    ofs<<"Computed: \n";
    ofs<<"tol: "<<tol;
    ofs<<"\nc-p: "<<V_c_sub_p;
    ofs<<"\nC-P: "<<V_C_sub_P;
    ofs<<"\nS*e~: "<<V_Se;
    ofs<<"\nK*e~: "<<e2<<endl;
    */
    ofs<<"Test--> 4.1: \n\n";
    if(abs((V_c_sub_p)-(V_Se-V_Ke)) <= tol){
        check_result = true;
    }else{
        check_result = false;
    }
    ofs<<"Result for passing equation 4.1 is "<<boolalpha<<check_result<<".\n\n";
    
    ofs<<"Test--> 4.2: \n\n";
    if((V_C_sub_P >= (V_Se-K-tol)) && (V_C_sub_P<=(S-V_Ke+tol))){
        check_result = true;
    }else{
        check_result = false;
    }
    ofs<<"Result for passing equation 4.2 is "<<check_result<<".\n\n";
    
    ofs<<"Test--> 4.3: \n\n";
    if((0<= FV_Eur_call && FV_Eur_call <= (V_Se+tol)) && (0<=FV_Am_call && FV_Am_call <= S+tol) &&
       (0<= FV_Eur_put && FV_Eur_put <= (V_Se+tol)) && (0<=FV_Am_put && FV_Am_put<= K+tol)){
        check_result = true;
    }
    else{
        check_result = false;
    }
    ofs<<"Result for passing equation 4.3 is "<<check_result<<".\n\n";
    
    ofs.close();
    return 0;
}

int main(int argc, const char * argv[]) {
    
    double yield;
    compute_Question1(yield);
    ofstream of;
    of.open("Solution_for_Q1.txt", ios_base::app);
    of<<"yield without rounding to 2 decimal is "<<yield<<"\n";
    of.close();
    
    Question4_test();
    return 0;
}
