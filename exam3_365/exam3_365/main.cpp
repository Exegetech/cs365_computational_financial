//
//  main.cpp
//  exam3_365
//
//  Created by XingSong  Lin on 11/20/17.
//  Copyright © 2017 XingSong  Lin. All rights reserved.
//

#include <iostream>
#include <cmath>
using namespace std;

class Binomial_model{
public:
    double S, K, r, q, sigma, T, t0, V, dt, df, growth, u, d, p_prob, q_prob;
    double* S_tmp;
    double* V_tmp;
    double** stock_nodes;
    double** option_nodes;
    bool call, American;
    int n;
    
    Binomial_model(){
    }
    
    int binomial_simple(double S, double K, double r, double q, double sigma, double T, double t0, bool call, bool American, int n, double & V){
        //7.2: Validation tests
        if(n<1 || S<=0 || T <= t0 || sigma <= 0.0){
            return 1;
        }
        
        //7.3: Parameters
        dt = (T-t0)/double(n);
        df = exp(-r*dt);
        growth = exp((r-q)*dt);
        u = exp(sigma*sqrt(dt));
        d = 1.0/u;
        p_prob = (growth - d)/(u-d);
        q_prob = 1.0 - p_prob;
        
        /*
        cout<<"dt = "<<dt<<endl;
        cout<<"df = "<<df<<endl;
        cout<<"growth = "<<growth<<endl;
        cout<<"u = "<<u<<endl;
        cout<<"d = "<<d<<endl;
        cout<<"p_prob = "<<p_prob<<endl;
        cout<<"q_prob = "<<q_prob<<endl;
        */
        
        if(p_prob<0.0 || p_prob>1.0){
            return 1;
        }
        
        return 0;
    }
    
    void binomial_simple_Array(){
        // allocate memory
        stock_nodes = new double*[n+1];
        option_nodes = new double*[n+1];
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
    }
    
    //7.5 setup stock price in nodes
    void binomial_simple_Setup(){
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
    }
    
    //7.6 terminal payoff
    void terminal_payoff(int layer){
        int i=layer;
        S_tmp = stock_nodes[i];
        V_tmp = option_nodes[i];
        for(int j=0; j<=layer; ++j){
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
    }
    
    //7.7: main valuation loop
    void main_loop(){
        
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
    }
    
    //7.8: Option Fair Value
    double option_fair_value(){
        int i=0;
        V_tmp = option_nodes[i];
        V=V_tmp[0];
        return V;
    }
    
    //7.10: Test
    void test_1(bool boolean1, bool boolean2){
        //double S, double K, double r, double q, double sigma, double T, double t0, bool call, bool American, int n, double & V
        call = boolean1; //call ->1, put->0
        American = boolean2; // European->0, American->1
        binomial_simple(S, K, r, q, sigma, T,  t0, call, American, n, V);
        
        binomial_simple_Array();
        
        //cout<<"7.5: Test function: Using lecture_17a.pdf data \ndefault: "<<endl;
        //print_array(stock_nodes);
        
        binomial_simple_Setup();

        //cout<<"After setup: "<<endl;
        //print_array(stock_nodes);
        
        //cout<<"\n===================================================="<<endl;
        //cout<<"7.6: Test function: terminal payoff"<<endl;
        
        terminal_payoff(n);
        
        //print_array(option_nodes);
        
        //cout<<"============================================"<<endl;
        //cout<<"7.7: Main validation Loop"<<endl;
        
        main_loop();
        
        //print_array(option_nodes);
    }
    
    void print_array(double** array){
        for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= n; ++j) {
                cout<<"at ("<<i<<", "<<j<<") with stock price "<<array[i][j]<<endl;
            }
            cout<<endl;
        }
    }
    
    //7.9
    ~Binomial_model(){
        for (int i = 0; i <= n; ++i) {
            delete [] stock_nodes[i];
            delete [] option_nodes[i];
        }
        delete [] stock_nodes;
        delete [] option_nodes;
    }
    
};

int main(int argc, const char * argv[]) {
    
    Binomial_model* model_check = new Binomial_model();
    
    // boolean “call” (true for a call, false for a put).
    // boolean “American” (true for an American option, false for European).
    cout<<"7.10: Test"<<endl;
    model_check->S = 100;
    model_check->K = 100;
    model_check->r = 0.1;
    model_check->q = 0.0;
    model_check->sigma = 0.5;
    model_check->T = 0.3;
    model_check->t0 = 0;
    model_check->n = 3;
    model_check->V = 100;
    model_check->test_1(true, true);
    cout<<"American Call: "<<endl;
    cout<<"The option fair value, C = "<<model_check->option_fair_value()<<endl;

    model_check->test_1(false, true);
    cout<<"American Put: "<<endl;
    cout<<"The option fair value, P = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(true, false);
    cout<<"European Call: "<<endl;
    cout<<"The option fair value, c = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(false, false);
    cout<<"European Put: "<<endl;
    cout<<"The option fair value, p = "<<model_check->option_fair_value()<<endl;
    
    cout<<"======================================="<<endl;
    //7.11
    model_check->q = 0.1;
    model_check->T = 1;
    model_check->n = 100;
    model_check->test_1(true, true);
    cout<<"Change q=0.1, T=1, n=100."<<endl;
    cout<<"American Call: "<<endl;
    cout<<"The option fair value, C = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(false, true);
    cout<<"American Put: "<<endl;
    cout<<"The option fair value, P = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(true, false);
    cout<<"European Call: "<<endl;
    cout<<"The option fair value, c = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(false, false);
    cout<<"European Put: "<<endl;
    cout<<"The option fair value, p = "<<model_check->option_fair_value()<<endl;
    
    //additional tests:
    cout<<"======================================="<<endl;
    model_check->S = 100;
    model_check->K = 100;
    model_check->r = 0.1;
    model_check->q = 0.1;
    model_check->sigma = 0.4;
    model_check->T = 0.4;
    model_check->t0 = 0;
    model_check->n = 100;
    model_check->V = 100;
    cout<<"Addition test: RESET to default value but sigma=0.4, T=0.4."<<endl;
    cout<<"American Call: "<<endl;
    cout<<"The option fair value, C = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(false, true);
    cout<<"American Put: "<<endl;
    cout<<"The option fair value, P = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(true, false);
    cout<<"European Call: "<<endl;
    cout<<"The option fair value, c = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(false, false);
    cout<<"European Put: "<<endl;
    cout<<"The option fair value, p = "<<model_check->option_fair_value()<<endl;
    
    cout<<"======================================="<<endl;
    model_check->S = 110;
    model_check->K = 110;
    model_check->r = 0.1;
    model_check->q = 0.1;
    model_check->sigma = 0.5;
    model_check->T = 1;
    model_check->t0 = 0;
    model_check->n = 100;
    model_check->V = 100;
    cout<<"Addition test: RESET to default value but S=K=110."<<endl;
    cout<<"American Call: "<<endl;
    cout<<"The option fair value, C = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(false, true);
    cout<<"American Put: "<<endl;
    cout<<"The option fair value, P = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(true, false);
    cout<<"European Call: "<<endl;
    cout<<"The option fair value, c = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(false, false);
    cout<<"European Put: "<<endl;
    cout<<"The option fair value, p = "<<model_check->option_fair_value()<<endl;
    cout<<"======================================="<<endl;
    model_check->S = 100;
    model_check->K = 100;
    model_check->r = 0.3;
    model_check->q = 0.3;
    model_check->sigma = 0.5;
    model_check->T = 1;
    model_check->t0 = 0;
    model_check->n = 100;
    model_check->V = 100;
    cout<<"Addition test: RESET to default value but r=q=0.3"<<endl;
    cout<<"American Call: "<<endl;
    cout<<"The option fair value, C = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(false, true);
    cout<<"American Put: "<<endl;
    cout<<"The option fair value, P = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(true, false);
    cout<<"European Call: "<<endl;
    cout<<"The option fair value, c = "<<model_check->option_fair_value()<<endl;
    
    model_check->test_1(false, false);
    cout<<"European Put: "<<endl;
    cout<<"The option fair value, p = "<<model_check->option_fair_value()<<endl;
    
    //7.9: Memory Dealocation.
    model_check->~Binomial_model();
    return 0;
}
