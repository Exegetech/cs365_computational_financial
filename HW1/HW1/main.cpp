//
//  main.cpp
//  HW1
//
//  Created by Eric Lin on 9/3/17.
//  Copyright © 2017 XingSong  Lin. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
using namespace std;

//1.1 Future value
double future_value(double F0, double t0, double t1, double r)
{
    double r_decimal = 0.01*r;
    double F1 = F0*exp(r_decimal*(t1-t0));
    return F1;
}

//1.2 Discount factor
double discount_factor(double F1, double t0, double t1, double r){
    double r_decimal = 0.01*r;
    double df = F1*exp((-1)*r_decimal*(t1-t0));
    return df;
}

int df_and_r(double F0, double F1, double t0, double t1, double & df, double & r)
{
    if (t1-t0 == 0.0) {
        df = 0;
        r = 0;
        return -1;
    }
    if ((F0 < 0.0) || (F1 < 0.0)) {
        // *** you figure it out ***
        df=0;
        r=0;
        return -2;
    }
    // *** you have to write the rest ***
    // calculate discount_factor value df
    df = discount_factor(F1, t0, t1, r);
    return 0;
}


//1.3 Bond Price from yield
void price_from_yield(double F, double c, double y, int n, double & B){
    B=0.0;
    double y_decimal = 0.01*y;
    for(int i=1; i<=n; i++){
        if(i!=F){
            B+= (0.5*c)/(pow(1+0.5*y_decimal, i));
        }else{
            B+=(F+0.5*c)/(pow(1+0.5*y_decimal, i));
        }
    }
}

string check_result(double a, double b){
    if(a==b){
        return "They are equal";
    }else{
        return "They are not equal, A is "+to_string(a)+" and B is "+to_string(b);
    }
}
//1.4 Yield from bond price
int yield_from_price(double F, double c, int n, double B_market, double tol, int max_iter, double & y)
{
    //step 1:
    int y_low=0;
    int y_heigh=100;
    return 0;
}


int main(int argc, const char * argv[]) {
    // insert code here...
    double F0, F1,t0, t1, r, df;
    double F, c, y, n, B;
    double temp;
    //set initial value random.
    F0=123;
    t0=23;
    t1=35;
    r=5;
    
    //call function, return value.
    F1=future_value(F0, t0, t1, r);
    cout<<"For 1.1: Future value, the F0, t0, t1. and r is given in code, reult F1 is "+to_string(F1)<<endl;
    cout<<"=====================================\n\n"<<endl;
    
    //calculated discount factor: do it in method.
    cout<<"For 1.2: discount factor, the check result is "+to_string(df_and_r(F0, F1, t0, t1, df, r))<<endl;
    cout<<"The value for df is "+to_string(df)<<endl;
    cout<<"=====================================\n\n"<<endl;
    
    //given initial value for B and y;
    B=0.0;
    F=100;
    y=10;
    c=10;
    n=100;
    //check B
    cout<<"For 1.3: "<<endl;
    price_from_yield(F, c, y, n, B);
    cout<<"Case 1: F=100, if y=c, result for B is "+to_string(B)<<endl;
    //set y=0 check it
    y=0;
    price_from_yield(F, c, y, n, B);
    cout<<"Case 2: checking if y=0"<<endl;
    temp=F+(n*c)/2;
    cout<<check_result(B, temp)<<endl;
    //set y=12 back, c=0
    y=10;
    c=0;
    price_from_yield(F, c, y, n, B);
    cout<<"Case 3: checking if c=0"<<endl;
    temp = F / pow((1+(0.5*y*0.01)), n);
    cout<<check_result(B, temp)<<endl;
    
    cout<<"======================================\n\n"<<endl;
    //set value of c back
    c=10;
    
    
    return 0;
}
