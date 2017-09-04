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
    double B, B_y_low, B_y_high, y_low, y_high;
    y_low=0; y_high=100;
    
    for(int i=0; i<max_iter; ++i){
        y=(y_low+y_high)/2.0;
        price_from_yield(F,c, y_low, n, B_y_low);
        price_from_yield(F, c, y_high, n, B_y_high);
        price_from_yield(F, c, y, n, B);
        
        if((y_high-y_low)<=tol){
            return 0;
        }
        
        if(abs(B-B_market)<=tol){
            return 0;
        }else if (B>B_market && B_market<= tol){
            y_low=y;
        }else{
            y_high=y;
        }
        
        if(abs(B_y_low-B_market) <= tol){
            y=y_low;
            return 0;
        }
        
        if(abs(B_y_high-B_market)<=tol){
            y=y_high;
            return 0;
        }
        
        if(B_y_low< B_market || B_y_high>B_market){
            y=0;
            return 1;
        }
        if(i==max_iter){
            y=0;
            return 1;
        }
    }
    return 0;
}


int main(int argc, const char * argv[]) {
    // insert code here...
    double F0, F1,t0, t1, r, df;
    double F, c, y, n, B;
    double temp;
    //set initial value random.
    F0=123; t0=23; t1=35; r=5;
    
    //call function, return value.
    F1=future_value(F0, t0, t1, r);
    cout<<"For 1.1: Future value, the F0, t0, t1. and r is given in code, result F1 is "+to_string(F1)<<endl;
    cout<<"=====================================\n\n"<<endl;
    
    //calculated discount factor: do it in method.
    cout<<"For 1.2: discount factor, the check result is "+to_string(df_and_r(F0, F1, t0, t1, df, r))<<endl;
    cout<<"The value for df is "+to_string(df)<<endl;
    cout<<"=====================================\n\n"<<endl;
    
    //given initial value for B and y;
    F=100; y=10; c=10; n=100;
    //check B
    cout<<"For 1.3: "<<endl;
    price_from_yield(F, c, y, n, B);
    cout<<"Case 1: F=100, if y=c, result for B is "+to_string(B)<<endl;
    //set y=0 check it
    y=0; B=0;
    price_from_yield(F, c, y, n, B);
    cout<<"Case 2: checking if y=0, result: "<<endl;
    temp=F+(n*c)/2;
    cout<<check_result(B, temp)<<endl;
    //set y=12 back, c=0
    y=10; c=0; B=0;
    price_from_yield(F, c, y, n, B);
    cout<<"Case 3: checking if c=0, result: "<<endl;
    temp = F / pow((1+(0.5*y*0.01)), n);
    cout<<check_result(B, temp)<<endl;
    
    cout<<"======================================\n\n"<<endl;
    cout<<"1.4: Yield from bond price"<<endl;
    //set value of c back
    c=10; F=100; y=10; n=100; F0=123; t0=0; t1=15; r=5;
    //continues steps in function yield_from_price, for checking
    //step 20: test case
    
    //step 21:
    double B_market =100; double tol =25; int max_iter=100;
    cout<<"test case for step 21: "<<endl;
    cout<<yield_from_price(F, c, n, B_market, tol, max_iter, y)<<endl;
    cout<<"y = "<<to_string(y)<<" c = "<<to_string(c)<<endl; //error
    
    B=0; y=10;
    price_from_yield(F, c, y, n, B);
    cout<<"check F=100, set y==c, B is "<<to_string(B)<<endl;
    
    //step 22:
    B_market=50;
    y=10;
    cout<<"step 22:\ncase 1: B_market <100 "<<endl;
    cout<<yield_from_price(F, c, n, B_market, tol, max_iter, y)<<endl;
    if(y>c){
        cout<<"Correct, y>c, where y="+ to_string(y)+" c="+to_string(c)<<endl;
    }else{
        cout<<"error, y<c where y="+to_string(y)+ " c="+to_string(c)<<endl;
    }
    //reset
    y=10;
    B_market=200;
    cout<<"case 2: B_market>100"<<endl;
    cout<<yield_from_price(F, c, n, B_market, tol, max_iter, y)<<endl;
    if(y<c){
        cout<<"correct, y<c, where y="+ to_string(y)+" c="+to_string(c)<<endl;
    }else{
        cout<<"error, y>c where y="+to_string(y)+ " c="+to_string(c)<<endl;
    }
    //step 23:
    y=10;
    B_market=150;
    c=0;
    temp = F / pow((1+(0.5*y*0.01)), n);
    cout<<"Step 23 test case:"<<endl;
    cout<<yield_from_price(F, c, n, B_market, tol, max_iter, y)<<endl;
    cout<<"Zero coupon bond case check: "<<endl;
    cout<<"y is "+to_string(y)<<endl;
    return 0;
}
