//
//  main.cpp
//  Hw2
//
//  Created by Eric Lin on 9/10/17.
//  Copyright © 2017 XingSong  Lin. All rights reserved.
//

#include <iostream>
#include <string>
#include <cmath>
using namespace std;

//2.1: Bond price from yield

void price_from_yield(double F, double c, double y, int n, double & B, double & D_mac, double & D_mod){
    double y_decimal = 0.01*y;
    B=0.0;
    D_mac=0.0;
    
    double temp1 = 1.0/(1+0.5*y_decimal);
    double temp2 = 1.0;
    
    for(int i=1; i<=n; i++){
        
        temp2*=temp1;
        
        if(i!=n){
            B+= (0.5*c)/temp2;
            D_mac+= (i/2)*((0.5*c)/temp2);
            D_mod=D_mac/(1+0.5*y);
        }else{
            B+=(F+0.5*c)/temp2;
            D_mac+= (i/2)*((F+0.5*c)/temp2);
            D_mac=D_mac/B;
            D_mod=D_mac/(1+0.5*y);
        }
    }
}


int main(int argc, const char * argv[]) {
    // insert code here...
    double F=100;
    double y=10;
    int n=100;
    double c=0;
    double B, D_mac, D_mod;
    //2.1 zero coupon boud test, c=0;
    price_from_yield(F, c, y, n, B, D_mac, D_mod);
    cout<<B<<endl;
    cout<<D_mac<<endl;
    cout<<D_mod<<endl;
    
    
    return 0;
}
