//
//  main.cpp
//  exam_bonus
//
//  Created by Eric Lin on 9/17/17.
//  Copyright © 2017 XingSong  Lin. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;

class Bond
{
    
private:
    double face; //F
    double maturity; //date
    vector<double> coupon;
    vector<double> coupon_date;
    Bond(const Bond &); // don’t write a copy constructor
    Bond & operator= (const Bond &); // don’t write an assignment operator
public:
    Bond(double f, double T){
        face = f;
        maturity = T;
    }
    
    ~Bond(){
        //no pointer created
    }
    
    void set_flows(int n, const double c[], const double d[]){
        for(int i=0; i<n; i++){
            coupon.push_back(c[i]);
            coupon_date.push_back(d[i]);
        }
        
    }
    
    void price(double yield, double t0, double & B) const{
        //initial
        B=0.0;
        yield = yield*0.01;
        double power, value1, value2;
        for(int i=0; i<coupon.size(); i++){
            power = 2*(coupon_date.at(i)-t0);
            value1 = (1+0.5*yield);
            value2 = 0.5*coupon.at(i);
            B+= value2 /pow(value1, power);
        }
        B+= face/pow((1+0.5*yield), (2*(maturity-t0)));
    }
    
    int yield(double target, double t0, double tol, int max_iter, double & y, int & num_iter) const{
        
        double B, B_y_low, B_y_high, y_low, y_high;
        y_low=0;
        price(y_low, t0, B_y_low);
        if(abs(B_y_low-target) <= tol){
            y=y_low;
            cout<<y_low<<endl;
            return y;
        }
        
        y_high=100;
        price(y_high, t0, B_y_high);
        if(abs(B_y_high-target)<=tol){
            y=y_high;
            cout<<y_high<<endl;
            return y;
        }
        cout<<B_y_low<<endl;
        cout<<B_y_high<<endl;
        
        if(B_y_low< target || B_y_high>target){
            cout<<"error"<<endl;
            return y;
        }
        
        //step 1-8 end.
        
        for(int i=num_iter; i<max_iter; ++i){
            y=(y_low+y_high)/2.0;
            price(y, t0, B);
            
            if(abs(B-target)<=tol){
                cout<<"1"<<endl;
                return y;
            }else if (B>target){
                y_low=y;
            }else{
                y_high=y;
            }
            
            if((y_high-y_low)<=tol){
                cout<<"2"<<endl;
                return y;
            }
            
            if(i==(max_iter-1)){
                cout<<"3"<<endl;
                return y;
            }
        }
        return y;
    }
    
};


int main(int argc, const char * argv[]) {
    
    
    Bond Bond(100,10);
    
    int n=3;
    double c[3] = {3,2,1};
    double d[3] = {4,6,7};
    double y=3;
    double t0=3;
    double B=0;
    int num_iterator =0;
    Bond.set_flows(n, c, d);
    const class Bond& object = Bond;
    object.price(y, t0, B);
    cout<<"Bond is "<<B<<endl;
    cout<<"++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    
    object.yield(102, t0, 0.01, 50, y, num_iterator);
    cout<<"yield is "<<y<<endl;
    
    
    return 0;
}
