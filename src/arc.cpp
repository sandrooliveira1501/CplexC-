#include <iostream>
#include "arc.h"
using namespace std;

Arc::Arc(){
    this->a = 0;
    this->b = 0;
}
Arc::Arc(int a,int b){
    this->a = a;
    this->b = b;
}

bool Arc::Equals(Arc arc){

    return arc.a == this->a && arc.b == this->b;
}
