#include<iostream>
#include<typeinfo>

int main(){
    int qrray[]{1,1,12,3,4,5}; 
    int v{5};
    int *test = new int[7];
    qrray[2] =5 ;
    char carray{'H'};
    int* ptr{qrray};
    int number{5};

    int *ptr1;
    ptr1 = &number;
    std:: cout<<&ptr1<<" "<<'\n';
    std:: cout<<typeid(ptr).name()<<'\n';
    std::cout<< "addres of the pointer "<< &ptr<< '\n';
    std::cout<< "element 0 holds "<< qrray[2]<<std::endl;
    std::cout<< "The array decays to a pointer holding the address "<< *qrray<<'\n';
    std::cout<< &carray<<'\n';
    return 0;
}