#include <iostream>
 
int main()
{
    int *value{5};
    //std::cout<< *value<<std::endl; 

    if (value)
        std::cout<<"non-null pointer";

    else
        std::cout<<"null pointer";
    
    return 0;
}