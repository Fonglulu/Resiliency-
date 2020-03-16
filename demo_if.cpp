#include <iostream>
int main()
{
    std::cout <<"Enter two numbers:"<< std::endl;
    int v1, v2;
    std::cin  >> v1>>v2;
    int lower, upper;
    if (v1 <= v2){
        lower = v1;
        upper = v2;
    } 
    else {lower = v2;
          upper = v1;
    }
     
    std:: cout << lower<< upper<< std::endl;
    return 0;
    
}