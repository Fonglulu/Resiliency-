#include <iostream>
int main()
{
    int sum =0, val =1;
    while (val <= 10) {
        sum += val;
        val +=1;
    }

    std::cout <<sum<<std::endl; 

    return sum;
}