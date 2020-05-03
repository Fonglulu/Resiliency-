#include<iostream>

int main()
{
    const int value{5};
    const int *ptr = &value;
    *ptr = 6;
    return 0;    
}