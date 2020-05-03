#include <iostream>
 
int main()
{
    char value{ 'a' };
    char *ptr{ &value };
 
    std::cout << ptr << '\n';
    std::cout << ptr+1 << '\n';
    std::cout << ptr+2 << '\n';
    std::cout << ptr+3 << '\n';
    std::cout<< *ptr<< " "<< *(ptr+3)<<'\n';
 
    return 0;
}