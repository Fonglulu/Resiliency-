#include <iostream>
 
// parameter ptr contains a copy of the array's address
void changeArray(int *ptr)
{
    *ptr = 23; // so changing an array element changes the _actual_ array
}
 
int main()
{
    int array[]{ 1, 1, 2, 3, 5, 8, 13, 21 };
    std::cout << "Element 0 has value: " << array[0] << '\n';
 
    changeArray(array);
 
    std::cout << "Element 0 has value: " << array[0] << '\n';
 
     return 0;
}