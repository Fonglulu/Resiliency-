#include <iostream>

    int *returnInteger()
{
  int v =5;
  int *ptr = &v;
  //delete ptr;
  return ptr; // can't return &v
}


// int *returnArray(){

//     int array[1]{1};

//     return array;
// }



int *returnArray_2(){
    int *array = new int[1]{1};

    return array;
}
 
int main()
{

    int *ptr{new int};

    *ptr = 7;

    std::cout<< *ptr<<'\n';
    //delete ptr;
    //ptr = 0;
    //std::cout<<*ptr<<'\n';
    int* a = returnArray_2();
    std::cout << a <<'\n';


 
    return 0;
}