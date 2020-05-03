#include <vector>
#include<iostream>

using namespace std;

void print(int a, int * ipa) {
    cout<< a<<" "<< (long)ipa<<" " <<(long)&a<<endl;
}

int main()
// {
//     vector<int> vect{10,20,30,40};

//     for (int x : vect){
//          x = x+5;

//         cout<< x <<endl;
// }
    

//      for (int x : vect)

//         cout<< x<<" ";

//     return 0;
// }
{
int a =46;
int  ipa = &a;
cout<< a<<" "<< (long)ipa<<" " <<(long)&a<<endl;
//print(a,ipa);

}

