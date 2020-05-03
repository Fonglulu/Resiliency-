#include<iostream>

class Simple
{
    private:
        int m_id; // private member data
    public:

        Simple(int id): m_id{1} 
        {

        } // ddaefacult constructor: to initialise a class with private member variables.
        

        void setid(int id){
            m_id =id;
        }

        int getid(){
            return m_id;
        }

};

using namespace std;

int main()
{
    Simple simple(3); // call Simple(int 1)
    std::cout <<simple.getid()<<endl;
    simple.setid(22);
    std::cout<<simple.getid()<<endl;

    return 0;
}

