#include<iostream>
#include<string>
class Employee{
    private:
       int m_id;
       std::string m_name;

    public:

        Employee(int id = 0, const std::string &name= ""):
            m_id{id}, m_name{name}
        {  
            std::cout<<"Employee "<<m_name<<"created.\n";    

        } 

        Employee(const std::string &name):Employee(0,""){};

};


int main(){
    Employee emp1('Jack');
    return 0;
}