#include<iostream>
#include<vector>

class Grid{
    private:
        unsigned int N;
        std::vector<double> nodes; //declare an empty vector

    public:


        Grid(){
            N=0
        }

        


        Grid(int n): N(n), nodes(N*N,0.0) //  member initialiser list constructor  
        {
          
        }

    
    void print(){
        std::cout<< N<< std::endl;
        //std::cout<<nodes<<std::endl;
    }
};

int main(){
    Grid grid1;
    grid1.print();
    Grid grid2{5};
    grid2.print();

    return 0;
}