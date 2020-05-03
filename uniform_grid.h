#ifndef _UNIFORM_GRID
#define _UNIFORM_GRID
#include<iostream>
#include<vector>
#include<utility>


class Grid{

    private:
        // nodes per dimension
        unsigned int _N;
        // nodes in the grid
        unsigned int _Ntot;
        //1st component stores the rhs, 2nd stores the approx. values
        std::vector<pair<double, double>> values ; 

    
    public:
        //intialise private data

        // Set 0 node
        Grid(unsigned n):_N{n}, _Ntot{n*n}, values(n*n, (0.0,0.0))
        {

        };




        // Grid(): Grid(0,0){};
        // Grid(unsigned int N) : Grid (N, 0.0) {};
        // Grid(unsigned int N, yini);
        
        // //Get the pointer of each node
        // pair<double, double> * get_node()
 
        // //Assign value in the grid
        // void set_value(std::vector<pair<double, double>> &nodes);;
        // void set_value(unsigned int i, pair<double, double> &value );

        // //Declare Grid_index -> Coordinate function

        // std::vector<unsigned int> index_list(unsigned index i);

        // // Get info from the grid

        // unsigned int get_N();
        // unsigned int get_Ntot();

        // // Convert coordinate  -> index
        // unsigned int get_index_2d(unsigned int ix, unsigned int iy);

        








     

























}