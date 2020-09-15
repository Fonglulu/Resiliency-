#ifndef UNIFORM_MG.H
#define UNIFORM_MG.H


#include<iostream>
#include<vector>
#include<utility>
#include <iterator>
#include<cmath>
#include<stdexcept>

class Grid{

        private:
        // nodes per dimension
        unsigned int _N;
        // nodes in the grid
        unsigned int _Ntot;
        //1st component stores the rhs, 2nd stores the approx. values
        std::vector<std::pair<double, double>> values ;

        std::vector<double> residual; 

         double spacing{1.0/(_N -1)};

    
    public:
        //Intialise private data

        // Constructor 
        Grid(unsigned n):_N{n}, _Ntot{n*n}, values{_Ntot, std::make_pair(0.0,0.0)}, residual(_Ntot, 0.0)
        {
        };
        // Deconstructor
        ~Grid()
        {
        }
        // print the residual and approximation of nodes in the grid
        void print_nodes();
        
        // print the residual and approximation of nodes in the interior grid
        void print_int();
        
        // print the residual of nodes in the grid
        void print_res();

        // print the residual and approximation of node from values vector
        void print_values(const std::vector<std::pair<double, double>> &valueref);

        // convert a node index in 1D vector to a tuple in 2D grid
        std::pair<unsigned int, unsigned int> grid_index(const unsigned int i) const ;

        // convert a node tuple index in 2D grid to a index in 1D vector
        unsigned int vector_index(const unsigned int ix, const unsigned int iy) const;

        // get the values on the grid
        std::vector<std::pair<double, double>> get_values() const;

        // 
        std::pair<double, double> get_node(const unsigned int i) const;

        // get the approximation for a given node
        double get_approx(const unsigned int i)const;

        // get the rhs stored on a given node
        double get_rhs(const unsigned int i)const;




        





}