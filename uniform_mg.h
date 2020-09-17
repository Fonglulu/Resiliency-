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
        double get_rhs(const unsigned int i) const;

        // get the residual stored on a given node
        double get_res(const unsigned int i) const;

        // get the size of the grid per dim
        unsigned int get_size() const;

        // get the spacinf of the grid
        double get_h() const;

        // set new approximation on a given node
        void set_approx( double & approx, const int &i);

        // set new rhs on a given node
        void set_rhs( double & rhs, const int & i);

        // set new residual on a given node
        void set_res(double & res, const int & i);
        
        // reset approx to zeros.
        void reset()

        // Given a node index in 1D vector tell if it is on bnd.
        bool boundary (const unsigned int i) const;

        // calculate the l2 norm of residual across the board
        double l2_res();


}


class MG{

    private:

        unsigned int Nlevel;

        Grid **ptrgrids; // array/vector of pointers to grids via doubler pointer

        std::vector<std::pair<double, double>> values;

    
    public:

        MG(unsigned int l) //l is number of levels
        {
            Nlevel = l;
            ptrgrids = new Grid*[l]; //array of ptrs
            unsigned int count = 0;
            for ( int i =l-1; i>-1; --i){
                std::cout<< "l "<< l<< '\n';
                std::cout<< "i "<<i<< '\n';
                unsigned int number = pow(2, i+1)+1; 
                ptrgrids[count] = new Grid{number}; // dynamically allocated grid size to memory
                std::cout<<"count "<< count<<'\n';
                ++count;
            }
        }

        ~MG(){

            for(int i=0; i< Nlevel ;i++)

            delete ptrgrids[i]; 

            delete [] ptrgrids;

            // delete grids first.
        }

        // print the total number of grids
        void print_level_numbe();
        
        // get the ith grid. lvl 0 is the coarest grid.
        Grid* get_grid(unsigned i);

        // full weight restriction on residual
        void restriction(unsigned int level);

        // full weight interpolation. 
        void interpolation(unsigned int level);
        

