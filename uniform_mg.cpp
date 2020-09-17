#include "uniform_grid.h"
#include<iostream>


Grid::Grid(unsigned n):_N{n}, _Ntot{n*n}, values{_Ntot, std::make_pair(0.0,0.0)}, residual(_Ntot, 0.0)
        {
        }

Grid::~Grid()
        {
        }

void Grid::print_nodes(){
    std::cout<< "print_nodes "<<'\n';
            //ref.values
            // print the component of values
            std::vector<std::pair<double, double>> &ref2{values};//refence to the variable nodes, can use auto?
            for (int i{0}; i<ref2.size(); i++){
                 std::cout <<ref2[i].first << "  " << ref2[i].second <<'\n';
                 }
}


void Grid::print_int(){
    std::vector<std::pair<double, double>> &ref2{values};

            for (int i{0}; i<ref2.size(); i++){

                if (boundary(i) == false){
                    
                    std::cout<<ref2[i].first<< "  "<< ref2[i].second<<'\n';
                    
                }
            }
}