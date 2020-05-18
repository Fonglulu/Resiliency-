
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
        Grid(unsigned n):_N{n}, _Ntot{n*n}, values{_Ntot, std::make_pair(1.0,1.0)}, residual(_Ntot, 0.0)
        {
        };
        // Deconstructor
        ~Grid()
        {
 
        }


        void print_nodes() //  
        {   
            //ref.values
            // print the component of values
            std::vector<std::pair<double, double>> &ref2{values};//refence to the variable nodes, can use auto?
            for (int i{0}; i<ref2.size(); i++){
                 std::cout <<ref2[i].first << ',' << ref2[i].second <<'\n';
                 }
        }

        // Convert vector index to grid index

        std::pair<unsigned int, unsigned int> grid_index(const unsigned int i) const{
            

            unsigned int ix{i/ _N};
            unsigned int iy{i% _N};
            //std::cout<<ix<<" ix "<< iy<< " iy "<<"\n";
            std::pair<unsigned int, unsigned int> indices;

            indices.first = ix;
            indices.second = iy;


            if (i>=  _Ntot){
                throw std::invalid_argument(" unbounded index");

            }

            return indices;
         }
         unsigned int vector_index(const unsigned int ix, const unsigned int iy) const {

             unsigned int index {_N*ix+iy};

             //std::cout<<index<<"index"<<'\n';

             return index;
         }




        // Get the values of the grid
         std::vector<std::pair<double, double>> get_values() const{

             return values;
         }


         void print_values(const std::vector<std::pair<double, double>> &valueref){

             for (int i=0; i<valueref.size(); i++){
                 std::cout <<valueref[i].first << ',' << valueref[i].second <<'\n';
             }
         }


        // fetch elemtent of index i

        std::pair<double, double> get_node(const unsigned int i) const{

            //  Get the data member values
            //std::vector<std::pair<double, double>> nodes = rgrid.values;

            // reference to the nodes
            const std::vector<std::pair<double, double>> &ref{values};
            

            // Get ith element
            std::pair<double, double> node_i = ref[i];

            return node_i;


        



        }  



        double get_approx(const unsigned int i)const{

            const std::vector<std::pair<double, double>> &ref{values};

            std::pair<double, double> node_i = ref[i];

            return node_i.first;

        }

        double get_rhs(const unsigned int i)const{

            const std::vector<std::pair<double, double>> &ref{values};

            std::pair<double, double> node_i = ref[i];

            return node_i.second;

        }

        unsigned int get_size() const {

        //std::cout<<_N<< " size"<<"\n";

        return _N;
        }


        double get_h() const{
            return spacing;
        }


        double get_res(const unsigned int i) const{


            const std::vector<double> &ref{residual};

            double res_i = ref[i];

            return res_i;

            


        }


         std::vector<double> get_residual() {

             return residual;
         }




        // Set values to nodes

        void set_approx( double& approx, const int & i){

            values[i].first = approx;


        }


        // Set rhs to nodes

        void set_rhs( double & rhs, const int & i){

            values[i].second =rhs;
        }


        void set_res(double & res, const int & i ){

            residual[i] = res;
        }

        //friend class MG;

        bool boundary(const unsigned int i) const {

            if (grid_index(i).first == 0 ||  grid_index(i).first   == _N-1  || grid_index(i).second == 0 || grid_index(i).second == _N-1)
            {
                return true; 
            }
            else
            {
                return false;
                
                
            }
            
        }


    
};

// MG is a class implementing relax, restrict, interpolate
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

            for (unsigned int i =0; i<l; i++){
                unsigned int number = pow(2, i+1)+1; 
                ptrgrids[i] = new Grid{number}; // dynamically allocated grid size to memory
            }
        }

        ~MG(){

            for(int i=0; i< Nlevel ;i++)

            delete ptrgrids[i]; 

            delete [] ptrgrids;

            // delete grids first.
        }




         void print_grid(unsigned int i){


             //boundcheck


             //double value1 = 2.0;
             //ptrgrids[i] ->set_approx(value1, 0);
             ptrgrids[i] ->print_nodes();

             //return ptrgrids[i]->get_node(i);
         }

        // Get ith grid from the Grid pointer array
         Grid* get_grid(unsigned i) {

             //const ptr/variable


             return ptrgrids[i];
         }

        // restriction on residuals
         void restriction(unsigned int level){

             // get current grid

             Grid* gridabove = get_grid(level);

             // get size of the gridabove

             unsigned int size_above{gridabove->get_size()};


             Grid* gridbelow = get_grid(level+1);

             unsigned int size_below{gridbelow -> get_size()};

             for (int i = 1; i< size_below-1; i++){

                for (int j = 1; j<size_below-1; j++){

                    

                    unsigned int index = gridabove->vector_index(i,j);

                    // Boundary check
                    if (gridbelow->boundary(index) == false){

                    // Find the (2x-1, 2y) node index in upper grid
                    int right_above = gridabove->vector_index(2* i-1, 2* j);

                    // Find the (2x+1, 2y) node index in upper grid
                    int left_above = gridabove->vector_index(2*i+1, 2* j);

                    // Find the (2x,2y-1) node index in upper grid
                    int down_above = gridabove->vector_index(2*i,   2* j);

                    // Find the (2x, 2y+1) node index in upper grid
                    int upper_above = gridabove->vector_index(2* i,   2* j+1);

                    // Find the (2x-1, 2y-1) node
                    int downleft_above = gridabove->vector_index(2*i-1, 2*j-1);
                
                    // Find the (2x-1, 2y+1) node
                    int upperleft_above = gridabove->vector_index(2*i-1, 2*j+1);

                    // Find the (2x+1, 2y-1) node
                    int downright_above = gridabove->vector_index(2*i+1, 2*j-1);

                    // Find the (2x+1, 2y+1) node
                    int upperright_above = gridabove->vector_index(2*i+1, 2*j+1);


                    double residual_below = gridabove->get_res(index)+ 0.5*(gridabove->get_res(right_above)+gridabove->get_res(left_above)\
                            +gridabove->get_res(upper_above) +gridabove->get_res(down_above)+)+0.25*(gridabove->get_res(upperleft_above)\
                            +gridabove->get_res(upperright_above)+gridabove->get_res(downleft_above)+gridabove->get_res(downright_above));


                    gridbelow->set_rhs(residual_below, gridbelow->vector_index(i,j));
                    }

                    else{
                        std::cout<<' Out of bound'<<'\n';
                    }


                    


                }
                
            }
            
        }// end of restriction


            void interpolation(unsigned int level){
                // Prolonge residual and correct the fine-grid approximation

                // Fetch the current grid, a coarse grid 
                Grid* gridbelow = get_grid(level);
                
                // Get the size of the corase grid
                unsigned int size_below{gridbelow->get_size()};

                // Fetch the above fine grid
                Grid* gridabove = get_grid(level+1);
                
                // Get the size of the fine grid
                unsigned int size_above{gridabove -> get_size()};

                for (int i = 0; i< size_below; i++){

                        for (int j = 0; j<size_below; j++){

                            unsigned int index = gridabove->vector_index(i,j);




                            
                            double error_above = gridbelow->get_res(index);


                            gridabove->set_res(error_above, gridabove->vector_index(2*i,2*j));
                    }

                }
                
                for (int i=0; i<size_above; i+=2){

                    for (int j = 1; i<size_above-1; i+=2){

                            // for k in range(depth):
                            //     for i in range(0, ynodes, 2):
                            //          for j in range(1, xnodes-1, 2):
                            //              uf[k,i,j]=0.5*(uf[k,i,j-1]+uf[k,i,j+1])

                            unsigned int index = gridabove->vector_index(i,j);

                            unsigned int index_below = gridabove->vector_index(i, j-1);

                            unsigned int index_above = gridabove-> vector_index(i,j+1);

                            double error_above = 0.5* (gridabove->get_res(index_above)+ gridabove ->get_res(index_below));

                            gridabove ->set_res(error_above, index);
                    }
                }

                for (int i =1; i<size_above-1; i+=2){
                    
                   for (int j = 0; i<size_above; i+=2){

                        unsigned int index = gridabove-> vector_index(i,j);

                        unsigned int index_left = gridabove ->vector_index(i-1, j);

                        unsigned int index_right = gridabove -> vector_index(i+1, j);

                        double error_above = 0.5* (gridabove->get_res(index_left) + gridabove -> get_res(index_right));

                        gridabove ->set_res(error_above, index);

                    }
                }

                for (int i =1; i< size_above -1; i+=2){
                    
                    for (int j= 1; i<size_above -1; i+=2){

                            unsigned int index = gridabove->vector_index(i,j);

                            unsigned int index_below = gridabove->vector_index(i, j-1);

                            unsigned int index_above = gridabove-> vector_index(i,j+1);

                            unsigned int index_left = gridabove ->vector_index(i-1, j);

                            unsigned int index_right = gridabove -> vector_index(i+1, j);

                            double error_above = 0.25*( gridabove ->get_res(index_below) +gridabove->get_res(index_above)/
                                                          +gridabove ->get_res(index_left) + gridabove->get_res(index_right));
                            gridabove ->set_res(error_above, index);
                    }
                }


            }



};

         





class MG_solver{

    private:
        

        unsigned int _N; // size per dim on finest 

        unsigned int _totN; // #nodes on finest

        unsigned int _Nlevel; // #levels 

        unsigned int _Nmin; // size per dim on coarsest // redundent info 


        // Data needed

        MG _solutions; 

        //solution

        // convergence criterion

        double _stopping_criterion{1e-5};

        // Residual information

        double _res;

        double _res_ini;

        //smoother parameter

        unsigned int _sweeps_coarse{2};

        unsigned int _sweeps_finest{2};

        unsigned int _max_Vcycle{10};

        //Book-keeping variables


        unsigned int _current_Vcycle;

        // Internal methods

        void _Gauss_Seidel(Grid * grid){  

            unsigned int xdim{grid->get_size()};

            double h ;
            h =  grid->get_h();

            for (int i = 0; i< xdim; i++){

                for (int j = 0; j<xdim; j++){

                    unsigned int index = grid->vector_index(i,j);

                    if (grid->boundary(index) == false){

                        unsigned int index_up = grid->vector_index(i,j+1);

                        unsigned int index_down = grid->vector_index(i,j-1);

                        unsigned int index_right = grid->vector_index(i+1,j);

                        unsigned int index_left = grid->vector_index(i-1, j);

                        double new_approx{(grid->get_approx(index_up)+grid->get_approx(index_down)+grid->get_approx(index_left)\
                                         +grid->get_approx(index_right)+h*h*grid->get_rhs(index))/int{4}};

                        // store a copy of solution
                        grid->set_approx(new_approx, index);


                    }
                }


            }



        }// End of _Gauss_Seidel


        void _calculate_new_res( const unsigned int level){

            Grid *grid = _solutions.get_grid(level);

            unsigned int xdim{grid->get_size()};

            double h{grid->get_h()};

            for (int i = 0; i< xdim; i++){

                for (int j = 0; j<xdim; j++){

                    unsigned int index = grid->vector_index(i,j);

                    if (grid->boundary(index) == false){


                        unsigned int index_up = grid->vector_index(i,j+1);

                        unsigned int index_down = grid->vector_index(i,j-1);

                        unsigned int index_right = grid->vector_index(i+1,j);

                        unsigned int index_left = grid->vector_index(i-1, j);

                        double new_res{grid->get_rhs(index) + (grid->get_approx(index_down) +grid->get_approx(index_up)\
                                      +grid->get_approx(index_left) + grid->get_approx(index_right) - 4*grid->get_approx(index))\
                                      /double(h*h)};

                        grid->set_res(new_res, index);

                        // store a copy of solution


                       }

                }
            }

        }

        void _solve_current_grid(const unsigned int level){

            unsigned int sweeps;

            Grid* grid_level = _solutions.get_grid(level); // Get the grid at current level

            if (level == 0){

                sweeps = _sweeps_finest;
            }

            else{
                sweeps = _sweeps_coarse;
            }

            for (int i = 0; i <sweeps; i++){

            _Gauss_Seidel(grid_level);// apply Gauss Seidel smoother on this level 

            }
        } // End of _solve_current_grid


        void _correction(const unsigned int level){

            // 
            
            
            Grid *grid_level = _solutions.get_grid(level); // get the grid pointer of given level

            unsigned int size = grid_level->get_size(); // get the size 


            // // get the residual on this grid

            // std::vector<double> interpolated_residual;  // get the residual vecctor 

            // interpolated_residual = grid_level->get_residual();

            for (int i=0; i<size; i++){

                double corrected_approx = grid_level->get_approx(i) +grid_level->get_res(i);

                grid_level->set_approx(corrected_approx, i);

            }

            


        }


        void _recursive_go_down(unsigned int level){

            // check if we are already at bottom 

            if (level >= _Nlevel){
                return; 
            } 

            //calculate new residual on grid below

             _solutions.restriction(level); // restriction on residual

            _solve_current_grid(level+1);

            _calculate_new_res(level+1);

            

            _recursive_go_down(level+1);

            
        }// End of _recursive_go_down


        void _recursive_go_up(unsigned int level){

            _solutions.interpolation(level); // prolongation on residual 

            // Add residual on cuurent level to the solution current level

            _correction(level-1);




            _solve_current_grid(level - 1);

            _recursive_go_down(level -1);

             
            
        }


        void _check_convergence(){


        }









        public:

        void solve(){

             _current_Vcycle = 0;



             _solve_current_grid(0);\


             if (_check_convergence()) return; 

             while(1) {
            
            // Goes down from the finest grid [0] to the coarest
             _recursive_go_down(0);

             _solve_current_grid(_Nlevel-1);


            // Goes up from the coarest grid  _Nlevel-1 
             _recursive_go_up(_Nlevel-1);


             if (_check_convergence()) break;}



             



        }





        
        
        







        
};













int main(){
    //Grid grid{5};
    //grid.get_size();
    //auto values = grid.get_values();
    //grid.print_values(values);
    //grid.grid_index(9);
    //grid.vector_index(1,4);
    //nodes = grid.values;
    //grid.print_nodes();
    //std::cout<<"Now assign values to node"<<'\n';
    // Assign 2.0 to the first node
    //double v11 = 2.0;
    //grid.set_approx(v11 , 0);
    //grid.print_nodes();
    MG Vcycle{3};

    std::cout<<"Print the 3*3 grid"<<"\n";

    Vcycle.print_grid(0);

    Grid* grid2 = Vcycle.get_grid(2);

    std::cout<<"The size of grid"<<"\n";

    grid2->get_size();

    //Vcycle.restriction(*grid2);

    Grid gridbelow{Vcycle.restriction(*grid2)}; 

    gridbelow.print_nodes();

    
    //Vcycle.get_values(2);
    return 0;
}
