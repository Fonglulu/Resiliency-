
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
        {   std::cout<< "print_nodes "<<'\n';
            //ref.values
            // print the component of values
            std::vector<std::pair<double, double>> &ref2{values};//refence to the variable nodes, can use auto?
            for (int i{0}; i<ref2.size(); i++){
                 std::cout <<ref2[i].first << "  " << ref2[i].second <<'\n';
                 }
        }

        void print_int()
        {

            std::vector<std::pair<double, double>> &ref2{values};

            for (int i{0}; i<ref2.size(); i++){

                if (boundary(i) == false){
                    
                    std::cout<<ref2[i].first<< "  "<< ref2[i].second<<'\n';
                    
                }
            }


        }

        void print_res()

        {
             std::vector<double>  & _ref_res{residual};

             for (int i = 0; i< _ref_res.size(); i++){

                 std:: cout<<_ref_res[i]<<'\n';
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

        std::cout<<_N<< " size"<<"\n";

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
            unsigned int count = 0;
            for ( int i =l-1; i>-1; --i){
                std::cout<<"l"<< l<< '\n';
                std::cout<< "i"<<i<< '\n';
                unsigned int number = pow(2, i+1)+1; 
                ptrgrids[count] = new Grid{number}; // dynamically allocated grid size to memory
                ++count;
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

                    std::cout<< gridbelow->get_size()<<"gridbelowsize"<<'\n';
                    
                    if (gridabove->boundary(index) == false){   //If it's not a boundary point on the finer grid.

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
                            +gridabove->get_res(upper_above) +gridabove->get_res(down_above))+0.25*(gridabove->get_res(upperleft_above)\
                            +gridabove->get_res(upperright_above)+gridabove->get_res(downleft_above)+gridabove->get_res(downright_above));


                    gridbelow->set_rhs(residual_below, gridbelow->vector_index(i,j));
                    }

                    else{
                        std::cout<<"Out of bound"<<'\n';
                    }


                    


                }
                
            }
            
        }// end of restriction


            void interpolation(unsigned int level){
                // Prolonge residual and correct the fine-grid approximation

                // Fetch the current grid, a coarse grid 
              
                Grid* gridbelow = get_grid(level);
               std::cout<<"get gridbelow"<<" "<<level<<'\n';
                // Get the size of the corase grid
                unsigned int size_below{gridbelow->get_size()};

                // Fetch the above fine grid
                Grid* gridabove = get_grid(level-1);

                std::cout<<"get gridabove"<<'\n';
                std::cout<<  level+1<< '\n';
                gridabove -> get_size();
                std::cout<< "interpolate"<<'\n';
                // Get the size of the fine grid
                unsigned int size_above = gridabove -> get_size();
                std::cout<< "interpolate"<<'\n';
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


        double _h{1.0/(_N -1)};



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
            // Get the size of current grid
            unsigned int xdim{grid->get_size()};

            double h ;
            h =  grid->get_h();

            // take a copy of grid

            Grid * _copy_grid = grid;

            std::cout<<*(&grid)<< *(&_copy_grid)<<"copy"<<'\n';

            
            // Loop through each node on board
            for (int i = 1; i< xdim-1; i++){

                for (int j = 1; j<xdim-1; j++){
                    
                    unsigned int index = grid->vector_index(i,j);
                    // For interior grid
                    if (grid->boundary(index) == false){

                        unsigned int index_up = grid->vector_index(i+1,j);
                        std::cout<<index<<" "<< index_up << " " <<grid->get_rhs(index_up)<< "node above"<< '\n';

                        unsigned int index_down = grid->vector_index(i-1,j);

                        unsigned int index_right = grid->vector_index(i,j+1);

                        unsigned int index_left = grid->vector_index(i, j-1);

                        std::cout<< "rhs "<< grid->get_rhs(index)<<'\n';

                        double new_approx = (_copy_grid->get_approx(index_up)+_copy_grid->get_approx(index_down)+_copy_grid->get_approx(index_left)\
                                         + _copy_grid->get_approx(index_right)+h*h*_copy_grid->get_rhs(index))/4 ;

                        // store a copy of solution

                        std::cout<< new_approx<< " new approx"<<'\n';
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
                        std::cout<< new_res<< " res"<< '\n';

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

            grid_level->print_int();



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

            if (level >= _Nlevel-1){
                return; 
            } 

            //calculate new residual on grid below

             _solutions.restriction(level); // restriction on residual

            _solve_current_grid(level+1);

            _calculate_new_res(level+1);         

            _recursive_go_down(level+1);

            
        }// End of _recursive_go_down


        void _recursive_go_up(unsigned int level){

            std::cout<<"starting going up"<<" level " << level<<'\n';


            _solutions.interpolation(level); // prolongation on residual 
        
            // Add residual on cuurent level to the solution current level

            _correction(level-1);
            std::cout<<"Back to finest "<< level -1<<'\n';
            _solve_current_grid(level - 1);


            if (level  >1){

                std::cout<<"level "<< level<<'\n';
                _recursive_go_up(level -1);
            } 
            else{
                return;}
            }


        bool _check_convergence(){

            double err = _res_ini != 0.0 ? _res/_res_ini : 1.0;

            bool converged = false;

            // if (_res <_stopping_criterion){

            // converged = true;}

            if (_current_Vcycle >= _max_Vcycle ){

                converged = true;
            }

            return converged;




        }










        public:


        

        // member initialiser list for Constructors
        MG_solver(unsigned int n):_N {n},_totN{n*n},_Nlevel{(unsigned int) log2(n-1)}, _solutions{_Nlevel}{}

        



        void solve(){

            std::cout<<"solve rountine"<<'\n';

             _current_Vcycle = 0;



             _solve_current_grid(0);

             _calculate_new_res(0);

            std::cout<< "Finised solving on finest grid"<<'\n';
             if (_check_convergence()) return; 

             while(true) {
            
            // Goes down from the finest grid [0] to the coarest
             _recursive_go_down(0);

             _solve_current_grid(_Nlevel-1);
             std::cout<< "solve on the bottom level"<< " "<< _Nlevel-2<<'\n';


            // Goes up from the coarest grid  _Nlevel-2

             _recursive_go_up(_Nlevel-1);

             _current_Vcycle++;


             if (_check_convergence()) break;}

             Grid* _finest = _solutions.get_grid(0);

            std::cout<< " now prints the solutions"<< '\n';
             _finest->print_nodes();
             _finest->print_res();



            


             



        }

            void _set_boundary(){

                Grid* _finest = _solutions.get_grid(0);

                _finest->get_size();

                for (int i = 0; i< _N; i++){

                    for (int j = 0; j<_N; j++){

                        unsigned int index = _finest->vector_index(i,j);

                        if(_finest->boundary(index) == true){

                            double _boundary = exp(i*_h)*exp(j*_h);
                            _finest ->set_approx(_boundary, index);

                        }
                    }
                }

            } // end of _set_boundary


            void _make_source(){

            // Get the rhs from finest grid in _solution

            Grid* _finest = _solutions.get_grid(0);

            _finest->get_size();

            for (int i= 0; i<_N;i++){

                for (int j=0; j<_N; j++){

                    double source = -2*exp(i*_h)*exp(j*_h);

                    unsigned int index = _finest->vector_index(i,j);

                    _finest->set_rhs(source, index); 

                   
                   


                }


            }

            //_finest->print_nodes();
        }











        
        
        







        
};


int main(){


    MG_solver Laplacian(5);
    Laplacian._make_source();
    Laplacian._set_boundary();
    Laplacian.solve();
    

}













// int main(){
//     //Grid grid{5};
//     //grid.get_size();
//     //auto values = grid.get_values();
//     //grid.print_values(values);
//     //grid.grid_index(9);
//     //grid.vector_index(1,4);
//     //nodes = grid.values;
//     //grid.print_nodes();
//     //std::cout<<"Now assign values to node"<<'\n';
//     // Assign 2.0 to the first node
//     //double v11 = 2.0;
//     //grid.set_approx(v11 , 0);
//     //grid.print_nodes();
//     MG Vcycle{3};

//     std::cout<<"Print the 3*3 grid"<<"\n";

//     Vcycle.print_grid(0);

//     Grid* grid2 = Vcycle.get_grid(2);

//     std::cout<<"The size of grid"<<"\n";

//     grid2->get_size();

//     //Vcycle.restriction(*grid2);

//     Grid gridbelow{Vcycle.restriction(*grid2)}; 

//     gridbelow.print_nodes();

    
//     //Vcycle.get_values(2);
//     return 0;
// }
