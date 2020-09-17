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

void Grid::print_res(){
    std::vector<double>  & _ref_res{residual};

             for (int i = 0; i< _ref_res.size(); i++){

                 std:: cout<<_ref_res[i]<<" residual at "<< i<<'\n';
             }
}

 std::pair<unsigned int, unsigned int> 
 Grid::grid_index(const unsigned int i) const{

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


unsigned int 
Grid::vector_index(const unsigned int ix, const unsigned int iy) const{

             unsigned int index {_N*ix+iy};

             return index;
         }


std::vector<std::pair<double, double>> 
Grid::get_values() const{

             return values;
         }

void Grid::print_values(const std::vector<std::pair<double, double>> &valueref){

             for (int i=0; i<valueref.size(); i++){
                 std::cout <<valueref[i].first << ',' << valueref[i].second <<'\n';
             }
         }


std::pair<double, double>  Grid::get_node(const unsigned int i) const{

            //  Get the data member values
            //std::vector<std::pair<double, double>> nodes = rgrid.values;

            // reference to the nodes
            const std::vector<std::pair<double, double>> &ref{values};
            

            // Get ith element
            std::pair<double, double> node_i = ref[i];

            return node_i;

        }  

double Grid::get_approx(const unsigned int i)const{

            const std::vector<std::pair<double, double>> &ref{values};

            std::pair<double, double> node_i = ref[i];

            return node_i.first;

        }

double Grid::get_rhs(const unsigned int i)const{

            const std::vector<std::pair<double, double>> &ref{values};

            std::pair<double, double> node_i = ref[i];

            return node_i.second;

        }

unsigned int Grid::get_size() const {

        //std::cout<<_N<< " size"<<"\n";

        return _N;
        }


double Grid::get_h() const{
            return spacing;
        }


double Grid::get_res(const unsigned int i) const{


            const std::vector<double> &ref{residual};

            double res_i = ref[i];

            return res_i;
        }



std::vector<double> Grid::get_residual() {

             return residual;
         }


void Grid::set_approx( double& approx, const int & i){

            values[i].first = approx;

        }

void Grid::reset(){

            for(int i =0; i<_Ntot;i++){

                double zero = 0.0;

                set_approx(zero, i);
            }
        }


void Grid::set_rhs( double & rhs, const int & i){

            values[i].second =rhs;
        }

void Grid::set_res(double & res, const int & i ){

            residual[i] = res;
        }

bool Grid::boundary(const unsigned int i) const {

            if (grid_index(i).first == 0 ||  grid_index(i).first   == _N-1  || grid_index(i).second == 0 || grid_index(i).second == _N-1)
            {
                return true; 
            }
            else
            {
                return false; 
            }
            
        }


double Grid::l2_res() {

            //int size = vector.size();

            double sum = 0;

            for(int i = 0; i<_Ntot; i++){


                sum += pow(residual[i],2);
            }
         
            return sqrt(sum);

            }


MG::MG(unsigned int l) //l is number of levels
        {
            Nlevel = l;
            ptrgrids = new Grid*[l]; //array of ptrs
            unsigned int count = 0;
            for ( int i =l-1; i>-1; --i){
                std::cout<<"l"<< l<< '\n';
                std::cout<< "i"<<i<< '\n';
                unsigned int number = pow(2, i+1)+1; 
                ptrgrids[count] = new Grid{number}; // dynamically allocated grid size to memory
                std::cout<<"count "<< count<<'\n';
                ++count;
            }
        }


MG::~MG(){

            for(int i=0; i< Nlevel ;i++)

            delete ptrgrids[i]; 

            delete [] ptrgrids;

            // delete grids first.
        }

void MG::print_level_number(){

             std::cout<<"levels of Vcycle " <<Nlevel<<'\n'; 
         }


Grid* MG::get_grid(unsigned i) {

             //const ptr/variable


             return ptrgrids[i];
         }




        // restriction on residuals
void MG::restriction(unsigned int level){

             // get current grid

             Grid* gridabove = get_grid(level);

             // get size of the gridabove

             unsigned int size_above{gridabove->get_size()};


             Grid* gridbelow = get_grid(level+1);

             unsigned int size_below{gridbelow -> get_size()};

             for (int i = 1; i< size_below-1; i++){

                for (int j = 1; j<size_below-1; j++){

                    

                    unsigned int index = gridabove->vector_index(2*i,2*j);

                

                    // Boundary check

                    //std::cout<< gridbelow->get_size()<<" gridbelowsize"<<'\n';
                    
                    if (gridabove->boundary(index) == false){   //If it's not a boundary point on the finer grid.

                    // Find the (2x-1, 2y) node index in upper grid
                    int right_above = gridabove->vector_index(2* i-1, 2* j);

                    // Find the (2x+1, 2y) node index in upper grid
                    int left_above = gridabove->vector_index(2*i+1, 2* j);

                    // Find the (2x,2y-1) node index in upper grid
                    int down_above = gridabove->vector_index(2*i,   2* j-1);

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


                    double residual_below =4* gridabove->get_res(index)+ 2*(gridabove->get_res(right_above)+gridabove->get_res(left_above)\
                            +gridabove->get_res(upper_above) +gridabove->get_res(down_above))+(gridabove->get_res(upperleft_above)\
                            +gridabove->get_res(upperright_above)+gridabove->get_res(downleft_above)+gridabove->get_res(downright_above));

                    double res = residual_below/16;

                    //std::cout<<res << " Rhs"<< " "<<index<< '\n';

                    gridbelow->set_rhs(res, gridbelow->vector_index(i,j));
                    }

                    else{
                        std::cout<<"Out of bound"<<'\n';
                    }

                }
                
            }
            
        }// end of restriction


 void MG::interpolation(unsigned int level){
                // Prolonge residual and correct the fine-grid approximation

                // Fetch the current grid, a coarse grid 
              
                Grid* gridbelow = get_grid(level);
               //std::cout<<"get gridbelow"<<" "<<level<<'\n';
                // Get the size of the corase grid
                unsigned int size_below{gridbelow->get_size()};

                // Fetch the above fine grid
                Grid* gridabove = get_grid(level-1);

                //std::cout<<"get gridabove"<< " "<< level-1<<'\n';

                gridabove -> get_size();
            
                // Get the size of the fine grid
                unsigned int size_above = gridabove -> get_size();
                //std::cout<< "interpolate"<<'\n';
                for (int i = 0; i< size_below; i++){

                        for (int j = 0; j<size_below; j++){

                            unsigned int index = gridbelow->vector_index(i,j);
                            //std::cout<< " approx gridbelow " << index<<'\n';
                            
                            double error_above = gridbelow->get_approx(index);

                            gridabove->set_res(error_above, gridabove->vector_index(2*i,2*j));
                            //std::cout<< "Residual on interpolation "<< error_above<<" at "<< gridabove->vector_index(2*i,2*j)<< '\n';
                    }

                }
                
                for (int i=0; i<size_above; i+=2){

                    for (int j = 1; j<size_above-1; j+=2){

                            // for k in range(depth):
                            //     for i in range(0, ynodes, 2):
                            //          for j in range(1, xnodes-1, 2):
                            //              uf[k,i,j]=0.5*(uf[k,i,j-1]+uf[k,i,j+1])

                            unsigned int index = gridabove->vector_index(i,j);

                            unsigned int index_left = gridabove->vector_index(i, j-1);

                            unsigned int index_right = gridabove-> vector_index(i,j+1);

                            double error_above = 0.5* (gridabove->get_res(index_left)+ gridabove ->get_res(index_right));

                            gridabove ->set_res(error_above, index);
                            //std::cout<< "Residual on interpolation "<< error_above<<" at " << index<< '\n';
                    }
                }

                for (int i =1; i<size_above-1; i+=2){
                    
                   for (int j = 0; j<size_above; j+=2){

                        unsigned int index = gridabove-> vector_index(i,j);

                        unsigned int index_down = gridabove ->vector_index(i-1, j);

                        unsigned int index_above = gridabove -> vector_index(i+1, j);

                        double error_above = 0.5* (gridabove->get_res(index_above) + gridabove -> get_res(index_down));

                        gridabove ->set_res(error_above, index);
                        //std::cout<< "Residual on interpolation "<< error_above<<" at "<< index<< '\n';

                    }
                }

                for (int i =1; i< size_above -1; i+=2){

                    //std::cout<<"odd i "<< i <<" "<< size_above -1<<'\n';
                
                    for (int j= 1; j<size_above -1; j+=2){


                            unsigned int index = gridabove->vector_index(i,j);

                            unsigned int index_below = gridabove->vector_index(i-1, j);

                            unsigned int index_above = gridabove-> vector_index(i+1,j);

                            unsigned int index_left = gridabove ->vector_index(i, j-1);

                            unsigned int index_right = gridabove -> vector_index(i, j+1);

                            double error_above = 0.25*( gridabove ->get_res(index_below) +gridabove->get_res(index_above)\
                                                          + gridabove ->get_res(index_left) + gridabove->get_res(index_right));
                            gridabove ->set_res(error_above, index);


                    }
                }          

            }



};