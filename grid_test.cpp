
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

         double spacing{1.0/(_N -1)};

    
    public:
        //Intialise private data

        // Constructor 
        Grid(unsigned n):_N{n}, _Ntot{n*n}, values{_Ntot, std::make_pair(1.0,1.0)}
        {
        };
        // Deconstructor
        ~Grid()
        {
            //delete []  values; //1. raise error, cannot delete expression of type.
                                 //2. maybe redundant since we don't have any dynamically allocated memory.
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

        std::pair<unsigned int, unsigned int> grid_index(const unsigned int i){
            

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
            const std::vector<std::pair<double, double>> &ref2{values};
            

            // Get ith element
            std::pair<double, double> node_i = ref2[i];

            return node_i;


        



        }  

        unsigned int get_size() const {

        //std::cout<<_N<< " size"<<"\n";

        return _N;
        }

        // Set values to nodes

        void set_approx( double& approx, const int & i){

            values[i].first = approx;


        }


        // Set rhs to nodes

        void set_rhs( double & rhs, const int & i){

            values[i].second =rhs;
        }

        //friend class MG;

        bool boundary(const unsigned int i){

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

            for (unsigned int i =0; i<l; i+=1){
                unsigned int number = pow(2, i+1)+1; 
                ptrgrids[i] = new Grid{number}; // dynamically allocated grid size to memory
            }
        }

        ~MG(){

            for(int i=0; i< Nlevel ;i++)

            delete ptrgrids[i]; 

            // delete grids first.
        }




         void print_grid(unsigned int i){


             //double value1 = 2.0;
             //ptrgrids[i] ->set_approx(value1, 0);
             ptrgrids[i] ->print_nodes();

             //return ptrgrids[i]->get_node(i);
         }

        // Get ith grid from the Grid pointer array
         Grid* get_grid(unsigned i){

             return ptrgrids[i];
         }


         auto restriction(const Grid &gridabove){


             // Get the size of current grid
             unsigned int sizeuabove{gridabove.get_size()};

             //std::cout<<sizeuabove<<"sizeuabove"<<'\n';

             // Get the value of current grid
             std::vector<std::pair<double, double>> valabove = gridabove.get_values();

             // Initialise valbelow

            unsigned int sizeubelow{(sizeuabove+1)/2};

            Grid gridbelow{sizeubelow};

            unsigned int totbelow = sizeubelow * sizeubelow;



            //std::cout<<sizeubelow<<"sizeubelow"<<'\n';

            //std::vector<std::pair<double, double>> ubelow{sizeubelow*sizeubelow, std::make_pair(0.0,0.0)};


            for (int i=0; i<totbelow; i++)
            {
                if (gridbelow.boundary(i) == false)
                {

                std::cout<< i<<" i"<<"\n";
                // Conver each node index in lower grid to a 2D grid index
                std::pair<unsigned int, unsigned int>  ixiy_below{gridbelow.grid_index(i)};

                // Find the (2x,2y) node index in upper grid
                int centre_above = gridabove.vector_index(2* ixiy_below.first, 2* ixiy_below.second);

                std::cout<< centre_above<<" centre"<<'\n';
                
                // Find the (2x-1, 2y) node index in upper grid
                int right_above = gridabove.vector_index(2* ixiy_below.first-1, 2* ixiy_below.second);

                // Find the (2x+1, 2y) node index in upper grid
                int left_above = gridabove.vector_index(2* ixiy_below.first+1, 2* ixiy_below.second);

                // Find the (2x,2y-1) node index in upper grid
                int down_above = gridabove.vector_index(2* ixiy_below.first,   2* ixiy_below.second-1);

                // Find the (2x, 2y+1) node index in upper grid
                int upper_above = gridabove.vector_index(2* ixiy_below.first,   2* ixiy_below.second+1);

                // Find the (2x-1, 2y-1) node
                int downleft_above = gridabove.vector_index(2*ixiy_below.first-1, 2*ixiy_below.second-1);
                
                // Find the (2x-1, 2y+1) node
                int upperleft_above = gridabove.vector_index(2*ixiy_below.first-1, 2*ixiy_below.second+1);

                // Find the (2x+1, 2y-1) node
                int downright_above = gridabove.vector_index(2*ixiy_below.first+1, 2*ixiy_below.second-1);

                // Find the (2x+1, 2y+1) node
                int upperright_above = gridabove.vector_index(2*ixiy_below.first+1, 2*ixiy_below.second+1);

                std::pair<double, double> centre = gridabove.get_node(centre_above);

                centre.first=centre.first*1.0;

                centre.second = centre.second*1.0;

                std::pair<double, double> right = gridabove.get_node(right_above);

                right.first = right.first*0.5;

                right.second = right.second *0.5;

                //std::cout<< node2.first<< "node.first"<<"\n";

                std::pair<double, double> left = gridabove.get_node(left_above);

                left.first = left.first*0.5;

                left.second = left.second*0.5;

                std::pair<double, double> down = gridabove.get_node(down_above);

                down.first = down.first*0.5;

                down.second = down.second*0.5;

                std::pair<double, double> upper = gridabove.get_node(upper_above);

                upper.first = upper.first*0.5;

                upper.second = upper.second*0.5;

                std::pair<double, double> downleft = gridabove.get_node(downleft_above);

                downleft.first = downleft.first* 0.25;

                downleft.second = downleft.second*0.25;

                std::pair<double, double> upperleft = gridabove.get_node(upperleft_above);

                upperleft.first = upperleft.first*0.25;

                upperleft.second = upperleft.second*0.25;

                std::pair<double, double> downright = gridabove.get_node(downright_above);

                downright.first  = downright.first*0.25;

                downright.second = downright.second*0.25;

                std::pair<double, double> upperright = gridabove.get_node(upperright_above);

                upperright.first = upperright.first *0.25;

                upperright.second = upperright.second *0.25;

                double new_values = centre.first+right.first+left.first +upper.first + down.first + upperleft.first + upperright.first+ downleft.first+ downright.first;

                double new_rhs    = centre.second+ right.second+left.second +upper.second + down.second + upperleft.second + upperright.second+ downleft.second+ downright.second;




                gridbelow.set_approx(new_values, i);
                gridbelow.set_rhs(new_rhs, i);

                }


            }


            return gridbelow;
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

    Grid* grid2 = Vcycle.get_grid(1);

    std::cout<<"The size of grid"<<"\n";

    grid2->get_size();

    //Vcycle.restriction(*grid2);

    Grid gridbelow{Vcycle.restriction(*grid2)}; 

    gridbelow.print_nodes();

    
    //Vcycle.get_values(2);
    return 0;
}
