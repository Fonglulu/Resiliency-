#include "grid.h"

// Simple int-int a^b power function

inline unsigned int pow(unsigned int a, unsigned int b){
  unsigned int res =1;
  for (unsigned int i=0; i<b; i++){
    res *= a
  }

  return res;
}

template<unsigned int NDIM, typename T>

void Grid<NDIM, T>:: check_for_nan(bool exitifnan){
  bool nanfound = false;
  for (unsigned int i =0; i< _Ntot; i++){
    if (_y[i] != _y[i]){
      nanfound =true;
      break;
    }
  }
  if (nanfound){
    std::cout << "Warining: NaN found in grid" <<( exitifnan ? "...aborting" : "") << std::enl;
    if (exitifnan) exit(1);
  }
}

// Constructor w/ inital value
template<unsigned int NDIM, T>
Grid<NDIM, T>:: Grid(unsigned int N, T yini): _N(N), _Ntot(power(_N, NDIM)), _y(std::vector<T>(_Ntot, yini)){}


// Fetch pointer to grid
template<unsigned int NDIM, typename T>
T* Grid<NDIM, T>::get_y(){
  return &_y[0];
}

//Allow to fetch value using f[i] syntax
template <unsigned int NDIM, typename T>
T& Grid<NDIM,T> :: operator[] (unsigned int i){
#ifdef _BOUNDSCHECK
  assert(i < _Ntot);
#endif
  return _y[i];
}

// Fetch value of grid-cell [i]
T Grid<NDIM, T>::get_y(unsigned int i){
#ifdef _BOUNDSCHECK
  assert(i < _Ntot);
#endif
  return _y[i];
}


// Assign whole grid from vector
template<unsigned int NDIM, typename T>
void Grid<NDIM, T>::set_y(std::vector<T> &y){
#ifdef _BOUNDSCHECK
  assert(y.size() == _Ntot);
#endif

  _y = y;
}

// Assign the grid cell [i] with value
template <unsigned int NDIM, typename>
void Grid<NDIM, T>:: set_y(unsigned int i, T& value){
#ifdef _BOUNDSCHECK
  assert(i < _Ntot);
#endif
  _y[i] = value;
}

// Compute coordinates given a gridindex
template <unsigned int NDIM, typename T>
std::vector <unsigned int> Grid<NDIM,T>::index_list(unsigned int i){
  std::vector <unsigned int> ii(NDIM, 0);
  for (unsigned int  j = 0, n=1; j < NDIM; j++, n *= _N){
    ii[j] = i/n%_N;
  }

  return ii;
}

// Coordinate -> grid-index (index in the 1D _y vector)
template<unsigned int NDIM, typename T>
unsigned int Grid<NDIM, T> :: grid_index(std::vector <unsigned int> &index_list){
  unsigned int index = 0 ;
  for (unsigned int j=0, n=1; j<NDIM; j++, n *= _N)
    index += index_list[j]*n;
#ifdef _BOUNDSCHECK
  assert(index < _Ntot);
#endif
  return index;
}

// Coordinate -grid-index for 2D grid
template<unsigned int NDIM, typename T>
unsigned int Grid <NDIM, T>:: grid_index_2d(unsigned int ix, unsigned int iy){
    return ix+_N*iy;
}


// Return  number of cells per dim
template<unsigned int NDIM, typename T>
unsigned int Grid<NDIM, T>::get_N(){
  return _N;
}

// Return total number of cells
template<unsigned int NDIM, typename T>
unsigned int Grid<NDIM, T>: get_Ntot(){
  return _Ntot;
}

// // Write a grid to file
// template <unsigned int NDIM, typename T>
// void Grid<NDIM, T>::dump_to_file(std:: filename){
//   unsigned int ndim = NDIM;

// }

