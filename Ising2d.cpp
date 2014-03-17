#include "Ising2d.h"


/* 
Constructor
parm   : m is the row or column size  MUST BE M*M !
         matrixType is 0-4 that creates a matrix using
         the createNewMatrix function.
method : sets the size of the state vector and creates a
		 matrix depending on the matrixType passed. 
*/	 
Ising2d::Ising2d(int m, int matrixType){
	m_m = m;
	m_state.resize(m_m*m_m);
	createNewMatrix(matrixType);
}


/*
fuction : createNewMatrix(int matrixType)
params  : takes value 0-4 each value creating a matrix
          of various spin alignments. 
method  : 0  fills matrix with -1
		  1  fills matrix with  1
		  2  fills matrix with random spins
		  3  fills matrix with alternating elements up/down
		  4  fiils matrix with alternating rows up/down
*/
void Ising2d::createNewMatrix(int matrixType){


	if(matrixType == 0)
		fill(m_state.begin() , m_state.end(), -1);
	
	if(matrixType == 1)
		fill(m_state.begin() , m_state.end(), 1);

	if(matrixType == 2)
		for (int i = 0; i <m_m*m_m; i++){
			m_state[i] = (rand() % 2) * 2 - 1;
		}

bool flipNext = false;
	if(matrixType == 3)
		for (int i = 0; i <m_m*m_m; i++){
			if(flipNext == false){
				m_state[i] = 1;
				flipNext = true;
			}else {
				if(flipNext ==true){
					m_state[i] = -1;
					flipNext = false;
				}
			}
		}	

int s = 1.0;	
	if(matrixType == 4)
		for (int i = 0 ; i <m_m*m_m; i++){
			m_state[i]= 1.0*s;
		if((i+1) %m_m == 0.0){
			s = -1.0*s;
		}
	}
}


/*
 fuction ; getStateEnergy()
 method  : calculates energy of a state configuration by summing 
           the product of a spin and its neighbors.

	       s is the current spin
	       n is the contrubution from nearest neighboring spins  
				   			   top [i-1,j]
				   				|
				   				|

		    [i, j-1] left ---  spin  --- right [i, j+1]
								
								|
								|
							  bottom [i+1,j]
			the sum is then divided by 2 to prevent double count
			-ing and also divided by the number of spin sites to 
			obtain the average energy per site. Valve is stored
			in the variable m_stateEne.
*/
void Ising2d::calcStatesEnergy(){
	m_stateEne = 0.0;
	int m = m_m - 1; //used in spin args to keep proper indexing
	for (int i = 0 ; i < m_m ; i++){
		for (int j = 0 ; j < m_m ; j++){
			int s = spin(i,j);
			int neighbors, top, bottom, left, right; 
 			
 			//testing if on top/bottom edge. 
 			if ( (i - 1) < 0 ){
 				top = spin( m  , j); 
 			}else{
 				top = spin(i - 1 , j);
 			};

 			if ( (i + 1) > m){
 				bottom = spin( 0 , j);
 			}else{
 				bottom = spin( (i + 1) , j);
 			};

			//testing if on right/left edge. 
 			if ( (j - 1) < 0 ){
 				left = spin( i, m ); 
 			}else{
 				left = spin( i, j - 1);
 			};

 			if ( (j + 1) > m){
 				right = spin( i , 0);
 			}else{
 				right = spin( i, (j + 1));
 			};

 			neighbors = right+left+top+bottom;
			m_stateEne += -double(s)*double(neighbors)*m_j;	
		}
	}
	 m_stateEne = m_stateEne/(2.0*m_m*m_m); //prevent double counting 

}

/*
fuction : getStatesEnergy()
returns : the energy of the current state m_stateEne.
method  : first calls calcStatesEnergy() to calculate the states
		  energy then returns the value. 
*/
double Ising2d::getStatesEnergy(){
	calcStatesEnergy();
	return m_stateEne;
}

/*
 fuction ; getChangeEnergy(int i, int j)
 params  : i is row index j is the column index. 
 returns : energy of the current states configuration.
 method  : calculates energy of a state configuration taking
           the product of a spin and the sum of its neighbors.

	       s is the current spin
	       n is the contrubution from nearest neighboring spins  
				   			   top [i-1,j]
				   				|
				   				|

		    [i, j-1] left ---  spin  --- right [i, j+1]
								
								|
								|
							  bottom [i+1,j]
			unlike getStatesEnergy this function does not loop
			through each spinsite summing their products.
*/
double Ising2d::getChangeEnergy(int i , int j){
	m_changeEne=0.0;
	int m = m_m - 1;
	int s = spin(i,j);
	int neighbors, top, bottom, left, right; 
	
	//testing if on top/bottom edge. 
	if ( (i - 1) < 0 ){
		top = spin( m  , j); 
	}else{
		top = spin(i - 1 , j);
	};

	if ( (i + 1) > m){
		bottom = spin( 0 , j);
	}else{
		bottom = spin( (i + 1) , j);
	};

	//testing if on right/left edge. 
	if ( (j - 1) < 0 ){
		left = spin( i, m ); 
	}else{
		left = spin( i, j - 1);
	};

	if ( (j + 1) > m){
		right = spin( i , 0);
	}else{
		right = spin( i, (j + 1));
	};

	neighbors = right+left+top+bottom;
	m_changeEne = 2*double(s)*double(neighbors)*m_j;

	return m_changeEne;
}

/*
fuction : getStatesMag()
returns : m_mag 
method : calcuates the magnetization of the state
*/
double Ising2d::getStatesMag(){
	double totSpin = 0;
	double m = double(m_m)*double(m_m);
	for (int i = 0; i < m_m*m_m; i++){
		totSpin = double(totSpin) + double(m_state[i]);
	}
	m_mag = totSpin/m;
	return m_mag;
}

/*
fuction : getStatesMag()
params  : row index i column index j 
returns : element i,j of the state matrix
method  : preserves 2d indexing of the flattend matrix
*/
int Ising2d::spin(int i, int j ){
	//cout << i*m_m +j << endl;
	return m_state[i*m_m + j]; //allows i,j indexing to be perserved 
}

/*
fuction : getStatesMag()
params  : row index i column index j and a value to set element to
returns : element i,j of the state matrix
method  : Overloaded function allowing the value of element i,j to 
          be changed. 
*/
void Ising2d::spin(int i , int j , int value){
	 m_state[i*m_m + j] = value;
}

/*
fuction : printMatrix()
method  : Prints the states configuration in the terminal. + for 1 
		  - for -1 . Output is colored to better distigush spin types.
*/
void Ising2d::printMatrix(){

	
	for (int i = 0; i < m_m*m_m; i++){
		
		if(m_state[i]==1.0){
			cout<< "\033[1;31m"<< "+ " <<"\033[0m";// \033 ASCII escape char #;#m color and other params
		}else{
			cout<< "\033[1;36m"<< "- " << "\033[0m" ;
		}
		if ( (i+1) % m_m == 0){
			cout << endl;
		}

	}
}		
	





