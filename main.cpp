#include <time.h>
#include <fstream>

#include "Ising2d.h"

const int maxStep = 500000;
const int warmUpScalar = maxStep/10;
const int latticeSize = 32;  //   MxM lattice 
int numOfBins = 10;
int binSize = maxStep/numOfBins;

int maxHeatingStep;

double T = 7.0;
const double finalTemp = 1.5;
const double tempStep = 0.1;
const double K_b =1.0; //actual value : 1.3806503 x 10^(-23)





vector<double> energyPerStepBin;
vector<double> avgEnergyPerStepBin; 

vector<double> eneSquaredBin;
vector<double> avgEneSquaredBin;


vector<double> magnetizationPerStepBin; 
vector<double> avgMagnetizationPerStepBin;

vector<double> magSquaredBin;
vector<double> avgMagSquaredBin;

//MEASUREMENTS
vector<double> tempValues;
vector<double> energyMeasurements;
vector<double> eneSquaredMeasurements;
vector<double> magnetizationMeasurements;
vector<double> magSquaredMeasurements;

vector<double>magSuscept;
vector<double>specficHeat;


void saveToCsv();
double average(vector<double>& v);

using namespace std;

int main(){

	

	Ising2d IsingSim (latticeSize,2);

	/* Prints Inital Config Note May overflow terminal Depending
	on the size of the Matrix */
	//cout <<"====== inital Configuration ======"<<endl;
	//IsingSim.printMatrix();

	double E = IsingSim.getStatesEnergy();
	double M = IsingSim.getStatesMag();

	cout << " "<< endl;
	cout << " binsize : " << binSize <<endl;
	cout << " Inital Temperature of system : " << T <<endl;
	cout << " inital energy of the state : " << E << endl;
	cout << " inital magnetization  : " << M << endl;
	cout <<" "<<endl;

	clock_t start, end;

	start = clock();

	while(T>finalTemp){

		int count = 0;

		//SET maxHeatingStep depending on T

		if (T<=6.0) maxHeatingStep = 4*warmUpScalar;
		if (T<=5.0) maxHeatingStep = 5*warmUpScalar;
		if (T<=4.0) maxHeatingStep = 6*warmUpScalar;
		if (T<=3.0) maxHeatingStep = 7*warmUpScalar;
		if (T<=2.0) maxHeatingStep = 8*warmUpScalar;
		if (T<=1.0) maxHeatingStep = 10*warmUpScalar;
		

	    
		cout <<" +++++++++ Performing Initial Sweeps +++++++++++++++ "<<endl;
		cout<< "at Temerature : " << T <<endl;
		cout<< "with "<< maxHeatingStep <<" iterations "<<endl;
		for(int step = 0 ; step < maxHeatingStep ; step++ ){
		 
		 	int i = rand()%latticeSize;
		 	int j = rand()%latticeSize;

		 	double delEne = IsingSim.getChangeEnergy(i, j);
			double boltzman = exp((-1.0*delEne)/(K_b*T));

			if(delEne <= 0.0){
				IsingSim.flipSpin(i,j);//flip random spin
			}else{

				double n =((double)rand()/(double)RAND_MAX); //some voodoo to make random double between 0 to 1
				
				if(n <= boltzman){
					IsingSim.flipSpin(i,j);		
				}
			}	
		}


		cout <<" +++++++++ MONTE CARLO START at Temp: " << T << " +++++++++++++++ "<<endl;

		int rejectedFlips = 0;
		int delEneFlip = 0;
		int boltzmanFlip = 0;

		for(int step = 0 ; step < maxStep ; step++ ){
		 
		 	int i = rand()%latticeSize;
		 	int j = rand()%latticeSize;

		 	double delEne = IsingSim.getChangeEnergy(i, j);
			double boltzman = exp((-1.0*delEne)/(K_b*T));
			
			if(delEne <= 0.0){

				IsingSim.flipSpin(i,j);//flip random spin
				energyPerStepBin.push_back(IsingSim.getStatesEnergy());
				eneSquaredBin.push_back(IsingSim.getStatesEnergy()*IsingSim.getStatesEnergy());
				delEneFlip++;
		 		
			}else{

				double n =((double)rand()/(double)RAND_MAX); //some voodoo to make random double between 0 to 1

				if(n <= boltzman){

					IsingSim.flipSpin(i,j);	
					energyPerStepBin.push_back(IsingSim.getStatesEnergy());
					eneSquaredBin.push_back(IsingSim.getStatesEnergy()*IsingSim.getStatesEnergy());
					boltzmanFlip++;
					
		 				
				}else{

				energyPerStepBin.push_back(IsingSim.getStatesEnergy());
				eneSquaredBin.push_back(IsingSim.getStatesEnergy()*IsingSim.getStatesEnergy());
				rejectedFlips++;
		
				}
			}

			magnetizationPerStepBin.push_back(IsingSim.getStatesMag());
			magSquaredBin.push_back(IsingSim.getStatesMag()*IsingSim.getStatesMag());

			while( energyPerStepBin.size() == binSize ){
				
				double E = average(energyPerStepBin);
				avgEnergyPerStepBin.push_back(E);

				double E2 = average(eneSquaredBin);
				avgEneSquaredBin.push_back(E2);

				double M = average(magnetizationPerStepBin);
				avgMagnetizationPerStepBin.push_back(M);

				double M2 = average(magSquaredBin);
				avgMagSquaredBin.push_back(M);

				energyPerStepBin.clear();
				eneSquaredBin.clear();
				magnetizationPerStepBin.clear();
				magSquaredBin.clear();
		
			}
		}

		//take Measurements;
		tempValues.push_back(T);
		double E = average(avgEnergyPerStepBin);
		energyMeasurements.push_back(E);

		double E2 = average(avgEneSquaredBin);
		eneSquaredMeasurements.push_back(E2);

		double c = (E2-(E*E))/(T*T);
		specficHeat.push_back(c);

		double M = average(avgMagnetizationPerStepBin);
		magnetizationMeasurements.push_back(M);
		double M2 = average(avgMagSquaredBin);
		magSquaredMeasurements.push_back(M2);
		double x = (M2 - (M*M))/(T);
		magSuscept.push_back(x);

		avgEnergyPerStepBin.clear();
		avgEneSquaredBin.clear();
		avgMagnetizationPerStepBin.clear();
		avgMagSquaredBin.clear();

		cout <<" +++++++++ MONTE CARLO FINISH at Temp : " << T << " +++++++++++++++ "<< endl;
		cout <<endl;
		cout<<" < E > : "<< energyMeasurements[energyMeasurements.size()-1] <<endl;
		cout<<" < E^2 > : " << eneSquaredMeasurements[eneSquaredMeasurements.size()-1] <<endl;
		cout<<" < M > : " << magnetizationMeasurements[magnetizationMeasurements.size()-1] << endl;
		cout<<" < M^2 > : " << magSquaredMeasurements[magSquaredMeasurements.size()-1] << endl;
		cout<<"Specfic Heat : " <<specficHeat[specficHeat.size()-1] <<endl;
		cout<<"Magnetic Susceptability : "<<magSuscept[magSuscept.size()-1] <<endl;
		cout<<" Number of Rejected Flips : " << rejectedFlips <<endl;
		cout<<" Number of delE <= 0.0 Flips : " << delEneFlip <<endl;
		cout<<" Number of random n <= boltzman Flips : " << boltzmanFlip<<endl;
		
		// uncomment the line below to print the matrix to the terminal 
		// Note this may overflow the line for large matrices

		//IsingSim.printMatrix();


		T = T - tempStep;
	}

	end = clock();

	cout<<" Total Simulation Runtime: "<< double((end - start)) / double(CLOCKS_PER_SEC)<<"s"<<endl;
	cout<<" Final Energy : "<< energyMeasurements[energyMeasurements.size()-1]<< endl;
	cout<<" Final Manetization : "  << magnetizationMeasurements[magnetizationMeasurements.size()-1] << endl;
	saveToCsv();

}


void saveToCsv(){
	ofstream file;
	file.open ("/Users/Eric/Github_Repos/Ising2D/Data.csv"); //change to directory that your repo is saved.
	file <<"step,T,E,M,X,C"<<endl;
	for (int i = 0; i <tempValues.size(); i++){
		file << i << "," << tempValues[i] << "," << energyMeasurements[i]<< "," << magnetizationMeasurements[i] << "," << magSuscept[i] << "," << specficHeat[i]<<endl; 
	}
	file.close();
}

double average(vector<double>& v){
	double sum;
	for(int i = 0; i < v.size(); i++){
		sum = sum + v[i];
	}
	return sum/v.size();
}