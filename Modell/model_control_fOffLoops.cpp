#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <deque>
#include <math.h>

#include <stdlib.h>
#include <stdio.h>

using namespace std;
void ErrorhandlerExit(const char* p1,const char* p2="")
{
  cerr << p1 << "  " << p2 << endl;
  std::exit(1);
}

int main(int argc,char **argv){

//--------------------VARIABLE DECLARATION---------------//
//**MODEL PARAMETERS**//
// procuction rates
double k0=1.0;					// TF mRNA transcription rate
double k1=1.0;					// TF protein translation rate
double k2=1.0; 					// target mRNA transcription rate
// degradation rates		
double d0=1.0;				    // TF mRNA degradation
double d1=1.0;					// TF protein degrataion rate
double d2=1.0; 				    // target mRNA degradation rate
// amounts
double mRNA=0.0;				// TF mRNA level
double ptf=0.0; 				// TF protein level
double target=0.0;				// target mRNA level
// frequency
double frequency;				// frequency
double f_iter[]{0.001, 0.002, 0.003,0.005, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.3, 0.5, 0.6, 0.9, 1, 1.4, 1.78, 9};	// frequency of k0 input signal
// offset
double offSet;
double offSet_iter[]{0, 0.1, 0.5, 1, 5, 10, 50}; //


// model running parameters
double t=0.0;					// elapsed time
double t_max = 1000;			// run time
double dt=0.0007;  				// times steps
long counter=0;					// model iteration counter; used to regulate output
double max_target=0.0;			// holds largest occured value; used to indicate equilibrium state


//**PROGRAMM PARAMETERS**//
long program_counter=0; 		// program counter; indicates programm progress
long max_program_counter;		// indicates # of iterations programm will run
int i, j;						// loop iterators

//**OUTPUT FILE**//
ofstream file_out;
string file_name, max_name, fstart="sine_";
stringstream concat_file_name;



//--------------------PROGRAM CODE---------------//

max_program_counter = (sizeof(f_iter)/sizeof(f_iter[0])) * (sizeof(offSet_iter)/sizeof(offSet_iter[0]));

// offset loop
for(i=0; i<(sizeof(offSet_iter)/sizeof(offSet_iter[0])); i++) 
{

	// select offset for current iteration
	offSet = offSet_iter[i];

	for(j=0; j<(sizeof(f_iter)/sizeof(f_iter[0])); j++)
	{
		
		// select frequency for current iteration
		frequency = f_iter[j];

		// prepare output file
		concat_file_name.str("");
		concat_file_name << fstart << "f-" << frequency << "_d0-" << d0 << "_k1-" << k1 << "_off-" << offSet;
		
		file_name = concat_file_name.str() + ".txt";
		file_out.open(file_name.c_str());
		file_out << std::setprecision(15);
		file_out << "t" << "\t" << "k0" << "\t" << "mRNA" <<  "\t" << "ptf" << "\t" << "target" << endl;
		
		// reset parameters for next round
		t = 0;
		mRNA = offset;
		ptf = 0;
		target = 0;
		max_target = 0;
		counter = 0;
		
		//model loop
		while(t<t_max)
		{

			// OEDs
			k0=sin(2*3.1415*frequency*t);
			mRNA+=(k0-mRNA*d0)*dt;
			ptf+=(mRNA*k1-d1*ptf)*dt;
			target+=(ptf*k2-d2*target)*dt;
	
			// keep track of max_target 
			if((target>max_target) && (t > 990))
			{
				max_target=target;
			}

			// write parameters to file
			file_out << t << "\t" << k0 << "\t" << mRNA <<  "\t" << ptf << "\t" << target << endl;
			
	
			counter++;
			t+=dt;
			
			
		

		} // model loop end

		file_out.close();
		
		program_counter++;
		cout << "Finished " << program_counter << " of " << max_program_counter << " combinations; f = " << frequency << ", offset = " << offSet << endl;
		//max_out.close();
		
	} // frequency loop end
	
} //offSet loop end

} // main end





