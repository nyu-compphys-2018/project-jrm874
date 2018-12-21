#include "itensor/all.h"
#include "tools.h"
#include <iostream>
#include <fstream>

using namespace itensor;
using namespace std;


int 
main(int argc, char* argv[])
    {
        //Imput form an imput file (so there's no need to recompile every time any variable is changed)
        if(argc < 2)
        {
            printfln("Usage: %s <input_file>", argv[0]);
            return 0;
        }
        auto infile = InputFile(argv[1]);
        auto input = InputGroup(infile,"input");
        //Get the variables from input file
        auto N = input.getInt("N"); //number of lattice points
        auto beta = input.getReal("beta"); // 1/temperature
        auto hsteps = input.getReal("hsteps"); //number of magnetic field steps
        auto hstep = input.getReal("hstep"); // magnetic field step size
        auto hinit = input.getReal("hinit"); // initial magnetic field value
        auto J1 = input.getReal("J1"); // Coupling constant
        auto J2 = input.getReal("J2"); // Coupling constant
        auto u = input.getReal("u"); // Measure magnetizarion between u and v
        auto v = input.getReal("v");
        auto tau = input.getReal("tau",0.1); // break beta in steps of size Tau for Trotter decomposition
        auto cutoff = input.getReal("cutoff"); //Cutoff parameter for regauging the MPS
        auto maxm = input.getInt("maxm",5000);
        auto MMETTS = input.getInt("MMETTS",50000); //Total number of metts generated
        auto Malias = input.getInt("Malias",5); //Number of METTS generated as aliasing
        
        ofstream mH("mH.txt");
        
        //Add a header with the parameters of the simulation in the output file
        mH<<"N = "<<N<<";  beta = "<<beta<<";  tau = "<<tau<<";  cutoff = "<<cutoff<<";  maxm = "<<maxm<<";  MMETTS = "<<MMETTS<<";  Malias = "<<Malias<<";  J1 = "<<J1<<";  J2 = "<<J2<<endl;
        mH<<endl;
        
        Real h = hinit;
        
        //Define the hamiltonian to study
        auto hamiltonian = "1D Heisenberg Second Nearest Neighbour";
        
        //Calls METTS and saves the magnetization for different values of magnetic field
        for (int i =0; i < hsteps; i++)
        {
            //Obtaining the magnetizartion
            cout<<"Computing magnetization for h =  "<<h<<endl;
            
            auto Mag = METTS_MAG(N, h, J1, J2, u, v, beta, tau, Malias, MMETTS, maxm, cutoff, {"Hamiltonian=", hamiltonian});
            
            cout<<"Done computing"<<endl;
            cout<<"Saving data"<<endl;
            
            //exporting the data
            mH<<h<<"  ";
            for (int j = 0; j < MMETTS; j++)
            {
                Real m = Mag.at(j);
                mH<<m<<"  ";
            }
            mH<<endl;
            cout<<"Done saving data"<<endl;
            
            h += hstep;
        }
        
    
    }



