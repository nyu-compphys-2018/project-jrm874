// In this script I will implement all the tools required to run the METTS algorithm
// Includes code to:
//      -Measure certain observables (Magnetization)
//      -Collapse state to a random cps with probability |<psi|cps>|^2
//      -Create a Neel initial CPS
//      -Generate and exponentiate the wanted hamiltonian
//      -Run METTS algorithm


// Including all the required libraries
#include <string>
#include <iostream>
#include "itensor/all.h"

using namespace itensor;
using namespace std;
//using namespace math;


//----COLLAPSE ANY MPS TO A CPS----
//Accepsts the MPS (METTS) and the direction to which it it will be collapsed
MPS collapse(MPS psi, SiteSet sites, Args const& args = Global::args())
{
    float sqrt2 = 1.41421356237;
    Complex isqrt2 = 1.41421356237_i;
    Real prob_up = 0.;
    
    auto N = psi.N(); // Size of the system
    auto CPS = vector<int> (N+1);
    
    auto direction = args.getString("Direction"); // Direction in which we project the state
    
    //collapsed state. Initialized to all up
    InitState Nstate(sites,"Up");
    

    //Collapse site by site
    for(int i = 1; N >= i; i++ )
    {
        // obtain index j in the MPS
        Index ind_i = sites(i);
        Index pind_i = prime(ind_i);

        // Defining the projector +, can be in z, x or y directions
        auto P_Up = ITensor(ind_i, pind_i);
        auto P_Dn = ITensor(ind_i, pind_i);
        
        //Defining the up and down states tensors
        auto up = ITensor(ind_i);
        auto dn = ITensor(ind_i);
        
        // Here we specify if it is z,x direction and compute projectors and states up and down on site i
        if(direction == "z")
        {
            P_Up.set(ind_i(1), pind_i(1), 1.0);
            P_Dn.set(ind_i(2), pind_i(2), 1.0);
        }
        
        else if(direction == "x")
        {
            P_Up.set(ind_i(1), pind_i(1), 0.5);
            P_Up.set(ind_i(1), pind_i(2), 0.5);
            P_Up.set(ind_i(2), pind_i(1), 0.5);
            P_Up.set(ind_i(2), pind_i(2), 0.5);
            
            P_Dn.set(ind_i(1), pind_i(1), 0.5);
            P_Dn.set(ind_i(1), pind_i(2), -0.5);
            P_Dn.set(ind_i(2), pind_i(1), -0.5);
            P_Dn.set(ind_i(2), pind_i(2), 0.5);
        }
        else Error("Direction '" + direction + "' not allowed");
        
        // Move ortogonality center to position i
        psi.position(i);
 
        // Compute the probability of finding site j on state up
        prob_up = ( dag(prime(psi.A(i),Site)) * P_Up * psi.A(i) ).real();
        Real R = Global::random();
        
        // generate a random number to cecide whether collapsing to up or down
        if(R < prob_up)
        {

            Nstate.set(i,"Up");
            CPS.at(i) = 1;
            
        }
        if(R >= prob_up)
        {

            Nstate.set(i,"Dn");
            CPS.at(i) = -1;
            
        }
    }
    
    auto Npsi = MPS(Nstate);
        
    for(int i = 1; i <= N; i++)
        {
            // obtain index j in the MPS
            Index ind_i = sites(i);
            Index pind_i = prime(ind_i);
            
            // Defining the projector +, can be in z, x or y directions
            auto P_Up = ITensor(ind_i, pind_i);
            auto P_Dn = ITensor(ind_i, pind_i);
            
            //Defining the up and down states tensors
            auto up = ITensor(ind_i);
            auto dn = ITensor(ind_i);
            
            // Here we specify if it is z,x direction and compute projectors and states up and down on site i
            if(direction == "z")
            {
                P_Up.set(ind_i(1), pind_i(1), 1.0);
                P_Dn.set(ind_i(2), pind_i(2), 1.0);
            }
            
            else if(direction == "x")
            {
                P_Up.set(ind_i(1), pind_i(1), 0.5);
                P_Up.set(ind_i(1), pind_i(2), 0.5);
                P_Up.set(ind_i(2), pind_i(1), 0.5);
                P_Up.set(ind_i(2), pind_i(2), 0.5);
                
                P_Dn.set(ind_i(1), pind_i(1), 0.5);
                P_Dn.set(ind_i(1), pind_i(2), -0.5);
                P_Dn.set(ind_i(2), pind_i(1), -0.5);
                P_Dn.set(ind_i(2), pind_i(2), 0.5);
            }
            else Error("Direction '" + direction + "' not allowed");
            
            if (CPS.at(i) == 1)
            {
                auto newA = P_Up * prime(Npsi.A(i),Site);
                auto nr = norm(newA);
                newA = newA/nr;
                Npsi.setA(i, newA);
            }
            
            if (CPS.at(i) == -1)
            {
                auto newA = P_Dn * prime(Npsi.A(i),Site);
                auto nr = norm(newA);
                newA = newA/nr;
                Npsi.setA(i, newA);
            }
            
        }

    return Npsi;
}


//----MEASUREMENTS----
//Accepsts the MPS (METTS), the SiteSet variable of the given MPS to measure and the physical quantities to be measured (<M> for magnetization, )
// u and v are the positions in between the observable will be measured. Both sites included. To avoit the boundaries v has to be greater than u
Real measure(MPS psiM, SiteSet sites, int u, int v,Args const& args = Global::args())
{
    float w = v-u;
    cout<<w<<endl;
    if (w<0)
    {
        Error("value of v has to be greater than the value of u");
    }
    auto N = psiM.N(); // Size of the system
    
    Real O = 0.; // Value to return
    auto observable = args.getString("observable"); // Observable to measure
    
    //Measure magnetization. Local observable can be measured site by site.
    if (observable == "<M>")
    {
        Real M_i = 0.;
        for(int i = u; v >= i; i++)
        {
            // obtain index j in the MPS
            Index ind_i = sites(i);
            Index pind_i = prime(ind_i);
            
            // Defining the operator Sz
            auto Sz_i = ITensor(ind_i, prime(ind_i));
            Sz_i.set(ind_i(1), pind_i(1), 1.0);
            Sz_i.set(ind_i(2), pind_i(2), -1.0);
            
            //Move orthoganality center to position i
            psiM.position(i);
            //Measure spectation value on site i
            M_i = ( dag(prime(psiM.A(i),Site)) * Sz_i * psiM.A(i) ).real();
            // Sume it up
            O += M_i/(w+1);
        }
    }
    return O;
    
}

//----INITIAL CPS----
// Genereate the first CPS of alternating up and down spins 1D chain
MPS initial(int N, SiteSet sites)
{
    auto s = InitState(sites,"Up");
    for(int j = 1; N >= j; j++)
    {
        auto st = (j%2 == 1 ? "Up" : "Dn");
        s.set(j,st);
    }
    auto psi = MPS(s);
    return psi;
}



//----GENERATES THE HAMILTONIAN FOR THE SPECCIFIED PHYSICAL SYSTEM----
//Inputs:
//  -h (value of the magnetic field)
//  -J1 and J2 are the exchange constants
//  -N number of spins
//  -sites is an object with the name of the index for each site
//  -Hamiltonian (string with the name og the hamiltonian to study, choices are ("1D Heisenberg Nearest Neighbour", "1D Heisenberg Second Nearest Neighbour", "Dimer chain", "Triangular lattice 3Ny"))
//
//Outputs:
//  -The hamiltonian exponentiated as an MPO
MPO HAMILTONIAN( Real h, Real J1, Real J2, int N, SiteSet sites, Real tau,Args const& args = Global::args())
{
    auto hamiltonian = args.getString("Hamiltonian"); // Hamiltonian to construct
    cout<<"Making MPO for Hamiltonian: "<<hamiltonian<<". For h = "<<h<<endl;
    
    //1D HEISENBERG CHAIN NEAREST NEIGHBOUR INTERACTIONS
    if (hamiltonian == "1D Heisenberg Nearest Neighbour")
    {
        auto ampo = AutoMPO(sites);
        for(int j = 1; N >= j; j++)
        {
            if (j < N)
            {
                ampo += -0.5*J1,"S+",j,"S-",j+1;
                ampo += -0.5*J1,"S-",j,"S+",j+1;
                ampo += -1.0*J1,"Sz",j,"Sz",j+1;
            }
            ampo += -h, "Sz",j;
        }
        auto H = MPO(ampo);
        cout<<"Done Making MPO for Hamiltonian: "<<endl;
        return toExpH<ITensor>(ampo,tau);
    }
    
    //1D HEISENBERG CHAIN SECOND NEAREST NEIGHBOUR INTERACTIONS
    if (hamiltonian == "1D Heisenberg Second Nearest Neighbour")
    {
        auto ampo = AutoMPO(sites);
        for(int j = 1; N >= j; j++)
        {
            if (j < N)
            {
                ampo += -0.5*J1,"S+",j,"S-",j+1;
                ampo += -0.5*J1,"S-",j,"S+",j+1;
                ampo += -1.0*J1,"Sz",j,"Sz",j+1;
            }
            
            if (j < N-1)
            {
                ampo += -0.5*J2,"S+",j,"S-",j+2;
                ampo += -0.5*J2,"S-",j,"S+",j+2;
                ampo += -1.0*J2,"Sz",j,"Sz",j+2;
            }
            ampo += -h, "Sz",j;
        }
        auto H = MPO(ampo);
        cout<<"Done Making MPO for Hamiltonian: "<<endl;
        return toExpH<ITensor>(ampo,tau);
    }
    
    //DIMER CHAIN
    if (hamiltonian == "Dimer chain")
    {
        if (N%3 > 0)
        {
            Error("Total number of spins must be a multiple of 4. Periodicity for the Dimer chain is 4");
        }
        
        auto ampo = AutoMPO(sites);
        //Horizontal interactions
        for(int j = 5; N-3 >= j; j = j+4)
        {
            ampo += -0.5*J1,"S+",j,"S-",j-1;
            ampo += -0.5*J1,"S-",j,"S+",j-1;
            ampo += -1.0*J1,"Sz",j,"Sz",j-1;
            
            ampo += -0.5*J2,"S+",j,"S-",j+1;
            ampo += -0.5*J2,"S-",j,"S+",j+1;
            ampo += -1.0*J2,"Sz",j,"Sz",j+1;
            
            ampo += -0.5*J2,"S+",j,"S-",j+2;
            ampo += -0.5*J2,"S-",j,"S+",j+2;
            ampo += -1.0*J2,"Sz",j,"Sz",j+2;
        }
        
        for(int j = 4; N >= j; j = j+4)
        {
            ampo += -0.5*J2,"S+",j,"S-",j-1;
            ampo += -0.5*J2,"S-",j,"S+",j-1;
            ampo += -1.0*J2,"Sz",j,"Sz",j-1;
            
            ampo += -0.5*J2,"S+",j,"S-",j-2;
            ampo += -0.5*J2,"S-",j,"S+",j-2;
            ampo += -1.0*J2,"Sz",j,"Sz",j-2;
        }
        
        for(int j = 2; N-2 >= j; j = j+4)
        {
            ampo += -0.5*J1,"S+",j,"S-",j+1;
            ampo += -0.5*J1,"S-",j,"S+",j+1;
            ampo += -1.0*J1,"Sz",j,"Sz",j+1;
        }
        //Left boundary condition
        ampo += -0.5*J2,"S+",1,"S-",2;
        ampo += -0.5*J2,"S-",1,"S+",2;
        ampo += -1.0*J2,"Sz",1,"Sz",2;
        
        ampo += -0.5*J2,"S+",1,"S-",3;
        ampo += -0.5*J2,"S-",1,"S+",3;
        ampo += -1.0*J2,"Sz",1,"Sz",3;
        
        //Magnetic field term
        for(int j = 1; j <= N; j++)
        {
            ampo += -h, "Sz",j;
        }
        
        auto H = MPO(ampo);
        cout<<"Done Making MPO for Hamiltonian: "<<endl;
        return toExpH<ITensor>(ampo,tau);
    }
    
    //Triangular lattice
    if (hamiltonian == "Triangular lattice 3Ny")
    {
        if (N%3 > 0)
        {
            Error("Total number of spins must be a multiple of 3. N/3 = Nx");
        }
        Real Nx = N/3.;
        auto ampo = AutoMPO(sites);
        
        for(int j = N-1; 5 <= j; j = j-3)
        {
            ampo += -0.5*J1,"S+",j,"S-",j-1;
            ampo += -0.5*J1,"S-",j,"S+",j-1;
            ampo += -1.0*J1,"Sz",j,"Sz",j-1;
            
            ampo += -0.5*J1,"S+",j,"S-",j+1;
            ampo += -0.5*J1,"S-",j,"S+",j+1;
            ampo += -1.0*J1,"Sz",j,"Sz",j+1;
            
            ampo += -0.5*J1,"S+",j,"S-",j-2;
            ampo += -0.5*J1,"S-",j,"S+",j-2;
            ampo += -1.0*J1,"Sz",j,"Sz",j-2;
            
            ampo += -0.5*J1,"S+",j,"S-",j-3;
            ampo += -0.5*J1,"S-",j,"S+",j-3;
            ampo += -1.0*J1,"Sz",j,"Sz",j-3;
        }
        
        for(int j = N-2; 4 <= j; j = j-3)
        {
            ampo += -0.5*J1,"S+",j,"S-",j-3;
            ampo += -0.5*J1,"S-",j,"S+",j-3;
            ampo += -1.0*J1,"Sz",j,"Sz",j-3;
            
            ampo += -0.5*J1,"S+",j,"S-",j-2;
            ampo += -0.5*J1,"S-",j,"S+",j-2;
            ampo += -1.0*J1,"Sz",j,"Sz",j-2;
        }
        
        for(int j = N; 6 <= j; j = j-3)
        {
            ampo += -0.5*J1,"S+",j,"S-",j-3;
            ampo += -0.5*J1,"S-",j,"S+",j-3;
            ampo += -1.0*J1,"Sz",j,"Sz",j-3;
            
            ampo += -0.5*J1,"S+",j,"S-",j-2;
            ampo += -0.5*J1,"S-",j,"S+",j-2;
            ampo += -1.0*J1,"Sz",j,"Sz",j-2;
            
            ampo += -0.5*J1,"S+",j,"S-",j-5;
            ampo += -0.5*J1,"S-",j,"S+",j-5;
            ampo += -1.0*J1,"Sz",j,"Sz",j-5;
        }
        
        ampo += -0.5*J1,"S+",1,"S-",2;
        ampo += -0.5*J1,"S-",1,"S+",2;
        ampo += -1.0*J1,"Sz",1,"Sz",2;
        
        ampo += -0.5*J1,"S+",3,"S-",1;
        ampo += -0.5*J1,"S-",3,"S+",1;
        ampo += -1.0*J1,"Sz",3,"Sz",1;
        
        ampo += -0.5*J1,"S+",2,"S-",3;
        ampo += -0.5*J1,"S-",2,"S+",3;
        ampo += -1.0*J1,"Sz",2,"Sz",3;
        
        
        for(int j = 1; j <= N; j++)
        {
            ampo += -h, "Sz",j;
        }
        auto H = MPO(ampo);
        cout<<"Done Making MPO for Hamiltonian: "<<endl;
        return toExpH<ITensor>(ampo,tau);
    }
    return 0;
}

    


//----GENERAL METTS IMPLEMENTATION----
//Inputs:
//  -N number of spins
//  -J1 and J2 are the exchange constants
//  -u and v are the positions in between we measure the magnetization
//  -Hamiltonian (string with the name og the hamiltonian to study, choices are ("1D Heisenberg Nearest Neighbour", "1D Heisenberg Second Nearest Neighbour", "Dimer chain"))
//  -h (value of the magnetic field)
//  -beta (1/T) in dimenssionless units
//  -tau (intervals in which we break betta into)
//  -Malias (METTS that will be run before computing any expectation value)
//  -MMETTS (number of METTS used to average over the desired quantities)
//  -maxm (max bond dimension for the exponentiation part)
//  -cutoff (for exponentiation method)
//
//Outputs:
//  -NMETTS magnetization values as a vector

vector<float> METTS_MAG(int N, Real h, Real J1, Real J2, int u, int v, Real beta, Real tau, int Malias, int MMETTS, Real maxm, Real cutoff, Args const& args = Global::args())
{
    
    auto M = MMETTS + Malias; //Total number of operations
    auto hamiltonian = args.getString("Hamiltonian"); // Hamiltonian to time evolve
    
    cout<<hamiltonian<<endl;
    //Imaginary time evolution total time
    auto ttotal = beta/2.;
    const int nt = int( ttotal/tau + 1e-9*(ttotal/tau) );
    if(fabs(nt*tau-ttotal) > 1E-9)
    {
        Error("Timestep not commensurate with total time");
    }
    
    //Output magnetization vector
    auto Mag = vector<float>(MMETTS);
    
    //Define the sites of the MPS
    auto sites = SpinHalf(N);
    //Set observable to measure and initial projection direction
    auto dir = "z";
    auto obs = "<M>";
    
    //Get the exponentiated hamiltonian to time evolve
    auto expH = HAMILTONIAN( h, J1, J2, N, sites, tau,{"Hamiltonian=", hamiltonian});
    
    
    //Generate the updown initial CPS
    auto psi0 = initial(N, sites);
    auto CPS = initial(N, sites);
    
    //Generate an ensemble of METTS
    int av = 0;
    for(int i = 1; M >= i; i++)
    {
        //Generate the METTS
        cout<<"Making METTS "<<i<<" out of "<<M<<endl;
        float mu = (nt/4.);
        for(int j = 1; j <= nt; ++j)
        {
            if(j > mu)
            {
                Real cutoffF = 10E-8;
                psi0 = applyMPO(expH,psi0,{ "Method=", "DensityMatrix", "Maxm=", maxm, "Cutoff=", cutoffF, "Normalize", true});
            }
            if(j <= mu)
            {
                psi0 = applyMPO(expH,psi0,{ "Method=", "DensityMatrix", "Maxm=", maxm, "Cutoff=", cutoff, "Normalize", true});
            }
            
                
        }
        cout<<"Done making METTS "<<i<<endl;
        
        
        if (i > Malias)
        {
            
            Mag.at(av) = measure(psi0, sites, u, v, { "observable=", obs});
            av++;
            
            
        }
        
        //Decide if collapse to z or x directions
        if (i%2 > 0)
        {
            dir = "z";
        }
        
        if (i%2 == 0)
        {
            dir = "x";
        }
        
        //Collapse to new CPS
        cout<<"Collapsing METTS "<<i<<endl;
        CPS = collapse( psi0, sites,{"Direction=",dir} );
        psi0 = CPS;
        cout<<"Done Collapsing METTS "<<i<<endl;
        
        
    }
    
    
    return Mag;
    
}

