/*Electrostatic solver -- Used for the purposes of computing the electric fields between plate capacitors when given a voltage to the cathode
The intent here is to implement a function that accepts voltage as input which then outputs an electric field that
can be used as input for the drift-diffusion perturbation equation problem*/
#include "/g/g11/frye11/mfem/mfem.hpp"
#include "/g/g11/frye11/mfem/fem/NSnonlininteg.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[]){
    const char *meshfile = "../../mfem/data/ref-segment.mesh";
    int order = 2;
    
}