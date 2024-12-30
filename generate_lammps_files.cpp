#include "declarations.h"
#include "functions.h"                                
#include "read_mol2.h"
#include "write_mol2.h"
#include "read_zmat.h"
#include "write_zmat.h"
#include "read_xyz.h"

int main(int argc, char** argv){

	std::string VMDCODE="/include/vmd-1.9.3/build/bin/./vmd -dispdev text ./Packed_mol2.mol2 < vmdscriptfile.tcl";
	std::string VMD=INSTDIR.append(VMDCODE);
	system(VMD.c_str());	
}	
