//................................................................................................
//Read (Edited!) Z-matrix file 
//Please edit the connection ofthe first three atoms, which are usually left blancks.
//e.g.,
//C 0 r1 0 r2 0 d2
//C 1 r2 0 a2 0 d2
//C 2 r3 1 a3 0 d3
//etc.
//Variables:
//r1= 0
//a1= 0
//d1= 0
//etc.
//................................................................................................

void read_zmat(std::ifstream& zmat_file, int& total_atoms, double& charge, double& multiplicity, Vec_s& element_names, 
	       Vec_i& bond_id, Vec_s& bond_name, Vec_d& bond_value,	   	
	       Vec_i& angle_id, Vec_s& angle_name, Vec_d& angle_value, 	   	
	       Vec_i& dihedral_id, Vec_s& dihedral_name, Vec_d& dihedral_value)
{ 	  
	std::string read_string;

	element_names.resize(total_atoms);
	bond_id.resize(total_atoms);
	bond_name.resize(total_atoms);
	bond_value.resize(total_atoms);
	angle_id.resize(total_atoms);
	angle_name.resize(total_atoms);
	angle_value.resize(total_atoms);
 	dihedral_id.resize(total_atoms);
	dihedral_name.resize(total_atoms);
	dihedral_value.resize(total_atoms);

        zmat_file >> read_string; //Remove everything after # in the gzmat file. If you want to keep it, edit this line. 
                                          //This is just so we can view it easily with avogadro. 
        zmat_file >> charge >> multiplicity;
            
        for(int i = 0; i < total_atoms; i++)
        	zmat_file >> element_names[i] 
               		  >> bond_id[i] >> bond_name[i] 
                	  >> angle_id[i] >> angle_name[i] 
                     	  >> dihedral_id[i] >> dihedral_name[i];
            
        zmat_file >> read_string;// "Variables:"
        for(int i = 0; i < total_atoms; i++){
	        zmat_file >> read_string >> bond_value[i];//Imp: no space between 'rxx' and '='. e.g., r1= 0.0
        	zmat_file >> read_string >> angle_value[i];
	        zmat_file >> read_string >> dihedral_value[i];
	}
            
        std::cout << " Z-matrix file read " << std::endl;		 
}
