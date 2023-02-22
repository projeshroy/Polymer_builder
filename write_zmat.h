void write_zmat(std::ofstream& zmat_file, int& total_atoms, double& charge, double& multiplicity, Vec_s& element_names, 
	        Vec_i& bond_id, Vec_s& bond_name, Vec_d& bond_value,	   	
	        Vec_i& angle_id, Vec_s& angle_name, Vec_d& angle_value, 	   	
	        Vec_i& dihedral_id, Vec_s& dihedral_name, Vec_d& dihedral_value)
{

	zmat_file << "#\n\n  \n" <<std::endl;//Header.
        zmat_file << charge << "  " << multiplicity << std::endl;//Charge and multplicity.  

	for(int i = 0; i < total_atoms; i++){
	zmat_file << element_names[i];
	if(bond_id[i] > 0)
	zmat_file << "   " << bond_id[i] << "   " << bond_name[i];
	if(angle_id[i] > 0)
	zmat_file << "   " << angle_id[i] << "   " << angle_name[i];
	if(dihedral_id[i] > 0)
	zmat_file << "   " << dihedral_id[i] << "   " << dihedral_name[i];
      	zmat_file << std::endl;
	}

	zmat_file << "Variables:"<< std::endl;

	for(int i = 0; i < total_atoms; i++){
		zmat_file << bond_name[i] << "= " << bond_value[i] << std::endl;
		zmat_file << angle_name[i] << "= " << angle_value[i] << std::endl;
		zmat_file << dihedral_name[i] << "= " << dihedral_value[i] << std::endl;
	}

       std::cout << " Z-matrix written " << std::endl;			 
}
 
