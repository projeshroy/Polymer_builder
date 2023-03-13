//................................................................................................
//Write mol2 file 
//................................................................................................
void write_mol2(std::ofstream& mol2_file, Vec_s& element_names, Mat_d& coordinates, 
	        Vec_s& atom_type_names, Vec_i& molecule_type_indices, Vec_s& molecule_names, 
	        Vec_d& atom_charges, Mat_i& bonds, Vec_s& bond_types)
{
	mol2_file << "@<TRIPOS>MOLECULE" << std::endl;
	mol2_file << "<Molecule name>" << std::endl;
	mol2_file << "  "  << element_names.rows() << "  " << bond_types.rows() << " 0  0  0"  << std::endl;
	mol2_file << "SMALL" << std::endl;
	mol2_file << "USER \n"  << std::endl;
	mol2_file << "@<TRIPOS>ATOM" << std::endl;	

	for(int i = 0 ; i < element_names.rows(); i++)
		mol2_file << " " << std::setfill(' ') << std::setw(10) << (i+1) 
			  << " " << std::setfill(' ') << std::setw(4) << element_names[i] 
			  << " " << std::setfill(' ') << std::setw(10) << coordinates(0, i) 
			  << " " << std::setfill(' ') << std::setw(10) << coordinates(1, i) 
			  << " " << std::setfill(' ') << std::setw(10) << coordinates(2, i) 
			  << " " << std::setfill(' ') << std::setw(15) << atom_type_names[i] 
			  << " " << std::setfill(' ') << std::setw(10) << molecule_type_indices[i] 
			  << " " << std::setfill(' ') << std::setw(10) << molecule_names[i] 
			  << " " << std::setfill(' ') << std::setw(10) << atom_charges[i] 
			  << std::endl;

	mol2_file << "@<TRIPOS>BOND" << std::endl;
	for(int i = 0 ; i < bond_types.rows(); i++)
		mol2_file << " " << std::setfill(' ') << std::setw(10) << (i+1) 
			  << " " << std::setfill(' ') << std::setw(10) << bonds(i, 0)
			  << " " << std::setfill(' ') << std::setw(10) << bonds(i, 1)
			  << " " << std::setfill(' ') << std::setw(10) << bond_types[i]
			  << std::endl;
	
	std::cout << " mol2 file writen " << std::endl;		
}
