void read_xyz(std::ifstream& xyz_file, Vec_s& element_names, Mat_d& coordinates)
{
	int total_atoms;
	std::string read_string;
	xyz_file >> total_atoms;
	element_names.resize(total_atoms);
	coordinates.resize(3, total_atoms);

	xyz_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//<Molecule name>
	xyz_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');//<Molecule name>

	for(int i = 0; i < total_atoms; i++)
		xyz_file >> element_names[i] >> coordinates(0, i) >> coordinates(1, i) >> coordinates(2, i);
}

