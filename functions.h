std::string getFileAddress(std::string& directory, std::string file_name)
{
	std::string temp_directory = directory;
	temp_directory.append(file_name);
	return (temp_directory);
}

double getRand(double range_ini, double range_fin)
{
	std::mt19937 generator = std::mt19937(std::random_device()());
	std::uniform_real_distribution<double> distribution(range_ini, range_fin);
	double _rand = distribution(generator);

	return (_rand);
}

                                
