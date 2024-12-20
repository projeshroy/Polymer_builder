#include "declarations.h"
#include "functions.h"                                
#include "read_mol2.h"
#include "write_mol2.h"
#include "read_zmat.h"
#include "write_zmat.h"
#include "read_xyz.h"
#include "variable_replace.h"

int main(int argc, char** argv){

        std::string directory = std::string("./");
            
        std::string input_file_address = getFileAddress(directory, std::string("input_polymer_builder.in"));
        std::ifstream input_file(input_file_address.c_str());    
       
        std::string Polymer_mol2_file_address = getFileAddress(directory, std::string("Polymer_mol2.mol2"));
        std::ofstream Polymer_mol2_file;
        Polymer_mol2_file.open(Polymer_mol2_file_address); 
 
        std::string Polymer_zmat_file_address = getFileAddress(directory, std::string("Polymer_zmat.gzmat"));
        std::ofstream Polymer_zmat_file;
        Polymer_zmat_file.open(Polymer_zmat_file_address); 

        std::string Polymer_charge_file_address = getFileAddress(directory, std::string("Polymer_charge.dat"));
        std::ofstream Polymer_charge_file;
        Polymer_charge_file.open(Polymer_charge_file_address); 
	
//................................................................................................
//Read input.in file 
//................................................................................................

        int Capped_Monomer_atom_count, Total_fragments, Polymer_length, Cap1_atom_count, Cap2_atom_count;
        Vec_i Cap1_atom_indices, Cap2_atom_indices;
     
        std::string read_string, Monomer_zmat_file_address, Monomer_mol2_file_address;
     
	input_file >> read_string >> Monomer_zmat_file_address;
	std::ifstream Monomer_zmat_file(Monomer_zmat_file_address.c_str());
	input_file >> read_string >> Monomer_mol2_file_address;
	std::ifstream Monomer_mol2_file(Monomer_mol2_file_address.c_str());

        input_file >> read_string >> Capped_Monomer_atom_count;
        input_file >> read_string >> Total_fragments;
        input_file >> read_string >> Polymer_length;

        std::cout << " Capped_Monomer_atom_count " << Capped_Monomer_atom_count << std::endl;
        std::cout << " Total_fragments " << Total_fragments << std::endl;
        std::cout << " Polymer_length " << Polymer_length << std::endl;
           
        input_file >> read_string >> Cap1_atom_count;
	std::cout << " Cap1_atom_count " << Cap1_atom_count << " Cap1_atom_indices ";
        Cap1_atom_indices.resize(Cap1_atom_count);
        input_file >> read_string;
        for(int i = 0; i < Cap1_atom_count; i++){
            input_file >> Cap1_atom_indices[i];	
        std::cout << Cap1_atom_indices[i] << " ";
     	}
	std::cout << std::endl;

        input_file >> read_string >> Cap2_atom_count;
        std::cout << " Cap2_atom_count " << Cap2_atom_count << " Cap2_atom_indices ";
        Cap2_atom_indices.resize(Cap2_atom_count);
        input_file >> read_string;
        for(int i = 0; i < Cap2_atom_count; i++){
            input_file >> Cap2_atom_indices[i];	
	std::cout << Cap2_atom_indices[i] << " ";
     	}
	std::cout << std::endl;

        int Polymer_atom_count = Capped_Monomer_atom_count * Polymer_length;
       	std::cout << " Polymer_atom_count " << Polymer_atom_count << std::endl;

//................................................................................................
//Read mol2 file 
//................................................................................................

        Vec_s Element_names, Monomer_atom_type_names, Monomer_fragment_names, Monomer_mol2_bond_types;
        Vec_d Monomer_atom_charges;
        Vec_i Monomer_fragment_type_indices;
	Mat_i Monomer_mol2_bonds; 
	Mat_d Monomer_coordinates;

	read_mol2(Monomer_mol2_file, Element_names, Monomer_coordinates, 
		  Monomer_atom_type_names, Monomer_fragment_type_indices, Monomer_fragment_names, 
		  Monomer_atom_charges, Monomer_mol2_bonds, Monomer_mol2_bond_types);

	int Capped_Monomer_bond_count = Monomer_mol2_bond_types.rows();
	int Polymer_bond_count = Capped_Monomer_bond_count * Polymer_length;

	std::cout << " Capped_Monomer_bond_count " << Capped_Monomer_bond_count << std::endl;
       	std::cout << " Polymer_bond_count " << Polymer_bond_count << std::endl;

//................................................................................................
//Read Zmat file
//Set the three "connection atoms" at the top of the coordinate file, 
//e.g. connection atom 1 (bond), connection atom 2 (angle), connection atom 3 (dihedral)
//Convert the coordinate file to Zmatrix using avogadro. 
//set missing connecting points (r1 to d3) and their variables properly in the Zmatrix file.
//................................................................................................

    
        Vec_s bond_connection_names, angle_connection_names, dihedral_connection_names;
        Vec_i bond_connections, angle_connections, dihedral_connections;
        Vec_d bond_connection_values, angle_connection_values, dihedral_connection_values;
	double Monomer_charge, Monomer_multiplicity;

	read_zmat(Monomer_zmat_file, Capped_Monomer_atom_count, Monomer_charge, Monomer_multiplicity, Element_names, 
	          bond_connections, bond_connection_names, bond_connection_values,
	          angle_connections, angle_connection_names, angle_connection_values,
   	          dihedral_connections, dihedral_connection_names, dihedral_connection_values);

//................................................................................................
//Change names of the capping atoms to remove them in the babel generated polymer xyz file. 
//Dont give any element/atom-type names 'X' in the input Z-matrix/mol2 file
//................................................................................................
    
        Vec_s Ini_Element_names = Element_names;
        Vec_s Fin_Element_names = Element_names;
        Vec_s Mid_Element_names = Element_names;
        Vec_s Ini_Monomer_atom_type_names = Monomer_atom_type_names;
        Vec_s Fin_Monomer_atom_type_names = Monomer_atom_type_names;
        Vec_s Mid_Monomer_atom_type_names = Monomer_atom_type_names;
        
        for(int c1 = 0; c1 < Cap1_atom_count; c1++){ 
        	Mid_Element_names[Cap1_atom_indices[c1]-1] = std::string("X");
	        Fin_Element_names[Cap1_atom_indices[c1]-1] = std::string("X");
	        Mid_Monomer_atom_type_names[Cap1_atom_indices[c1]-1] = std::string("X");
	        Fin_Monomer_atom_type_names[Cap1_atom_indices[c1]-1] = std::string("X");
        }
        for(int c2 = 0; c2 < Cap2_atom_count; c2++){
        	Mid_Element_names[Cap2_atom_indices[c2]-1] = std::string("X");
	        Ini_Element_names[Cap2_atom_indices[c2]-1] = std::string("X");
	        Mid_Monomer_atom_type_names[Cap2_atom_indices[c2]-1] = std::string("X");
	        Ini_Monomer_atom_type_names[Cap2_atom_indices[c2]-1] = std::string("X");
	}

//................................................................................................
//Connection points.
//Make sure connecting atoms are at the top of the Zmatrix list.
//................................................................................................

	int x = 1;
	for(int i = 0; i < x; i++)
		bond_connections[i] -= Capped_Monomer_atom_count;
	x++;
	for(int i = 0; i < x; i++)
		angle_connections[i] -= Capped_Monomer_atom_count;
	x++;
	for(int i = 0; i < x; i++)
		dihedral_connections[i] -= Capped_Monomer_atom_count;

//................................................................................................
//Convert everything to matrix form to facilitate variable replace
//................................................................................................

        Mat_s Element_names_matrix, Mid_Element_names_matrix, Ini_Element_names_matrix, Fin_Element_names_matrix;
        Element_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        Mid_Element_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        Ini_Element_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        Fin_Element_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
     
        Mat_s bond_connection_names_matrix, angle_connection_names_matrix, dihedral_connection_names_matrix;
        bond_connection_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        angle_connection_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        dihedral_connection_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count); 	
     
        Mat_i bond_connections_matrix, angle_connections_matrix, dihedral_connections_matrix;
        bond_connections_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        angle_connections_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        dihedral_connections_matrix.resize(Polymer_length, Capped_Monomer_atom_count);       
     
        Mat_d bond_connection_values_matrix, angle_connection_values_matrix, dihedral_connection_values_matrix;
        bond_connection_values_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        angle_connection_values_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        dihedral_connection_values_matrix.resize(Polymer_length, Capped_Monomer_atom_count);       
     
        Mat_s Monomer_atom_type_names_matrix, Mid_Monomer_atom_type_names_matrix, Ini_Monomer_atom_type_names_matrix, Fin_Monomer_atom_type_names_matrix;
        Monomer_atom_type_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        Mid_Monomer_atom_type_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        Ini_Monomer_atom_type_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        Fin_Monomer_atom_type_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
     
        Mat_s Monomer_fragment_names_matrix;
        Mat_d Monomer_atom_charges_matrix;
        Mat_i Monomer_fragment_type_indices_matrix;
        Monomer_fragment_names_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        Monomer_atom_charges_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
        Monomer_fragment_type_indices_matrix.resize(Polymer_length, Capped_Monomer_atom_count);
     
        for(int i = 0; i < Polymer_length; i++){
        	int monomer_id = i+1;
     
	        for(int j = 0; j < Capped_Monomer_atom_count; j++){ 
        		Element_names_matrix(i, j) = Element_names[j];
	        	Mid_Element_names_matrix(i, j) = Mid_Element_names[j];
	            	Ini_Element_names_matrix(i, j) = Ini_Element_names[j];
        	    	Fin_Element_names_matrix(i, j) = Fin_Element_names[j];
     
			std::string new_bond_name = bond_connection_names[j];
	        	bond_connection_names_matrix(i, j) = new_bond_name.append(std::string("_pol")).append(std::to_string(monomer_id));;
			std::string new_angle_name = angle_connection_names[j];
	      	        angle_connection_names_matrix(i, j) = new_angle_name.append(std::string("_pol")).append(std::to_string(monomer_id));;
			std::string new_dihedral_name = dihedral_connection_names[j];
        	        dihedral_connection_names_matrix(i, j) = new_dihedral_name.append(std::string("_pol")).append(std::to_string(monomer_id));;
     
	                bond_connections_matrix(i, j) = bond_connections[j];
        	        angle_connections_matrix(i, j) = angle_connections[j];
                	dihedral_connections_matrix(i, j) = dihedral_connections[j];
                
	                bond_connection_values_matrix(i, j) = bond_connection_values[j];
        	        angle_connection_values_matrix(i, j) = angle_connection_values[j];
                	dihedral_connection_values_matrix(i, j) = dihedral_connection_values[j];
                
	                Monomer_atom_type_names_matrix(i, j) = Monomer_atom_type_names[j];
        	        Mid_Monomer_atom_type_names_matrix(i, j) = Mid_Monomer_atom_type_names[j];
                	Ini_Monomer_atom_type_names_matrix(i, j) = Ini_Monomer_atom_type_names[j];
	                Fin_Monomer_atom_type_names_matrix(i, j) = Fin_Monomer_atom_type_names[j];
                
        	        Monomer_fragment_names_matrix(i, j) = Monomer_fragment_names[j];
	                Monomer_atom_charges_matrix(i, j) = Monomer_atom_charges[j];
	                Monomer_fragment_type_indices_matrix(i, j) = Monomer_fragment_type_indices[j];	
        	}
        }
     
//................................................................................................
//Read variable replace section in input.in file
//................................................................................................

        input_file >> read_string >> read_string;

        if(read_string == std::string("yes")){
	variable_replace(input_file,
		Polymer_length,
		Element_names_matrix, 
		Mid_Element_names_matrix, 
		Ini_Element_names_matrix, 
		Fin_Element_names_matrix, 
		bond_connection_values_matrix, 
		angle_connection_values_matrix, 
		dihedral_connection_values_matrix,
            	Monomer_atom_type_names_matrix,		
            	Mid_Monomer_atom_type_names_matrix,	
            	Ini_Monomer_atom_type_names_matrix,	
            	Fin_Monomer_atom_type_names_matrix,	
            	Monomer_fragment_names_matrix,		
            	Monomer_atom_charges_matrix,		
            	Monomer_fragment_type_indices_matrix);	
	}

//................................................................................................
//Build Polymer
//................................................................................................
    
        std::cout << " Building polymer ..... " << std::endl;			
     	double Polymer_charge = Monomer_charge * Polymer_length;//Change if required. 
	double Polymer_multiplicity = 1;//Change if required.

        Vec_s Polymer_atom_type_names, Polymer_fragment_names;
        Vec_d Polymer_atom_charges;
        Vec_i Polymer_fragment_type_indices;
        Polymer_atom_type_names.resize(Polymer_atom_count);
        Polymer_fragment_names.resize(Polymer_atom_count);
        Polymer_atom_charges.resize(Polymer_atom_count);
        Polymer_fragment_type_indices.resize(Polymer_atom_count);

        Vec_s Polymer_Element_names;
	Polymer_Element_names.resize(Polymer_atom_count);	

	Vec_s Polymer_bond_connection_names, Polymer_angle_connection_names, Polymer_dihedral_connection_names;
        Polymer_bond_connection_names.resize(Polymer_atom_count);
        Polymer_angle_connection_names.resize(Polymer_atom_count);
        Polymer_dihedral_connection_names.resize(Polymer_atom_count); 	
     
        Vec_i Polymer_bond_connections, Polymer_angle_connections, Polymer_dihedral_connections;
        Polymer_bond_connections.resize(Polymer_atom_count);
        Polymer_angle_connections.resize(Polymer_atom_count);
        Polymer_dihedral_connections.resize(Polymer_atom_count);       
     
        Vec_d Polymer_bond_connection_values, Polymer_angle_connection_values, Polymer_dihedral_connection_values;
        Polymer_bond_connection_values.resize(Polymer_atom_count);
        Polymer_angle_connection_values.resize(Polymer_atom_count);
        Polymer_dihedral_connection_values.resize(Polymer_atom_count);       

	Vec_i Polymer_extra_atom_shift;
        Polymer_extra_atom_shift.resize(Polymer_atom_count);
	Polymer_extra_atom_shift.setZero();

	Vec_s Polymer_mol2_bond_types;
	Mat_i Polymer_mol2_bonds;
	Polymer_mol2_bond_types.resize(Polymer_atom_count*2);
	Polymer_mol2_bonds.resize(Polymer_atom_count*2, 2); 
	

        int Polymer_atom_counter = -1;
	int Polymer_mol2_atom_counter = -1;
	int Polymer_mol2_bond_counter = -1;
	
        for(int i = 0; i < Polymer_length; i++){
		Mat_s new_Element_names_matrix = Mid_Element_names_matrix;
		Mat_s new_Monomer_atom_type_names_matrix = Mid_Monomer_atom_type_names_matrix;

        	if(i == 0){
		new_Element_names_matrix = Ini_Element_names_matrix;
        	new_Monomer_atom_type_names_matrix = Ini_Monomer_atom_type_names_matrix;
		}
	        if(i == (Polymer_length - 1)){
		new_Element_names_matrix = Fin_Element_names_matrix;
		new_Monomer_atom_type_names_matrix = Fin_Monomer_atom_type_names_matrix;
		}

		int atom_shift = Capped_Monomer_atom_count * i;

        	for(int j = 0; j < Capped_Monomer_atom_count; j++){ 
  			
			//With 'X'...
	            	Polymer_atom_counter++;
			Polymer_Element_names[Polymer_atom_counter] = new_Element_names_matrix(i, j);
			Polymer_bond_connections[Polymer_atom_counter] = bond_connections_matrix(i, j) + atom_shift;
			Polymer_angle_connections[Polymer_atom_counter] = angle_connections_matrix(i, j) + atom_shift;
			Polymer_dihedral_connections[Polymer_atom_counter] = dihedral_connections_matrix(i, j) + atom_shift;
			Polymer_bond_connection_names[Polymer_atom_counter] = bond_connection_names_matrix(i, j);
			Polymer_angle_connection_names[Polymer_atom_counter] = angle_connection_names_matrix(i, j);
			Polymer_dihedral_connection_names[Polymer_atom_counter] = dihedral_connection_names_matrix(i, j);
			Polymer_bond_connection_values[Polymer_atom_counter] = bond_connection_values_matrix(i, j);
			Polymer_angle_connection_values[Polymer_atom_counter] = angle_connection_values_matrix(i, j);
			Polymer_dihedral_connection_values[Polymer_atom_counter] = dihedral_connection_values_matrix(i, j);

			//Without 'X'...
			if(new_Element_names_matrix(i, j) != std::string("X")){
			Polymer_mol2_atom_counter++;

			Polymer_extra_atom_shift[Polymer_atom_counter] = Polymer_mol2_atom_counter - Polymer_atom_counter;
		        Polymer_atom_type_names[Polymer_mol2_atom_counter] = new_Monomer_atom_type_names_matrix(i, j);
	            	Polymer_fragment_names[Polymer_mol2_atom_counter] = Monomer_fragment_names_matrix(i, j);
        	    	Polymer_fragment_type_indices[Polymer_mol2_atom_counter] = Monomer_fragment_type_indices_matrix(i, j) + Total_fragments*i;
            		Polymer_atom_charges[Polymer_mol2_atom_counter] = Monomer_atom_charges_matrix(i, j);
			}
		}

		//Edit bond connections...
		for(int j = 0; j < Capped_Monomer_bond_count; j++){
			if((new_Element_names_matrix(i, (Monomer_mol2_bonds(j, 0)-1)) != std::string("X")) &&
			   (new_Element_names_matrix(i, (Monomer_mol2_bonds(j, 1)-1)) != std::string("X"))){
			Polymer_mol2_bond_counter++;
			Polymer_mol2_bonds(Polymer_mol2_bond_counter, 0) = Monomer_mol2_bonds(j, 0) + atom_shift;
			Polymer_mol2_bonds(Polymer_mol2_bond_counter, 0) += Polymer_extra_atom_shift[Polymer_mol2_bonds(Polymer_mol2_bond_counter, 0)-1];
			Polymer_mol2_bonds(Polymer_mol2_bond_counter, 1) = Monomer_mol2_bonds(j, 1) + atom_shift;
			Polymer_mol2_bonds(Polymer_mol2_bond_counter, 1) += Polymer_extra_atom_shift[Polymer_mol2_bonds(Polymer_mol2_bond_counter, 1)-1];
			Polymer_mol2_bond_types[Polymer_mol2_bond_counter] = Monomer_mol2_bond_types[j];
			}
		}

		//Extra bond for polymer connection.
		if(i > 0){
		Polymer_mol2_bond_counter++;		
		Polymer_mol2_bonds(Polymer_mol2_bond_counter, 0) = 1 + atom_shift;
		Polymer_mol2_bonds(Polymer_mol2_bond_counter, 0) += Polymer_extra_atom_shift[Polymer_mol2_bonds(Polymer_mol2_bond_counter, 0)-1];
		Polymer_mol2_bonds(Polymer_mol2_bond_counter, 1) = bond_connections_matrix(i, 0) + atom_shift;
		Polymer_mol2_bonds(Polymer_mol2_bond_counter, 1) += Polymer_extra_atom_shift[Polymer_mol2_bonds(Polymer_mol2_bond_counter, 1)-1];
		Polymer_mol2_bond_types[Polymer_mol2_bond_counter] = std::string("1");
		}
	}

//................................................................................................
//Write Z-matrix file
//................................................................................................

	write_zmat(Polymer_zmat_file, Polymer_atom_count, Polymer_charge, Polymer_multiplicity, Polymer_Element_names, 
	           Polymer_bond_connections, Polymer_bond_connection_names, Polymer_bond_connection_values,
	           Polymer_angle_connections, Polymer_angle_connection_names, Polymer_angle_connection_values,
   	           Polymer_dihedral_connections, Polymer_dihedral_connection_names, Polymer_dihedral_connection_values);

//................................................................................................
//Convert to xyz format with obabel. 
//Please install babel from http://openbabel.org/wiki/Category:Installation
//................................................................................................
	
	Mat_d Polymer_coordinates;	
	system("obabel Polymer_zmat.gzmat -O Polymer_xyz.xyz");  	

        std::string Polymer_xyz_file_address = getFileAddress(directory, std::string("Polymer_xyz.xyz"));
        std::ifstream Polymer_xyz_file(Polymer_xyz_file_address.c_str());    
	Vec_s new_Polymer_Element_names;
	Mat_d new_Polymer_coordinates;

	read_xyz(Polymer_xyz_file, new_Polymer_Element_names, new_Polymer_coordinates);

//................................................................................................
//Write Polymer mol2 file
//................................................................................................

	Vec_s new_Polymer_atom_type_names = Polymer_atom_type_names.topRows(Polymer_mol2_atom_counter+1);
	Vec_i new_Polymer_fragment_type_indices = Polymer_fragment_type_indices.topRows(Polymer_mol2_atom_counter+1);
	Vec_s new_Polymer_fragment_names = Polymer_fragment_names.topRows(Polymer_mol2_atom_counter+1);
	Vec_d new_Polymer_atom_charges = Polymer_atom_charges.topRows(Polymer_mol2_atom_counter+1);
	Mat_i new_Polymer_mol2_bonds = Polymer_mol2_bonds.topRows(Polymer_mol2_bond_counter+1);
	Vec_s new_Polymer_mol2_bond_types = Polymer_mol2_bond_types.topRows(Polymer_mol2_bond_counter+1);

	write_mol2(Polymer_mol2_file, new_Polymer_Element_names, new_Polymer_coordinates, 
	       	   new_Polymer_atom_type_names, new_Polymer_fragment_type_indices, new_Polymer_fragment_names, 
		   new_Polymer_atom_charges, new_Polymer_mol2_bonds, new_Polymer_mol2_bond_types);

//................................................................................................
//Write charge file
//................................................................................................
	
	for(int i = 0; i < new_Polymer_atom_charges.size(); i++)
		Polymer_charge_file << new_Polymer_atom_charges[i] << std::endl;

	std::cout << "Total charge on the polymer = " << new_Polymer_atom_charges.sum() << std::endl;

}
