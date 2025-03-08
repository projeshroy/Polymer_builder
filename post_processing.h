#include "replace_type.h"

void post_processing(std::ifstream& input_file,
		int& Polymer_length,
		Mat_s& Element_names_matrix, 
		Mat_s& Mid_Element_names_matrix, 
		Mat_s& Ini_Element_names_matrix, 
		Mat_s& Fin_Element_names_matrix, 
		Mat_d& bond_connection_values_matrix, 
		Mat_d& angle_connection_values_matrix, 
		Mat_d& dihedral_connection_values_matrix,
            	Mat_s& Monomer_atom_type_names_matrix,		
            	Mat_s& Mid_Monomer_atom_type_names_matrix,	
            	Mat_s& Ini_Monomer_atom_type_names_matrix,	
            	Mat_s& Fin_Monomer_atom_type_names_matrix,	
            	Mat_s& Monomer_fragment_names_matrix,		
            	Mat_d& Monomer_atom_charges_matrix,		
            	Mat_i& Monomer_fragment_type_indices_matrix){	

	std::string read_string;
        int post_processing_group_count;
        input_file >> read_string >> post_processing_group_count;

	Vec_i post_processing_count, monomer_count;
	post_processing_count.resize(post_processing_group_count);
	monomer_count.resize(post_processing_group_count);
	Vec_s group_name, replace_type;	
	group_name.resize(post_processing_group_count);
	replace_type.resize(post_processing_group_count);
       	Mat_i monomer_id, complementary_monomer_id; 
	monomer_id.resize(post_processing_group_count, Polymer_length);
	complementary_monomer_id.resize(post_processing_group_count, Polymer_length);
       	monomer_id.setZero(); complementary_monomer_id.setZero();
   
	for(int i = 0; i < post_processing_group_count; i++){	
		input_file >> read_string >> group_name[i] 
			   >> read_string >> post_processing_count[i] 
			   >> read_string >> monomer_count[i] 
			   >> read_string >> replace_type[i]; 

		select_monomer_id(i, input_file, Polymer_length, monomer_count, 
				  group_name, replace_type, monomer_id, complementary_monomer_id);	

		if((Polymer_length-monomer_count[i]) > 0){
		int count = -1;
		for(int c1 = 0; c1 < Polymer_length; c1++){
			bool ignore = false;
			for(int c2 = 0 ; c2 < monomer_count[i]; c2++){
				if((c1+1) == monomer_id(i, c2)){
				ignore = true;
				break;
				}
			}
			if(!ignore){
			count++;
			complementary_monomer_id(i, count) = c1+1;
			}
		}}
		
		std::cout << " Variable replace group name " << group_name[i] 
			  << " Atom count " << post_processing_count[i]
			  << " monomer_count " << monomer_count[i] 
			  << " replace type " << replace_type[i] << std::endl;
		std::cout << " monomer_id " << monomer_id.row(i).leftCols(monomer_count[i]) << std::endl;
		std::cout << " complementary_monomer_id " << complementary_monomer_id.row(i).leftCols(Polymer_length-monomer_count[i]) << std::endl;

		//..........................................................
		input_file >> read_string;//atom_wise
		if(post_processing_count[i] > 0){
		for(int j = 0; j < post_processing_count[i]; j++){
            		int atom_id, fragment_type_index;
			double bond, angle, dihedral, charge;
            		std::string name, atom_type_name, fragment_name;

			input_file >> atom_id >> name >> bond >> angle >> dihedral;
	            	input_file >> atom_type_name >> charge >> fragment_name >> fragment_type_index;
		
	            	std::cout << " atom_id " << atom_id  
            		  << " name " << name << " bond " << bond << " angle " << " dihedral " << dihedral
            		  << " atom_type_name " << atom_type_name << " charge " << charge 
            		  << " fragment_name " << fragment_name << " fragment_type_index " << fragment_type_index << std::endl;
		
            		for(int k = 0; k < monomer_count[i]; k++){
            	    		Element_names_matrix(monomer_id(i, k)-1, atom_id-1) = name;
            	    		Mid_Element_names_matrix(monomer_id(i, k)-1, atom_id-1) = name;	
            	    		Ini_Element_names_matrix(monomer_id(i, k)-1, atom_id-1) = name;
                    		Fin_Element_names_matrix(monomer_id(i, k)-1, atom_id-1) = name;
     
            	    		bond_connection_values_matrix(monomer_id(i, k)-1, atom_id-1) = bond;
            	    		angle_connection_values_matrix(monomer_id(i, k)-1, atom_id-1) = angle;
            	    		dihedral_connection_values_matrix(monomer_id(i, k)-1, atom_id-1) = dihedral;
            	
            	    		Monomer_atom_type_names_matrix(monomer_id(i, k)-1, atom_id-1) = atom_type_name;
            	    		Mid_Monomer_atom_type_names_matrix(monomer_id(i, k)-1, atom_id-1) = atom_type_name;
            	    		Ini_Monomer_atom_type_names_matrix(monomer_id(i, k)-1, atom_id-1) = atom_type_name;
            	    		Fin_Monomer_atom_type_names_matrix(monomer_id(i, k)-1, atom_id-1) = atom_type_name;
     
            	    		Monomer_fragment_names_matrix(monomer_id(i, k)-1, atom_id-1) = fragment_name;
            	    		Monomer_atom_charges_matrix(monomer_id(i, k)-1, atom_id-1) = charge;
            	    		Monomer_fragment_type_indices_matrix(monomer_id(i, k)-1, atom_id-1) = fragment_type_index;	
            		}		
        	}}

		int atom_id_ini, atom_id_fin, fragment_type_index;
		std::string fragment_name;
		input_file >> read_string >> atom_id_ini >> atom_id_fin >> fragment_name >> fragment_type_index;//fragment_wise; fin_index > start_index;

		for(int k = 0; k < monomer_count[i]; k++){
			for(int l = 0; l < (atom_id_fin-atom_id_ini+1); l++){
				Monomer_fragment_names_matrix(monomer_id(i, k)-1, atom_id_ini+l-1) = fragment_name;	
				Monomer_fragment_type_indices_matrix(monomer_id(i, k)-1, atom_id_ini+l-1) = fragment_type_index;
				}
		}

		}
//....
}
