void residue_replace(std::ifstream& input_file,
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
        int variable_replace_group_count;
        input_file >> read_string >> variable_replace_group_count;

	Vec_i variable_replace_count, monomer_count;
	variable_replace_count.resize(variable_replace_group_count);
	monomer_count.resize(variable_replace_group_count);
	Vec_s group_name, replace_type;	
	group_name.resize(variable_replace_group_count);
	replace_type.resize(variable_replace_group_count);
       	Mat_i monomer_id, complementary_monomer_id; 
	monomer_id.resize(variable_replace_group_count, Polymer_length);
	complementary_monomer_id.resize(variable_replace_group_count, Polymer_length);
       	monomer_id.setZero(); complementary_monomer_id.setZero();
   
	for(int i = 0; i < variable_replace_group_count; i++){	
		input_file >> group_name[i] >> variable_replace_count[i] >> monomer_count[i] >> replace_type[i]; 

		//Selecting Monomer id for the replacement group.........
		// replace_type = manual, all, from, sameas, complto, random.
		//monomer_count variable must be equal to the monomer counts in the following entries. 

		if(replace_type[i] == std::string("manual")){//manual entry
		for(int j = 0; j < monomer_count[i]; j++)
			input_file >> monomer_id(i, j);
		}
		else if (replace_type[i] == std::string("all")){//throughout the polymer chain
		for(int j = 0; j < Polymer_length; j++)
                        monomer_id(i, j) = j+1;
		}
		else if (replace_type[i] == std::string("from")){//from an integer to another. Both ini and fin must be less than Polymer_length
		int ini, fin;
		input_file >> ini >> read_string >> fin;
		for(int j = 0; j < (fin-ini+1); j++)
			monomer_id(i, j) = ini+j;
		}
		else if (replace_type[i] == std::string("sameas")){//same as a previous group which is already given!
		input_file >> read_string;
		for(int g = 0; g < i; g++){
			if(read_string == group_name[g]){
			monomer_id.row(i) = monomer_id.row(g);
			break;
			}
		}}
		else if (replace_type[i] == std::string("complto")){
		input_file >> read_string;

                for(int g = 0; g < i; g++){
                        if(read_string == group_name[g]){
                        monomer_id.row(i) = complementary_monomer_id.row(g);
			break;
			}
                }}
            	else if(replace_type[i] == std::string("random")){//random selection
            	for(int j = 0; j < monomer_count[i]; j++){
            		int rand = 0;
            		bool check = true;
     
            		while(check){
            		rand = std::floor(getRand(1, Polymer_length));
			check = false;
            		for(int k = 0; k < monomer_count[i]; k++){
            		    	if(rand == monomer_id(i, k)){
            		    	check = true;
            		    	break;
            			}
			}
            		if(!check)
            		monomer_id(i, j) = rand;
            	     }
            	}}
		else if(replace_type[i] == std::string("random_exclude")){//random selection but excluding some previously mentioned group.
		int exclude_group_count;
		input_file >> exclude_group_count;

		Vec_i exclude_group_id;
		exclude_group_id.resize(exclude_group_count); exclude_group_id.setZero();
		Mat_i exclude_monomer_id;
		exclude_monomer_id.resize(exclude_group_count, Polymer_length);	exclude_monomer_id.setZero();	

		for(int g1 = 0; g1 < exclude_group_count; g1++){
			input_file >> read_string;

			for(int g2 = 0; g2 < i; g2++){
				if(read_string == group_name[g2]){
				exclude_group_id[g1] = g2;
				exclude_monomer_id.row(g1) = monomer_id.row(g2);		
				break;
				}
			}
		}

		for(int j = 0; j < monomer_count[i]; j++){
			int rand = 0;
			bool check = true;
        
			while(check){
			rand = std::floor(getRand(1, Polymer_length));
		    	check = false;
        
			for(int k = 0; k < monomer_count[i]; k++){
			    if(rand == monomer_id(i, k)){
			    check = true;
			    break;
			    }
			}

			if(!check){
			for(int g1 = 0; g1 < exclude_group_count; g1++){
				for(int g2 = 0; g2 < monomer_count[exclude_group_id[g1]]; g2++){
					if(rand == exclude_monomer_id(g1, g2)){
					check = true;
                            		break;
                            		} 
				}
			}}

            		if(!check)
            		monomer_id(i, j) = rand;
            	     	}}
		}

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
			  << " monomer_count " << monomer_count[i] 
			  << " replace type " << replace_type[i] << std::endl;
		std::cout << " monomer_id " << monomer_id.row(i).leftCols(monomer_count[i]) << std::endl;
		std::cout << " complementary_monomer_id " << complementary_monomer_id.row(i).leftCols(Polymer_length-monomer_count[i]) << std::endl;
		//..........................................................
		for(int j = 0; j < variable_replace_count[i]; j++){
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
        	}		
	}
//....
}
