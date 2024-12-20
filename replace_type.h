void select_monomer_id(int& i,
		std::ifstream& input_file,
		int& Polymer_length,
		Vec_i& monomer_count,
		Vec_s& group_name,
		Vec_s& replace_type,
		Mat_i& monomer_id, 
		Mat_i& complementary_monomer_id){	

		//Selecting Monomer id for the replacement group.........
		// replace_type = manual, all, from, sameas, complto, random.
		//monomer_count variable must be equal to the monomer counts in the following entries. 
		std::string read_string;

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
		else if(replace_type[i] == std::string("shift")){//Shift the monomer index w.r.t previous group

		std::string shift_group_name;		
		int shift;
		input_file >> shift_group_name >> shift;
std::cout << shift_group_name << "  " << shift;
		for(int g = 0; g < i; g++){
			if(shift_group_name == group_name[g]){
			for(int j = 0; j < Polymer_length; j++)
				monomer_id(i, j) = monomer_id(g, j) + shift;		
			break;
			}
		}
		}
	
//....
}
