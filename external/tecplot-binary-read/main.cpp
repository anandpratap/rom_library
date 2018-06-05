/*
  Copyright (C) 2015 Philippe Miron
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "tecplotread.hpp"
#include "../armadillo/include/armadillo"
#include "../cxxopts/include/cxxopts.hpp"

// this is a really simple main file

void compile_snapshots(int start, int end, std::string directory, bool residual){
	int zone_id = 0;
	
	arma::mat x(64000, end-start+1);
	arma::mat x_residual(64000, end-start+1);

	
	std::vector<std::string> v = {"Static_Pressure","U","V","Temperature","CH4_mf","O2_mf","H2O_mf","CO2_mf"};
	std::vector<std::string> v_residual = {"Eq_1_residual", "Eq_2_residual", "Eq_3_residual", "Eq_4_residual", "Eq_5_residual", "Eq_6_residual", "Eq_7_residual", "Eq_8_residual"};

	for(int time=start; time<= end; time++){
		std::cout<<"Reading "<< (time-start) <<" of "<< (end-start+1)<< "\r";
		int idx = time - start;
		std::string filename = directory + "/gemsma_cmb_"+std::to_string(time)+".plt";
		shared_ptr<tecplotread> tpobj = make_shared<tecplotread>(filename);
		std::vector<int> v_map(8);
		std::vector<int> v_map_residual(8);
		for (size_t j(0); j<tpobj->zones[zone_id]->data_double.size(); j++){
			for(int k=0; k<8; k++){
				if(tpobj->variable_names[tpobj->zones[zone_id]->variable_index[1][j]] == v[k]){
					v_map[k] = j;
				}
				if(residual){
					if(tpobj->variable_names[tpobj->zones[zone_id]->variable_index[1][j]] == v_residual[k]){
						v_map_residual[k] = j;
					}
				}

			}
		}
		
		int n = tpobj->zones[zone_id]->data_double[v_map[0]].size();
		for(int i=0; i<n; i++){
			for(int k=0; k<8; k++){
				x(i*8+k, idx) = tpobj->zones[zone_id]->data_double[v_map[k]][i];
				if(residual){
					x_residual(i*8+k, idx) = tpobj->zones[zone_id]->data_double[v_map_residual[k]][i];
				}
			}
		}
	}

	std::string output_filename = "snapshots_solution.bin";
	x.save(output_filename, arma::arma_binary);
	if(residual){
		x_residual.save(output_filename+"_residual", arma::arma_binary);
	}
	std::cout<<"done!"<< std::endl;
		
}

int main(int argc, char** argv)
{
	cxxopts::Options options("MyProgram", "One line description of MyProgram");
	options.add_options()("residual", "Also save residual");
	options.add_options()("directory", "Directory of plt files", cxxopts::value<std::string>()->default_value(""));
	options.add_options()("start", "start", cxxopts::value<int>()->default_value("0"));
	options.add_options()("end", "end", cxxopts::value<int>()->default_value("0"));
	options.add_options()("help", "Print help");

	auto parser = options.parse(argc, argv);

	int start = parser["start"].as<int>();
	int end = parser["end"].as<int>();
	std::string directory = parser["directory"].as<std::string>();
	compile_snapshots(start, end, directory, parser["residual"].as<bool>());
	return 0;
}
