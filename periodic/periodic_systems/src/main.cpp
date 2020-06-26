// c++ xyz_periodic_handler.cpp -o main && ./main gnr_7_periodic -p X 4.26 -t 1.8
// ./build/Application gnr_7_periodic.xyz -p X 4.26 -t 1.8

#include <algorithm> //std::sort()
#include <complex>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip> //std::setw()
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <string>
#include <typeinfo>
#include <vector>

#include "TBTK/Array.h"
#include "TBTK/BrillouinZone.h"
#include "TBTK/Functions.h"
#include "TBTK/Model.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Vector3d.h"

using namespace TBTK;
using comperator_t = std::function<bool(double, double)>;
using bonds_t = std::vector<std::pair<int, int>>;
using std::string;
using std::size_t;
using std::vector;
using std::complex;
using std::ofstream;

const double k_eps { 0.0001 };
const double atomic_radii_C { 0.68 };
const double threshold { 1e-4 };
const int k_size_brillouin_zone{ 500 };
size_t k_num_atoms_unit_cell { 0 };
std::complex<double> i(0, 1);
Array<double> spinAndSiteResolvedDensity;
double U { 0.0 };
double spin_and_site_resolved_density_tol { 1e-5 };
double DENSITY_TOLLERANCE { 1e-7 };
double TARGET_DENSITY_PER_SITE { 1.0 };
double MIXING_PARAMETER { 0.5 };
Model model;

static int print_helper_01{ 0 };
static int print_helper_02{ 0 };


#define PRINTVAR(x) std::cout << #x << " = " << (x) << std::endl;
void print_help(bool full);

struct XYZCoordinate{
    float x;
    float y;
    float z;

    friend std::ostream& operator<<(std::ostream& out, XYZCoordinate &xyz) {
        out << std::setw(7) << xyz.x << "\ty: " << std::setw(7)
            << xyz.y << "\tz: " << std::setw(7) << xyz.z << "\n";
        return out;
    }
};


struct XYZFileLine{
    XYZCoordinate xyz;
    std::string name;

    friend std::istream& operator>>(std::istream& in, XYZFileLine &xyz_file) {
        in >> xyz_file.name >> xyz_file.xyz.x
           >> xyz_file.xyz.y >> xyz_file.xyz.z;
        return in;
    }
    friend std::ostream& operator<<(std::ostream& out, XYZFileLine &xyz_file) {
        out << "Atom: " << std::setw(2) << xyz_file.name << "\tx: "
            << std::setw(7) << xyz_file.xyz.x << "\ty: " << std::setw(7)
            << xyz_file.xyz.y << "\tz: " << std::setw(7) << xyz_file.xyz.z;
        return out;
    }
};


/** One atom consists of XYZFileLine (atom name and coordiantes)
 *  and an unique id to distingush and access the atoms
**/
class Atom{
public:
    Atom() : id_(counter_++) {};

    friend std::ostream& operator<<(std::ostream& out, Atom &atom) {
        out << atom.name_coords_ << "\n";
        return out;
    }

    friend std::istream& operator>>(std::istream& in, Atom &atom) {
        in >> atom.name_coords_.name >> atom.name_coords_.xyz.x
           >> atom.name_coords_.xyz.y >> atom.name_coords_.xyz.z;
        return in;
    }

    std::ostream& print_atom(std::ostream& out){
        out << "Id: " << std::setw(3) << id_ << "\tAtom: " << std::setw(7)
            << name_coords_ << std::setw(7);
        return out;
    }

    XYZFileLine name_coords_;

private:
    friend class Molecule;
    const int id_;
    static int counter_;
};
int Atom::counter_ = 0;

class Molecule{
public:
    Molecule(std::istream& in){
        std::string n_atom;
        getline(in, n_atom);
        std::cout << "Number of atoms in file: " << n_atom << '\n';

        std::string comment_line;
        getline(in, comment_line);

        Atom atom;
        while (in >> atom) {
            if(atom.name_coords_.name != "C") continue;
            atoms_.push_back(atom);
            val_x.push_back(atom.name_coords_.xyz.x);
            val_y.push_back(atom.name_coords_.xyz.y);
        }
        set_min_max();
    }

    Molecule(Molecule other_mol,
            std::vector<int> indices, double periodicity_distance){
        std::cout << "Mol fed to func:" << '\n'<< other_mol << '\n';
        for(auto index : indices){
            atoms_.push_back(other_mol.atoms_[index]);
            val_x.push_back(other_mol.atoms_[index].name_coords_.xyz.x);
            val_y.push_back(other_mol.atoms_[index].name_coords_.xyz.y);
            std::cout << "Printing atom " << atoms_.size()-1 << " " << atoms_.at(atoms_.size()-1);
        }
        for(auto index : indices){
            atoms_.push_back(other_mol.atoms_[index]);
            atoms_[atoms_.size()-1].name_coords_.xyz.x += periodicity_distance;
            val_x.push_back(other_mol.atoms_[atoms_.size()-1].name_coords_.xyz.x);
            val_y.push_back(other_mol.atoms_[index].name_coords_.xyz.y);
            std::cout << "Printing atom " << atoms_.size()-1 << " " << atoms_.at(atoms_.size()-1);
        }
        set_min_max();
    }

    void set_min_max(){
        std::vector<double> temp_x = val_x;
        comperator_t compFunctor = [](double elem1, double elem2){
				return elem1 < elem2;
		};
        std::sort(temp_x.begin(), temp_x.end(), compFunctor);

        min_x = temp_x[0];
        max_x = temp_x[temp_x.size()-1];

        std::vector<double> temp_y = val_y;
        std::sort(temp_y.begin(), temp_y.end(), compFunctor);

        min_y = temp_y[0];
        max_y = temp_y[temp_y.size()-1];
    }

    size_t size(){
        return atoms_.size();
    }

    friend std::ostream& operator<<(std::ostream& out, const Molecule& atoms){
        for(auto i : atoms.atoms_){
            out << i;
        }
        return out;
    }

    double get_x_coords(int index){
        return atoms_[index].name_coords_.xyz.x;
    }
    double get_y_coords(int index){
        return atoms_[index].name_coords_.xyz.y;
    }
    double get_x_max(){
        return max_x;
    }
    double get_x_min(){
        return min_x;
    }
    double get_y_max(){
        return max_y;
    }
    double get_y_min(){
        return min_y;
    }

private:
    std::vector<Atom> atoms_;
    std::vector<double> val_x;
    std::vector<double> val_y;
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    friend class Bonds;
    friend class Atom;
};

class Bonds{
public:

    void print_bonds(){
        int bonds_counter = 0;
        std::cout << "-------------------------------------------------------\n"
                  << "Bonds added: " << '\n';
        for(auto i : bonds_){
            std::cout << "{" << i.first << "," << i.second << "}  ";
            ++bonds_counter;
        }
        std::cout << "\nTotal bonds created: " << bonds_counter
                  << "\n-------------------------------------------------------\n";
    }

    Bonds(Molecule mol, double threshold, complex<double> t){
        bond_threshold_ = threshold;
        t_ = t;
        std::cout << "Size of the considered molecule (no H): " << mol.size() << "\n";
        for(size_t i = 0; i < mol.size(); ++i){
            double x_1 = mol.atoms_[i].name_coords_.xyz.x;
            double y_1 = mol.atoms_[i].name_coords_.xyz.y;
            double z_1 = mol.atoms_[i].name_coords_.xyz.z;

            for (size_t j = i; j < mol.size(); j++){
                if(i == j) continue;
                double x_2 = mol.atoms_[j].name_coords_.xyz.x;
                double y_2 = mol.atoms_[j].name_coords_.xyz.y;
                double z_2 = mol.atoms_[j].name_coords_.xyz.z;
                double x_diff = x_1 > x_2 ? x_1 - x_2 : x_2 - x_1;
                double y_diff = y_1 > y_2 ? y_1 - y_2 : y_2 - y_1;
                double z_diff = z_1 > z_2 ? z_1 - z_2 : z_2 - z_1;
                double result = std::sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);

                // if two atoms are inside the sum of their atomic radii
                // plus a threshold that is provided by the user
                // there is a bond between them
                if(result < bond_threshold_){
                    // if(i % 2 == 1 && j % 2 == 1){
                    //     continue;
                    // }
                    bonds_.push_back(std::make_pair(i, j));
                }
            }
        }
    }

    void print_hoppings(){
        int hoppings_counter = 0;
        for(auto i : bonds_){
            std::cout << "model << HoppingAmplitude(-t, {" << std::setw(3)
                      << i.first << ",s}, {" << std::setw(3)
                      << i.second << ",s})+HC;" << std::endl;
            ++hoppings_counter;
        }
        std::cout << "Amount of hoppings created: " << hoppings_counter << '\n';
    }

    bonds_t get_bonds(){
        return bonds_;
    }

    size_t size(){
        return bonds_.size();
    }

    bonds_t bonds_;
private:
    double bond_threshold_;
    complex<double> t_;
    friend class Molecule;
};


bool bond_is_in_unit_cell(bonds_t bonds_in_unit_cell,
        int first_atom, int second_atom){
    for(auto bond : bonds_in_unit_cell){
        if(bond.first == first_atom && bond.second == second_atom){
            return true;
        }
    }
    return false;
}


enum h_state_enum {h_unit_cell_enum,
    h_border_crossing_enum,
    h_both_enum,
    undefined_enum
};

std::ostream& operator<<(std::ostream& os, h_state_enum c)
{
    switch(c)
    {
        case h_both_enum: os << "Both Unit Cell & Border Crossing"; break;
        case h_unit_cell_enum: os << "Unit Cell"; break;
        case h_border_crossing_enum : os << "Border Crossing"; break;
        case undefined_enum : os << "Undefined"; break;
        default : os.setstate(std::ios_base::failbit);
    }
    return os;
}


//Initialize the spin and site resolved density with random numbers between 0
//and 1.
void initSpinAndSiteResolvedDensity(Array<double>& spinAndSiteResolvedDensity){
	srand(time(nullptr));
	for(unsigned int spin = 0; spin < 2; spin++){
		for(unsigned int site = 0; site < k_num_atoms_unit_cell; site++){
			spinAndSiteResolvedDensity[
				{spin, site}
			] = (rand()%100)/100.0; //1.0;
		}
	}
}

//Callback function responsible for returning the current value of H_U for the
//given indices. Since H_U is supposed to be diagonal, toIndex and fromIndex is
//assumed to be equal and only the fromIndex is used to determine the spin and
//site. Compare to the model specification where H_C only is passed to the
//model with both the 'to' and 'from' indices equal.
complex<double> H_U(const Index &toIndex, const Index &fromIndex){
	unsigned int spin = fromIndex[3];
	int site = fromIndex[4];

    // std::cout << "Printing in H_U: " << '\n';
    // PRINTVAR((spin + 1)%2);
    // PRINTVAR(site);
    //PRINTVAR(spinAndSiteResolvedDensity[{(spin + 1)%2, site}]);

	return U*spinAndSiteResolvedDensity[{(spin + 1)%2, site}];
}

//Adjusts the models chemical potential to fix the density. The algorithm is
//implemented as a binary search. The density is calculated, and depending on
//whether the density is too small or too large, the chemical potential is
//increased or decreased, respectively. This is iterated with a doubled step
//length until the calculation overshoots, at which point the procedure
//continues as before, but now cutting the step in two each iteration. The
//procedure stops once a density is found that differes from the
//TARGET_DENSITY_PER_SITE by at most DENSITY_TOLLERANCE.
void fixDensity(PropertyExtractor::BlockDiagonalizer &propertyExtractor){
    double stepLength = 1;

	//Get the eigenvalues.
	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();

    std::cout << "----------- Reached position 1 -----------" << '\n';
	//Perform binary search. The flags stepDirection and hasOvershot
	//indicates in which direction the previous step was taken and whether
	//the calculation has yet overshot. The initial step can be taken in
	//either direction, therefore stepDirection is initialized to zero.
	//Once the step direction has been set in the first iteration of the
	//loop, the first change in step direction will indicate that the step
	//has overshot the target value. At this point it is time to start
	//halving the step size since we are in the vicinity of the target
	//density. Before the overshot occurs, the step size is instead doubled
	//each step to rapidly cause the overshot independently of the initial
	//stepLength.
	int stepDirection = 0;
	bool hasOvershot = false;

    const int number_of_sites = model.getBasisSize()/(2 * k_size_brillouin_zone);
    std::cout << "number of sites: " << number_of_sites << '\n';
	while(true){
		//Calculate the density per unit cell.
		double total_num_el = 0;
		for(int i_spin = 0; i_spin < 2; i_spin++){
    		for(int i_state = 0; i_state < number_of_sites; i_state++){
                for(int i_k = 0; i_k < k_size_brillouin_zone; i_k++){
                    // We are looping through both spin channels
                    // and we add the occupation of each orbital
                    // giving us the total number of electron


                    total_num_el
                    		+= Functions::fermiDiracDistribution(
                    			propertyExtractor.getEigenValue({i_k, 0, 0, i_spin}, i_state),
                    			model.getChemicalPotential(),
                    			model.getTemperature()
                    		);
                    double occ = Functions::fermiDiracDistribution(
                        propertyExtractor.getEigenValue({i_k, 0, 0, i_spin}, i_state),
                        model.getChemicalPotential(),
                        model.getTemperature()
                    );

                    // if(print_helper_01++ < 25){
                    //     PRINTVAR(i_spin);
                    //     PRINTVAR(i_k);
                    //     PRINTVAR(i_state);
                    //     PRINTVAR(propertyExtractor.getEigenValue({i_k, 0, 0, i_spin}, i_state));
                    //     PRINTVAR(occ);
                    // }
                    // std::cout << i_spin << " " << i_k << " " << i_state << " "
                    //     << eigenValues({i_k, 0, 0, i_spin}, i_state) << " "
                    //     << std::endl;

                }
            }
		}


        double el_per_site_per_k_point = total_num_el / (number_of_sites * k_size_brillouin_zone);

        // if(print_helper_02++ < 5){
        std::cout << "El per site " << el_per_site_per_k_point << std::endl;
        std::cout << "Chem pot " << model.getChemicalPotential() << std::endl;
        // }

		//Exit the loop if the target density is met within the given
		//tolerance.
        PRINTVAR(el_per_site_per_k_point);
        PRINTVAR(TARGET_DENSITY_PER_SITE);
        PRINTVAR(abs(el_per_site_per_k_point - TARGET_DENSITY_PER_SITE));
        PRINTVAR(DENSITY_TOLLERANCE);

		if(
			abs(el_per_site_per_k_point - TARGET_DENSITY_PER_SITE)
			< DENSITY_TOLLERANCE
		){
            //std::cout << "Density per site: " << densityPerUnitCell << '\n';
            //std::cout << "Final chemical potetntial: " << model.getChemicalPotential() << '\n';
			break;
		}


		//Determine whether an overshot has occured and step the chemical
		//potential.
		if(el_per_site_per_k_point < TARGET_DENSITY_PER_SITE){
			if(stepDirection == -1)
				hasOvershot = true;

			stepDirection = 1;
		}
		else{
			if(stepDirection == 1)
				hasOvershot = true;

			stepDirection = -1;
		}
		model.setChemicalPotential(
			model.getChemicalPotential() + stepDirection*stepLength
		);


		//Scale the stepLength depending on whether the overshot has
		//occurred or not.
		if(hasOvershot)
			stepLength /= 2.0;
		else
			stepLength *= 2.0;
	}
}



//Callback function that is to be called each time the model Hamiltonian has
//been diagonalized. The function first fixes the density to the target
//density, then calculates the spin and site resolved density. The new spin and
//site resolved density is mixed with the previous values to stabilize the
//self-consistent calculation.
bool selfConsistencyCallbackPeriodic(Solver::BlockDiagonalizer &solver){
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);

    // std::cout << "Printing initial Eigenvalues in selfConsistencyCallbackPeriodic:" << '\n';
    // for (int k = 0; k < k_size_brillouin_zone; k++) {
    //     std::cout << "Eigenvalues for k = " << k << '\n';
    //     for (int spin = 0; spin < 2; spin++) {
    //
    //         //for(size_t i = 0; i < 2; ++i){
    //             std::cout << propertyExtractor.getEigenValue({k, 0, 0, spin}, 0) << "\n";
    //         //}
    //     }
    // }

	//Fix the density.
	fixDensity(propertyExtractor);

	//Save the old result for later comparison.
	Array<double> oldSpinAndSiteResolvedDensity
		= spinAndSiteResolvedDensity;

	//Calculate the spin and site resolved density. Note that the k-indices
	//are summed over using the IDX_SUM_ALL specifier, while the density is
	//stored separately for each spin and site because of the IDX_ALL
	//specifier.

    Property::Density density = propertyExtractor.calculateDensity({
	    {IDX_SUM_ALL, IDX_SUM_ALL, IDX_SUM_ALL, IDX_ALL, IDX_ALL}
    });

    const int number_of_sites = model.getBasisSize()/(2 * k_size_brillouin_zone);
    PRINTVAR(number_of_sites);

	//Update the spin and site resolved density. Mix with the previous
	//value to stabilize the self-consistent calculation.
    for(int spin = 0; spin < 2; spin++){
		for(int site = 0; site < k_num_atoms_unit_cell; site++){
			spinAndSiteResolvedDensity[{spin, site}]
				= MIXING_PARAMETER*spinAndSiteResolvedDensity[
					{spin, site}
				] + (1 - MIXING_PARAMETER)*density({
                    IDX_SUM_ALL,
                    IDX_SUM_ALL,
                    IDX_SUM_ALL,
					(int)spin,
                    (int)site
				}) / k_size_brillouin_zone; ///(model.getBasisSize()/k_num_atoms_unit_cell);

            PRINTVAR((spinAndSiteResolvedDensity[spin, site]));
		}
	}

    // std::cout << "Printing Eigenvalues in selfConsistencyCallbackPeriodic:" << '\n';
    // for (int k = 15; k < 20; k++) {
    //     std::cout << "Eigenvalues for k = " << k << '\n';
    //     for (int spin = 0; spin < 2; spin++) {
    //
    //         //for(size_t i = 0; i < 2; ++i){
    //             std::cout << propertyExtractor.getEigenValue({k, 0, 0, spin}, 0) << "\n";
    //         //}
    //     }
    // }

	//Calculate the maximum difference between the new and old spin and
	//site resolved density.
    int max_spin {0}, max_site{0};
	double maxDifference = 0;
	for(unsigned int spin = 0; spin < 2; spin++){
		for(int site = 0; site < k_num_atoms_unit_cell; site++){
			double difference = abs(
				spinAndSiteResolvedDensity[{spin, site}]
				- oldSpinAndSiteResolvedDensity[{spin, site}]
			);
			if(difference > maxDifference){
				maxDifference = difference;
                max_spin = spin; max_site = site;
            }

		}
	}

    static int iteration = 0;
    std::cout << "Iteration number: " << ++iteration << '\n';
    std::cout << "Conv. difference:"
                << std::abs(maxDifference - spin_and_site_resolved_density_tol)
                << "\tabs. value: "
                << spinAndSiteResolvedDensity[{max_spin, max_site}]
                << "\tlocation: " << max_spin << " " << max_site << '\n';

	//Return whether the spin and site resolved density has converged. The
	//self-consistent loop will stop once true is returned.
	if(maxDifference > spin_and_site_resolved_density_tol)
		return false;
	else
		return true;
}

std::vector<int> get_atoms_in_unit_cell(Molecule molecule,
        double x_min, double x_max){

    std::vector<int> atoms_in_unit_cell;
    for(size_t i = 0; i < molecule.size(); ++i){
        double mol_x_coord = molecule.get_x_coords(i);
        std::cout << "Testing for value " << mol_x_coord << '\n';
        if(mol_x_coord + k_eps < x_max && mol_x_coord >= x_min - k_eps){
            atoms_in_unit_cell.push_back(i);
            std::cout << "Successful, lies within " << x_min << " and " << x_max << '\n';
        }
    }
    return atoms_in_unit_cell;
}

bool atom_in_unit_cell(int index){
    return index < k_num_atoms_unit_cell;
}


bool atoms_in_unit_cell(int index_1, int index_2){
    return (index_1 < k_num_atoms_unit_cell && index_2 < k_num_atoms_unit_cell);
}


void print_atom_indices(std::vector<int> v){
    for(size_t j = 0; j < v.size(); ++j){
        std::cout << v[j] << " ";
    }
    std::cout << '\n';
}


int main(int argc, char **argv) {
    string periodicity_direction { "" };
    double periodicity_distance { 0.0 };
    complex<double> t { 1.0 };
    double threshold { 1.7 };
    double hubbard { 0.0 };
    double temperature { 0.01 };

    if(argc == 1){
        std::cout << "You need to provide at least a file." << '\n';
        print_help(false); exit(0);
    }
    string file = argv[1];
    std::ifstream in(file);

    for(int i = 0; i < argc; ++i){
        try{
            if(!strcmp(argv[i], "-p") || !strcmp(argv[i], "--periodic")){
                periodicity_direction = argv[++i];
                periodicity_distance = std::stod(argv[++i]);
            }
            if(!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bond_threshold")){
                threshold = std::stod(argv[++i]);
            }
            if(!strcmp(argv[i], "-t") || !strcmp(argv[i], "--hopping_amplitude")){
                t = std::stod(argv[++i]);
            }
            if(!strcmp(argv[i], "-H") || !strcmp(argv[i], "--Hubbard")){
                hubbard = std::stod(argv[++i]);
                U = hubbard;
            }
            if(!strcmp(argv[i], "-T") || !strcmp(argv[i], "--temperature")){
                temperature = std::stod(argv[++i]);
            }
            if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
                print_help(true); exit(0);
            }
        }
        catch(...){
            std::cout << "Invalid amount of parameters. Don't forget to provide a value for the flags you're setting!" << '\n';
            print_help(false); exit(0);
        }
    }

    std::cout << "Parameter values: " << '\n';
    PRINTVAR(periodicity_direction); PRINTVAR(periodicity_distance);
    PRINTVAR(t); PRINTVAR(threshold); PRINTVAR(hubbard); PRINTVAR(temperature);

    //Set the natural units for this calculation.
	UnitHandler::setScales({"1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});

	const int k_points_per_path = k_size_brillouin_zone; // / 5;
	vector<unsigned int> numMeshPoints = {
		k_points_per_path,
        1,
        1
	};

    Molecule whole_molecule(in);

    double x_min = whole_molecule.get_x_min();
    double x_max = x_min + periodicity_distance;

    std::vector<int> atoms_in_unit_cell = get_atoms_in_unit_cell(
            whole_molecule, x_min, x_max);
    k_num_atoms_unit_cell = atoms_in_unit_cell.size();

    std::cout << "Selected atom in original unit cell: " << '\n';
    print_atom_indices(atoms_in_unit_cell);

    double y_min = whole_molecule.get_y_min();
    double y_max = whole_molecule.get_y_max();

    std::cout << "The unit cell spans from (" << x_min << ", " << y_min
              << ") to (" << x_max << ", " << y_max << ")." << '\n';

    Molecule molecule(whole_molecule,
            atoms_in_unit_cell, periodicity_distance);

    std::vector<int> new_atoms_in_unit_cell;
    new_atoms_in_unit_cell = get_atoms_in_unit_cell(molecule, x_min, x_max);
    std::cout << "Selected atom in constructed molecule: " << '\n';
    print_atom_indices(new_atoms_in_unit_cell);

    std::cout << "Minimal molecule: " << '\n';
    std::cout << molecule << std::endl;

    Bonds bonds(molecule, threshold, t);
    bonds.print_bonds();

    bonds_t bonds_in_unit_cell = bonds.bonds_;
    std::cout << "Bonds in molecule: " << '\n';
    for(auto el : bonds_in_unit_cell){
        std::cout << "(" << el.first << "," << el.second << "), ";
    }
    std::cout << '\n';
    std::cout << std::boolalpha;


    // Setup lattice vector. By using three three-dimensional vectors
	// instead of two two-dimensional vectors, the cross product expression
	// for the reciprocal vectors can be expressed in terms of cross
	// products.
    Vector3d r[3];
    r[0] = Vector3d({periodicity_distance, 0, 0});
    r[1] = Vector3d({0,	periodicity_distance, 0});
    r[2] = Vector3d({0, 0, periodicity_distance});

    PRINTVAR(r[0]); PRINTVAR(r[1]); PRINTVAR(r[2]);
    Vector3d k_orig[3];
	for(unsigned int n = 0; n < 3; n++){
		k_orig[n] = 2*M_PI*r[(n+1)%3]*r[(n+2)%3]/(
			Vector3d::dotProduct(r[n], r[(n+1)%3]*r[(n+2)%3])
		);
	}

    PRINTVAR(k_orig[0]); PRINTVAR(k_orig[1]); PRINTVAR(k_orig[2]);

    BrillouinZone brillouinZone(
		{
			{k_orig[0].x,k_orig[0].y,k_orig[0].z},
			{k_orig[1].x,k_orig[1].y,k_orig[1].z},
			{k_orig[2].x,k_orig[2].y,k_orig[2].z},
		},
		SpacePartition::MeshType::Nodal
	);

	vector<vector<double>> mesh = brillouinZone.getMinorMesh(
		numMeshPoints
	);

    Array<double> dummy_array({2, molecule.size()});
    spinAndSiteResolvedDensity = dummy_array;

    initSpinAndSiteResolvedDensity(spinAndSiteResolvedDensity);

    std::cout << "mesh size: " << mesh.size() << '\n' << '\n' << '\n';
    std::cout << "Printing debug messages for m = 3:" << '\n';
    static int printer = 0;
    const int comp_val = 3;
    for(unsigned int m = 0; m < mesh.size(); m++){
        Index kIndex = brillouinZone.getMinorCellIndex(
            mesh[m],
            numMeshPoints
        );

        Vector3d kmesh({mesh[m][0], mesh[m][1], mesh[m][2]});

        const complex<double> one(1, 0);

        h_state_enum h_state = undefined_enum;

        for(int s = 0; s < 2; s++){
            for(auto bond : bonds_in_unit_cell){
                bool first_in_unit_cell = atom_in_unit_cell(bond.first);
                bool second_in_unit_cell = atom_in_unit_cell(bond.second);
                complex<double> h { -t * one };

                if(!first_in_unit_cell){
                    assert(("Error in bonds: Second atom in unit cell but \
                             not the first one." , !second_in_unit_cell));
                    if(printer == comp_val){
                        std::cout << "Skipping hopping for bond (" << bond.first
                              << ", " << bond.second << ")." << '\n';
                    }
                    continue;
                }

                if(!second_in_unit_cell){
                    if(bond_is_in_unit_cell(bonds_in_unit_cell,
                            bond.first-k_num_atoms_unit_cell, bond.second) ||
                            bond_is_in_unit_cell(bonds_in_unit_cell,
                            bond.first, bond.second-k_num_atoms_unit_cell)){
                        h = -t * (one + exp(-i*Vector3d::dotProduct(kmesh, r[0])));
                        h_state = h_both_enum;
                    }
                    else{
                        h = -t * exp(-i*Vector3d::dotProduct(kmesh, r[0]));
                        h_state = h_border_crossing_enum;
                    }
                }
                else{
                    if(bond_is_in_unit_cell(bonds_in_unit_cell,
                            bond.first+k_num_atoms_unit_cell, bond.second) ||
                            bond_is_in_unit_cell(bonds_in_unit_cell,
                            bond.first, bond.second+k_num_atoms_unit_cell)){
                        if(printer == comp_val){
                            std::cout << "Skipping hopping for bond (" << bond.first
                                  << ", " << bond.second << ")." << '\n';
                        }
                        continue;
                    }
                    h_state = h_unit_cell_enum;
                }

                int first_atom { bond.first % k_num_atoms_unit_cell };
                int second_atom { bond.second % k_num_atoms_unit_cell };
                model << HoppingAmplitude(
        			h,
        			{kIndex[0], kIndex[1], kIndex[2], s, first_atom},
                    {kIndex[0], kIndex[1], kIndex[2], s, second_atom}
        		) + HC;
                if(printer == comp_val){
                    std::cout << "Added hoppings from " << first_atom
                            << " to " << second_atom
                            << " as " << h_state << " with h value " << h << '\n';
                }


                // if(h_state == h_border_crossing_enum){
                //     h = -t * exp(-i*Vector3d::dotProduct(kmesh, r[0]));
                //     model << HoppingAmplitude(
            	// 		h,
            	// 		{kIndex[0], kIndex[1], kIndex[2], s, (first_atom % k_num_atoms_unit_cell)},
                //         {kIndex[0], kIndex[1], kIndex[2], s, (second_atom % k_num_atoms_unit_cell))}
            	// 	) + HC;
                //     if(printer == comp_val){
                //         std::cout << "Added hoppings from " << (second_atom - (1-first_is_odd))
                //                 << " to " << (first_atom - first_is_odd)
                //                 << " as " << h_state << " with h value " << h << '\n';
                //     }
                //     continue;
                // }
                // if(first_atom % 2 == 1){
                //     first_atom -= 1;
                //
                //     model << HoppingAmplitude(
            	// 		h,
            	// 		{kIndex[0], kIndex[1], kIndex[2], s, second_atom},
            	// 		{kIndex[0], kIndex[1], kIndex[2], s, first_atom}
            	// 	) + HC;
                //
                //     if(printer == comp_val){
                //         std::cout << "Added hoppings from " << second_atom
                //                 << " to " << first_atom
                //                 << " as " << h_state << " with h value " << h << '\n';
                //     }
                // }
                // else{
                //     bool second_atom_odd { second_atom % 2 == 1 };
                //     model << HoppingAmplitude(
            	// 		h,
            	// 		{kIndex[0], kIndex[1], kIndex[2], s, first_atom},
            	// 		{kIndex[0], kIndex[1], kIndex[2], s, second_atom - second_atom_odd}
            	// 	) + HC;
                //     if(printer == comp_val){
                //         std::cout << "Added hoppings from " << first_atom
                //                 << " to " << second_atom
                //                 << " as " << h_state << " with h value " << h << '\n';
                //     }
                // }

                h_state = undefined_enum;
            }
        }

        if(hubbard){
            for(int atom = 0; atom < new_atoms_in_unit_cell.size(); ++atom){
                for (int s = 0; s < 2; s++){
                    model << HoppingAmplitude(
                        H_U,
                        {kIndex[0], kIndex[1], kIndex[2], s, atom},
                        {kIndex[0], kIndex[1], kIndex[2], s, atom}
                    );
                }
            }
        }


        printer++;
    }
    std::cout << '\n';

    model.construct();


    //Get the HoppingAmplitudeSet from the Model
    //and extract the basis size.
    const HoppingAmplitudeSet &hoppingAmplitudeSet
        = model.getHoppingAmplitudeSet();
    unsigned int basisSize
        = hoppingAmplitudeSet.getBasisSize();
    //Initialize the Hamiltonian on a format most
    //suitable for the algorithm at hand.
    Array<complex<double>> hamiltonian(
        {basisSize, basisSize},
        0.
    );
    //Iterate over the HoppingAmplitudes.
    for(
        HoppingAmplitudeSet::ConstIterator iterator
            = hoppingAmplitudeSet.cbegin();
        iterator != hoppingAmplitudeSet.cend();
        ++iterator
    ){
        //Extract the amplitude and physical
        //indices from the HoppingAmplitude.
        complex<double> amplitude
            = (*iterator).getAmplitude();
        const Index &toIndex
            = (*iterator).getToIndex();
        const Index &fromIndex
            = (*iterator).getFromIndex();

        //Convert the physical indices to linear
        //indices.
        unsigned int row
            = hoppingAmplitudeSet.getBasisIndex(
                toIndex
            );
        unsigned int column
            = hoppingAmplitudeSet.getBasisIndex(
                fromIndex
            );

        //Write the amplitude to the Hamiltonian
        //that will be used in this algorithm.
        hamiltonian[{row, column}] += amplitude;
    }

    // Print the Hamiltonian
    std::cout << "Hamiltonian for k point 0:" << '\n';
    for(unsigned int row = 0; row < atoms_in_unit_cell.size()*2; row++){
        for(unsigned int column = 0; column < atoms_in_unit_cell.size()*2; column++){
            Streams::out << real(hamiltonian[{row, column}])
                << "\t";
        }
        Streams::out << "\n";
    }
    std::cout << '\n';
    std::cout << "Hamiltonian for k point 1:" << '\n';
    unsigned int twice_uc = 2 * atoms_in_unit_cell.size() * 1;
    for(unsigned int row = twice_uc; row < twice_uc + atoms_in_unit_cell.size() * 2; row++){
        for(unsigned int column = twice_uc; column < twice_uc + atoms_in_unit_cell.size() * 2; column++){
            Streams::out << real(hamiltonian[{row, column}])
                << "\t";
        }
        Streams::out << "\n";
    }

	Solver::BlockDiagonalizer solver;
	solver.setModel(model);

    if(hubbard){
        solver.setSelfConsistencyCallback(selfConsistencyCallbackPeriodic);
	    solver.setMaxIterations(1000);
    }

	solver.run();

    if(hubbard){
        std::cout << "Final chemical potential: "
                  << model.getChemicalPotential() << '\n';
    }

	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);

    // std::cout << "Printing Eigenvalues in Main Function 1:" << '\n';
    // for (int k = 15; k < 20; k++) {
    //     std::cout << "Eigenvalues for k = " << k << '\n';
    //     for (int spin = 0; spin < 2; spin++) {
    //
    //         //for(size_t i = 0; i < 2; ++i){
    //             std::cout << propertyExtractor.getEigenValue({k, 0, 0, spin}, 0) << "\n";
    //         //}
    //     }
    // }

	Property::DOS dos = propertyExtractor.calculateDOS();

    Vector3d lowest({-M_PI / periodicity_distance, 0, 0});
    Vector3d highest({M_PI / periodicity_distance, 0, 0});


    double x_diff_test = whole_molecule.get_x_coords(0)
            - whole_molecule.get_x_coords(1);
    double y_diff_test = whole_molecule.get_y_coords(0)
            - whole_molecule.get_y_coords(1);
    Vector3d dist_vec_test({x_diff_test, y_diff_test, 0});
    std::cout << "Distance between atoms: " << dist_vec_test << '\n';
    std::cout << "K Space: " << lowest.x << " to " << highest.x << '\n';
    std::cout << "Periodicity of molecule: " << periodicity_distance << '\n';

    vector<vector<Vector3d>> paths = {
        {lowest, highest}
    };


	Range interpolator(0, 1, k_points_per_path);

    ofstream myfile;
	myfile.open("results.txt");
    PRINTVAR(new_atoms_in_unit_cell.size());
	std::cout << "Writing results to file" << '\n';
    //for (size_t i = 0; i < atoms_in_unit_cell.size(); i++) {
    for (size_t i = 0; i < new_atoms_in_unit_cell.size(); i++) {
        myfile << "value" << i;
        if(i != new_atoms_in_unit_cell.size()-1){
            myfile << ", ";
        }
        else{
            myfile << std::endl;
        }
    }

	for(unsigned int p = 0; p < 1; p++){
		Vector3d startPoint = paths[p][0];
		Vector3d endPoint = paths[p][1];

        std::cout << "startPoint: " << startPoint << '\n';
        std::cout << "endPoint: " << endPoint << '\n';
		for(int n = 0; n < k_points_per_path; n++){
			Vector3d k = (
				interpolator[n]*endPoint
				+ (1 - interpolator[n])*startPoint
			);

			Index kIndex = brillouinZone.getMinorCellIndex(
				{k.x, k.y, k.z},
                numMeshPoints
			);

            //std::cout << "Printing eigenvalues: " << '\n';
            // std::cout << "Eigenvalues for k = " << kIndex[0] << '\n';
            // PRINTVAR(kIndex[0]);
            // PRINTVAR(kIndex[1]);
            // PRINTVAR(kIndex[2]);
            // PRINTVAR(n);
            // PRINTVAR(k);
            // for (int spin = 0; spin < 2; spin++) {
            //     std::cout << propertyExtractor.getEigenValue({kIndex[0], 0, 0, spin}, 0) << "\n";
            // }


            // std::cout << propertyExtractor.getEigenValue({5,0,0,1}, 1) << '\n';
            // std::cout << propertyExtractor.getEigenValue({5,0,0,0}, 1) << '\n';
            // std::cout << "Printing Eigenvalues in Main Function 2:" << '\n';
            // for (int k = 15; k < 20; k++) {
            //     std::cout << "Eigenvalues for k = " << k << '\n';
            //     for (int spin = 0; spin < 2; spin++) {
            //
            //         //for(size_t i = 0; i < 2; ++i){
            //             std::cout << propertyExtractor.getEigenValue({k, 0, 0, spin}, 0) << "\n";
            //         //}
            //     }
            // }

            // std::cout << "New atoms in unit cell:" << '\n';
            // for(size_t i = 0; i < new_atoms_in_unit_cell.size(); ++i){
            //     std::cout << new_atoms_in_unit_cell[i] << " ";
            // }


            for(size_t i = 0; i < new_atoms_in_unit_cell.size(); ++i){
                myfile << propertyExtractor.getEigenValue({kIndex[0], kIndex[1], kIndex[2], 0}, new_atoms_in_unit_cell[i]); //Old: new_atoms_in_unit_cell[i]
                if(i != new_atoms_in_unit_cell.size()-1)
                    myfile << ", ";
                else{
                    myfile << std::endl;
                }
            }
		}
	}

	std::cout << "Writing done." << '\n';
	myfile.close();

    return 0;
}

void print_help(bool full){
    if(full){
        printf("/**********************************************************************/\n");
        printf("// Tight-Binding and Mean-Field Hubbard approximation                 //\n");
        printf("// Bachelor's thesis, Spring Semester 2019, ETH Zurich                //\n");
        printf("// Authors: Robin Worreby                                             //\n");
        printf("// License: Use if you like, but give me credit.                      //\n");
        printf("/**********************************************************************/\n");
        printf("\n");
    }
    printf("Usage: \n./Application molecule.xyz [parameters]\nParameters:\n");
    printf("  -p or --periodic \t\t- sets the periodicity direction and distance, two parameters needed [X, Y, Z] and [distance]\n");
    printf("  -t or --hopping_amplitude \t- sets the Hamiltonian parameter t value (Hopping amplitude), default is 1.0\n");
    printf("  -b or --bond_threshold \t- sets the bond threshold in Ångström\n");
    printf("  -H or --Hubbard \t\t- sets the Hamiltonian parameter U, default is 0.0 \n");
    printf("  -T or --temperature \t\t- sets the temperature for the system in K, default is 0.01 K\n");
    printf("  -h or --help \t\t\t- prints this help info.\n");
}
