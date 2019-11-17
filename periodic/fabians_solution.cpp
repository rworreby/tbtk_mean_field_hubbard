/* Copyright 2019 Kristofer Björnson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
​
#include <TBTK/BrillouinZone.h>
#include <TBTK/Model.h>
#include <TBTK/Property/DOS.h>
#include <TBTK/PropertyExtractor/BlockDiagonalizer.h>
#include <TBTK/Range.h>
#include <TBTK/Solver/BlockDiagonalizer.h>
#include <TBTK/Streams.h>
#include <TBTK/UnitHandler.h>
#include <TBTK/Vector3d.h>
​
using namespace std;
using namespace TBTK;
//using namespace Plot;
​
complex<double> i(0, 1);
​
int main(int argc, char **argv){
	//Set the natural units for this calculation.
	UnitHandler::setScales({"1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});
​
	//Define parameters.
	double t = 3;	//eV
	double a = 2.5;	//Ångström
	unsigned int BRILLOUIN_ZONE_RESOLUTION = 1000;
	vector<unsigned int> numMeshPoints = {
		BRILLOUIN_ZONE_RESOLUTION,
    // NOTE: Resolution is just one in the not-existing directions!
		1,
    1
	};
	const int K_POINTS_PER_PATH = 100;
	const double ENERGY_LOWER_BOUND = -10;
	const double ENERGY_UPPER_BOUND = 10;
	const int ENERGY_RESOLUTION = 1000;
​
	//Setup lattice vector. By using three three-dimensional vectors
	//instead of two two-dimensional vectors, the cross product expression
	//for the reciprocal vectors can be expressed in terms of cross
	//products.
	Vector3d r[3];
	r[0] = Vector3d({a,	0,		0});
  // NOTE: dummy variable to build r_AB, set to (0,a,0) later
	r[1] = Vector3d({-a/2,	a*sqrt(3)/2,	0});
	r[2] = Vector3d({0,	0,		a});
​
	//Define nearest neighbor vectors for A site.
	Vector3d r_AB[3];
	r_AB[0] = (r[0] + 2*r[1])/3.;
	r_AB[1] = -r[1] + r_AB[0];
	r_AB[2] = -r[0] - r[1] + r_AB[0];
	r[1] = Vector3d({0,	a,	0});
​
	//Calculate the reciprocal lattice vectors.
	Vector3d k[3];
	for(unsigned int n = 0; n < 3; n++){
		k[n] = 2*M_PI*r[(n+1)%3]*r[(n+2)%3]/(
			Vector3d::dotProduct(r[n], r[(n+1)%3]*r[(n+2)%3])
		);
	}
​
  std::cout << k[0];
  std::cout << k[1];
  std::cout << k[2];
​
	//Setup the BrillouinZone.
  // NOTE: just a regular old cube!
	BrillouinZone brillouinZone(
		{
			{k[0].x,k[0].y,k[0].z},
			{k[1].x,k[1].y,k[1].z},
			{k[2].x,k[2].y,k[2].z},
		},
		SpacePartition::MeshType::Nodal
	);
​
	//Create mesh.
	vector<vector<double>> mesh = brillouinZone.getMinorMesh(
		numMeshPoints
	);
​
	//Setup model.
	Model model;
	for(unsigned int n = 0; n < mesh.size(); n++){
		//Get the Index representation of the current k-point.
		Index kIndex = brillouinZone.getMinorCellIndex(
			mesh[n],
			numMeshPoints
		);
​
		//Calculate the matrix element.
    // NOTE: k = (val, 0, 0) -> x-projection only!
		Vector3d k({mesh[n][0], mesh[n][1], mesh[n][2]});
​
    /* Unit cell with hoppings:
       o             -> 7
        \     /
           o         -> 6
           |
           o         -> 4
        /     \
       o             -> 5
       |
       o             -> 3
        \     /
           o         -> 2
           |
           o         -> 0
        /     \
       o             -> 1
    */

		complex<double> h_01 = -t*(
			 exp(-i*Vector3d::dotProduct(k, r_AB[1]))
			+ exp(-i*Vector3d::dotProduct(k, r_AB[2]))
		);
		complex<double> h_02 = -t*(
			exp(-i*Vector3d::dotProduct(k, r_AB[0]))
		);
		complex<double> h_23 = -t*(
			 exp(-i*Vector3d::dotProduct(k, r_AB[1]))
			+ exp(-i*Vector3d::dotProduct(k, r_AB[2]))
		);
		complex<double> h_35 = -t*(
			exp(-i*Vector3d::dotProduct(k, r_AB[0]))
		);
		complex<double> h_45 = -t*(
			 exp(-i*Vector3d::dotProduct(k, r_AB[1]))
			+ exp(-i*Vector3d::dotProduct(k, r_AB[2]))
		);
		complex<double> h_46 = -t*(
			exp(-i*Vector3d::dotProduct(k, r_AB[0]))
		);
		complex<double> h_67 = -t*(
			 exp(-i*Vector3d::dotProduct(k, r_AB[1]))
			+ exp(-i*Vector3d::dotProduct(k, r_AB[2]))
		);
​
		//Add the matrix element to the model.
    // NOTE: kIndex = [indexVal, 0, 0]
		model << HoppingAmplitude(
			h_67,
			{kIndex[0], kIndex[1], kIndex[2], 6},
			{kIndex[0], kIndex[1], kIndex[2], 7}
		) + HC;
		model << HoppingAmplitude(
			h_46,
			{kIndex[0], kIndex[1], kIndex[2], 4},
			{kIndex[0], kIndex[1], kIndex[2], 6}
		) + HC;
		model << HoppingAmplitude(
			h_45,
			{kIndex[0], kIndex[1], kIndex[2], 4},
			{kIndex[0], kIndex[1], kIndex[2], 5}
		) + HC;
		model << HoppingAmplitude(
			h_35,
			{kIndex[0], kIndex[1], kIndex[2], 3},
			{kIndex[0], kIndex[1], kIndex[2], 5}
		) + HC;
		model << HoppingAmplitude(
			h_23,
			{kIndex[0], kIndex[1], kIndex[2], 2},
			{kIndex[0], kIndex[1], kIndex[2], 3}
		) + HC;
		model << HoppingAmplitude(
			h_02,
			{kIndex[0], kIndex[1], kIndex[2], 0},
			{kIndex[0], kIndex[1], kIndex[2], 2}
		) + HC;
		model << HoppingAmplitude(
			h_01,
			{kIndex[0], kIndex[1], kIndex[2], 0},
			{kIndex[0], kIndex[1], kIndex[2], 1}
		) + HC;
	}
	model.construct();
​
	//Setup the solver.
	Solver::BlockDiagonalizer solver;
	solver.setModel(model);
	solver.run();
​
	//Setup the property extractor.
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
	propertyExtractor.setEnergyWindow(
		ENERGY_LOWER_BOUND,
		ENERGY_UPPER_BOUND,
		ENERGY_RESOLUTION
	);
​
	//Calculate the density of states.
  /*
	Property::DOS dos = propertyExtractor.calculateDOS();
	Plotter plotter;
	plotter.setLabelX("Energy");
	plotter.setLabelY("DOS");
	plotter.plot(dos, 0.03);
	plotter.save("figures/DOS.png");
  */
​
	//Define high symmetry points.
	Vector3d Gamma({-M_PI/a,		0,			0});
	Vector3d M({M_PI/a, 0,	0});
​
	//Define paths between high symmetry points.
	vector<vector<Vector3d>> paths = {
		{Gamma, M}
	};
​
	//Calculate the band structure along the path Gamma -> M
//	Array<double> bandStructure({2, K_POINTS_PER_PATH}, 0);
	Range interpolator(0, 1, K_POINTS_PER_PATH);
	for(unsigned int p = 0; p < 1; p++){
		//Select the start and end points for the current path.
		Vector3d startPoint = paths[p][0];
		Vector3d endPoint = paths[p][1];
​
		//Loop over a single path.
		for(unsigned int n = 0; n < K_POINTS_PER_PATH; n++){
			//Interpolate between the paths start and end point.
			Vector3d k = (
				interpolator[n]*endPoint
				+ (1 - interpolator[n])*startPoint
			);
​
			//Get the Index representation of the current k-point.
			Index kIndex = brillouinZone.getMinorCellIndex(
				{k.x, k.y, k.z},
				numMeshPoints
			);
​
			//Extract the eigenvalues for the current k-point.
      // NOTE: for some reason you can increase the eigenvalue number
      //       but you'll just get repeated values...
			std::cerr << propertyExtractor.getEigenValue(kIndex, 0) << " ";
			std::cerr << propertyExtractor.getEigenValue(kIndex, 1) << " ";
			std::cerr << propertyExtractor.getEigenValue(kIndex, 2) << " ";
			std::cerr << propertyExtractor.getEigenValue(kIndex, 3) << " ";
			std::cerr << propertyExtractor.getEigenValue(kIndex, 4) << " ";
			std::cerr << propertyExtractor.getEigenValue(kIndex, 5) << " ";
			std::cerr << propertyExtractor.getEigenValue(kIndex, 6) << " ";
			std::cerr << propertyExtractor.getEigenValue(kIndex, 7) << std::endl;
		}
	}
​
/*
	//Find max and min value for the band structure.
	double min = bandStructure[{0, 0}];
	double max = bandStructure[{1, 0}];
	for(unsigned int n = 0; n < 1*K_POINTS_PER_PATH; n++){
		if(min > bandStructure[{0, n}])
			min = bandStructure[{0, n}];
		if(max < bandStructure[{1, n}])
			max = bandStructure[{1, n}];
	}
*/
/*
	//Plot the band structure.
	plotter.clear();
	plotter.setHold(true);
	plotter.setLabelX("k");
	plotter.setLabelY("Energy");
	plotter.plot(bandStructure.getSlice({0, _a_}));
	plotter.plot(bandStructure.getSlice({1, _a_}));
	plotter.plot({K_POINTS_PER_PATH, K_POINTS_PER_PATH}, {min, max});
	plotter.plot({2*K_POINTS_PER_PATH, 2*K_POINTS_PER_PATH}, {min, max});
	plotter.save("figures/BandStructure.png");
*/
​
	return 0;
}
