// Translate this file with
//
// g++ -O3 assignment-code.cpp -o assignment-code
//
// Run it with
//
// ./demo-code
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2020 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>



#include <cmath>


double t			= 0;
double tFinal			= 0;
double tPlot			= 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;
int InitNumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double	 timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double	 maxV = 0.0;

/**
 * Minimum distance between two elements.
 */
double	 minDx = std::numeric_limits<double>::max();


/**
 * Set up scenario from the command line.
 *
 * If you need additional helper data structures, you can
 * initialise them here. Alternatively, you can introduce a
 * totally new function to initialise additional data fields and
 * call this new function from main after setUp(). Either way is
 * fine.
 *
 * This operation's semantics is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
	NumberOfBodies = (argc-4) / 7;
	InitNumberOfBodies = NumberOfBodies;
	x		 = new double*[NumberOfBodies];
	v		 = new double*[NumberOfBodies];
	mass = new double [NumberOfBodies];

	int readArgument = 1;

	tPlotDelta	 = std::stof(argv[readArgument]); readArgument++;
	tFinal		 = std::stof(argv[readArgument]); readArgument++;
	timeStepSize = std::stof(argv[readArgument]); readArgument++;

	for (int i=0; i<NumberOfBodies; i++) {
		x[i] = new double[3];
		v[i] = new double[3];

		x[i][0] = std::stof(argv[readArgument]); readArgument++;
		x[i][1] = std::stof(argv[readArgument]); readArgument++;
		x[i][2] = std::stof(argv[readArgument]); readArgument++;

		v[i][0] = std::stof(argv[readArgument]); readArgument++;
		v[i][1] = std::stof(argv[readArgument]); readArgument++;
		v[i][2] = std::stof(argv[readArgument]); readArgument++;

		mass[i] = std::stof(argv[readArgument]); readArgument++;
		if (mass[i]<=0.0) {
			std::cerr << "invalid mass for body " << i << std::endl;
			exit(-2);
		}
	}
	std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
	
	if (tPlotDelta<=0.0) {
		std::cout << "plotting switched off" << std::endl;
		tPlot = tFinal + 1.0;
	}
	else {
		std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
		tPlot = 0.0;
	}
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
	videoFile.open( "result.pvd" );
	videoFile << "<?xml version=\"1.0\"?>" << std::endl
						<< "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
						<< "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
	videoFile << "</Collection>"
						<< "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
	static int counter = -1;
	counter++;
	std::stringstream filename;
	filename << "result-" << counter <<  ".vtp";
	std::ofstream out( filename.str().c_str() );
	out << "<VTKFile type=\"PolyData\" >" << std::endl
			<< "<PolyData>" << std::endl
			<< " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
			<< "	<Points>" << std::endl
			<< "	 <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//			<< "	 <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

	for (int i=0; i<NumberOfBodies; i++) {
		out << x[i][0]
				<< " "
				<< x[i][1]
				<< " "
				<< x[i][2]
				<< " ";
	}

	out << "	 </DataArray>" << std::endl
			<< "	</Points>" << std::endl
			<< " </Piece>" << std::endl
			<< "</PolyData>" << std::endl
			<< "</VTKFile>"  << std::endl;

	videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}


void update_collisions()
{
	for (int pri_part=0; pri_part<NumberOfBodies; pri_part++) {
		for (int sec_part=pri_part+1; sec_part<NumberOfBodies; sec_part++) {
			const double distance = sqrt(
				(x[pri_part][0]-x[sec_part][0]) * (x[pri_part][0]-x[sec_part][0]) +
				(x[pri_part][1]-x[sec_part][1]) * (x[pri_part][1]-x[sec_part][1]) +
				(x[pri_part][2]-x[sec_part][2]) * (x[pri_part][2]-x[sec_part][2])
			);
			if (distance <= ((pow(10, -2)/InitNumberOfBodies) * (mass[pri_part] + mass[sec_part])) + 2e-10 && pri_part != sec_part)
			{
				int merge_par = pri_part > sec_part ? sec_part : pri_part;
				int rem = pri_part > sec_part ? pri_part : sec_part;

				v[merge_par][0] = ((mass[pri_part] * v[pri_part][0]) + (mass[sec_part] * v[sec_part][0]))/(mass[pri_part] + mass[sec_part]);
				v[merge_par][1] = ((mass[pri_part] * v[pri_part][1]) + (mass[sec_part] * v[sec_part][1]))/(mass[pri_part] + mass[sec_part]);
				v[merge_par][2] = ((mass[pri_part] * v[pri_part][2]) + (mass[sec_part] * v[sec_part][2]))/(mass[pri_part] + mass[sec_part]);

				x[merge_par][0] = ((mass[pri_part] * x[pri_part][0]) + (mass[sec_part] * x[sec_part][0]))/(mass[pri_part] + mass[sec_part]);
				x[merge_par][1] = ((mass[pri_part] * x[pri_part][1]) + (mass[sec_part] * x[sec_part][1]))/(mass[pri_part] + mass[sec_part]);
				x[merge_par][2] = ((mass[pri_part] * x[pri_part][2]) + (mass[sec_part] * x[sec_part][2]))/(mass[pri_part] + mass[sec_part]);

				mass[merge_par] = mass[pri_part] + mass[sec_part];

				x[rem] = x[NumberOfBodies-1];
				v[rem] = v[NumberOfBodies-1];
				mass[rem] = mass[NumberOfBodies-1];
				pri_part = merge_par;
				sec_part = -1;
				NumberOfBodies--;
				minDx = std::min(minDx, distance);
			}
		}
	}
}
			
//Aim here is to be as accurate as possible	
void update_acceleration(double f[][3])
{
	double comp_force[3];
	double distance = 0;
	for (int pri_part=0; pri_part<NumberOfBodies; pri_part++) {
		#pragma omp simd reduction(+: f[pri_part:NumberOfBodies][:3]) reduction(min: minDx) private(distance, comp_force)
		for (int sec_part=pri_part+1; sec_part<NumberOfBodies; sec_part++) {
			distance = sqrt(
				(x[pri_part][0]-x[sec_part][0]) * (x[pri_part][0]-x[sec_part][0]) +
				(x[pri_part][1]-x[sec_part][1]) * (x[pri_part][1]-x[sec_part][1]) +
				(x[pri_part][2]-x[sec_part][2]) * (x[pri_part][2]-x[sec_part][2])
			);
			comp_force[0] = (x[sec_part][0]-x[pri_part][0]) * mass[pri_part] * mass[sec_part] / distance / distance / distance;
                        comp_force[1] = (x[sec_part][1]-x[pri_part][1]) * mass[pri_part] * mass[sec_part] / distance / distance / distance;
                        comp_force[2] = (x[sec_part][2]-x[pri_part][2]) * mass[pri_part] * mass[sec_part] / distance / distance / distance;

			f[pri_part][0] += comp_force[0];
			f[pri_part][1] += comp_force[1];
			f[pri_part][2] += comp_force[2];

			f[sec_part][0] += -comp_force[0];
			f[sec_part][1] += -comp_force[1];
			f[sec_part][2] += -comp_force[2];
			
			minDx = std::min(minDx, distance);
		}
	}
}

/**
 * This is the main operation you should change in the assignment. You might
 * want to add a few more variables or helper functions, but this is where the
 * magic happens.
*/
void updateBody() {
	/* Given state t-1 are then updating to reflect change in state over interval. Change will be different if collision
	   occured at end of previous state so must update for these first before updating for acceleration. Could update for 
	   collisions at end of function, however, in state t=0 particles may already be close enough to collide, therefore, must 
	   update collision and then motion in interval*/
	double f[NumberOfBodies][3];
	for(int i = 0; i < NumberOfBodies; ++i)
	{
		for (int t = 0; t < 3; t++)
		{
			f[i][t] = 0;
		}
	}
	update_collisions();
	update_acceleration(f);
	
	//Moves first particle accoriding to forces from other paticles just calculated
	//x is updated using old v -> can do this without force but won't have same opertion
	//If duplicate v can then update x whenever -> update new v as you go 
	for (int part=0; part<NumberOfBodies; part++) {
		x[part][0] = x[part][0] + timeStepSize * v[part][0];
		x[part][1] = x[part][1] + timeStepSize * v[part][1];
		x[part][2] = x[part][2] + timeStepSize * v[part][2];

		//Updates velocity of first particle according to forces calculated
		v[part][0] = v[part][0] + (timeStepSize * f[part][0] / mass[part]);
		v[part][1] = v[part][1] + (timeStepSize * f[part][1] / mass[part]);
		v[part][2] = v[part][2] + (timeStepSize * f[part][2] / mass[part]);

		double v_s = v[part][0]*v[part][0] + v[part][1]*v[part][1] + v[part][2]*v[part][2];

		const double speed = std::sqrt(v_s);
		maxV = std::max(maxV, speed);
	}
	t += timeStepSize;
}


/**
 * Main routine.
 *
 * No major changes in assignment. You can add a few initialisation
 * or stuff if you feel the need to do so. But keep in mind that you
 * may not alter what the program plots to the terminal.
 */
int main(int argc, char** argv) {
	if (argc==1) {
		std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
							<< "	snapshot				interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
							<< "	final-time			simulated time (greater 0)" << std::endl
							<< "	dt							time step size (greater 0)" << std::endl
							<< std::endl
							<< "Examples:" << std::endl
							<< "0.01	100.0  0.001		0.0 0.0 0.0  1.0 0.0 0.0	1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
							<< "0.01	100.0  0.001		0.0 0.0 0.0  1.0 0.0 0.0	1.0			0.0 1.0 0.0  1.0 0.0 0.0	1.0  \t One spiralling around the other one" << std::endl
							<< "0.01	100.0  0.001		3.0 0.0 0.0  0.0 1.0 0.0	0.4			0.0 0.0 0.0  0.0 0.0 0.0	0.2			2.0 0.0 0.0  0.0 0.0 0.0	1.0 \t Three body setup from first lecture" << std::endl
							<< "0.01	100.0  0.001		3.0 0.0 0.0  0.0 1.0 0.0	0.4			0.0 0.0 0.0  0.0 0.0 0.0	0.2			2.0 0.0 0.0  0.0 0.0 0.0	1.0			2.0 1.0 0.0  0.0 0.0 0.0	1.0			2.0 0.0 1.0  0.0 0.0 0.0	1.0 \t Five body setup" << std::endl
							<< std::endl
							<< "In this naive code, only the first body moves" << std::endl;

		return -1;
	}
	else if ( (argc-4)%7!=0 ) {
		std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
		std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
		std::cerr << "run without arguments for usage instruction" << std::endl;
		return -2;
	}

	std::cout << std::setprecision(15);

	setUp(argc,argv);

	openParaviewVideoFile();

	int snapshotCounter = 0;
	if (t > tPlot) {
		printParaviewSnapshot();
		std::cout << "plotted initial setup" << std::endl;
		tPlot = tPlotDelta;
	}

	int timeStepCounter = 0;
	while (t<=tFinal) {
		updateBody();
		timeStepCounter++;
		if (t >= tPlot) {
			printParaviewSnapshot();
			std::cout << "plot next snapshot"
						<< ",\t time step=" << timeStepCounter
						<< ",\t t="					<< t
				<< ",\t dt="				<< timeStepSize
				<< ",\t v_max="			<< maxV
				<< ",\t dx_min="		<< minDx
				<< std::endl;

			tPlot += tPlotDelta;
		}
	}

	std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
	std::cout << "Position of first remaining object: " << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;

	closeParaviewVideoFile();

	return 0;
}
