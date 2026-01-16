/**
 * \file JavaNativeAccessWrappers.cpp
 * \author Tom Eichlersmith <eichl008@umn.edu>
 *
 * Wrappers around construction, getters, setters, and destruction
 * functions of Gbl classes.
 *
 * These wrappers are put into an `extern "C"` block so that there
 * names are _not_ mangled and therefore can be accessed by name
 * within a Java Native Access (JNA) class accessing the GBL library.
 *
 * \note This file is large however it is not intended to do anything besides
 * interface between JNA and GBL. If any logic in this file affects the
 * calculations done by GBL, that is considered a bug with the implementation.
 *
 * Information on using these wrappers to access the GBL C++ library
 * from with java using JNA can be found in on the \ref jnausage
 * page.
 */

/**
 * \page jnausage JNA Usage
 * \brief how to use wrapper functions within java
 *
 * More information about JNA can be found at its
 * [GitHub repository](https://github.com/java-native-access/jna).
 *
 * In general, the most difficult part of implementing the JNA-GBL
 * interaction is getting the transfer of data between the java objects
 * and the C++/GBL objects done well and safely. The JNA documentation
 * apart of its repository contains more information about how to do
 * this well - here I will simply walk through an example accessing a
 * single class from GBL within java. In reality, a working solution
 * would require a similar setup for almost all of the different GBL
 * classes one would be interested in using.
 *
 * **It is suggested to enable the JNA_DEBUG build option when
 * first developing a JNA-based java usage of GBL. This will
 * help you make sure that you aren't leaking memory and avoid
 * a program crash.**
 *
 * ### Wrapping MilleBinary
 * It is easiest to simply look at the java source code for this
 * class so one can see how to interface with it. The basic idea
 * is to define a singleton that represents the GBL library as
 * loaded by JNA and then have classes whose job is to wrap
 * these functions in simpler, object-oriented forms.
 * ```java
 * import com.sun.jna.Library;
 * import com.sun.jna.Native;
 * import com.sun.jna.Pointer;
 *
 * // extension of JNA Library that represents the GBL dynamic library
 * public interface GBLInterface extends Library {
 *     GBLInterface INSTANCE = (GBLInterface) Native.loadLibrary("GBL", GBLInterface.class);
 *     Pointer MilleBinaryCtor(String filename, int filenamesize, int doublePrec, int keepZeros, int aSize);
 *     void MilleBinary_close(Pointer self);
 *     // would add more functions whose signatures "match" the ones defined in this file
 *     // (more on what "match" means below)
 * }
 *
 * // java class that represents a MilleBinary file
 * public class MilleBinary {
 *     // hold onto the object dynamically allocated from C++
 *     private Pointer self;
 *
 *     // provide a constructor that is easier for a user than the raw C-wrapped function
 *     // make parameters easier, for example doing the boolean->integer conversion for the user
 *     public MilleBinary(String fileName, boolean doublePrec, boolean keepZeros, int aSize) {
 *         self = GBLInterface.INSTANCE.MilleBinaryCtor(fileName, fileName.length(),
 *             doubePrec ? 0 : 1, keepZeros ? 0 : 1, aSize);
 *     }
 *
 *     // dynamically allocated memory done by the native library is not cleaned up
 *     // by the JVM so one must manually clean it up
 *     public void close() {
 *         GBLInterface.INSTANCE.MilleBinary_close(self);
 *     }
 * };
 * ```
 *
 * ### Matching Function Signatures
 * For function signatures to "match" between the functions defined in the JNA
 * library extensions (`GBLInterface` above) and the ones defined here, the
 * return value type, the name, and the argument types need to match. The name is easy,
 * but the types are slightly more complicated since the typename between java
 * and C++ are different.
 *
 * - The simple types (`int` and `double`) are the same.
 * - A java `String` is converted to a C-style string `char *`.
 * - A pointer to any structure in C is represented by Pointer on the JNA side.
 * - If the C side needs a variable passed by reference, one needs to use
 *   `IntByReference` (or `DoubleByReference`) on the java side and a pointer
 *   on the C side.
 * - If an array is of a known length, one can allocate the array in `java`
 *   and then simply pass the pointer to the array to the C side.
 *   - e.g. A length-3 array is common to represent position. Both "sides"
 *      could use the `double position[3]` syntax and one just has to make
 *      sure to allocate the correct size on the java side and the C function
 *      will simply write to those addresses.
 * - If the array length _must_ be determined by the C side, then one
 *   must use `PointerByReference`.
 *
 * ### Memory Handling
 * When using JNA to call functions from a native library, the memory
 * allocated by those functions _is not_ monitored and cleaned up by
 * java's garbage collector. Effectively, this means you _need_ to
 * have a `delete` call for every `Ctor` call you make.
 *
 * Since this can get complicated very quickly, it is recommended
 * to develop your java program with JNA_MEMORY_MONITOR enabled
 * in the GBL C++ library. This will print out a summary of the 
 * GBL structures still allocated at the end of running allowing
 * you to make sure that the GBL structures you created are also
 * deleted.
 * ```
 * cmake -B build -S . -DJNA_MEMORY_MONITOR=ON
 * ```
 */

#include "Mille/MilleFactory.h"
#include "GblPoint.h"
#include "GblTrajectory.h"
#include "GblUtilities.h" 

#include "Eigen/Core"

const int NROW = 5;
const int NCOL = 5;

using namespace gbl;
using namespace Eigen;

/**
 * Do the memory monitoring if either JNA_DEBUG or JNA_MEMORY_MONITOR
 * is defined.
 */
#if (defined (JNA_DEBUG) || defined (JNA_MEMORY_MONITOR))
#define JNA_DO_MONITOR 1
#else
#define JNA_DO_MONITOR 0
#endif

#if JNA_DO_MONITOR
#include <iostream>

long int num_gbl_point = 0;
long int num_gbl_traj  = 0;
long int num_mille_bin = 0;
long int num_gbl_det_layer = 0;
long int num_gbl_simple_helix = 0;
long int num_gbl_helix_prediction = 0;

/**
 * Print the status of the running counts
 *
 * This function is only compiled when either JNA_DEBUG or JNA_MEMORY_MONITOR
 * are enabled. We include a GCC attribute so that the function is run when
 * the GBL library is offloaded during execution so we see the status at the
 * end of program running.
 *
 * https://gcc.gnu.org/onlinedocs/gcc/Common-Function-Attributes.html
 *
 * There are similar attributes for other compilers, but I leave that
 * for future contributors.
 */
__attribute__((destructor))
void print_status() {
	std::cout
		<< "  GBL Structures Left\n"
		<< "GBL Points:       " << num_gbl_point << "\n"
		<< "GBL Trajectories: " << num_gbl_traj << "\n"
		<< "Mille Binaries:   " << num_mille_bin << "\n"
		<< "GBL Det Layers:   " << num_gbl_det_layer << "\n"
		<< "GBL Simple Helix: " << num_gbl_simple_helix << "\n"
		<< "GBL Helix Pred:   " << num_gbl_helix_prediction << "\n"
		<< std::flush;
}
#endif

/**
 * \brief convert the pointer array of gbl::GblPoint into a vector holding the objects
 *
 * This is a helper function and should not be bound to a function by JNA,
 * so we are keeping it _outside_ the `extern "C"` block so its name is mangled.
 *
 * \note We *copy* the data pointed to into the vector, so the points
 * input into this function **still need to be deleted**.
 *
 * \param [in] points array of pointers to gblGblPoint to put into vector
 * \param [in] npoints number of points (size of array)
 * \return vector of GblPoints with same content as array
 */
std::vector<GblPoint> ptr_array_to_vector(GblPoint* points[], int npoints) {
	std::vector<GblPoint> points_vec;
	// since we already know the size of the vector, we 'reserve' the size
	// so that the vector doesn't need to waste time copying/moving the GblPoints
	// around as it grows in size
	// we do *not* use 'resize' since that would involve default-constructing
	// all the GblPoints
	points_vec.reserve(npoints);

	for (int i{0}; i < npoints; ++i) {
		// get the pointer
		GblPoint* gblpoint = points[i];
		// COPY the data into the vector,
		//  this copy-constructs a /new/ gbl point
		points_vec.emplace_back(*(gblpoint));
#if JNA_DO_MONITOR
		++num_gbl_point;
#endif
#ifdef JNA_DEBUG
		std::cout << "COPY GblPoint " << gblpoint << " -> " << &(points_vec.back()) << std::endl;
#endif
	}

	return points_vec;
}

extern "C" { 
/**
 * \brief Dynamically create new MilleBinary file.
 *
 * We return the raw pointer to this object so that JNA can effectively access it.
 *
 * Unfortunately the translation of booleans to/from java is not very stable, so it
 * is safer to simply pass integers and check if they are 0 (for false) or nonzero (for true).
 *
 * \param [in] fileName name of file in a C-style string
 * \param [in] filenamesize length of filename for accessing C-style string
 * \param [in] doublePrecision non-zero if should use doublePrecision, zero if should use float precision
 * \param [in] keepZeros non-zero to keep zeros, zero to not keep them
 * \param [in] aSize size of buffer to keep in memory
 */
MilleRecord* MilleBinaryCtor(const char* fileName, int filenamesize, int doublePrecision, int keepZeros, int aSize) {
#ifdef JNA_DEBUG
	std::cout << "MilleBinaryCtor(" << fileName << ", " << filenamesize << ", " 
		<< doublePrecision << ", " << keepZeros << ", " << aSize << ")" << std::endl;
#endif
	std::string binName(fileName,filenamesize);
	MilleRecord* mb = spawnMilleRecord(binName,doublePrecision!=0,keepZeros!=0,aSize).release();
#ifdef JNA_DEBUG
	std::cout << "MilleBinary created at " << mb << std::endl;
#endif
#if JNA_DO_MONITOR
	++num_mille_bin;
#endif
	return mb;
}

/**
 * \brief Closing a gbl::MilleBinary file is the same as destructing it
 *
 * Since the destructor of the gbl::MilleBinary class is what handles
 * performing the final write operations, we simply `delete` the 
 * object pointed to by the passed pointer.
 *
 * @note This means the object on the JNA side will be invalid
 * and will cause a program crash if it is accessed after using this
 * function on it!
 *
 * \param [in] self gbl::MilleBinary to delete
 */
void MilleRecord_close(MilleRecord* self) {
#ifdef JNA_DEBUG
	std::cout << "MilleBinary_close(" << self << ")" << std::endl;
#endif
#if JNA_DO_MONITOR
	--num_mille_bin;
#endif
	if (self) delete self;
}

/**
 * \brief create new gbl::GblPoint from a jacobian
 *
 * \param [in] matrixArray C-style double array listing the jacobian row-by-row.
 * \return dynamically created GblPoint
 */
GblPoint* GblPointCtor(double matrixArray[NROW*NCOL]) {
	Map<Matrix5d> jacobian(matrixArray,5,5);
	GblPoint* self = new GblPoint(jacobian);
#if JNA_DO_MONITOR
	++num_gbl_point;
#endif
#ifdef JNA_DEBUG
	std::cout << "GblPointCtor at " << self << " " << num_gbl_point << std::endl;
#endif
	return self;
}

/**
 * \brief delete gbl::GblPoint
 *
 * provided so that users of JNA can clean up the memory that is not handled
 * automatically by java itself
 *
 * \param [in] self gbl::GblPoint to delete
 */
void GblPoint_delete(GblPoint* self) {
#if JNA_DO_MONITOR
	--num_gbl_point;
#endif
#ifdef JNA_DEBUG
	std::cout << "GblPoint_delete(" << self << ") " << num_gbl_point << std::endl;
#endif
	if (self) delete self;
}


/**
 * \brief call gbl::GblPoint::printPoint on self
 */
void GblPoint_printPoint(const GblPoint* self, unsigned int level) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_printPoint(" << self << ", " << level << ")" << std::endl;
#endif
	self->printPoint(level);
}

/**
 * \brief calculate number of measurements in a gbl::GblPoint
 *
 * Since java struggles to handle the C++ iterators, we calculate
 * the difference between the iterators here so that the java side
 * can access how many measurements a GblPoint has.
 *
 * \param [in] self gbl::GblPoint to operate on
 * \return number of measurements
 */
int GblPoint_getNumMeasurements(GblPoint* self) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_getNumMeasurements(" << self << ")" << std::endl;
#endif
	return (self->getMeasEnd() - self->getMeasBegin());
}

/**
 * \brief add a 2D measurement
 *
 * \note We are only supporting a 2D position residual
 * and a 2x2 projection matrix!
 *
 * \see gbl::GblPoint::addMeasurement
 *
 * \param [in] self gbl::GblPoint to operate on
 * \param [in] projArray length-4 array holding the entries in the 2x2 proj matrix
 * \param [in] resArray length-2 array holding residuals
 * \param [in] precArray length-2 array holding the precisions
 * \param [in] minPrecision Minimal precision to accept measurement
 */
void GblPoint_addMeasurement2D(GblPoint* self, 
								 double *projArray,
								 double *resArray,
								 double *precArray, 
								 double minPrecision) { 
#ifdef JNA_DEBUG
	std::cout << "GblPoint_addMeasurement2D("
		<< self << ", " << projArray << ", " << resArray << ", " << precArray << ", " << minPrecision
		<< ")" << std::endl;
#endif
	
	Map<Matrix2d> aProjection(projArray,2,2);
	Map<Vector2d> aResiduals(resArray, 2);
	Map<Vector2d> aPrecision(precArray,2);
	
	self->addMeasurement(aProjection, aResiduals, aPrecision, minPrecision);
}

/**
 * \brief add a vector-precision scatterer
 * \note only supporting vector-precision at this time
 * \see gbl::GblPoint::addScatterer
 * \param [in] self gbl::GblPoint to operate on
 * \param [in] resArray length-2 residuals array
 * \param [in] precArray length-2 precision array
 */
void GblPoint_addScatterer(GblPoint* self, double *resArray, double *precArray) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_addScatterer(" << self << ", " 
		<< resArray << ", " << precArray << ")" << std::endl;
#endif
	// chose to do the Vector2d addScatterer since
	// PF's original comment "only support vector precision"
	Eigen::Vector2d aResiduals(resArray);
	Eigen::Vector2d aPrecision(precArray);
	
	self->addScatterer(aResiduals,aPrecision);
}

/**
 * \brief add global derivatives to the point
 * \see gbl::GblPoint::addGlobals
 * \param [in] self gbl::GblPoint to operate on
 * \param [in] labels int array of labels nlabels long
 * \param [in] nlabels number of labels
 * \param [in] derArray double array of derivatives nlabels long
 */
void GblPoint_addGlobals(GblPoint* self, int *labels, int nlabels, double* derArray) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_addGlobals(" << self
		<< ", " << labels << ", " << nlabels << ", " << derArray << ")" << std::endl;
#endif
	std::vector<int> aLabels;
	for (int i=0; i<nlabels; i++) {
		aLabels.push_back(labels[i]);
	}
	Map<Eigen::MatrixXd> derivatives(derArray,1,nlabels);
	self->addGlobals(aLabels, derivatives);
}

/**
 * \brief get global derivatives and their labels from the point
 * \see gbl::GblPoint::getGlobalLabelsAndDerivatives
 * \param [in] self GblPoint to operate on
 * \param [out] nlabels pointer to int where number of labels will be stored
 * \param [out] labels pointer to array where labels will be stored
 * \param [out] ders pointer to array where derivatives will be stored
 */
void GblPoint_getGlobalLabelsAndDerivatives(GblPoint* self, int* nlabels, int** labels, double** ders) {
#ifdef JNA_DEBUG
	std::cout << "GblPoint_getGlobalLabelsAndDerivatives("
		<< self << ", " << labels << ", " << ders
		<< ")" << std::endl;
#endif
	std::vector<int> glabels;
	std::vector<double> gders;

#ifdef JNA_DEBUG
	std::cout << "  GblPoint has " << std::flush
		<< self->getMeasEnd() - self->getMeasBegin()
		<< " measurements." << std::endl;
#endif

	//Should I add the number of derivatives? -  Row/Col? CHECK CHECK CHECK
	self->getGlobalLabelsAndDerivatives(
			0 /* aMeas */, 0 /* aRow	*/,
			glabels, gders);

#ifdef JNA_DEBUG
	std::cout <<"  Aquired " << glabels.size() << " labels and " << gders.size() << " derivatives." <<std::endl;
#endif
	
	*nlabels = glabels.size();
	*labels = new int[*nlabels];
	*ders   = new double[*nlabels];

#ifdef JNA_DEBUG
	std::cout <<"  Allocated return arrays." << std::endl;
#endif
	
	
	for (std::size_t il{0}; il < glabels.size(); ++il) {
		(*labels)[il] = glabels.at(il);
		//std::cout<<glabels.at(il)<<std::endl;
	}

	//std::cout<<"GblPointWrapper::gders"<<std::endl;
	
	for (std::size_t id{0}; id < glabels.size(); ++id) {
		(*ders)[id] = gders.at(id);
		//std::cout<<gders.at(il)<<std::endl;
	}

	// set array using Eigen
	//Map<MatrixXd>(ders,1,gders.size()) = gders;
#ifdef JNA_DEBUG
	std::cout << "  Done with assignment and leaving." << std::endl;
#endif
}

/**
 * \brief simple construction of a new gbl::GblTrajectory
 *
 * \see gbl::GblTrajectory::GblTrajectory
 *
 * Just constructing a trajectory from a set of points. No seeding matrix.
 * 
 * \param [in] points array of pointers to gbl::GblPoint to put into trajectory
 * \param [in] npoints number of points in array
 * \param [in] flagCurv use q/p - non-zero for true, zero for false
 * \param [in] flagU1dir use in u1 direction - non-zero for true, zero for false
 * \param [in] flagU2dir use in u2 direction - non-zero for true, zero for false
 * \return dynamically allocated trajectory wrapping input points
 */
GblTrajectory* GblTrajectoryCtorPtrArray(GblPoint* points[], int npoints, 
										 int flagCurv, int flagU1dir, int flagU2dir) {
#ifdef JNA_DEBUG
	std::cout << "GblTracjectoryCtorPtrArray("
		<< points << ", " << npoints << ", "
		<< flagCurv << ", " << flagU1dir << ", " << flagU2dir
		<< ")" << std::endl;
#endif
#if JNA_DO_MONITOR
	++num_gbl_traj;
#endif
	
	return new GblTrajectory(ptr_array_to_vector(points, npoints), 
			flagCurv!=0, flagU1dir!=0, flagU2dir!=0);
}

/**
 * \brief construct new gbl::GblTrajectory with a seed matrix
 * \see gbl::GblTrajectory::GblTrajectory
 * \param [in] points array of pointers to gbl::GblPoint to put into trajectory
 * \param [in] npoints number of points in array
 * \param [in] aLabel integer label for seed
 * \param [in] seedArray double-array of length 25 listing seed matrix elements row-wise
 * \param [in] flagCurv use q/p - non-zero for true, zero for false
 * \param [in] flagU1dir use in u1 direction - non-zero for true, zero for false
 * \param [in] flagU2dir use in u2 direction - non-zero for true, zero for false
 * \return dynamically allocated trajectory wrapping input points
 */
GblTrajectory* GblTrajectoryCtorPtrArraySeed(GblPoint* points[], int npoints,
											 int aLabel, double seedArray[],
											 int flagCurv, int flagU1dir, int flagU2dir) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectoryCtorPtrArraySeed("
		<< points << ", " << npoints << ", "
		<< aLabel << ", " << seedArray << ", "
		<< flagCurv << ", " << flagU1dir << ", " << flagU2dir
		<< ")" << std::endl;
#endif
#if JNA_DO_MONITOR
	++num_gbl_traj;
#endif
	
	Map<Matrix5d> seed(seedArray,5,5);
	
	return new GblTrajectory(ptr_array_to_vector(points, npoints), aLabel, seed, 
			flagCurv!=0, flagU1dir!=0, flagU2dir!=0);
}

/**
 * \brief compose trajectory for a 2-body decay
 * \see gbl::GblTrajectory::GblTrajectory
 * \param [in] points_1 array of gbl::GblPoint for one track
 * \param [in] npoints_1 number of points
 * \param [in] trafo_1 double array listing elements of 2x3 track matrix
 * \param [in] points_2 array of gbl::GblPoint for one track
 * \param [in] npoints_2 number of points
 * \param [in] trafo_2 double array listing elements of 2x3 track matrix
 * \return dynamically created trajectory composed of two tracks
 */
GblTrajectory* GblTrajectoryCtorPtrComposed(GblPoint* points_1[], int npoints_1, double trafo_1[],
											GblPoint* points_2[], int npoints_2, double trafo_2[]) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectoryCtorPtrComposed("
		<< points_1 << ", " << npoints_1 << ", " << trafo_1 << ", "
		<< points_2 << ", " << npoints_2 << ", " << trafo_2 << ")"
		<< std::endl;
#endif
#if JNA_DO_MONITOR
	++num_gbl_traj;
#endif
	
	
	MatrixXd inner_1(2,3);
	inner_1(0,0)=trafo_1[0];
	inner_1(0,1)=trafo_1[1];
	inner_1(0,2)=trafo_1[2];

	inner_1(1,0)=trafo_1[3];
	inner_1(1,1)=trafo_1[4];
	inner_1(1,2)=trafo_1[5];

	std::pair<std::vector<GblPoint>, MatrixXd> track_trafo_1 = std::make_pair(ptr_array_to_vector(points_1, npoints_1), inner_1);
	
	//second track

	MatrixXd inner_2(2,3);
	inner_2(0,0)=trafo_2[0];
	inner_2(0,1)=trafo_2[1];
	inner_2(0,2)=trafo_2[2];

	inner_2(1,0)=trafo_2[3];
	inner_2(1,1)=trafo_2[4];
	inner_2(1,2)=trafo_2[5];

	std::pair<std::vector<GblPoint>, MatrixXd> track_trafo_2 = std::make_pair(ptr_array_to_vector(points_2, npoints_2), inner_2);
	
	return new GblTrajectory({track_trafo_1, track_trafo_2});

}

/**
 * \brief call gbl::GblTrajectory::fit on self
 * \param [in] self gbl::GblTrajectory to operate on
 * \param [out] Chi2 pointer to double where Chi2 result will be stored
 * \param [out] Ndf pointer to double where Ndf result will be stored
 * \param [out] lostWeight pointer to double where lostWeight result will be stored
 * \param [in] c_optionList C-style string listing options
 * \param [in] aLabel integer label for fit
 */
void GblTrajectory_fit(GblTrajectory* self, double* Chi2, int* Ndf, double* lostWeight, char* c_optionList, unsigned int aLabel) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_fit("
		<< self << ", " << Chi2 << ", " << Ndf << ", " << lostWeight << ", " << c_optionList << ", " << aLabel << ")"
		<< std::endl;
#endif
	
	std::string optionList(c_optionList);
	self->fit(*Chi2, *Ndf, *lostWeight, optionList,aLabel);
}

/**
 * \brief delete self
 *
 * Cleanup for JNA, handles cleaning up all GblPoints within it!
 */
void GblTrajectory_delete(GblTrajectory* self) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_delete(" << self << ")" << std::endl;
#endif
#if JNA_DO_MONITOR
	if (self) {
		--num_gbl_traj;
		num_gbl_point -= self->getNumPoints();
	}
#endif
	if (self) delete self;
}

/**
 * \brief call gbl::GblTrajectory::isValid
 * \param [in] self gbl::GblTrajectory to operate on
 * \return 0 if true, 1 if false
 */
int GblTrajectory_isValid(GblTrajectory* self) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_isValid(" << self << ")" << std::endl;
#endif
	return self->isValid() ? 0 : 1;
}

/**
 * \brief call gbl::GblTrajectory::getNumPoints
 * \param [in] self gbl::GblTrajectory to operate on
 * \return number of points
 */
int GblTrajectory_getNumPoints(GblTrajectory* self) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_getNumPoints(" << self << ")" << std::endl;
#endif
	return (int) self->getNumPoints();
}

/**
 * \brief call gbl::GblTrajectory::printTrajectory
 * \param [in] self gbl::GblTrajectory to operate on
 * \param [in] level integer level for printing
 */
void GblTrajectory_printTrajectory(GblTrajectory* self, int level) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_printTrajectory(" << self << ", " << level << ")" << std::endl;
#endif
	return self->printTrajectory(level);
}

/**
 * \brief call gbl::GblTrajectory::printData
 * \param [in] self gbl::GblTrajectory to operate on
 */
void GblTrajectory_printData(GblTrajectory* self) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_printData(" << self << ")" << std::endl;
#endif
	return self->printData();
}

/**
 * \brief call gbl::GblTrajectory::printPoints
 * \param [in] self gbl::GblTrajectory to operate on
 * \param [in] level integer level for printing
 */
void GblTrajectory_printPoints(GblTrajectory* self, int level) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_printPoints(" << self << ", " << level << ")" << std::endl;
#endif
	return self->printPoints(level);
}

/**
 * \brief get trajectory results
 * \note only supports 5-vector and 5x5 covariance matrix!
 * \param [in] self gbl::GblTrajectory to operate on
 * \param [in] aSignedLabel integer label for results
 * \param [out] localPar length 5 double array where local results will be stored
 * \param [out] nLocalPar integer length of double array (always set to 5)
 * \param [out] localCov length 25 double array where local covariance results will be stored
 * \param [out] sizeLocalCov integer dimension of covariance matrix (always set to 5)
 */
void GblTrajectory_getResults(GblTrajectory* self, int aSignedLabel, double* localPar, int* nLocalPar,
								double * localCov, int* sizeLocalCov) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_getResults(" << self 
		<< ", " << aSignedLabel << ", " << localPar << ", " << nLocalPar
		<< ", " << localCov << ", " << sizeLocalCov << ")" << std::endl;
#endif

	Eigen::VectorXd e_localPar(5);
	Eigen::MatrixXd e_localCov(5,5);
	
	self->getResults(aSignedLabel, e_localPar, e_localCov);

	//std::cout<<"gblTrajectoryWrapper::getResults"<<std::endl;
	//std::cout<<e_localPar<<std::endl;
	//std::cout<<e_localCov<<std::endl;
	
	Map<Vector5d>(localPar, 5) = e_localPar;
	Map<Matrix5d>(localCov, 5, 5) = e_localCov; 
	*nLocalPar = 5;
	*sizeLocalCov = 5;
}

/**
 * \brief get measurement results from trajectory
 * \see gbl::GblTrajectory::getMeasResults
 * \param [in] self gbl::GblTrajectory to operate on
 * \param [in] aLabel integer label of results
 * \param [out] numData length of resulting arrays (expect 2)
 * \param [out] aResiduals double array of residuals
 * \param [out] aMeasErrors double array of measurement errors
 * \param [out] aResErrors double array of residual errors
 * \param [out] aDownWeights double array of down weights
 */
int GblTrajectory_getMeasResults(GblTrajectory* self, int aLabel, int* numData, 
									double* aResiduals, double* aMeasErrors, double* aResErrors, 
									double* aDownWeights) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_getMeasResults(" << self 
		<< ", " << aLabel << ", " << numData << ", " << aResiduals
		<< ", " << aMeasErrors << ", " << aResErrors << ", " << aDownWeights << ")" << std::endl;
#endif
	
	Eigen::VectorXd e_aResiduals(2);
	Eigen::VectorXd e_aMeasErrors(2);
	Eigen::VectorXd e_aResErrors(2);
	Eigen::VectorXd e_aDownWeights(2);
	unsigned int num_data = 0;
	
	unsigned int out = self->getMeasResults(aLabel, num_data, e_aResiduals, e_aMeasErrors,
											e_aResErrors, e_aDownWeights);
#ifdef JNA_DEBUG
	std::cout << "  getMeasResults returned the status code " << out << std::endl;
#endif
	
	*numData = num_data;
	
	for (unsigned int i = 0; i < num_data; i++) {
		aResiduals[i] = e_aResiduals(i);
		aMeasErrors[i] = e_aMeasErrors(i);
		aResErrors[i] = e_aResErrors(i);
		aDownWeights[i] = e_aDownWeights(i);
	}

	return out;
}

/**
 * \brief write a trajectory to a gbl::MilleBinary file
 * \see gbl::GblTrajectory::milleOut
 * \param [in] self gbl::GblTrajectory to write out
 * \param [in] millebinary gbl::MilleBinary to write to
 */
void GblTrajectory_milleOut(GblTrajectory* self, MilleRecord* millebinary) {
#ifdef JNA_DEBUG
	std::cout << "GblTrajectory_milleOut(" << self 
		<< ", " << millebinary << ")" << std::endl;
#endif
	self->milleOut(millebinary);
}

/**
 * \brief construct a detector layer
 * \see gbl::GblDetectorLayer::GblDetectorLayer
 * \param [in] aName C-style string name
 * \param [in] aLayer integer ID for layer
 * \param [in] aDim dimension of layer
 * \param [in] thickness thickness of alyer
 * \param [in] aCenter double array of length 3 defining center of layer
 * \param [in] aResolution double array of length 2
 * \param [in] aPrecision double array of length 2
 * \param [in] measTrafo double array defining 3x3 matrix row wise
 * \param [in] alignTrafo double array defining 3x3 matrix row wise
 * \return newly constructed detector layer
 */
GblDetectorLayer* GblDetectorLayerCtor(const char* aName, int aLayer, int aDim, double thickness,
																			 double aCenter[], double aResolution[], double aPrecision[],
																			 double measTrafo[], double alignTrafo[]) {
#ifdef JNA_DEBUG
	std::cout << "GblDetectorLayerCtor(" << aName << ", " << aLayer
		<< ", " << aDim << ", " << thickness << ", " << aCenter << ", " << aResolution
		<< ", " << aPrecision << ", " << measTrafo << ", " << alignTrafo << ")" << std::endl;
#endif
#if JNA_DO_MONITOR
	++num_gbl_det_layer;
#endif

	// uses eigen's Map structure to decompose an array into our type of Vector/Matrix
	// need to then copy that object into a local variable since Map's only allow const
	// reference access
	Vector3d e_aCenter(Map<Vector3d>(aCenter,3));
	Vector2d e_aResolution(Map<Vector2d>(aResolution,2));
	Vector2d e_aPrecision(Map<Vector2d>(aPrecision,2));
	Matrix3d e_measTrafo(Map<Matrix3d>(measTrafo,3,3));
	Matrix3d e_alignTrafo(Map<Matrix3d>(alignTrafo,3,3));
	
	GblDetectorLayer* layer = new GblDetectorLayer(aName, aLayer, aDim, thickness, e_aCenter, e_aResolution, e_aPrecision, e_measTrafo, e_alignTrafo);
	return layer;
}

/**
 * \brief delete gbl::GblDetectorLayer
 */
void GblDetectorLayer_delete(GblDetectorLayer* self) {
#ifdef JNA_DEBUG
	std::cout << "GblDetectorLayer_delete(" << self << ")" << std::endl;
#endif
#if JNA_DO_MONITOR
	--num_gbl_det_layer;
#endif
	if (self) delete self;
}

/**
 * \brief call gbl::GblDetectorLayer::print
 * \param [in] self gbl::GblDetectorLayer to operate on
 */
void GblDetectorLayer_print(GblDetectorLayer* self) {
#ifdef JNA_DEBUG
	std::cout << "GblDetectorLayer_print(" << self << ")" << std::endl;
#endif
	self->print();
}

/**
 * \brief call gbl::GblDetectorLayer::getRadiationLength
 * \param [in] self gbl::GblDetectorLayer to operate on
 */
double GblDetectorLayer_getRadiationLength(GblDetectorLayer* self) {
#ifdef JNA_DEBUG
	std::cout << "GblDetectorLayer_getRadiationLength(" << self << ")" << std::endl;
#endif
	return self->getRadiationLength();
}

//TODO Implement wrappers for
//getResolution
//getPrecision
//getCenter
//getMeasSystemDirs
//getAlignSystemDirs


//Helix prediction on layer

/**
 * \brief construct a helix prediction from a layer
 * \param [in] self gbl::GblDetectorLayer to operate on
 * \param [in] hlx gbl::GblSimpleHelix to do helical calculations
 * \return newly constructed gbl::GblHelixPrediction
 */
GblHelixPrediction* GblDetectorLayer_intersectWithHelix(GblDetectorLayer* self, GblSimpleHelix* hlx) {
#ifdef JNA_DEBUG
	std::cout << "GblDetectorLayer_intersectWithHelix(" << self << ", " << hlx << ")" << std::endl;
#endif
#if JNA_DO_MONITOR
	++num_gbl_helix_prediction;
#endif
	
	Vector3d center = self->getCenter();
	Vector3d udir		= (self->getMeasSystemDirs()).row(0);
	Vector3d vdir		= (self->getMeasSystemDirs()).row(1);
	
	return new GblHelixPrediction(hlx->getPrediction(center, udir, vdir));
}

/**
 * \brief construct a new simple helix
 *
 * The input parameters are all the same as the C++ constructor
 * gbl::GblSimpleHelix::GblSimpleHelix.
 *
 * \return newly constructed gbl::GblSimpleHelix
 */
GblSimpleHelix* GblSimpleHelixCtor(double aRinv, double aPhi0, double aDca, double aDzds, double aZ0) {
#if JNA_DO_MONITOR
	++num_gbl_simple_helix;
#endif
	return new GblSimpleHelix(aRinv, aPhi0, aDca, aDzds, aZ0);
}

/**
 * \brief delete gbl::GblSimpleHelix
 */
void GblSimpleHelix_delete(GblSimpleHelix* self) {
#if JNA_DO_MONITOR
	--num_gbl_simple_helix;
#endif
	if (self) delete self;
}

/**
 * \brief call gbl::GblSimpleHelix::getPhi
 * \param [in] self gbl::GblSimpleHelix to operate on
 * \param [in] aRadius radius to calculate phi at
 * \return phi calculation
 */
double GblSimpleHelix_getPhi(GblSimpleHelix* self, double aRadius) {
	return self->getPhi(aRadius);
}

/**
 * \brief call gbl::GblSimpleHelix::getArcLengthR
 * \param [in] self gbl::GblSimpleHelix to operate on
 * \param [in] aRadius radius to calculate arc length of
 * \return arc length calculation
 */
double GblSimpleHelix_getArcLengthR(GblSimpleHelix* self, double aRadius) {
	return self->getArcLengthR(aRadius);
}

/**
 * \brief call gbl::GblSimpleHelix::getArcLengthXY
 * \param [in] self gbl::GblSimpleHelix to operate on
 * \param [in] xPos x-position 
 * \param [in] yPos y-position
 * \return arc length calculation
 */
double GblSimpleHelix_getArcLengthXY(GblSimpleHelix* self, double xPos, double yPos) {
	return self->getArcLengthXY(xPos, yPos);
}

/**
 * \brief call gbl::GblSimpleHelix::moveToXY
 * \param [in] self gbl::GblSimpleHelix to operate on
 * \param [in] xPos x-position 
 * \param [in] yPos y-position
 * \param [out] newPhi0 double address to store resulting phi
 * \param [out] newDca double address to store resulting Dca
 * \param [out] newZ0 double address to store resulting Z0
 */
void GblSimpleHelix_moveToXY(GblSimpleHelix* self, double xPos, double yPos,
														 double* newPhi0, double* newDca, double* newZ0) {
	self->moveToXY(xPos, yPos,
								 *newPhi0, *newDca, *newZ0);
}

/**
 * \brief get a helical prediction from a reference coordinate system
 * \see gbl::GblSimpleHelix::getPrediction
 * \param [in] self gbl::GblSimpleHelix to operate on
 * \param [in] refPos double array of length three containing reference position
 * \param [in] uDir 3-length double array defining u direction
 * \param [in] vDir 3-length double array defining v direction
 * \return new gbl::GblHelixPrediction from this coordinate system
 */
GblHelixPrediction* GblSimpleHelix_getPrediction(GblSimpleHelix* self, double refPos[], double uDir[], double vDir[]) {
	
	Map<Vector3d> e_refPos(refPos,3);
	Map<Vector3d> e_uDir(uDir,3);
	Map<Vector3d> e_vDir(vDir,3);
	GblHelixPrediction prediction = self->getPrediction(e_refPos,e_uDir,e_vDir);
	
	/*std::cout<<"Cross Check GBL predicted position!"<<std::endl;
	std::cout<<prediction->getPosition()<<std::endl;
	
	std::cout<<"Cross Check GBL meas predicted position!"<<std::endl;
	std::cout<<prediction->getMeasPred()<<std::endl;*/
	
#if JNA_DO_MONITOR
	++num_gbl_helix_prediction;
#endif
	// JNA only deals with pointers so we need to dynamically create a new copy
	return new GblHelixPrediction(prediction);
}

/**
 * \brief create a new helix prediction manually
 * \param [in] sArc arc length
 * \param [in] aPred length-2 double array predicted measurement
 * \param [in] tDir length-3 double array defining t direction
 * \param [in] uDir length-3 double array defining u direction
 * \param [in] vDir length-3 double array defining v direction
 * \param [in] nDir length-3 double array defining n direction
 * \param [in] aPos length-3 double array defining position
 */
GblHelixPrediction* GblHelixPredictionCtor(double sArc, double aPred[], double tDir[], double uDir[], double vDir[],
																				 double nDir[], double aPos[]) {
	Map<Vector2d> e_aPred(aPred,2);
	Map<Vector3d> e_tDir(tDir,3);
	Map<Vector3d> e_uDir(uDir,3);
	Map<Vector3d> e_vDir(vDir,3);
	Map<Vector3d> e_nDir(nDir,3);
	Map<Vector3d> e_aPos(aPos,3);
	
	
#if JNA_DO_MONITOR
	++num_gbl_helix_prediction;
#endif
	return new GblHelixPrediction(sArc, e_aPred, e_tDir, e_uDir, e_vDir, 
																e_nDir, e_aPos);
}

/**
 * \brief delete a gbl::GblHelixPrediction
 * \param [in] self gbl::GblHelixPrediction to delete
 */
void GblHelixPrediction_delete(GblHelixPrediction* self) {
#if JNA_DO_MONITOR
	--num_gbl_helix_prediction;
#endif
	if (self) delete self;
}

/**
 * \brief get the arc length of a helix prediction
 * \see gbl::GblHelixPrediction::getArcLength
 * \param [in] self gbl::GblHelixPrediction to operate on
 * \return arc length
 */
double GblHelixPrediction_getArcLength(GblHelixPrediction* self) {
	return self->getArcLength();
}

/**
 * \brief get the predicted measurement
 * \see gbl::GblHelixPrediction::getMeadPred
 * \param [in] self gbl::GblHelixPrediction to operate on
 * \param [out] prediction length-2 double array that will hold predicted measurement
 */
void GblHelixPrediction_getMeasPred(GblHelixPrediction* self, double* prediction) {
	Vector2d e_pred = self->getMeasPred();
	
	prediction[0] = e_pred(0);
	prediction[1] = e_pred(1);
}

/**
 * \brief get the position
 * \see gbl::GblHelixPrediction::getPosition
 * \param [in] self gbl::GblHelixPrediction to operate on
 * \param [out] position length-3 double array that will hold position
 */
void GblHelixPrediction_getPosition(GblHelixPrediction* self, double* position) {
	Vector3d e_pos = self->getPosition();
	
	position[0] = e_pos(0);
	position[1] = e_pos(1);
	position[2] = e_pos(2);
}

/**
 * \brief get the direction
 * \see gbl::GblHelixPrediction::getDirection
 * \param [in] self gbl::GblHelixPrediction to operate on
 * \param [out] direction length-3 double array that will hold direction
 */
void GblHelixPrediction_getDirection(GblHelixPrediction* self, double direction[]) {
	Vector3d e_dir = self->getDirection();
	
	direction[0] = e_dir(0);
	direction[1] = e_dir(1);
	direction[2] = e_dir(2);
}

/**
 * \brief get the cosine incidence
 * \see gbl::GblHelixPrediction::getCosIncidence
 * \param [in] self gbl::GblHelixPrediction to operate on
 * \return value of cosine incidence
 */
double GblHelixPrediction_getCosIncidence(GblHelixPrediction* self) {
	return self->getCosIncidence();
}

/**
 * \brief get the curvilinear directions
 * \see gbl::GblHelixPrediction::getCurvilinearDirs
 * \param [in] self gbl::GblHelixPrediction to operate on
 * \param [out] curvilinear length-6 double array that will hold the two direction vectors
 */
void GblHelixPrediction_getCurvilinearDirs(GblHelixPrediction* self, double curvilinear[]) {
	Matrix<double,2,3> curDirs = self->getCurvilinearDirs();
	
	curvilinear[0] = curDirs(0,0);
	curvilinear[1] = curDirs(0,1);
	curvilinear[2] = curDirs(0,2);
	curvilinear[3] = curDirs(1,0);
	curvilinear[4] = curDirs(1,1);
	curvilinear[5] = curDirs(1,2);
}
}

