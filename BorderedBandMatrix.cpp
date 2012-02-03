/*
 * BorderedBandMatrix.cpp
 *
 *  Created on: Aug 14, 2011
 *      Author: kleinwrt
 */

#include "BorderedBandMatrix.h"

/// Create bordered band matrix.
BorderedBandMatrix::BorderedBandMatrix() {

}

BorderedBandMatrix::~BorderedBandMatrix() {
	// TODO Auto-generated destructor stub
}

/// Resize bordered band matrix.
/**
 * \param nSize [in] Size of matrix
 * \param nBorder [in] Size of border (=1 for q/p + additional local parameters)
 * \param nBand [in] Band width (usually = 5, for simplified jacobians = 4)
 */
void BorderedBandMatrix::resize(unsigned int nSize, unsigned int nBorder,
		unsigned int nBand) {
	numSize = nSize;
	numBorder = nBorder;
	numCol = nSize - nBorder;
	numBand = 0;
	theBorder.ResizeTo(nBorder, nBorder);
	theBorder.Zero();
	theMixed.ResizeTo(nBorder, nSize - nBorder);
	theMixed.Zero();
	theBand.ResizeTo(nBand + 1, nSize - nBorder);
	theBand.Zero();
}

/// Add symmetric block matrix.
/**
 * Add (extended) block matrix defined by 'aVector * aWeight * aVector.T'
 * to bordered band matrix:
 * BBmatrix(anIndex(i),anIndex(j)) += aVector(i) * aWeight * aVector(j).
 * \param aWeight [in] Weight
 * \param anIndex [in] List of rows/colums to be used
 * \param aVector [in] Vector
 */
void BorderedBandMatrix::addBlockMatrix(double aWeight,
		const std::vector<unsigned int>* anIndex,
		const std::vector<double>* aVector) {
	int nBorder = numBorder;
	for (unsigned int i = 0; i < anIndex->size(); i++) {
		int iIndex = (*anIndex)[i] - 1; // anIndex has to be sorted
		for (unsigned int j = 0; j <= i; j++) {
			int jIndex = (*anIndex)[j] - 1;
			if (iIndex < nBorder) {
				theBorder(iIndex, jIndex) += (*aVector)[i] * aWeight
						* (*aVector)[j];
			} else if (jIndex < nBorder) {
				theMixed(jIndex, iIndex - nBorder) += (*aVector)[i] * aWeight
						* (*aVector)[j];
			} else {
				unsigned int nBand = iIndex - jIndex;
				theBand(nBand, jIndex - nBorder) += (*aVector)[i] * aWeight
						* (*aVector)[j];
				numBand = std::max(numBand, nBand); // update band width
			}
		}
	}
}

/// Retrieve symmetric block matrix.
/**
 * Get (compressed) block from bordered band matrix: aMatrix(i,j) = BBmatrix(anIndex(i),anIndex(j)).
 * \param anIndex [in] List of rows/colums to be used
 */
TMatrixDSym BorderedBandMatrix::getBlockMatrix(
		std::vector<unsigned int> anIndex) {

	TMatrixDSym aMatrix(anIndex.size());
	int nBorder = numBorder;
	for (unsigned int i = 0; i < anIndex.size(); i++) {
		int iIndex = anIndex[i] - 1; // anIndex has to be sorted
		for (unsigned int j = 0; j <= i; j++) {
			int jIndex = anIndex[j] - 1;
			if (iIndex < nBorder) {
				aMatrix(i, j) = theBorder(iIndex, jIndex);
			} else if (jIndex < nBorder) {
				aMatrix(i, j) = theMixed(jIndex, iIndex - nBorder);
			} else {
				unsigned int nBand = iIndex - jIndex;
				aMatrix(i, j) = theBand(nBand, jIndex - nBorder);
			}
		}
	}
	return aMatrix;
}

/// Solve linear equation system, partially calculate inverse.
/**
 * Solve linear equation A*x=b system with bordered band matrix A,
 * calculate bordered band part of inverse of A. Use decomposition
 * in border and band part for block matrix algebra.
 * \param [in] aRightHandSide Right hand side (vector) 'b' of A*x=b
 * \param [out] aSolution Solution (vector) x of A*x=b
 */
void BorderedBandMatrix::solveAndInvertBorderedBand(
		const TVectorD &aRightHandSide, TVectorD &aSolution) {

	int nBorder = numBorder;
	int nSize = numSize;
	// decompose band
	decomposeBand(); // TODO: check for positive definiteness
	// invert band
	TMatrixD inverseBand(invertBand());
	if (nBorder > 0) {
		// solve for mixed part
		TMatrixD auxMat = solveBand(theMixed);
		// solve for border part
		TMatrixDSym inverseBorder(nBorder);
		TMatrixD auxBorder(nBorder, nBorder), auxBorderT(nBorder, nBorder);
		TVectorD auxVec = aRightHandSide.GetSub(0, nBorder - 1)
				- auxMat * aRightHandSide.GetSub(nBorder, nSize - 1);
		auxBorder = (theBorder - theMixed * auxMat.T()) * 0.5;
		auxBorderT.Transpose(auxBorder);
		inverseBorder.SetSub(0, auxBorder + auxBorderT);
		inverseBorder.Invert();
		TVectorD borderSolution = inverseBorder * auxVec;
		// solve for band part
		TVectorD bandSolution = solveBand(
				aRightHandSide.GetSub(nBorder, nSize - 1));
		aSolution.SetSub(0, borderSolution);
		aSolution.SetSub(nBorder, bandSolution - auxMat * borderSolution);
		// parts of inverse
		theBorder = inverseBorder;
		theMixed = (inverseBorder * auxMat.T()) * -1.;
		theBand.SetSub(0, 0,
				inverseBand + bandOfAVAT(auxMat.T(), inverseBorder));
	} else {
		aSolution = solveBand(aRightHandSide);
		theBand.SetSub(0, 0, inverseBand);
	}
}

/// Print bordered band matrix.
void BorderedBandMatrix::printMatrix() {
	std::cout << "Border part: " << theBorder.GetNrows() << std::endl;
	theBorder.Print();
	std::cout << "Mixed  part: " << theMixed.GetNrows() << ","
			<< theMixed.GetNcols() << std::endl;
	theMixed.Print();
	std::cout << "Band   part: " << theBand.GetNrows() << ","
			<< theBand.GetNcols() << std::endl;
	theBand.Print();
}

/*============================================================================
 from Dbandmatrix.F (MillePede-II by V. Blobel, Univ. Hamburg)
 ============================================================================*/
/// (root free) Cholesky decomposition of band part: C=LDL^T
/**
 * Decompose band matrix into diagonal matrix D and lower triangular band matrix
 * L (diagonal=1). Overwrite band matrix with D and off-diagonal part of L.
 */
void BorderedBandMatrix::decomposeBand() {

	int nRow = numBand + 1;
	int nCol = numCol;
	TVectorD auxVec(nCol);
	for (int i = 0; i < nCol; i++) {
		auxVec(i) = theBand(0, i) * 16.0; // save diagonal elements
	}
	for (int i = 0; i < nCol; i++) {
		if ((theBand(0, i) + auxVec(i)) != theBand(0, i)) {
			theBand(0, i) = 1.0 / theBand(0, i);
		} else {
			theBand(0, i) = 0.0;
		}
		for (int j = 1; j < std::min(nRow, nCol - i); j++) {
			double rxw = theBand(j, i) * theBand(0, i);
			for (int k = 0; k < std::min(nRow, nCol - i) - j; k++) {
				theBand(k, i + j) -= theBand(k + j, i) * rxw;
			}
			theBand(j, i) = rxw;
		}
	}
}

/// Solve for band part.
/**
 * Solve C*x=b for band part using decomposition C=LDL^T
 * and forward (L*z=b) and backward substitution (L^T*x=D^-1*z).
 * \param [in] aRightHandSide Right hand side (vector) 'b' of C*x=b
 * \return Solution (vector) 'x' of C*x=b
 */
TVectorD BorderedBandMatrix::solveBand(const TVectorD &aRightHandSide) {

	int nRow = numBand + 1;
	int nCol = numCol;
	TVectorD aSolution(aRightHandSide);
	for (int i = 0; i < nCol; i++) // forward substitution
			{
		for (int j = 1; j < std::min(nRow, nCol - i); j++) {
			aSolution(j + i) -= theBand(j, i) * aSolution(i);
		}
	}
	for (int i = nCol - 1; i >= 0; i--) // backward substitution
			{
		double rxw = theBand(0, i) * aSolution(i);
		for (int j = 1; j < std::min(nRow, nCol - i); j++) {
			rxw -= theBand(j, i) * aSolution(j + i);
		}
		aSolution(i) = rxw;
	}
	return aSolution;
}

/// solve band part for mixed part (border rows).
/**
 * Solve C*X=B for mixed part using decomposition C=LDL^T
 * and forward and backward substitution.
 * \param [in] aRightHandSide Right hand side (matrix) 'B' of C*X=B
 * \return Solution (matrix) 'X' of C*X=B
 */
TMatrixD BorderedBandMatrix::solveBand(const TMatrixD &aRightHandSide) {

	int nRow = numBand + 1;
	int nCol = numCol;
	TMatrixD aSolution(aRightHandSide);
	for (unsigned int iBorder = 0; iBorder < numBorder; iBorder++) {
		for (int i = 0; i < nCol; i++) // forward substitution
				{
			for (int j = 1; j < std::min(nRow, nCol - i); j++) {
				aSolution(iBorder, j + i) -= theBand(j, i)
						* aSolution(iBorder, i);
			}
		}
		for (int i = nCol - 1; i >= 0; i--) // backward substitution
				{
			double rxw = theBand(0, i) * aSolution(iBorder, i);
			for (int j = 1; j < std::min(nRow, nCol - i); j++) {
				rxw -= theBand(j, i) * aSolution(iBorder, j + i);
			}
			aSolution(iBorder, i) = rxw;
		}
	}
	return aSolution;
}

/// Invert band part.
/**
 * \return Inverted band
 */
TMatrixD BorderedBandMatrix::invertBand() {

	int nRow = numBand + 1;
	int nCol = numCol;
	TMatrixD inverseBand(nRow, nCol);
	inverseBand.Zero();

	for (int i = nCol - 1; i >= 0; i--) {
		double rxw = theBand(0, i);
		for (int j = i; j >= std::max(0, i - nRow + 1); j--) {
			for (int k = j + 1; k < std::min(nCol, j + nRow); k++) {
				rxw -= inverseBand(abs(i - k), std::min(i, k))
						* theBand(k - j, j);
			}
			inverseBand(i - j, j) = rxw;
			rxw = 0.;
		}
	}
	return inverseBand;
}

/// Calculate band part of: 'anArray * aSymArray * anArray.T'.
/**
 * \return Band part of product
 */
TMatrixD BorderedBandMatrix::bandOfAVAT(const TMatrixD &anArray,
		const TMatrixDSym &aSymArray)

		{
	int nBand = numBand;
	int nCol = numCol;
	int nBorder = numBorder;
	double sum;
	TMatrixD aBand(nBand + 1, nCol);
	aBand.Zero();
	for (int i = 0; i < nCol; i++) {
		for (int j = std::max(0, i - nBand); j <= i; j++) {
			sum = 0.;
			for (int l = 0; l < nBorder; l++) {
				for (int k = 0; k < nBorder; k++) {
					sum += anArray[i][l] * aSymArray[l][k] * anArray[j][k];
				}
			}
			aBand[i - j][j] = sum;
			/*			aBand[i - j][j] =
			 ((anArray.GetSub(i, i, 0, nBorder - 1) * aSymArray)
			 * anArray.GetSub(j, j, 0, nBorder - 1))[0][0];*/
		}
	}
	return aBand;
}

