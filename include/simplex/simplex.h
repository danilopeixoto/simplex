// Copyright (c) 2019, Danilo Peixoto. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef SIMPLEX_SIMPLEX_H
#define SIMPLEX_SIMPLEX_H

#define SIMPLEX_NAME "Simplex"
#define SIMPLEX_VERSION "1.0.0"
#define SIMPLEX_AUTHOR "Danilo Peixoto"
#define SIMPLEX_COPYRIGHT "Copyright (c) 2019, Danilo Peixoto. All rights reserved."
#define SIMPLEX_LICENSE "BSD-3-Clause"

#define SIMPLEX_VERSION_MAJOR 1
#define SIMPLEX_VERSION_MINOR 0
#define SIMPLEX_VERSION_PATCH 0

#define SIMPLEX_NAMESPACE_BEGIN namespace simplex {
#define SIMPLEX_NAMESPACE_END };
#define SIMPLEX_NAMESPACE_USING using namespace simplex;

#include <vector>
#include <limits>
#include <ostream>

SIMPLEX_NAMESPACE_BEGIN

template<typename T>
using Array1D = std::vector<T>;
template<typename T>
using Array2D = Array1D<Array1D<T>>;

template<typename T>
std::ostream & operator <<(std::ostream & ostream, const Array1D<T> & array) {
    ostream << '[';
    
    size_t count = array.size();
    
    for (size_t i = 0; i < count; i++) {
        ostream << array[i];
        
        if (i < count - 1)
           ostream << ", ";
    }
    
    return ostream << ']';
}

template<typename T>
std::ostream & operator <<(std::ostream & ostream, const Array2D<T> & array) {
    ostream << '[';
    
    size_t count = array.size();
    
    for (size_t i = 0; i < count; i++) {
        ostream << array[i];
        
        if (i < count - 1)
            ostream << ", " << std::endl;
    }
    
    return ostream << ']';
}

template<typename T>
class Simplex {
public:
    typedef T Type;
    typedef Array1D<T> Array1DType;
    typedef Array2D<T> Array2DType;
    
    enum OptimizationType {
        Minimize = 0,
        Maximize
    };
    
    Simplex(const Array1DType & objectiveCoefficients,
            T objectiveConstant,
            const Array2DType & restrictionCoefficients,
            const Array1DType & restrictionLimits,
            OptimizationType optimizationType,
            size_t maximumIteration = std::numeric_limits<size_t>::max()) {
        this->objectiveCoefficients = objectiveCoefficients;
        this->objectiveConstant = objectiveConstant;
        this->restrictionCoefficients = restrictionCoefficients;
        this->restrictionLimits = restrictionLimits;
        this->optimizationType = optimizationType;
        this->maximumIteration = maximumIteration;
        
        reset();
    }
    ~Simplex() {}
    
    const Array1DType & getObjectiveCoefficients() const {
        return objectiveCoefficients;
    }
    const Type getObjectiveConstant() const {
        return objectiveConstant;
    }
    const Array2DType & getRestrictionCoefficients() const {
        return restrictionCoefficients;
    }
    const Array1DType & getRestrictionLimits() const {
        return restrictionLimits;
    }
    OptimizationType getOptimizationType() const {
        return optimizationType;
    }
    size_t getMaximumIteration() const {
        return maximumIteration;
    }
    
    size_t getIteration() const {
        return iteration;
    }
    size_t getRowCount() const {
        return rowCount;
    }
    size_t getColumnCount() const {
        return columnCount;
    }
    size_t getVariableCount() const {
        return variableCount;
    }
    size_t getRestrictionCount() const {
        return restrictionCount;
    }
    bool isOptimized() const {
        for (size_t i = 0; i < columnCount; i++) {
            if (optimizedMatrix[0][i] < (Type)0.0)
               return false;
        }
        
        return true;
    }
    const Array2DType & getOptimizedMatrix() const {
        return optimizedMatrix;
    }
    
    Simplex & reset() {
        iteration = 0;
        
        size_t variableCount = objectiveCoefficients.size();
        size_t restrictionCount = restrictionCoefficients.size();
        
        rowCount = restrictionCount + 1;
        columnCount = variableCount + rowCount + 1;
        
        size_t lastColumnIndex = columnCount - 1;
        size_t restrictionBeginIndex = variableCount + 1;
        
        optimizedMatrix.resize(rowCount);
        
        for (size_t i = 0; i < rowCount; i++) {
            Array1DType & row = optimizedMatrix[i];
            row.resize(columnCount);
        }
        
        for (size_t i = 0; i < columnCount; i++) {
            if (i == 0)
                optimizedMatrix[0][i] = (Type)1.0;
            else if (i > 0 && i <= variableCount)
                optimizedMatrix[0][i] = -objectiveCoefficients[i - 1];
            else if (i > variableCount && i < lastColumnIndex)
                optimizedMatrix[0][i] = (Type)0.0;
            else
                optimizedMatrix[0][i] = objectiveConstant;
        }
        
        for (size_t i = 1; i <= restrictionCount; i++) {
            optimizedMatrix[i][0] = (Type)0.0;
            optimizedMatrix[i][lastColumnIndex] = restrictionLimits[i - 1];
        }
        
        for (size_t i = 1; i <= restrictionCount; i++) {
            for (size_t j = 1; j <= variableCount; j++)
                optimizedMatrix[i][j] = restrictionCoefficients[i - 1][j - 1];
        }
        
        for (size_t i = 1; i <= restrictionCount; i++) {
            for (size_t j = restrictionBeginIndex; j < lastColumnIndex; j++) {
                if (i - 1 == j - restrictionBeginIndex)
                    optimizedMatrix[i][j] = (Type)1.0;
                else
                    optimizedMatrix[i][j] = (Type)0.0;
            }
        }
        
        return *this;
    }
    const Array2DType & next() {
        iteration++;
    }
    const Array2DType & solve() {
        reset();
        
        while (iteration < maximumIteration && !isOptimized())
            next();
        
        return optimizedMatrix;
    }
    
private:
    Array1DType objectiveCoefficients;
    Type objectiveConstant;
    Array2DType restrictionCoefficients;
    Array1DType restrictionLimits;
    OptimizationType optimizationType;
    size_t maximumIteration;
    
    size_t iteration;
    size_t rowCount;
    size_t columnCount;
    size_t variableCount;
    size_t restrictionCount;
    Array2DType optimizedMatrix;
};

SIMPLEX_NAMESPACE_END

#endif
