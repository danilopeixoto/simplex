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

#include <simplex/simplex.h>

#include <iostream>

SIMPLEX_NAMESPACE_USING

typedef Simplex<double> DoubleSimplex;

int main(int argc, char ** argv) {
    DoubleSimplex::Array1DType objectiveCoefficients = { 5.0, 7.0 };
    DoubleSimplex::Array2DType restrictionCoefficients = { { 3.0,  0.0 },
                                                           { 0.0,  1.5 },
                                                           { 0.25, 0.5 } };
    DoubleSimplex::Array1DType restrictionLimits = { 250.0,
                                                     100.0,
                                                     50.0 };
    
    DoubleSimplex simplex(objectiveCoefficients,
            0,
            restrictionCoefficients,
            restrictionLimits,
            DoubleSimplex::Maximize,
            10);
    
    const DoubleSimplex::Array2DType & optimizedMatrix = simplex.solve();
    
    std::cout << "Optimized Matrix:" << std::endl;
    std::cout << optimizedMatrix << std::endl;
    
    return 0;
}
