/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressZChecker.cpp
 *
 *  Created on: Jul 9, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "CompressZChecker.h"

namespace adios2
{
namespace core
{
namespace compress
{

CompressZChecker::CompressZChecker(const Params &parameters,
                                   const bool debugMode, Operator &op)
: Operator("z-checker", parameters, debugMode), m_Operator(op)
{
    printf("z-checker ... CompressZChecker\n");
}

size_t CompressZChecker::Compress(const void *dataIn, const Dims &dimensions,
                                  const size_t elementSize,
                                  const std::string type, void *bufferOut,
                                  const Params &parameters, Params &info) const
{
    // Init Z-Checker parameters (maybe a separate function) similar as in SZ

    printf("z-checker ... Compress\n");

    bool isFinished = false;
    while (isFinished == false)
    {
        // SZ,  compression
        const size_t sizeOut = m_Operator.Compress(
            dataIn, dimensions, elementSize, type, bufferOut, parameters, info);

        // deflated buffer

        //        m_Operator.Decompress()
        //
        //            const void *bufferIn,
        //            const size_t sizeIn, void *dataOut, const Dims
        //            &dimensions, const std::string type, const Params
        //            &parameters

        // do some checks or make decisions based on parameters
        //        if (counter == 0)
        //        {
        //            isFinished = true;
        //        }
    }

    return 0;
}

}
}
}
