/*
 * Copyright (c) 2011, The University of Oxford
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the University of Oxford nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */


#include "utility/oskar_work_init.h"
#include "utility/oskar_mem_init.h"
#include "stdlib.h"

#ifdef __cplusplus
extern "C" {
#endif

int oskar_work_init(oskar_Work* work, const int type, const int location)
{
    int error = 0;

    if (work == NULL)
        return OSKAR_ERR_INVALID_ARGUMENT;

    if (type != OSKAR_SINGLE && type != OSKAR_DOUBLE)
        return OSKAR_ERR_BAD_DATA_TYPE;

    if (location != OSKAR_LOCATION_CPU && location != OSKAR_LOCATION_GPU)
        return OSKAR_ERR_BAD_LOCATION;

    error = oskar_mem_init(&work->real, type, location, 0, 1);
    if (error) return error;
    error = oskar_mem_init(&work->complex, (type | OSKAR_COMPLEX), location, 0, 1);
    if (error) return error;
    error = oskar_mem_init(&work->matrix, (type | OSKAR_COMPLEX | OSKAR_MATRIX),
            location, 0, 1);
    if (error) return error;

    return error;
}


#ifdef __cplusplus
}
#endif
