/*
 * Copyright (c) 2012, The University of Oxford
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

#ifndef OSKAR_STATION_MODEL_SAVE_CONFIGURATION_H_
#define OSKAR_STATION_MODEL_SAVE_CONFIGURATION_H_

/**
 * @file oskar_station_model_save_configuration.h
 */

#include "oskar_global.h"
#include "station/oskar_StationModel.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Writes OSKAR station model data to an ASCII file.
 *
 * @details
 * Writes data from the station model to an ASCII file
 * consisting of a simple header describing the number of (antenna) elements
 * in the station, the station longitude, latitude and altitude followed by
 * a CSV list of the local horizontal x,y,z positions of the elements
 * (in metres) and remaining station data.
 *
 * Note:
 * - The oskar_Mem pointers holding the coordinates must reside on host (CPU).
 * - The coordinates of the station file must be in metres.
 *
 * @param[in] filename Pathname of file to write.
 * @param[in] station  Station model to write.
 *
 * @return An OSKAR error code.
 */
OSKAR_EXPORT
int oskar_station_model_save_config(const char* filename,
        const oskar_StationModel* station);

#ifdef __cplusplus
}
#endif

#endif /* OSKAR_STATION_MODEL_SAVE_CONFIGURATION_H_ */